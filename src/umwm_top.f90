module umwm_top

  implicit none

contains

  subroutine umwm_initialize(config, spectrum, grid, forcing)

    use umwm_config, only: config_type
    use umwm_forcing, only: forcing_type
    use umwm_grid, only: grid_type
    use umwm_init,  only: alloc, init
    use umwm_io,    only: output_grid
    use umwm_spectrum, only: spectrum_type
    use umwm_stokes,only: stokes_drift

    type(config_type), intent(in) :: config
    type(spectrum_type), intent(in) :: spectrum
    type(grid_type), intent(in) :: grid
    type(forcing_type), intent(inout) :: forcing

    call alloc(grid, spectrum) ! allocate unrolled arrays
    call output_grid(config, grid)    ! output a grid file
    call forcing % load(config, config % starttimestr, grid) ! read initial fields
    call init(config, spectrum, grid, forcing) ! initialize model variables
    call stokes_drift(spectrum, config, grid, 'init') ! initialize stokes drift arrays

  end subroutine umwm_initialize


  subroutine umwm_run(config, spectrum, grid, forcing)

#ifdef MPI
    use umwm_mpi, only: exchange_halo
    use umwm_module, only: ierr
    use mpi
#endif

    use umwm_config, only: config_type
    use umwm_module, only: cd, currenttime, dts, e, ef, f, first, &
                           firstdtg, iip, nproc, &
                           nproc_plot, oc, starttime, stoptime, sumt
    use umwm_forcing, only: forcing_type
    use umwm_grid, only: grid_type
    use umwm_physics, only: source, diag
    use umwm_advection,only: propagation, refraction
    use umwm_io, only: output_grid_nc, output_spectrum_nc
    use umwm_restart, only: restart_read, restart_write
    use umwm_spectrum, only: spectrum_type
    use umwm_stokes, only: stokes_drift
    use umwm_util, only: sigwaveheight, meanwaveperiod
    use umwm_stress, only: stress

    use umwm_source_functions, only: sin_d12, sds_d12, snl_d12, s_ice

    use datetime_module

    type(config_type), intent(in) :: config
    type(spectrum_type), intent(in) :: spectrum
    type(grid_type), intent(in) :: grid
    type(forcing_type), intent(inout) :: forcing

    character(19) :: currenttimestr
    logical :: fullhour

    ! convert start and stop time strings to datetime objects:
    starttime = strptime(config % starttimestr,'%Y-%m-%d %H:%M:%S')
    stoptime  = strptime(config % stoptimestr, '%Y-%m-%d %H:%M:%S')

    currenttime = starttime
    currenttimestr = trim(currenttime % strftime('%Y-%m-%d_%H:%M:%S'))

    ! read wave spectrum field from a restart file if necessary:
    if (first .and. config % restart) call restart_read(config % starttimestr, spectrum, grid)

    do while (currenttime < stoptime) ! outer time loop

      ! report current and next checkpoint time:
      if(nproc == 0)write(*,fmt='(a)')&
      'umwm: solver: current time is:     '//currenttimestr

      currenttime = currenttime + timedelta(seconds=nint(config % dtg))
      currenttimestr = trim(currenttime % strftime('%Y-%m-%d_%H:%M:%S'))

      if(nproc == 0)write(*,fmt='(a)')&
        'umwm: solver: integrating to time: '//currenttimestr

      ! read atmosphere and ocean input fields from file.
      ! in case of esmf coupling, fields are assumed to be updated
      ! externally, and this call is not used.
#ifndef ESMF
      call forcing % update(config, currenttimestr, grid)
#endif

      ! inner time loop: global time step
      sumt = 0
      do while (sumt < config % dtg)

#ifndef ESMF
        call forcing % interpolate(config, grid) ! interpolate force fields in time
#endif

        call sin_d12(config, spectrum, grid, forcing) ! compute source input term Sin
        call sds_d12(config, spectrum, grid) ! compute source dissipation term Sds
        call snl_d12(config, spectrum, grid, forcing) ! compute non-linear source term Snl
        call s_ice(config, spectrum, grid, forcing)   ! compute sea ice attenuation term Sice
        call source(config, spectrum, grid)  ! integrate source functions

#ifdef MPI
        call exchange_halo(config, spectrum, grid) ! exchange halo points
#endif

        call propagation(config, spectrum, grid, forcing) ! compute advection and integrate

#ifdef ESMF
        e(:,:,grid % istart:grid % iend) = ef(:,:,grid % istart:grid % iend) ! update
#endif

        call refraction(config, spectrum, grid, forcing)    ! compute refraction and integrate

        e(:,:,grid % istart:grid % iend) = ef(:,:,grid % istart:grid % iend) ! update

        call stress(config, 'atm', spectrum, grid, forcing) ! compute wind stress and drag coefficient

#ifdef ESMF
        call stress(config, 'ocn', spectrum, grid, forcing) ! compute stress into ocean top and bottom
#endif

        if (first) then

          ! diagnostic calculations before output
          call diag(spectrum, grid, forcing)

          if (config % outgrid > 0) call output_grid_nc(config, config % starttimestr, spectrum, grid, forcing)
          if (config % outspec > 0) call output_spectrum_nc(config, config % starttimestr, spectrum, grid, forcing)

#ifdef MPI
          call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

          if(nproc == 0)then
            ! diagnostic output header on screen
            write(*,fmt='(a)')'dtg frac   tstep [s]  wspd [m/s] wdir [rad]   '&
                            //'swh [m]    mwp [s]    cd*10^3   fprog [hz]'
          end if

          first = .false.

        end if

        ! diagnostic output on screen
        if(nproc == nproc_plot)then
          write(*,fmt=100)sumt/config % dtg,dts,forcing % wspd(iip),forcing % wdir(iip),      &
                          sigwaveheight(iip, spectrum),meanwaveperiod(iip, spectrum),&
                          cd(iip)*1e3,f(oc(iip))
        end if

      end do ! end while(sumt<dtg) loop

      if (firstdtg) firstdtg = .false.
      if (config % stokes) call stokes_drift(spectrum, config, grid)

      call diag(spectrum, grid, forcing) ! model diagnostics for output

      fullhour = currenttime % getminute() == 0 &
           .and. currenttime % getsecond() == 0

      ! gridded output
      if(config % outgrid > 0)then
        if(mod(currenttime % gethour(),config % outgrid) == 0 .and. fullhour)then
#ifndef ESMF
          call stress(config, 'ocn', spectrum, grid, forcing)
#endif
          call output_grid_nc(config, currenttimestr, spectrum, grid, forcing)
        end if
      elseif(config % outgrid == -1)then
#ifndef ESMF
        call stress(config, 'ocn', spectrum, grid, forcing)
#endif
        call output_grid_nc(config, currenttimestr, spectrum, grid, forcing)
      end if

      ! spectrum output
      if (config % outspec > 0) then
        if (mod(currenttime % gethour(),config % outspec) == 0 .and. fullhour) then
          call output_spectrum_nc(config, currenttimestr, spectrum, grid, forcing)
        end if
      else if (config % outspec == -1) then
        call output_spectrum_nc(config, currenttimestr, spectrum, grid, forcing)
      end if

      ! restart output
      if (config % outrst > 0) then
        if (mod(currenttime % gethour(), config % outrst) == 0 .and. fullhour)&
        call restart_write(currenttimestr, spectrum, grid)
      else if (config % outspec == -1) then
        call restart_write(currenttimestr, spectrum, grid)
      end if

    end do ! end outer loop

    100 format(2x, f4.2, 2x, f8.3, 6(1x, f9.6))

  end subroutine umwm_run


  subroutine umwm_finalize(grid, forcing)
    use umwm_forcing, only: forcing_type
    use umwm_grid, only: grid_type
    use umwm_util, only: dealloc
    use umwm_env, only: env_stop
    type(grid_type), intent(inout) :: grid
    type(forcing_type), intent(inout) :: forcing
    call dealloc() ! deallocate arrays
    call forcing % finalize()
    call grid % finalize()
    call env_stop()
  end subroutine umwm_finalize

end module umwm_top
