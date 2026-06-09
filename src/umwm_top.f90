module umwm_top

  implicit none

contains

  subroutine umwm_initialize(spectrum)

    use umwm_env, only: env_init
    use umwm_module,only: starttimestr => starttimestr_nml
    use umwm_init,  only: nmlread, initialize_spectrum, alloc, grid, masks, partition, remap, init
    use umwm_io,    only: input_nc, output_grid
    use umwm_spectrum, only: spectrum_type
    use umwm_stokes,only: stokes_drift

    type(spectrum_type), intent(out) :: spectrum

    call env_init()             ! initialize the environment
    call nmlread()              ! read the namelist
    call initialize_spectrum(spectrum)
    call alloc(1)               ! allocate 2-d arrays
    call grid()                 ! define model grid
    call masks()                ! define seamasks
    call partition()            ! domain partitioning
    call alloc(2, spectrum)     ! allocate unrolled arrays
    call remap()                ! remap 2-d arrays to 1-d
    call output_grid()          ! output a grid file
    call input_nc(starttimestr) ! read initial fields
    call init(spectrum)         ! initialize model variables
    call stokes_drift(spectrum, 'init') ! initialize stokes drift arrays

  end subroutine umwm_initialize


  subroutine umwm_run(starttimestr, stoptimestr, spectrum)

#ifdef MPI
    use umwm_mpi, only: exchange_halo
    use umwm_module, only: ierr
    use mpi
#endif

    use umwm_module, only: cd, currenttime, dtg, dts, e, ef, f, first, &
                           firstdtg, iend, iip, istart, nproc, &
                           nproc_plot, oc, outgrid, outspec, outrst, &
                           restart, starttime, stokes, stoptime, sumt, &
                           wdir, wspd
    use umwm_physics, only: source, diag
    use umwm_advection,only: propagation, refraction
    use umwm_forcing, only: forcinginput, forcinginterpolate
    use umwm_io, only: output_grid_nc, output_spectrum_nc
    use umwm_restart, only: restart_read, restart_write
    use umwm_spectrum, only: spectrum_type
    use umwm_stokes, only: stokes_drift
    use umwm_util, only: sigwaveheight, meanwaveperiod
    use umwm_stress, only: stress

    use umwm_source_functions, only: sin_d12, sds_d12, snl_d12, s_ice

    use datetime_module

    character(19), intent(in) :: starttimestr, stoptimestr
    type(spectrum_type), intent(in) :: spectrum

    character(19) :: currenttimestr
    logical :: fullhour

    ! convert start and stop time strings to datetime objects:
    starttime = strptime(starttimestr,'%Y-%m-%d %H:%M:%S')
    stoptime  = strptime(stoptimestr, '%Y-%m-%d %H:%M:%S')

    currenttime = starttime
    currenttimestr = trim(currenttime % strftime('%Y-%m-%d_%H:%M:%S'))

    ! read wave spectrum field from a restart file if necessary:
    if (first .and. restart) call restart_read(starttimestr, spectrum)

    do while (currenttime < stoptime) ! outer time loop

      ! report current and next checkpoint time:
      if(nproc == 0)write(*,fmt='(a)')&
      'umwm: solver: current time is:     '//currenttimestr

      currenttime = currenttime + timedelta(seconds=nint(dtg))
      currenttimestr = trim(currenttime % strftime('%Y-%m-%d_%H:%M:%S'))

      if(nproc == 0)write(*,fmt='(a)')&
        'umwm: solver: integrating to time: '//currenttimestr

      ! read atmosphere and ocean input fields from file.
      ! in case of esmf coupling, fields are assumed to be updated
      ! externally, and this call is not used.
#ifndef ESMF
      call forcinginput(currenttimestr)
#endif

      ! inner time loop: global time step
      sumt = 0
      do while (sumt < dtg)

#ifndef ESMF
        call forcinginterpolate() ! interpolate force fields in time
#endif

        call sin_d12(spectrum) ! compute source input term Sin
        call sds_d12(spectrum) ! compute source dissipation term Sds
        call snl_d12(spectrum) ! compute non-linear source term Snl
        call s_ice(spectrum)   ! compute sea ice attenuation term Sice
        call source(spectrum)  ! integrate source functions

#ifdef MPI
        call exchange_halo() ! exchange halo points
#endif

        call propagation(spectrum) ! compute advection and integrate

#ifdef ESMF
        e(:,:,istart:iend) = ef(:,:,istart:iend) ! update
#endif

        call refraction(spectrum)    ! compute refraction and integrate

        e(:,:,istart:iend) = ef(:,:,istart:iend) ! update

        call stress('atm', spectrum) ! compute wind stress and drag coefficient

#ifdef ESMF
        call stress('ocn', spectrum) ! compute stress into ocean top and bottom
#endif

        if (first) then

          ! diagnostic calculations before output
          call diag(spectrum)

          if (outgrid > 0) call output_grid_nc(starttimestr, spectrum)
          if (outspec > 0) call output_spectrum_nc(starttimestr, spectrum)

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
          write(*,fmt=100)sumt/dtg,dts,wspd(iip),wdir(iip),      &
                          sigwaveheight(iip, spectrum),meanwaveperiod(iip, spectrum),&
                          cd(iip)*1e3,f(oc(iip))
        end if

      end do ! end while(sumt<dtg) loop

      if (firstdtg) firstdtg = .false.
      if (stokes) call stokes_drift(spectrum)

      call diag(spectrum) ! model diagnostics for output

      fullhour = currenttime % getminute() == 0 &
           .and. currenttime % getsecond() == 0

      ! gridded output
      if(outgrid > 0)then
        if(mod(currenttime % gethour(),outgrid) == 0 .and. fullhour)then
#ifndef ESMF
          call stress('ocn', spectrum)
#endif
          call output_grid_nc(currenttimestr, spectrum)
        end if
      elseif(outgrid == -1)then
#ifndef ESMF
        call stress('ocn', spectrum)
#endif
        call output_grid_nc(currenttimestr, spectrum)
      end if

      ! spectrum output
      if (outspec > 0) then
        if (mod(currenttime % gethour(),outspec) == 0 .and. fullhour) then
          call output_spectrum_nc(currenttimestr, spectrum)
        end if
      else if (outspec == -1) then
        call output_spectrum_nc(currenttimestr, spectrum)
      end if

      ! restart output
      if (outrst > 0) then
        if (mod(currenttime % gethour(), outrst) == 0 .and. fullhour)&
        call restart_write(currenttimestr, spectrum)
      else if (outspec == -1) then
        call restart_write(currenttimestr, spectrum)
      end if

    end do ! end outer loop

    100 format(2x, f4.2, 2x, f8.3, 6(1x, f9.6))

  end subroutine umwm_run


  subroutine umwm_finalize
    use umwm_util, only: dealloc
    use umwm_env, only: env_stop
    call dealloc() ! deallocate arrays
    call env_stop()
  end subroutine umwm_finalize

end module umwm_top
