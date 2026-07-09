module umwm_restart
  ! Provides read and write subroutines for UMWM restart files
  use, intrinsic :: iso_fortran_env, only: real64
  use umwm_grid, only: grid_type
  use umwm_io, only: put_time_metadata, seconds_since_reference
  use umwm_module, only: e, f, ierr, k, mpisize, nproc, ustar
  use netcdf
  use umwm_spectrum, only: spectrum_type
  use umwm_util, only: raiseexception

#ifdef MPI
  use mpi
  use umwm_mpi
#endif 

  implicit none

  private
  public :: restart_read, restart_write

contains

  subroutine restart_read(timestr, spectrum, grid)
    character(19), intent(in) :: timestr
    type(spectrum_type), intent(in) :: spectrum
    type(grid_type), intent(in) :: grid
    character(19) :: timestrnew
    character(9999) :: filename
    integer :: stat, ncid, ustid, specid

    if (nproc == 0) &
      write(*, '(a)')'umwm: restart_read: reading restart file for ' // timestr

    timestrnew = timestr
    timestrnew(11:11) = '_'

    filename = 'restart/umwmrst_' // timestrnew // '.nc'

    stat = nf90_open(trim(filename), nf90_share, ncid)

    if (stat /= 0 .and. nproc == 0) then
      call raiseexception('abort', 'restart_read', nf90_strerror(stat))
      stop
    end if

    !TODO we should read the frequency and direction dimensions 
    !TODO and make sure they're consistent with the values in memory
    !TODO we should error-handle the values of NetCDF statuses

    stat = nf90_inq_varid(ncid, 'F', specid)
    stat = nf90_inq_varid(ncid, 'ust', ustid)
    stat = nf90_get_var(ncid, specid, e(:,:,grid % istart:grid % iend), &
                        start=[1, 1, grid % istart], &
                        count=[spectrum % num_frequencies, spectrum % num_directions, &
                               grid % iend - grid % istart + 1])
    stat = nf90_get_var(ncid, ustid, ustar(grid % istart:grid % iend), &
                        start=[grid % istart], count=[grid % iend - grid % istart + 1])
    stat = nf90_close(ncid)

#ifdef MPI
    call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

  end subroutine restart_read


  subroutine restart_write(timestr, reftimestr, spectrum, grid)
    character(19), intent(in) :: timestr
    character(*), intent(in) :: reftimestr
    type(spectrum_type), intent(in) :: spectrum
    type(grid_type), intent(in) :: grid
  
    integer :: i, nn
    integer :: stat, ncid, xdimid, fdimid, thdimid, tdimid
    integer :: kid, lonid, latid, freqid, thetaid, timeid, ustid, specid
    real :: lon_tmp(grid % im), lat_tmp(grid % im)
    real(real64) :: time_value(1)

    if (nproc == 0) then

      !TODO we should error-handle the values of NetCDF statuses
      stat = nf90_create('restart/umwmrst_' // timestr // '.nc', NF90_CLOBBER, ncid)

      stat = nf90_def_dim(ncid, 'x', grid % im, xdimid)
      stat = nf90_def_dim(ncid, 'f', spectrum % num_frequencies, fdimid)
      stat = nf90_def_dim(ncid, 'th', spectrum % num_directions, thdimid)
      stat = nf90_def_dim(ncid, 'time', 1, tdimid)

      stat = nf90_def_var(ncid, 'time', nf90_double, [tdimid], timeid)
      call put_time_metadata(ncid, timeid, reftimestr)

      stat = nf90_def_var(ncid, 'lon', nf90_float, [xdimid], lonid)
      stat = nf90_put_att(ncid, lonid, name='description', values='longitude')
      stat = nf90_put_att(ncid, lonid, name='units', values='degrees east')

      stat = nf90_def_var(ncid, 'lat', nf90_float, [xdimid], latid)
      stat = nf90_put_att(ncid, latid, name='description', values='latitude')
      stat = nf90_put_att(ncid, latid, name='units', values='degrees north')

      stat = nf90_def_var(ncid, 'frequency', nf90_float, [fdimid], freqid)
      stat = nf90_put_att(ncid, freqid, name='description', values='frequency')
      stat = nf90_put_att(ncid, freqid, name='units', values='hz')

      stat = nf90_def_var(ncid, 'theta', nf90_float, [thdimid], thetaid)
      stat = nf90_put_att(ncid, thetaid, name='description', values='directions')
      stat = nf90_put_att(ncid, thetaid, name='units', values='rad')

      stat = nf90_def_var(ncid, 'ust', nf90_float, [xdimid], ustid)
      stat = nf90_put_att(ncid, ustid, name='description', values='friction velocity')
      stat = nf90_put_att(ncid, ustid, name='units', values='m s^-1')

      stat = nf90_def_var(ncid, 'wavenumber', nf90_float, [fdimid, xdimid], kid)
      stat = nf90_put_att(ncid, kid, name='description', values='wavenumber')
      stat = nf90_put_att(ncid, kid, name='units', values='rad m^-1')

      stat = nf90_def_var(ncid, 'F', nf90_float, [fdimid, thdimid, xdimid], specid)
      stat = nf90_put_att(ncid, specid, name='description', values='wave energy spectrum')
      stat = nf90_put_att(ncid, specid, name='units', values='m^4 rad^-1')

      stat = nf90_enddef(ncid)

      ! fill in lon and lat arrays
      do i = 1, grid % im
        lon_tmp(i) = grid % lon(grid % mi(i), grid % ni(i))
        lat_tmp(i) = grid % lat(grid % mi(i), grid % ni(i))
      end do

      time_value(1) = seconds_since_reference(timestr, reftimestr)
      stat = nf90_put_var(ncid, timeid, time_value, start=[1], count=[1])
      stat = nf90_put_var(ncid, lonid, lon_tmp)
      stat = nf90_put_var(ncid, latid, lat_tmp)
      stat = nf90_put_var(ncid, freqid, spectrum % frequency)
      stat = nf90_put_var(ncid, thetaid, spectrum % direction)

      stat = nf90_close(ncid)

    end if

    ! loop over processes in order
    do nn = 0, mpisize - 1

      ! write to file if it is my turn
      if (nproc == nn) then

        stat = nf90_open('restart/umwmrst_' // timestr // '.nc', NF90_WRITE, ncid)
        stat = nf90_inq_dimid(ncid, 'x', xdimid)
        stat = nf90_inq_dimid(ncid, 'f', fdimid)
        stat = nf90_inq_dimid(ncid, 'th', thdimid)
        stat = nf90_inq_varid(ncid, 'F', specid)
        stat = nf90_inq_varid(ncid, 'wavenumber', kid)
        stat = nf90_inq_varid(ncid, 'ust', ustid)
        stat = nf90_put_var(ncid, specid, e(:,:,grid % istart:grid % iend), &
                            start=[1, 1, grid % istart], &
                            count=[spectrum % num_frequencies, spectrum % num_directions, &
                                   grid % iend - grid % istart + 1])
        stat = nf90_put_var(ncid, kid, k(:,grid % istart:grid % iend), &
                            start=[1, grid % istart], &
                            count=[spectrum % num_frequencies, grid % iend-grid % istart+1])
        stat = nf90_put_var(ncid, ustid, ustar(grid % istart:grid % iend), &
                            start=[grid % istart], count=[grid % iend-grid % istart+1])
        stat = nf90_close(ncid)

      end if

#ifdef MPI
      ! this call ensures that all processes wait for each other
      call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

    end do

    if (nproc == 0) &
      write(*, '(a)') 'umwm: restart_write: restart written to restart/umwmrst_' // timestr // '.nc'

  end subroutine restart_write

end module umwm_restart
