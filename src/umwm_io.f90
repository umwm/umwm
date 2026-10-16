module umwm_io
! NetCDF input/output helpers for the wave model
use, intrinsic :: iso_fortran_env, only: real64
use datetime_module, only: date2num, datetime, strptime
use netcdf
use umwm_env, only: nproc
use umwm_grid, only: grid_type

implicit none

private
public :: gather_field, nc_check, output_grid, put_description, &
          put_time_metadata, seconds_since_reference

contains


subroutine output_grid(grid)
! Outputs grid related fields into a netcdf file.
type(grid_type), intent(in) :: grid
integer :: ncid
integer :: xdimid,ydimid
integer :: lonid,latid,dlonid,dlatid,dxid,dyid,arid,maskid,did,nprocid
integer :: xid,yid,curvid

if(nproc == 0)then

  call nc_check(nf90_create('output/umwmout.grid',nf90_clobber,ncid))
  call nc_check(nf90_def_dim(ncid,'x',grid % mm,xdimid))
  call nc_check(nf90_def_dim(ncid,'y',grid % nm,ydimid))
  call nc_check(nf90_def_var(ncid,'lon',NF90_FLOAT,[xdimid,ydimid],lonid))
  call nc_check(nf90_def_var(ncid,'lat',NF90_FLOAT,[xdimid,ydimid],latid))
  call nc_check(nf90_def_var(ncid,'xx',NF90_FLOAT,[xdimid,ydimid],xid))
  call nc_check(nf90_def_var(ncid,'yy',NF90_FLOAT,[xdimid,ydimid],yid))
  call nc_check(nf90_def_var(ncid,'dlon',NF90_FLOAT,[xdimid,ydimid],dlonid))
  call nc_check(nf90_def_var(ncid,'dlat',NF90_FLOAT,[xdimid,ydimid],dlatid))
  call nc_check(nf90_def_var(ncid,'dx',NF90_FLOAT,[xdimid,ydimid],dxid))
  call nc_check(nf90_def_var(ncid,'dy',NF90_FLOAT,[xdimid,ydimid],dyid))
  call nc_check(nf90_def_var(ncid,'curvature',NF90_FLOAT,[xdimid,ydimid],curvid))
  call nc_check(nf90_def_var(ncid,'area',NF90_FLOAT,[xdimid,ydimid],arid))
  call nc_check(nf90_def_var(ncid,'depth',NF90_FLOAT,[xdimid,ydimid],did))
  call nc_check(nf90_def_var(ncid,'seamask',nf90_int,[xdimid,ydimid],maskid))
  call nc_check(nf90_def_var(ncid,'nproc',nf90_int,[xdimid,ydimid],nprocid))
  call nc_check(nf90_enddef(ncid))
  call nc_check(nf90_put_var(ncid,lonid,grid % lon))
  call nc_check(nf90_put_var(ncid,latid,grid % lat))
  call nc_check(nf90_put_var(ncid,xid,grid % x))
  call nc_check(nf90_put_var(ncid,yid,grid % y))
  call nc_check(nf90_put_var(ncid,dlonid,grid % dlon))
  call nc_check(nf90_put_var(ncid,dlatid,grid % dlat))
  call nc_check(nf90_put_var(ncid,dxid,grid % dx_2d))
  call nc_check(nf90_put_var(ncid,dyid,grid % dy_2d))
  call nc_check(nf90_put_var(ncid,curvid,grid % curv))
  call nc_check(nf90_put_var(ncid,arid,grid % ar_2d))
  call nc_check(nf90_put_var(ncid,maskid,grid % mask))
  call nc_check(nf90_put_var(ncid,did,grid % d_2d))
  call nc_check(nf90_put_var(ncid,nprocid,grid % nproc_out))
  call nc_check(nf90_close(ncid))

  write(*,'(a)')'umwm: output_grid: grid definition written in output/umwmout.grid'

end if

end subroutine output_grid


subroutine put_description(ncid, varid, description)
  integer, intent(in) :: ncid
  integer, intent(in) :: varid
  character(*), intent(in) :: description
  integer :: stat

  stat = nf90_put_att(ncid, varid, name='description', values=description)
  stat = nf90_put_att(ncid, varid, name='long_name', values=description)

end subroutine put_description


subroutine put_time_metadata(ncid, timeid, reftimestr)
  integer, intent(in) :: ncid
  integer, intent(in) :: timeid
  character(*), intent(in) :: reftimestr
  integer :: stat

  call put_description(ncid, timeid, 'time')
  stat = nf90_put_att(ncid, timeid, name='cartesian_axis', values='T')
  stat = nf90_put_att(ncid, timeid, name='units', values='seconds since '//trim(reftimestr))
  stat = nf90_put_att(ncid, timeid, name='time_origin', values=trim(reftimestr))

end subroutine put_time_metadata


real(real64) function seconds_since_reference(timestr, reftimestr) result(seconds)
  character(*), intent(in) :: timestr
  character(*), intent(in) :: reftimestr
  type(datetime) :: output_time
  type(datetime) :: reference_time

  output_time = strptime(normalized_timestamp(timestr), '%Y-%m-%d %H:%M:%S', tz=0._real64)
  reference_time = strptime(normalized_timestamp(reftimestr), '%Y-%m-%d %H:%M:%S', tz=0._real64)
  seconds = anint((date2num(output_time) - date2num(reference_time)) * 86400._real64)

end function seconds_since_reference


character(19) function normalized_timestamp(timestr) result(normalized)
  character(*), intent(in) :: timestr

  normalized = timestr
  if (normalized(11:11) == '_') normalized(11:11) = ' '

end function normalized_timestamp


subroutine gather_field(field, field_mn, grid)
! Gathers a tile field on the root process and remaps it to (m,n).
! Land points are set to NaN.
use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
#ifdef MPI
use umwm_mpi, only: gather_array
#endif

type(grid_type), intent(in) :: grid
real, intent(in) :: field(grid % istart:grid % iend)
real, intent(out) :: field_mn(grid % mm,grid % nm)

real :: field_ii(grid % imm)

field_ii = ieee_value(field_ii, ieee_quiet_nan)

#ifdef MPI
call gather_array(field, field_ii(1:grid % im), grid)
#endif

if(nproc == 0)then
  field_ii(grid % istart:grid % iend) = field
  field_mn = grid % remap_i2mn(field_ii)
end if

end subroutine gather_field


subroutine nc_check(stat)
! Checks for netcdf errors and if any, print and abort.

integer,intent(in) :: stat

if(stat /= nf90_noerr)then
  write(*,*)'error in netcdf i/o'
  write(*,*)trim(nf90_strerror(stat))
  stop
end if

end subroutine nc_check

end module umwm_io
