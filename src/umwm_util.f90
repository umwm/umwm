module umwm_util
  ! Utility procedures.

  implicit none

  private
  public :: raiseexception, distance_haversine, input_size_problem

contains

  subroutine raiseexception(exceptiontype, routinename, message, flag)
    ! Prints a message to stdout and sets the flag to .false. if present.
    character(len=*), intent(in) :: exceptiontype
    character(len=*), intent(in) :: routinename
    character(len=*), intent(in) :: message
    logical, intent(inout), optional :: flag

    if (present(flag)) flag = .false.

    write(unit=*, fmt='(a)') 'umwm: ' // trim(routinename) // ': ' &
                             // trim(exceptiontype) // ': ' // trim(message)

  end subroutine raiseexception


  pure function distance_haversine(lon1, lon2, lat1, lat2) result(distance)
    ! Shortest distance between two points on the unit sphere
    ! using the Haversine formula; input angles in radians.
    real, intent(in) :: lon1, lon2, lat1, lat2
    real :: distance
    real :: dlon, dlat

    dlon = abs(lon2 - lon1)
    dlat = abs(lat2 - lat1)

    distance = 2 * asin(sqrt((sin(0.5 * dlat))**2 + &
                             cos(lat1) * cos(lat2) * (sin(0.5 * dlon))**2))

  end function distance_haversine


  function input_size_problem(ncid, name, mm, nm) result(problem)
    ! Why the field name of the open netCDF file ncid does not fit a grid of
    ! mm x nm cells: empty if its first two dimensions, x and y, are mm and
    ! nm, or if the file has no such field, which its reader reports. A file
    ! for a larger grid would otherwise be read in part, without an error.
    use netcdf, only: nf90_inq_varid, nf90_inquire_variable, nf90_inquire_dimension, &
                      nf90_max_var_dims, nf90_noerr
    integer, intent(in) :: ncid, mm, nm
    character(*), intent(in) :: name
    character(:), allocatable :: problem
    integer :: varid, ndims, dimids(nf90_max_var_dims), nx, ny
    character(64) :: sizes

    problem = ''
    if (nf90_inq_varid(ncid, name, varid) /= nf90_noerr) return
    if (nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids) /= nf90_noerr) return
    nx = 1
    ny = 1
    if (ndims >= 1) then
      if (nf90_inquire_dimension(ncid, dimids(1), len=nx) /= nf90_noerr) return
    end if
    if (ndims >= 2) then
      if (nf90_inquire_dimension(ncid, dimids(2), len=ny) /= nf90_noerr) return
    end if
    if (nx == mm .and. ny == nm) return
    write(sizes, '(i0,a,i0,a,i0,a,i0)') nx, ' x ', ny, ', but the grid is ', mm, ' x ', nm
    problem = name // ' is ' // trim(sizes) // ' (mm x nm)'

  end function input_size_problem

end module umwm_util
