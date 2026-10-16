module umwm_util
  ! Utility procedures.

  implicit none

  private
  public :: raiseexception, distance_haversine

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

end module umwm_util
