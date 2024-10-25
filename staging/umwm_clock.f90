module umwm_clock
  use datetime_module, only: datetime, timedelta
  implicit none
    
  private
  public :: clock_type
    
  type :: clock_type
    type(datetime) :: start, stop, current
    type(timedelta) :: interval
  contains
    procedure :: tick
  end type clock_type

  interface clock_type
    module procedure clock_type_cons
  end interface clock_type

contains

  type(clock_type) elemental function clock_type_cons(start, stop, interval) result(res)
      type(datetime), intent(in) :: start, stop
      type(timedelta), intent(in) :: interval
      res % start = start
      res % stop = stop
      res % current = start
      res % interval = interval
  end function clock_type_cons
    
  elemental subroutine tick(self)
    class(clock_type), intent(inout) :: self
    self % current = self % current + self % interval
  end subroutine tick

end module umwm_clock