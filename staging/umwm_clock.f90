module umwm_clock
  use datetime_module, only: datetime, timedelta
  implicit none
    
  private
  public :: clock_type
    
  type :: clock_type
    type(datetime) :: start, stop, current
    type(timedelta) :: step
  contains
    procedure :: tick
  end type clock_type

  interface clock_type
    module procedure clock_type_cons
  end interface clock_type

contains

  type(clock_type) elemental function clock_type_cons(start, stop) result(res)
      type(datetime), intent(in) :: start, stop
      !type(timedelta), intent(in) :: step
      res % start = start
      res % stop = stop
      res % current = start
      !res % step = step
  end function clock_type_cons
    
  impure elemental subroutine tick(self)
    class(clock_type), intent(inout) :: self
    self % current = self % current + self % step
  end subroutine tick

end module umwm_clock