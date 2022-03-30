module umwm_domain

  !! domain_type is the main UMWM derived type.
  !! Its components are an instance of grid_type, which defines its geographical
  !! grid, and an instance of spectrum_type, which defines its spectral space.
  !!
  !! +-----------------------+
  !! | domain                |
  !! | +------+ +----------+ |
  !! | | grid | | spectrum | |
  !! | +------+ +----------+ |
  !! +-----------------------+

  use datetime_module, only: datetime, timedelta
  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type
  use umwm_state, only: state_type

  implicit none

  private
  public :: domain_type

  type :: domain_type
    type(datetime) :: start_time
    type(datetime) :: stop_time
    type(datetime) :: current_time
    type(grid_type) :: grid
    type(spectrum_type) :: spectrum
    type(state_type) :: state
  contains
    procedure :: run
    procedure :: step
  end type domain_type

  interface domain_type
    module procedure :: domain_type_cons
  end interface domain_type

contains

  type(domain_type) elemental function domain_type_cons( &
    start_time, stop_time, grid, spectrum) result(res)

    type(datetime), intent(in) :: start_time
    type(datetime), intent(in) :: stop_time
    type(grid_type), intent(in) :: grid
    type(spectrum_type), intent(in) :: spectrum
    integer :: stat

    res % start_time = start_time
    res % stop_time = stop_time
    res % current_time = res % start_time
    res % grid = grid
    res % spectrum = spectrum
    res % state = state_type(grid, spectrum)

  end function domain_type_cons


  impure elemental subroutine run(self)
    class(domain_type), intent(inout) :: self
    do while (self % current_time < self % stop_time)
      print *, self % current_time % strftime('%Y-%m-%d %H:%M:%S')
      call self % step()
    end do
  end subroutine run


  elemental subroutine step(self)
    class(domain_type), intent(inout) :: self
    self % current_time = self % current_time + timedelta(hours=1)
  end subroutine step

end module umwm_domain
