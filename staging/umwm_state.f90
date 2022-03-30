module umwm_state

  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type

  implicit none

  private
  public :: state_type

  type :: state_type
    real, allocatable :: variance(:,:,:,:)
  contains
  end type state_type

  interface state_type
    module procedure :: state_type_cons
  end interface state_type

contains

  type(state_type) elemental function state_type_cons(grid, spectrum) result(res)

    type(grid_type), intent(in) :: grid
    type(spectrum_type), intent(in) :: spectrum
    integer :: stat

    allocate(res % variance(spectrum % num_frequencies, &
                            spectrum % num_directions, &
                            grid % size_x, &
                            grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % variance.'
    res % variance = 0

  end function state_type_cons

end module umwm_state
