module umwm_domain

  !! domain_type is the main UMWM derived type.
  !! Its components are an instance of grid_type, which defines its geographical
  !! grid, and an instance of spectrum_type, which defines its spectral space.

  use datetime_module, only: datetime
  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type

  implicit none

  private
  public :: domain_type

  type :: domain_type
    type(datetime) :: start_time
    type(datetime) :: stop_time
    type(grid_type) :: grid
    type(spectrum_type) :: spectrum
    real, allocatable :: variance(:,:,:,:)
  end type domain_type

  interface domain_type
    module procedure :: domain_type_cons
  end interface domain_type

contains

  type(domain_type) pure function domain_type_cons(start_time, stop_time, grid, spectrum) result(res)
    type(datetime), intent(in) :: start_time
    type(datetime), intent(in) :: stop_time
    type(grid_type), intent(in) :: grid
    type(spectrum_type), intent(in) :: spectrum
    integer :: stat

    res % start_time = start_time
    res % stop_time = stop_time
    res % grid = grid
    res % spectrum = spectrum

    allocate(res % variance(spectrum % num_frequencies, &
                            spectrum % num_directions, &
                            grid % size_x, &
                            grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % variance.'

  end function domain_type_cons

end module umwm_domain
