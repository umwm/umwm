module umwm_spectrum

  use umwm_config, only: config_type

  implicit none

  private
  public :: spectrum_type

  type :: spectrum_type
    integer :: num_frequencies
    integer :: num_directions
    real, allocatable :: frequency(:)
    real, allocatable :: direction(:)
  end type spectrum_type

  interface spectrum_type
    module procedure :: spectrum_type_cons
  end interface spectrum_type

contains

  type(spectrum_type) pure function spectrum_type_cons(config) result(res)
    type(config_type), intent(in) :: config
    integer :: stat

    res % num_frequencies = config % num_frequencies
    res % num_directions = config % num_directions

    allocate(res % frequency(res % num_frequencies), stat=stat)
    if (stat /= 0) error stop 'Error allocating spectrum_type % frequency.'
    res % frequency = 0

    allocate(res % direction(res % num_directions), stat=stat)
    if (stat /= 0) error stop 'Error allocating spectrum_type % direction.'
    res % direction = 0

  end function spectrum_type_cons

end module umwm_spectrum
