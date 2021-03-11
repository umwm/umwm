module umwm_spectrum

  use umwm_config, only: config_type

  implicit none

  private
  public :: spectrum_type

  type :: spectrum_type
    integer :: num_frequencies
    integer :: num_directions
    real :: frequency_min
    real :: frequency_max
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
    res % frequency_min = config % frequency_min
    res % frequency_max = config % frequency_max

    !allocate(res % frequency(res % num_frequencies), stat=stat)
    !if (stat /= 0) error stop 'Error allocating spectrum_type % frequency.'
    res % frequency = frequency_logspace(res % frequency_min, &
                                         res % frequency_max, &
                                         res % num_frequencies)

    allocate(res % direction(res % num_directions), stat=stat)
    if (stat /= 0) error stop 'Error allocating spectrum_type % direction.'
    res % direction = 0

  end function spectrum_type_cons


  pure function frequency_logspace(frequency_min, frequency_max, num_frequencies) result(res)
    !! Computes the frequency array in the range [frequency_min, frequency_max]
    !! in a logarithmic space with a total of num_frequencies bins.
    real, intent(in) :: frequency_min
      !! Minimum frequency [Hz]
    real, intent(in) :: frequency_max
      !! Maximum frequency [Hz]
    integer, intent(in) :: num_frequencies
      !! Total number of frequency bins
    real :: res(num_frequencies)
    real :: frequency_spacing
    integer :: n

    frequency_spacing = (log(frequency_max) - log(frequency_min)) &
                      / real(num_frequencies - 1)

    do concurrent(n = 1:num_frequencies)
      res(n) = exp(log(frequency_min) + (n - 1) * frequency_spacing)
    end do

  end function frequency_logspace

end module umwm_spectrum
