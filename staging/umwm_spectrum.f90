module umwm_spectrum

  use iso_fortran_env, only: real64
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

    res % frequency = frequency_logspace(res % frequency_min, &
                                         res % frequency_max, &
                                         res % num_frequencies)

    res % direction = direction(res % num_directions)

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
                      / (num_frequencies - 1)

    res = [(exp(log(frequency_min) + (n - 1) * frequency_spacing), &
            n = 1, num_frequencies)]

  end function frequency_logspace

  pure function direction(num_directions) result(res)
    integer, intent(in) :: num_directions
    real :: res(num_directions)
    real, parameter :: PI = 4 * atan(1._real64)
    real :: directional_spacing
    integer :: n

    directional_spacing = 2 * PI / num_directions

    res = [((n - 0.5 * (num_directions + 1)) * directional_spacing, n = 1, num_directions)]

  end function direction

end module umwm_spectrum
