module umwm_spectrum

  use umwm_constants, only: rk, twopi

  implicit none

  private
  public :: spectrum_type

  type :: spectrum_type
    integer :: num_frequencies
    integer :: num_directions
    real(rk) :: frequency_min
    real(rk) :: frequency_max
    real(rk), allocatable :: frequency(:)
    real(rk), allocatable :: direction(:)
    real(rk) :: dlnf
    real(rk) :: dth
  end type spectrum_type

  interface spectrum_type
    module procedure :: spectrum_type_cons
  end interface spectrum_type

contains

  pure type(spectrum_type) function spectrum_type_cons( &
    num_frequencies, num_directions, frequency_min, frequency_max &
  ) result(res)

    integer, intent(in) :: num_frequencies
    integer, intent(in) :: num_directions
    real(rk), intent(in) :: frequency_min
    real(rk), intent(in) :: frequency_max

    res % num_frequencies = num_frequencies
    res % num_directions = num_directions
    res % frequency_min = frequency_min
    res % frequency_max = frequency_max

    allocate(res % frequency(res % num_frequencies))
    allocate(res % direction(res % num_directions))

    res % dlnf = (log(res % frequency_max) - log(res % frequency_min)) &
      / float(res % num_frequencies - 1)
    res % dth = twopi / float(res % num_directions)
    res % frequency = frequency_logspace(res % frequency_min, &
                                         res % frequency_max, &
                                         res % num_frequencies)
    res % direction = direction(res % num_directions)

  end function spectrum_type_cons


  pure function frequency_logspace(frequency_min, frequency_max, num_frequencies) result(res)
    !! Computes the frequency array in the range [frequency_min, frequency_max]
    !! in a logarithmic space with a total of num_frequencies bins.
    real(rk), intent(in) :: frequency_min
      !! Minimum frequency [Hz]
    real(rk), intent(in) :: frequency_max
      !! Maximum frequency [Hz]
    integer, intent(in) :: num_frequencies
      !! Total number of frequency bins
    real(rk) :: res(num_frequencies)
    real(rk) :: frequency_spacing
    integer :: n

    frequency_spacing = (log(frequency_max) - log(frequency_min)) &
                      / float(num_frequencies - 1)

    do n = 1, num_frequencies
      res(n) = exp(log(frequency_min) + (n - 1) * frequency_spacing)
    end do

  end function frequency_logspace


  pure function direction(num_directions) result(res)
    integer, intent(in) :: num_directions
    real(rk) :: res(num_directions)
    real(rk) :: directional_spacing
    integer :: n

    directional_spacing = twopi / float(num_directions)

    do n = 1, num_directions
      res(n) = (n - 0.5 * (num_directions + 1)) * directional_spacing
    end do

  end function direction

end module umwm_spectrum
