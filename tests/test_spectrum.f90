program test_spectrum
  use tuff, only: test, test_result
  use umwm_constants, only: rk, stderr, twopi
  use umwm_spectrum, only: spectrum_type

  implicit none

  integer, parameter :: num_frequencies = 33
  integer, parameter :: num_directions = 16
  real(rk), parameter :: frequency_min = 0.04_rk
  real(rk), parameter :: frequency_max = 0.5_rk

  type(test_result) :: suite

  suite = test('test_spectrum', [ &
    test(constructor_metadata), &
    test(endpoint_frequencies), &
    test(legacy_logspace_equivalence), &
    test(direction_grid_equivalence), &
    test(spacing_values) &
  ])

  if (.not. suite % ok) then
    error stop 1
  end if

contains

  function constructor_metadata() result(res)
    type(test_result) :: res
    type(spectrum_type) :: spectrum

    spectrum = spectrum_type(num_frequencies, num_directions, frequency_min, frequency_max)

    res = test('constructor_metadata', &
      spectrum % num_frequencies == num_frequencies .and. &
      spectrum % num_directions == num_directions .and. &
      spectrum % frequency_min == frequency_min .and. &
      spectrum % frequency_max == frequency_max .and. &
      size(spectrum % frequency) == num_frequencies .and. &
      size(spectrum % direction) == num_directions)
  end function constructor_metadata


  function endpoint_frequencies() result(res)
    type(test_result) :: res
    type(spectrum_type) :: spectrum

    spectrum = spectrum_type(num_frequencies, num_directions, frequency_min, frequency_max)

    res = test('endpoint_frequencies', &
      nearly_equal(spectrum % frequency(1), frequency_min) .and. &
      nearly_equal(spectrum % frequency(num_frequencies), frequency_max))
    if (.not. res % ok) then
      write(stderr, '(a)') 'test_spectrum: endpoint frequencies do not match inputs.'
    end if
  end function endpoint_frequencies


  function legacy_logspace_equivalence() result(res)
    type(test_result) :: res
    type(spectrum_type) :: spectrum
    real(rk) :: legacy_frequency(num_frequencies)
    real(rk) :: dlnf
    integer :: o

    spectrum = spectrum_type(num_frequencies, num_directions, frequency_min, frequency_max)

    dlnf = (log(frequency_max) - log(frequency_min)) / float(num_frequencies - 1)
    do o = 1, num_frequencies
      legacy_frequency(o) = exp(log(frequency_min) + (o - 1) * dlnf)
    end do

    res = test('legacy_logspace_equivalence', all_nearly_equal(spectrum % frequency, legacy_frequency))
    if (.not. res % ok) then
      write(stderr, '(a)') 'test_spectrum: log-spaced frequencies differ from legacy formula.'
    end if
  end function legacy_logspace_equivalence


  function direction_grid_equivalence() result(res)
    type(test_result) :: res
    type(spectrum_type) :: spectrum
    real(rk) :: legacy_direction(num_directions)
    real(rk) :: dth
    integer :: p

    spectrum = spectrum_type(num_frequencies, num_directions, frequency_min, frequency_max)

    dth = twopi / float(num_directions)
    do p = 1, num_directions
      legacy_direction(p) = (p - 0.5 * (num_directions + 1)) * dth
    end do

    res = test('direction_grid_equivalence', all_nearly_equal(spectrum % direction, legacy_direction))
    if (.not. res % ok) then
      write(stderr, '(a)') 'test_spectrum: direction grid differs from legacy formula.'
    end if
  end function direction_grid_equivalence


  function spacing_values() result(res)
    type(test_result) :: res
    type(spectrum_type) :: spectrum
    real(rk) :: expected_dlnf, expected_dth

    spectrum = spectrum_type(num_frequencies, num_directions, frequency_min, frequency_max)

    expected_dlnf = (log(frequency_max) - log(frequency_min)) / float(num_frequencies - 1)
    expected_dth = twopi / float(num_directions)

    res = test('spacing_values', nearly_equal(spectrum % dlnf, expected_dlnf) .and. &
                                 nearly_equal(spectrum % dth, expected_dth))
    if (.not. res % ok) then
      write(stderr, '(a)') 'test_spectrum: spacing values differ from legacy formulas.'
    end if
  end function spacing_values


  pure elemental logical function nearly_equal(a, b) result(res)
    real(rk), intent(in) :: a, b

    res = abs(a - b) <= 10._rk * epsilon(1._rk) * max(1._rk, abs(a), abs(b))
  end function nearly_equal


  pure logical function all_nearly_equal(a, b) result(res)
    real(rk), intent(in) :: a(:), b(:)

    res = size(a) == size(b)
    if (.not. res) return
    res = all(nearly_equal(a, b))
  end function all_nearly_equal


end program test_spectrum
