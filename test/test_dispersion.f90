program test_dispersion
  use, intrinsic :: iso_fortran_env, only: stdout => output_unit, stderr => error_unit
  use umwm_dispersion, only: angular_frequency, group_speed, wavenumber
  use umwm_spectrum, only: spectrum_type

  implicit none

  integer :: num_frequencies
  real :: frequency_min, frequency_max, depth
  type(spectrum_type) :: spectrum
  real, allocatable :: f(:), k(:), reldiff(:)
  logical :: ok = .true.

  real, parameter :: twopi = 2 * acos(-1.d0)

  num_frequencies = 1000
  frequency_min = 0.01
  frequency_max = 100
  depth = 1e2

  spectrum = spectrum_type(num_frequencies, 1, frequency_min, frequency_max)

  k = wavenumber(spectrum % frequency, depth, 1e3, 9.8, 0.074)
  f = angular_frequency(k, depth, 1e3, 9.8, 0.074) / twopi
  reldiff = (f - spectrum % frequency) / spectrum % frequency * 100

  if (maxval(abs(reldiff)) > 1e-4) then
    write(stderr, '(a)') 'test_dispersion: Maximum relative error < 0.01%.. failed.'
    ok = .false.
  end if

  if (sum(abs(reldiff)) / num_frequencies > 1e-5) then
    write(stderr, '(a)') 'test_dispersion: Mean relative error < 0.001%.. failed'
    ok = .false.
  end if

  if (ok) then
    write(stdout, '(a)') 'test_dispersion: All tests passed.'
  else
    write(stderr, '(a)') 'test_dispersion: One or more tests failed.'
    error stop 1
  end if

end program test_dispersion