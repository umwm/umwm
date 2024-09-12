program test_dispersion
  use, intrinsic :: iso_fortran_env, only: stdout => output_unit, stderr => error_unit
  use umwm_dispersion, only: wavenumber
  use umwm_spectrum, only: spectrum_type

  implicit none

  integer :: num_frequencies
  real :: frequency_min, frequency_max
  type(spectrum_type) :: spectrum
  real, allocatable :: k(:)
  logical :: ok = .true.

  num_frequencies = 100
  frequency_min = 0.02
  frequency_max = 50

  spectrum = spectrum_type(num_frequencies, 1, frequency_min, frequency_max)

  k = wavenumber(spectrum % frequency, 1e3, 1e3, 9.8, 0.074)

  print *, 'k = ', k

  if (ok) then
    write(stdout, '(a)') 'test_dispersion: All tests passed.'
  else
    write(stderr, '(a)') 'test_dispersion: One or more tests failed.'
    error stop 1
  end if

end program test_dispersion