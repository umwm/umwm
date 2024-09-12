program test_spectrum
  use, intrinsic :: iso_fortran_env, only: stdout => output_unit, stderr => error_unit
  use umwm_spectrum, only: spectrum_type

  implicit none

  type(spectrum_type) :: spectrum
  integer :: num_frequencies, num_directions
  real :: frequency_min, frequency_max
  logical :: ok = .true.

  num_frequencies = 50
  num_directions = 36
  frequency_min = 0.04
  frequency_max = 2

  spectrum = spectrum_type( &
    num_frequencies, num_directions, frequency_min, frequency_max &
  )

  if (spectrum % num_frequencies /= num_frequencies) then
    write(stderr, '(a)') 'spectrum % num_frequencies is expected.. failed'
    ok = .false.
  end if

  if (spectrum % num_directions /= num_directions) then
    write(stderr, '(a)') 'spectrum % num_directions is expected.. failed'
    ok = .false.
  end if

  if (spectrum % frequency_min /= frequency_min) then
    write(stderr, '(a)') 'spectrum % frequency_min is expected.. failed'
    ok = .false.
  end if
  
  if (spectrum % frequency_max /= frequency_max) then
    write(stderr, '(a)') 'spectrum % frequency_min is expected.. failed'
    ok = .false.
  end if

  if (ok) then
    write(stdout, '(a)') 'test_spectrum: All tests passed.'
  else
    write(stderr, '(a)') 'test_spectrum: One or more tests failed.'
    error stop 1
  end if

end program test_spectrum