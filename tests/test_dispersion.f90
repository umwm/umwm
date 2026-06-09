program test_dispersion
  use umwm_constants, only: rk, stderr, stdout, twopi
  use umwm_dispersion, only: angular_frequency, group_speed, wavenumber

  implicit none

  integer, parameter :: num_frequencies = 1000
  real(rk), parameter :: frequency_min = 0.01_rk
  real(rk), parameter :: frequency_max = 100._rk
  real(rk), parameter :: depth = 100._rk
  real(rk), parameter :: water_density = 1000._rk
  real(rk), parameter :: gravity = 9.8_rk
  real(rk), parameter :: surface_tension = 0.074_rk
  real(rk), parameter :: max_reldiff_pct_allowed = 0.01_rk
  real(rk), parameter :: mean_reldiff_pct_allowed = 0.001_rk

  integer :: o
  real(rk) :: dlnf, max_reldiff, mean_reldiff
  real(rk), allocatable :: frequency(:), frequency_roundtrip(:)
  real(rk), allocatable :: cg(:), k(:), reldiff(:)
  logical :: ok = .true.

  allocate(frequency(num_frequencies))

  dlnf = (log(frequency_max) - log(frequency_min)) / real(num_frequencies - 1, rk)
  do o = 1, num_frequencies
    frequency(o) = exp(log(frequency_min) + real(o - 1, rk) * dlnf)
  end do

  k = wavenumber(frequency, depth, water_density, gravity, surface_tension)
  frequency_roundtrip = angular_frequency(k, depth, water_density, gravity, surface_tension) / twopi
  cg = group_speed(k, depth, water_density, gravity, surface_tension)
  reldiff = abs((frequency_roundtrip - frequency) / frequency * 100._rk)

  max_reldiff = maxval(reldiff)
  mean_reldiff = sum(reldiff) / real(num_frequencies, rk)

  if (max_reldiff > max_reldiff_pct_allowed) then
    write(stderr, '(a,es12.4,a)') &
      'test_dispersion: maximum relative error was ', max_reldiff, '%.'
    ok = .false.
  end if

  if (mean_reldiff > mean_reldiff_pct_allowed) then
    write(stderr, '(a,es12.4,a)') &
      'test_dispersion: mean relative error was ', mean_reldiff, '%.'
    ok = .false.
  end if

  if (any(.not. (cg > 0._rk)) .or. any(cg /= cg)) then
    write(stderr, '(a)') &
      'test_dispersion: group speed must be finite and positive.'
    ok = .false.
  end if

  if (ok) then
    write(stdout, '(a)') 'test_dispersion: All tests passed.'
  else
    write(stderr, '(a)') 'test_dispersion: One or more tests failed.'
    error stop 1
  end if

end program test_dispersion
