program test_dispersion
  use tuff, only: test, test_result
  use umwm_constants, only: rk, stderr, twopi
  use umwm_dispersion, only: angular_frequency, group_speed, wavenumber
#ifdef MPI
  use umwm_env, only: env_init, env_stop
#endif

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

  type(test_result) :: suite

#ifdef MPI
  call env_init()
#endif
  suite = test('test_dispersion', [ &
    test(roundtrip_accuracy), &
    test(group_speed_is_finite_and_positive) &
  ])

#ifdef MPI
  call env_stop()
#endif
  if (.not. suite % ok) then
    error stop 1
  end if

contains

  function roundtrip_accuracy() result(res)
    type(test_result) :: res

    real(rk) :: max_reldiff, mean_reldiff
    real(rk), allocatable :: frequency(:), frequency_roundtrip(:)
    real(rk), allocatable :: k(:), reldiff(:)

    res % name = 'roundtrip_accuracy'
    res % ok = .true.

    frequency = logarithmic_frequencies()
    k = wavenumber(frequency, depth, water_density, gravity, surface_tension)
    frequency_roundtrip = angular_frequency(k, depth, water_density, gravity, surface_tension) / twopi
    reldiff = abs((frequency_roundtrip - frequency) / frequency * 100._rk)

    max_reldiff = maxval(reldiff)
    mean_reldiff = sum(reldiff) / real(num_frequencies, rk)

    if (max_reldiff > max_reldiff_pct_allowed) then
      write(stderr, '(a,es12.4,a)') &
        'test_dispersion: maximum relative error was ', max_reldiff, '%.'
      res % ok = .false.
    end if

    if (mean_reldiff > mean_reldiff_pct_allowed) then
      write(stderr, '(a,es12.4,a)') &
        'test_dispersion: mean relative error was ', mean_reldiff, '%.'
      res % ok = .false.
    end if
  end function roundtrip_accuracy


  function group_speed_is_finite_and_positive() result(res)
    type(test_result) :: res

    real(rk), allocatable :: cg(:), frequency(:), k(:)

    res % name = 'group_speed_is_finite_and_positive'

    frequency = logarithmic_frequencies()
    k = wavenumber(frequency, depth, water_density, gravity, surface_tension)
    cg = group_speed(k, depth, water_density, gravity, surface_tension)

    res = test(res % name, all(cg > 0._rk) .and. all(cg == cg))
    if (.not. res % ok) then
      write(stderr, '(a)') &
        'test_dispersion: group speed must be finite and positive.'
    end if
  end function group_speed_is_finite_and_positive


  function logarithmic_frequencies() result(frequency)
    real(rk) :: frequency(num_frequencies)

    integer :: o
    real(rk) :: dlnf

    dlnf = (log(frequency_max) - log(frequency_min)) / real(num_frequencies - 1, rk)
    do o = 1, num_frequencies
      frequency(o) = exp(log(frequency_min) + real(o - 1, rk) * dlnf)
    end do
  end function logarithmic_frequencies

end program test_dispersion
