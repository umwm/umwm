program test_source_functions
  use tuff, only: test, test_result
  use umwm_constants, only: twopi
#ifdef MPI
  use umwm_env, only: env_init, env_stop
#endif
  use umwm_source_functions, only: mss_moments, power, sds_coefficient, sin_d12

  implicit none

  integer, parameter :: nf = 37, nd = 36

  type(test_result) :: suite

#ifdef MPI
  call env_init()
#endif
  suite = test('test_source_functions', [ &
    test(mss_moments_match_direct_sum), &
    test(power_matches_intrinsic), &
    test(wind_input_matches_reference) &
  ])

#ifdef MPI
  call env_stop()
#endif
  if (.not. suite % ok) error stop 1

contains

  function mss_moments_match_direct_sum() result(res)
    ! The breaking dissipation coefficient from the mean-squared slope
    ! moments matches the direct sum over all pairs of directions.
    type(test_result) :: res
    real :: e(nf,nd), k3dk(nf), f(nf), th(nd), moments(nf,3)
    real :: coef(nf,nd), slope, coef_direct
    real, parameter :: mss_fac = 360, twopisds_fac = twopi * 42
    real :: dth, max_reldiff
    integer :: o, op, p, pp

    dth = twopi / nd
    th = [((p - 0.5 * (nd + 1)) * dth, p = 1, nd)]
    f = [(0.0313 * exp((o - 1) * 0.1155), o = 1, nf)]
    k3dk = [(1e-3 * o**2, o = 1, nf)]
    call random_number(e)
    e = 1e-2 * e**4

    call mss_moments(nf, nd, e, k3dk, cos(2 * th), sin(2 * th), dth, moments)
    call sds_coefficient(nf, nd, 1, nf, f, cos(2 * th), sin(2 * th), mss_fac, twopisds_fac, &
                         moments, coef)

    max_reldiff = 0
    do p = 1, nd
      do o = 1, nf
        slope = 0
        do op = 1, o - 1
          do pp = 1, nd
            slope = slope + e(op,pp) * cos(dth * (pp - p))**2 * dth * k3dk(op)
          end do
        end do
        coef_direct = twopisds_fac * f(o) * (1 + mss_fac * slope)**2
        max_reldiff = max(max_reldiff, abs(coef(o,p) - coef_direct) / coef_direct)
      end do
    end do

    res = test('mss_moments_match_direct_sum', max_reldiff < 1e-5)
  end function mss_moments_match_direct_sum


  function power_matches_intrinsic() result(res)
    ! power agrees with x**a, and underflowing results are zero.
    type(test_result) :: res
    integer, parameter :: n = 1000
    real :: x(n), y(n), a
    logical :: ok
    integer :: j

    ok = .true.
    do j = 1, 2
      a = merge(2.4, 1 / 2.4, j == 1)
      call random_number(x)
      x = 10.**(-15 + 15 * x)
      y = x
      call power(n, y, a)
      ok = ok .and. all(abs(y - x**a) <= 4 * spacing(x**a))
    end do

    x = [0., tiny(1.), 1e-20, (1., j = 4, n)]
    call power(n, x, 2.4)
    ok = ok .and. all(x(1:3) == 0) .and. all(x(4:) == 1)

    res = test('power_matches_intrinsic', ok)
  end function power_matches_intrinsic


  function wind_input_matches_reference() result(res)
    ! sin_d12 matches a step-by-step evaluation of the wind input function.
    type(test_result) :: res
    real :: wind(nf), cosd(nd), cp(nf), cth(nd), sth(nd), fkovg(nf), th(nd)
    real :: ssin(nf,nd), s, reference, terms
    real, parameter :: uc = 0.3, vc = -0.2, shelt = 0.09, sin_fac = 0.11
    real, parameter :: fieldscale1 = 0.1 / 0.11, fieldscale2 = 0.001 / 0.1, scale = 0.007
    integer, parameter :: oc = 25
    logical :: ok
    integer :: o, p

    th = [((p - 0.5 * (nd + 1)) * twopi / nd, p = 1, nd)]
    cth = cos(th)
    sth = sin(th)
    cosd = cos(0.3 - th)
    cp = [(9.81 / (twopi * 0.0313 * exp((o - 1) * 0.1155)), o = 1, nf)]
    wind = [(12 + 0.1 * o, o = 1, nf)]
    fkovg = [(1e-3 * o, o = 1, nf)]

    call sin_d12(nf, nd, nf, oc, wind, cosd, cp, uc, vc, cth, sth, shelt, sin_fac, &
                 fieldscale1, fieldscale2, scale, fkovg, ssin)

    ok = .true.
    do p = 1, nd
      do o = 1, nf
        s = wind(o) * cosd(p) - cp(o) - uc * cth(p) - vc * sth(p)
        ! the relative wind s cancels where the wind matches the phase speed, so
        ! rounding error is relative to the magnitude of its terms, not to s
        terms = abs(wind(o) * cosd(p)) + cp(o) + abs(uc * cth(p)) + abs(vc * sth(p))
        s = sin_fac * abs(s) * s
        if (s > 0) s = s * shelt / sin_fac
        if (s < 0) s = s * fieldscale1
        if (s < 0 .and. cosd(p) > 0) s = s * fieldscale2
        reference = scale * s * fkovg(o)
        if (o > oc) reference = max(reference, 0.)
        ok = ok .and. abs(ssin(o,p) - reference) <= 1e-6 * (abs(reference) + scale * sin_fac * terms**2 * fkovg(o))
      end do
    end do

    ok = ok .and. any(ssin(1:oc,:) < 0) .and. all(ssin(oc+1:,:) >= 0)

    res = test('wind_input_matches_reference', ok)
  end function wind_input_matches_reference

end program test_source_functions
