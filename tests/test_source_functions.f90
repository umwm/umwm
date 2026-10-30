program test_source_functions
  use tuff, only: test, test_result
  use umwm_constants, only: twopi
#ifdef MPI
  use umwm_env, only: env_init, env_stop
#endif
  use umwm_source_functions, only: mss_moments, power, sds_coefficient, sin_d12, snl_d12, &
                                   snl_shift

  implicit none

  integer, parameter :: nf = 37, nd = 36

  type(test_result) :: suite

#ifdef MPI
  call env_init()
#endif
  suite = test('test_source_functions', [ &
    test(mss_moments_match_direct_sum), &
    test(power_matches_intrinsic), &
    test(wind_input_matches_reference), &
    test(downshifting_conserves_energy), &
    test(downshifting_shift_is_resolution_independent), &
    test(downshifting_matches_two_bin_weights) &
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


  function downshifting_conserves_energy() result(res)
    ! At coarse to fine resolutions, and for any k dk of the bins, the
    ! downshifting conserves the energy e k dk of the prognostic range, also
    ! where it would move energy below the lowest bin, apart from the energy
    ! that it moves in from the diagnostic bins.
    type(test_result) :: res
    real, parameter :: snl_fac = 5, dlnf(4) = [0.2, 0.11548, 0.05, 0.02]
    real :: e(nf,nd), sds(nf,nd), kdk(nf), snl(nf,nd), shift, frac, moved, inflow, net
    integer :: n, m, oc, o, p
    logical :: ok

    ok = .true.
    do n = 1, size(dlnf)
      shift = snl_shift / dlnf(n)
      m = int(shift)
      frac = shift - m
      oc = nf - m - 4
      call random_number(e)
      call random_number(sds)
      call random_number(kdk)
      kdk = kdk + 0.1
      call snl_d12(nf, nd, oc, e, kdk, snl_fac, m, frac, sds, snl)
      ok = ok .and. all(snl(oc+1:,:) == 0)
      do p = 1, nd
        moved = snl_fac * sum(sds(1:oc,p) * e(1:oc,p) * kdk(1:oc))
        inflow = snl_fac * (sum(sds(oc+1:oc+m,p) * e(oc+1:oc+m,p) * kdk(oc+1:oc+m)) &
                            + frac * sds(oc+m+1,p) * e(oc+m+1,p) * kdk(oc+m+1))
        net = 0
        do o = 1, oc
          net = net + snl(o,p) * kdk(o)
        end do
        ok = ok .and. abs(net - inflow) <= 1e-5 * (moved + inflow)
      end do
    end do

    res = test('downshifting_conserves_energy', ok)
  end function downshifting_conserves_energy


  function downshifting_shift_is_resolution_independent() result(res)
    ! The energy that the downshifting moves within the prognostic range
    ! moves on average by snl_shift in ln(f), at any resolution, so that the
    ! mean ln(f) of the spectrum changes at the same rate.
    type(test_result) :: res
    real, parameter :: snl_fac = 5, dlnf(5) = [0.2, 0.11548, 0.0832, 0.05, 0.02]
    integer, parameter :: nmax = 400
    real :: e(nmax,1), sds(nmax,1), kdk(nmax), snl(nmax,1), lnf(nmax)
    real :: shift, frac, moved, dmean
    integer :: n, m, oc, nb, o
    logical :: ok

    ok = .true.
    do n = 1, size(dlnf)
      shift = snl_shift / dlnf(n)
      m = int(shift)
      frac = shift - m
      ! bins from 0.03 to 3 Hz, a spectrum peaked at 0.3 Hz, and the
      ! dissipation of a saturated deep-water tail, none above oc
      nb = int(log(100.) / dlnf(n)) + 1
      oc = nb - m - 1
      lnf(:nb) = [(log(0.03) + (o - 1) * dlnf(n), o = 1, nb)]
      kdk(:nb) = exp(4 * lnf(:nb)) * dlnf(n)
      e(:nb,1) = exp(-(lnf(:nb) - log(0.3))**2 / 0.1) / kdk(:nb)
      sds(:nb,1) = exp(lnf(:nb))
      sds(oc+1:nb,1) = 0
      call snl_d12(nb, 1, oc, e(:nb,:), kdk(:nb), snl_fac, m, frac, sds(:nb,:), snl(:nb,:))
      moved = snl_fac * sum(sds(:oc,1) * e(:oc,1) * kdk(:oc))
      dmean = sum(snl(:oc,1) * kdk(:oc) * lnf(:oc))
      ok = ok .and. abs(dmean / moved + snl_shift) < 1e-4 * snl_shift
    end do

    res = test('downshifting_shift_is_resolution_independent', ok)
  end function downshifting_shift_is_resolution_independent


  function downshifting_matches_two_bin_weights() result(res)
    ! At the resolution for which it was calibrated, dlnf = ln(2 / 0.0313) / 36,
    ! the downshifting moves the energy to the next two lower bins with the
    ! weights of Donelan et al. (2012), b1 and b2 normalized to a sum of 1.
    type(test_result) :: res
    real :: dlnf, b1, b2, shift
    dlnf = log(2 / 0.0313) / 36
    b1 = exp(-16 * dlnf**2)
    b2 = exp(-64 * dlnf**2)
    shift = snl_shift / dlnf
    res = test('downshifting_matches_two_bin_weights', &
               int(shift) == 1 .and. abs((shift - 1) - b2 / (b1 + b2)) < 1e-5)
  end function downshifting_matches_two_bin_weights

end program test_source_functions
