program test_modulation
  use tuff, only: test, test_result
  use umwm_constants, only: stderr, twopi
#ifdef MPI
  use umwm_env, only: env_init, env_stop
#endif
  use umwm_modulation, only: breaking_modulation, breaking_steepness, frequency_ratio, &
                             input_steepness, longwave_moments, modulated_group_speed, &
                             modulation_factor, modulation_fit, num_moments, num_sums, &
                             phase_means, relative_speed_amplitude, separation_bins, sin_mod, &
                             sin_mod_bin, slope_variances, speed_ratio, speed_steepness, &
                             speed_variance, wavenumber_ratio
  use umwm_source_functions, only: mss_moments, sin_d12_bin

  implicit none

  integer, parameter :: nf = 20, nd = 36
  real, parameter :: dth = real(twopi) / nd, dlnf = 0.1155, g = 9.81
  real, parameter :: sds_power = 2.4

  type(test_result) :: suite

#ifdef MPI
  call env_init()
#endif
  suite = test('test_modulation', [ &
    test(crest_wavenumber_matches_paper), &
    test(frequency_is_invariant_at_first_order), &
    test(phase_means_match_reference), &
    test(polynomials_fit_phase_means), &
    test(column_input_matches_bins), &
    test(separation_bins_match_wavenumber_ratio), &
    test(moments_match_direct_sums), &
    test(group_speed_mean_matches_direct_mean), &
    test(single_long_wave_gives_its_steepness), &
    test(relative_speed_amplitude_of_single_long_wave), &
    test(phase_averaged_input_matches_quadrature), &
    test(breaking_modulation_is_one_without_longer_waves) &
  ])
#ifdef MPI
  call env_stop()
#endif
  if (.not. suite % ok) error stop 1

contains

  function crest_wavenumber_matches_paper() result(res)
    ! The wavenumber modulation on the long-wave crests exceeds the linear
    ! solution, eps, by 16.9, 38.5, 66.4 and 104 % at eps = 0.1, 0.2, 0.3
    ! and 0.4 (Curcic 2025, section 3.1).
    type(test_result) :: res
    real, parameter :: eps(4) = [0.1, 0.2, 0.3, 0.4], excess(4) = [0.169, 0.385, 0.664, 1.04]
    logical :: ok
    ok = all(abs((wavenumber_ratio(eps, 1.) - 1) / eps - 1 - excess) < 0.005)
    res = test('crest_wavenumber_matches_paper', ok)
    if (.not. ok) write(stderr, '(a,4f8.4)') 'test_modulation: crest excess over linear: ', &
      (wavenumber_ratio(eps, 1.) - 1) / eps - 1
  end function crest_wavenumber_matches_paper


  function frequency_is_invariant_at_first_order() result(res)
    ! The first-order wavenumber and gravity modulations cancel in the
    ! intrinsic frequency, which is modulated at second order, while the
    ! phase speed is modulated by -eps cos(psi) at first order.
    type(test_result) :: res
    real, parameter :: eps = 0.05
    real :: cospsi, freq_dev, speed_harmonic
    integer :: j
    freq_dev = 0
    do j = 0, 100
      cospsi = -1 + 0.02 * j
      freq_dev = max(freq_dev, abs(frequency_ratio(eps, cospsi) - 1))
    end do
    speed_harmonic = 0.5 * (speed_ratio(eps, 1.) - speed_ratio(eps, -1.))
    res = test('frequency_is_invariant_at_first_order', &
               freq_dev < eps**2 .and. abs(speed_harmonic / eps + 1) < 0.05)
    if (.not. res % ok) write(stderr, '(a,2es12.4)') &
      'test_modulation: frequency deviation, speed harmonic / eps: ', freq_dev, speed_harmonic / eps
  end function frequency_is_invariant_at_first_order


  function phase_means_match_reference() result(res)
    ! The phase means of the speed ratio, of k sigma, and of the breaking
    ! enhancement at eps = 0.1 and 0.2 match a double-precision quadrature,
    ! and 16 phases are converged to 1e-5 against 64.
    type(test_result) :: res
    real :: fc, fks, fds, fc64, fks64, fds64
    logical :: ok
    call phase_means(0.1, 16, sds_power, fc, fks, fds)
    ok = abs(fc - 0.998790) < 1e-4 .and. abs(fks - 1.008925) < 1e-4 .and. abs(fds - 1.136211) < 1e-4
    call phase_means(0.2, 16, sds_power, fc, fks, fds)
    ok = ok .and. abs(fc - 0.995637) < 1e-4 .and. abs(fks - 1.037975) < 1e-4 &
         .and. abs(fds - 1.562879) < 2e-4
    call phase_means(0.2, 64, sds_power, fc64, fks64, fds64)
    ok = ok .and. abs(fc - fc64) < 1e-5 .and. abs(fks - fks64) < 1e-5 .and. abs(fds - fds64) < 1e-5
    res = test('phase_means_match_reference', ok)
    if (.not. ok) write(stderr, '(a,6f10.6)') 'test_modulation: phase means at 0.2, 16 and 64: ', &
      fc, fks, fds, fc64, fks64, fds64
  end function phase_means_match_reference


  function polynomials_fit_phase_means() result(res)
    ! The polynomials in the squared steepness fit the phase means to
    ! better than 1e-3 over the range, are 1 at zero steepness, and are
    ! limited to the range.
    type(test_result) :: res
    real, parameter :: eps_max = 0.3
    real :: coef_c(4), coef_ks(4), coef_ds(4), fc, fks, fds, eps, worst
    integer :: j
    logical :: ok
    call modulation_fit(eps_max, 16, sds_power, coef_c, coef_ks, coef_ds)
    worst = 0
    do j = 0, 30
      eps = eps_max * (j + 0.5) / 31
      call phase_means(eps, 16, sds_power, fc, fks, fds)
      worst = max(worst, abs(fit_factor(coef_c, eps_max, eps**2) / fc - 1), &
                  abs(fit_factor(coef_ks, eps_max, eps**2) / fks - 1), abs(fit_factor(coef_ds, eps_max, eps**2) / fds - 1))
    end do
    call phase_means(eps_max, 16, sds_power, fc, fks, fds)
    ok = worst < 1e-3 .and. fit_factor(coef_c, eps_max, 0.) == 1 .and. fit_factor(coef_ds, eps_max, 0.) == 1 &
         .and. fit_factor(coef_ds, eps_max, -1.) == 1 .and. abs(fit_factor(coef_ds, eps_max, 1.) / fds - 1) < 1e-3 &
         .and. fds > 2 .and. fc < 1 .and. fks > 1
    res = test('polynomials_fit_phase_means', ok)
    if (.not. ok) write(stderr, '(a,es12.4)') 'test_modulation: worst relative fit error: ', worst
  end function polynomials_fit_phase_means


  real function fit_factor(coef, eps_max, x)
    ! A phase mean from its polynomial of degree 4 at the squared steepness x.
    real, intent(in) :: coef(4), eps_max, x
    fit_factor = modulation_factor(coef(1), coef(2), coef(3), coef(4), eps_max**2, x)
  end function fit_factor


  function column_input_matches_bins() result(res)
    ! The two-pass evaluation of the phase-averaged input over a column
    ! matches the evaluation bin by bin, with the modulation of the phase
    ! speed and of k sigma of each bin from the moments, for bins with and
    ! without a sign change of the relative speed over the phase, and the
    ! diagnostic bins not below zero.
    type(test_result) :: res
    integer, parameter :: oc = nf - 4
    real, parameter :: eps_max = 0.3, shelt = 0.09, sin_fac = 0.11, fieldscale1 = 0.1 / 0.11, &
                       fieldscale2 = 0.01, scale = 7e-3, wspd = 4.
    real :: e(nf,nd), k(nf), kdk(nf), k3dk(nf), worb(nf), wstk(nf), th(nd), moments(nf,num_moments)
    real :: mss(nf,3), esum(nf,num_sums)
    real :: coef_c(4), coef_ks(4), coef_ds(4), cp(nf), wind(nf), fkovg(nf), cosd(nd), ssin(nf,nd)
    real :: mu_g, mu_gk, mu_k, b, fc, fks, direct, worst
    integer :: o, p, changes
    call wave_grid(k, kdk, k3dk, worb, wstk, th)
    call modulation_fit(eps_max, 16, sds_power, coef_c, coef_ks, coef_ds)
    call random_number(e)
    e = 5e-3 * e**2 / spread(kdk, 2, nd) * 0.2 ! a sea of moderate steepness
    call longwave_moments(nf, nd, 1, e, k3dk, worb, wstk, cos(th), sin(th), cos(2 * th), &
                          sin(2 * th), cos(3 * th), sin(3 * th), cos(4 * th), sin(4 * th), dth, &
                          moments, mss, esum)
    cp = sqrt(g / k)
    wind = wspd
    fkovg = k * sqrt(g * k) / g
    cosd = cos(0.3 - th)
    call sin_mod(nf, nd, nf, oc, wind, cosd, cp, 0.1, -0.2, cos(th), sin(th), cos(2 * th), &
                 sin(2 * th), cos(3 * th), sin(3 * th), cos(4 * th), sin(4 * th), shelt, sin_fac, &
                 fieldscale1, fieldscale2, scale, fkovg, moments, coef_c, coef_ks, eps_max**2, &
                 ssin)
    worst = 0
    changes = 0
    do p = 1, nd
      do o = 1, nf
        call slope_variances(moments(o,1), moments(o,2), moments(o,3), moments(o,4), &
                             moments(o,5), cos(2 * th(p)), sin(2 * th(p)), cos(4 * th(p)), &
                             sin(4 * th(p)), mu_g, mu_gk, mu_k)
        fc = modulation_factor(coef_c(1), coef_c(2), coef_c(3), coef_c(4), eps_max**2, &
                               speed_steepness(mu_g, mu_gk, mu_k)**2)
        fks = modulation_factor(coef_ks(1), coef_ks(2), coef_ks(3), coef_ks(4), eps_max**2, &
                                input_steepness(mu_g, mu_gk, mu_k)**2)
        b = relative_speed_amplitude(mu_g, mu_gk, mu_k, moments(o,6), moments(o,7), &
                                     moments(o,8), moments(o,9), moments(o,10), moments(o,11), &
                                     moments(o,12), cos(th(p)), sin(th(p)), cos(2 * th(p)), &
                                     sin(2 * th(p)), cos(3 * th(p)), sin(3 * th(p)), cp(o))
        if (b > abs(wind(o) * cosd(p) - cp(o) * fc - 0.1 * cos(th(p)) + 0.2 * sin(th(p)))) &
          changes = changes + 1
        direct = sin_mod_bin(wind(o), cosd(p), cp(o) * fc, 0.1, -0.2, cos(th(p)), sin(th(p)), &
                             shelt, sin_fac, fieldscale1, fieldscale2, scale, fkovg(o), fks, &
                             merge(0., -huge(0.), o > oc), b)
        worst = max(worst, abs(ssin(o,p) - direct) / max(abs(direct), 1e-6))
      end do
    end do
    res = test('column_input_matches_bins', worst < 1e-5 .and. changes > 0 &
               .and. all(ssin(oc+1:,:) >= 0) .and. any(ssin(1:oc,:) < 0))
    if (.not. res % ok) write(stderr, '(a,es12.4,i6)') &
      'test_modulation: worst relative difference, bins with a sign change: ', worst, changes
  end function column_input_matches_bins


  function separation_bins_match_wavenumber_ratio() result(res)
    ! At the default resolution, a wavenumber ratio of 4 is 6 bins, a
    ! factor 2 in frequency; a ratio of 1 is the next lower bin.
    type(test_result) :: res
    real :: dlnf0
    dlnf0 = log(2 / 0.0313) / 36
    res = test('separation_bins_match_wavenumber_ratio', &
               separation_bins(4., dlnf0) == 6 .and. separation_bins(1., dlnf0) == 1 &
               .and. separation_bins(16., dlnf0) == 12 .and. separation_bins(1.5, dlnf0) == 2 &
               .and. separation_bins(4., log(14 / 0.04) / 51) == 6)
  end function separation_bins_match_wavenumber_ratio


  function moments_match_direct_sums() result(res)
    ! The slope variances, the orbital velocity variance and the Stokes
    ! projections from the moments match the direct sums over the longer
    ! waves and all pairs of directions.
    type(test_result) :: res
    integer, parameter :: sep = 2
    real :: e(nf,nd), k3dk(nf), worb(nf), wstk(nf), th(nd), moments(nf,num_moments)
    real :: mss(nf,3), mss_direct(nf,3), esum(nf,num_sums)
    real :: mu_g, mu_gk, mu_k, var_u, sk1, sk3, d_g, d_gk, d_k, d_u, d_sk1, d_sk3, delta, worst
    integer :: o, p, op, pp
    call random_number(e)
    e = 1e-2 * e**4
    th = [((p - 0.5 * (nd + 1)) * dth, p = 1, nd)]
    k3dk = [(1e-3 * o**2, o = 1, nf)]
    worb = [(0.5 * o, o = 1, nf)]
    wstk = [(0.1 * o**1.5, o = 1, nf)]
    call longwave_moments(nf, nd, sep, e, k3dk, worb, wstk, cos(th), sin(th), cos(2 * th), &
                          sin(2 * th), cos(3 * th), sin(3 * th), cos(4 * th), sin(4 * th), dth, &
                          moments, mss, esum)
    call mss_moments(nf, nd, e, k3dk, cos(2 * th), sin(2 * th), dth, mss_direct)
    worst = maxval(abs(mss - mss_direct)) / maxval(mss_direct)
    do o = 1, nf
      worst = max(worst, abs(esum(o,1) - sum(e(o,:))) / sum(e(o,:)), &
                  abs(esum(o,4) - sum(e(o,:) * cos(4 * th))) / sum(e(o,:)))
    end do
    do p = 1, nd
      do o = 1, nf
        call slope_variances(moments(o,1), moments(o,2), moments(o,3), moments(o,4), &
                             moments(o,5), cos(2 * th(p)), sin(2 * th(p)), cos(4 * th(p)), &
                             sin(4 * th(p)), mu_g, mu_gk, mu_k)
        var_u = 0.5 * (moments(o,6) + cos(2 * th(p)) * moments(o,7) + sin(2 * th(p)) * moments(o,8))
        sk1 = cos(th(p)) * moments(o,9) + sin(th(p)) * moments(o,10)
        sk3 = cos(3 * th(p)) * moments(o,11) + sin(3 * th(p)) * moments(o,12)
        d_g = 0
        d_gk = 0
        d_k = 0
        d_u = 0
        d_sk1 = 0
        d_sk3 = 0
        do op = 1, o - sep
          do pp = 1, nd
            delta = th(pp) - th(p)
            d_g = d_g + e(op,pp) * k3dk(op) * dth
            d_gk = d_gk + e(op,pp) * k3dk(op) * dth * cos(delta)**2
            d_k = d_k + e(op,pp) * k3dk(op) * dth * cos(delta)**4
            d_u = d_u + e(op,pp) * worb(op) * dth * cos(delta)**2
            d_sk1 = d_sk1 + e(op,pp) * wstk(op) * dth * cos(delta)
            d_sk3 = d_sk3 + e(op,pp) * wstk(op) * dth * cos(3 * delta)
          end do
        end do
        worst = max(worst, abs(mu_g - d_g), abs(mu_gk - d_gk), abs(mu_k - d_k), &
                    abs(var_u - d_u) / maxval(worb), abs(sk1 - d_sk1) / maxval(wstk), &
                    abs(sk3 - d_sk3) / maxval(wstk))
      end do
    end do
    worst = worst / maxval(moments(:,1))
    res = test('moments_match_direct_sums', worst < 1e-4 .and. all(moments(1:sep,:) == 0))
    if (.not. res % ok) write(stderr, '(a,es12.4)') 'test_modulation: worst relative error: ', worst
  end function moments_match_direct_sums


  function single_long_wave_gives_its_steepness() result(res)
    ! A single long wave of steepness eps gives the aligned short waves the
    ! equivalent steepness eps for the speed, the breaking and the input
    ! modulation, and the short waves across it eps/2, 3/8 eps and eps/2.
    type(test_result) :: res
    integer, parameter :: ol = 5, pl = 10
    real, parameter :: a = 2.
    real :: e(nf,nd), k(nf), kdk(nf), k3dk(nf), worb(nf), wstk(nf), th(nd), moments(nf,num_moments)
    real :: mss(nf,3), esum(nf,num_sums)
    real :: eps, mu_g, mu_gk, mu_k
    integer :: o, p
    logical :: ok
    call wave_grid(k, kdk, k3dk, worb, wstk, th)
    e = 0
    e(ol,pl) = a**2 / (2 * kdk(ol) * dth)
    eps = a * k(ol)
    call longwave_moments(nf, nd, 1, e, k3dk, worb, wstk, cos(th), sin(th), cos(2 * th), &
                          sin(2 * th), cos(3 * th), sin(3 * th), cos(4 * th), sin(4 * th), dth, &
                          moments, mss, esum)
    o = ol + 3
    p = pl
    call slope_variances(moments(o,1), moments(o,2), moments(o,3), moments(o,4), moments(o,5), &
                         cos(2 * th(p)), sin(2 * th(p)), cos(4 * th(p)), sin(4 * th(p)), &
                         mu_g, mu_gk, mu_k)
    ok = abs(speed_steepness(mu_g, mu_gk, mu_k) / eps - 1) < 1e-4 &
         .and. abs(breaking_steepness(mu_g, mu_gk, mu_k) / eps - 1) < 1e-4 &
         .and. abs(input_steepness(mu_g, mu_gk, mu_k) / eps - 1) < 1e-4
    p = pl + nd / 4
    call slope_variances(moments(o,1), moments(o,2), moments(o,3), moments(o,4), moments(o,5), &
                         cos(2 * th(p)), sin(2 * th(p)), cos(4 * th(p)), sin(4 * th(p)), &
                         mu_g, mu_gk, mu_k)
    ok = ok .and. abs(speed_steepness(mu_g, mu_gk, mu_k) / eps - 0.5) < 1e-4 &
         .and. abs(breaking_steepness(mu_g, mu_gk, mu_k) / eps - 0.375) < 1e-4 &
         .and. abs(input_steepness(mu_g, mu_gk, mu_k) / eps - 0.5) < 1e-4
    ! the long wave does not modulate itself or the bins within the separation
    ok = ok .and. all(moments(1:ol,:) == 0)
    res = test('single_long_wave_gives_its_steepness', ok)
  end function single_long_wave_gives_its_steepness


  function relative_speed_amplitude_of_single_long_wave() result(res)
    ! The orbital velocity of a single long wave of phase speed C_L, less
    ! the modulated phase speed cp, modulates the wind-wave relative speed
    ! of the aligned short waves by eps (C_L - cp) cos(psi), and that of
    ! the short waves across it by eps cp/2 cos(psi), from their phase speed
    ! alone.
    type(test_result) :: res
    integer, parameter :: ol = 5, pl = 10
    real, parameter :: a = 2.
    real :: e(nf,nd), k(nf), kdk(nf), k3dk(nf), worb(nf), wstk(nf), th(nd), moments(nf,num_moments)
    real :: mss(nf,3), esum(nf,num_sums)
    real :: eps, cl, cp, mu_g, mu_gk, mu_k, b_along, b_across
    integer :: o, p
    call wave_grid(k, kdk, k3dk, worb, wstk, th)
    e = 0
    e(ol,pl) = a**2 / (2 * kdk(ol) * dth)
    eps = a * k(ol)
    cl = sqrt(g / k(ol))
    call longwave_moments(nf, nd, 1, e, k3dk, worb, wstk, cos(th), sin(th), cos(2 * th), &
                          sin(2 * th), cos(3 * th), sin(3 * th), cos(4 * th), sin(4 * th), dth, &
                          moments, mss, esum)
    o = ol + 4
    cp = sqrt(g / k(o))
    p = pl
    call slope_variances(moments(o,1), moments(o,2), moments(o,3), moments(o,4), moments(o,5), &
                         cos(2 * th(p)), sin(2 * th(p)), cos(4 * th(p)), sin(4 * th(p)), &
                         mu_g, mu_gk, mu_k)
    b_along = relative_speed_amplitude(mu_g, mu_gk, mu_k, moments(o,6), moments(o,7), &
                                       moments(o,8), moments(o,9), moments(o,10), moments(o,11), &
                                       moments(o,12), cos(th(p)), sin(th(p)), cos(2 * th(p)), &
                                       sin(2 * th(p)), cos(3 * th(p)), sin(3 * th(p)), cp)
    p = pl + nd / 4
    call slope_variances(moments(o,1), moments(o,2), moments(o,3), moments(o,4), moments(o,5), &
                         cos(2 * th(p)), sin(2 * th(p)), cos(4 * th(p)), sin(4 * th(p)), &
                         mu_g, mu_gk, mu_k)
    b_across = relative_speed_amplitude(mu_g, mu_gk, mu_k, moments(o,6), moments(o,7), &
                                        moments(o,8), moments(o,9), moments(o,10), moments(o,11), &
                                        moments(o,12), cos(th(p)), sin(th(p)), cos(2 * th(p)), &
                                        sin(2 * th(p)), cos(3 * th(p)), sin(3 * th(p)), cp)
    res = test('relative_speed_amplitude_of_single_long_wave', &
               abs(b_along / (eps * (cl - cp)) - 1) < 1e-3 .and. abs(b_across / (0.5 * eps * cp) - 1) < 1e-3)
    if (.not. res % ok) write(stderr, '(a,4f10.5)') &
      'test_modulation: b along and across, expected: ', b_along, b_across, eps * (cl - cp), 0.5 * eps * cp
  end function relative_speed_amplitude_of_single_long_wave


  function phase_averaged_input_matches_quadrature() result(res)
    ! The closed-form mean of the wind input over the long-wave phase
    ! matches a quadrature of the input of Donelan et al. (2012) with the
    ! phase speed modulated by b cos(psi), whether the relative speed keeps
    ! its sign over the phase or not, for following and opposing winds and
    ! for swell overrunning the wind.
    type(test_result) :: res
    integer, parameter :: nq = 4000
    real, parameter :: cp = 5., shelt = 0.09, sin_fac = 0.11, fieldscale1 = 0.1 / 0.11, &
                       fieldscale2 = 0.01, scale = 7e-3, fkovg = 0.8
    real, parameter :: a(6) = [10., -3., 1., -1., 0., -2.], b(6) = [2., 1., 3., 3., 2., 2.5]
    real, parameter :: cosd(6) = [1., 1., 1., -1., 1., -1.]
    real :: wind, direct, mean, worst, cospsi
    integer :: n, j
    worst = 0
    do n = 1, 6
      wind = (a(n) + cp) / cosd(n)
      direct = sin_mod_bin(wind, cosd(n), cp, 0., 0., 1., 0., shelt, sin_fac, fieldscale1, &
                           fieldscale2, scale, fkovg, 1., -huge(0.), b(n))
      mean = 0
      do j = 1, nq
        cospsi = cos((j - 0.5) * real(twopi) / nq)
        mean = mean + sin_d12_bin(wind, cosd(n), cp + b(n) * cospsi, 0., 0., 1., 0., shelt, &
                                  sin_fac, fieldscale1, fieldscale2, scale, fkovg, -huge(0.))
      end do
      mean = mean / nq
      worst = max(worst, abs(direct - mean) / (scale * fkovg * (a(n)**2 + 0.5 * b(n)**2)))
      ! without modulation, the input of Donelan et al. (2012)
      worst = max(worst, abs(sin_mod_bin(wind, cosd(n), cp, 0., 0., 1., 0., shelt, sin_fac, &
                                         fieldscale1, fieldscale2, scale, fkovg, 1., -huge(0.), 0.) &
                             - sin_d12_bin(wind, cosd(n), cp, 0., 0., 1., 0., shelt, sin_fac, &
                                           fieldscale1, fieldscale2, scale, fkovg, -huge(0.))) &
                         / (scale * fkovg * max(a(n)**2, 1.)))
    end do
    res = test('phase_averaged_input_matches_quadrature', worst < 1e-4)
    if (.not. res % ok) write(stderr, '(a,es12.4)') 'test_modulation: worst relative error: ', worst
  end function phase_averaged_input_matches_quadrature


  function breaking_modulation_is_one_without_longer_waves() result(res)
    ! Without longer waves the breaking enhancement is 1, and it exceeds 1
    ! with them.
    type(test_result) :: res
    real, parameter :: eps_max = 0.3
    real :: coef_c(4), coef_ks(4), coef_ds(4), moments(nf,num_moments), f(nf,nd), th(nd)
    integer :: p
    th = [((p - 0.5 * (nd + 1)) * dth, p = 1, nd)]
    call modulation_fit(eps_max, 16, sds_power, coef_c, coef_ks, coef_ds)
    moments = 0
    f = 0
    call breaking_modulation(nf, nd, 1, nf, moments, cos(2 * th), sin(2 * th), cos(4 * th), &
                             sin(4 * th), coef_ds, eps_max**2, f)
    res = test('breaking_modulation_is_one_without_longer_waves', all(f == 1))
    moments(:,1) = 0.02 ! an omnidirectional mean-squared slope of 0.02 of the longer waves
    call breaking_modulation(nf, nd, 1, nf, moments, cos(2 * th), sin(2 * th), cos(4 * th), &
                             sin(4 * th), coef_ds, eps_max**2, f)
    ! the enhancement of the bins written out matches the phase mean at the
    ! breaking steepness of the slope variances, mu_g = 0.02, mu_gk = 0.01
    ! and mu_k = 0.0075 for an omnidirectional slope
    res = test('breaking_modulation_is_one_without_longer_waves', res % ok .and. all(f > 1) &
               .and. abs(f(nf,1) / fit_factor(coef_ds, eps_max, &
                                               breaking_steepness(0.02, 0.01, 0.0075)**2) - 1) < 1e-5)
  end function breaking_modulation_is_one_without_longer_waves


  function group_speed_mean_matches_direct_mean() result(res)
    ! The modulated group speed from the direction sums matches the
    ! polynomial at the mean over directions, weighted by the spectrum, of
    ! the squared equivalent steepness of the speed modulation, evaluated
    ! directly over the directions.
    type(test_result) :: res
    real, parameter :: eps_max = 0.3
    real :: e(nf,nd), k(nf), kdk(nf), k3dk(nf), worb(nf), wstk(nf), th(nd), moments(nf,num_moments)
    real :: mss(nf,3), esum(nf,num_sums), coef_c(4), coef_ks(4), coef_ds(4), cg0(nf), cgm(nf)
    real :: mu_g, mu_gk, mu_k, x, w, worst
    integer :: o, p
    call wave_grid(k, kdk, k3dk, worb, wstk, th)
    call modulation_fit(eps_max, 16, sds_power, coef_c, coef_ks, coef_ds)
    call random_number(e)
    e = 5e-3 * e**2 / spread(kdk, 2, nd) * 0.2
    call longwave_moments(nf, nd, 1, e, k3dk, worb, wstk, cos(th), sin(th), cos(2 * th), &
                          sin(2 * th), cos(3 * th), sin(3 * th), cos(4 * th), sin(4 * th), dth, &
                          moments, mss, esum)
    cg0 = 0.5 * sqrt(g / k)
    call modulated_group_speed(nf, esum, moments, coef_c, eps_max**2, cg0, cgm)
    worst = 0
    do o = 1, nf
      x = 0
      w = 0
      do p = 1, nd
        call slope_variances(moments(o,1), moments(o,2), moments(o,3), moments(o,4), &
                             moments(o,5), cos(2 * th(p)), sin(2 * th(p)), cos(4 * th(p)), &
                             sin(4 * th(p)), mu_g, mu_gk, mu_k)
        x = x + e(o,p) * speed_variance(mu_g, mu_gk, mu_k)
        w = w + e(o,p)
      end do
      worst = max(worst, abs(cgm(o) / (cg0(o) * fit_factor(coef_c, eps_max, x / w)) - 1))
    end do
    res = test('group_speed_mean_matches_direct_mean', worst < 1e-5 .and. all(cgm(8:) < cg0(8:)) &
               .and. all(cgm(1:1) == cg0(1:1)))
    if (.not. res % ok) write(stderr, '(a,es12.4)') 'test_modulation: worst relative error: ', worst
  end function group_speed_mean_matches_direct_mean


  subroutine wave_grid(k, kdk, k3dk, worb, wstk, th)
    ! Deep-water wavenumbers and the weights of the moments on a
    ! logarithmic frequency grid, as in the model.
    real, intent(out) :: k(nf), kdk(nf), k3dk(nf), worb(nf), wstk(nf), th(nd)
    real :: sigma
    integer :: o, p
    do o = 1, nf
      sigma = real(twopi) * 0.1 * exp((o - 1) * dlnf)
      k(o) = sigma**2 / g
      kdk(o) = 2 * k(o)**2 * dlnf ! k dk, dk = dsigma / cg
      k3dk(o) = k(o)**2 * kdk(o)
      worb(o) = sigma**2 * kdk(o)
      wstk(o) = 2 * sigma * k(o) * kdk(o)
    end do
    th = [((p - 0.5 * (nd + 1)) * dth, p = 1, nd)]
  end subroutine wave_grid

end program test_modulation
