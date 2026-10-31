module umwm_modulation
  ! Hydrodynamic modulation of short waves by longer waves, after Curcic
  ! (2025, J. Fluid Mech. 1015, A20): the steady solutions of the crest and
  ! action conservation of a short wave riding on a long wave of steepness
  ! eps, as functions of the long-wave phase psi, their means over the phase,
  ! the moments of the longer-wave spectrum from which each spectral bin gets
  ! the steepness of the equivalent single long wave, and the wind input
  ! averaged over the long-wave phase.
  !
  ! The model carries the phase-averaged spectrum on a grid of intrinsic
  ! frequencies, which the modulation leaves invariant at first order for
  ! short waves aligned with the long waves: the modulation lives within the
  ! bins. For a long wave at the angle dth to a short wave, the first-order
  ! modulations over the long-wave phase are
  !
  !   dk/k         =  eps cos(psi) cos^2(dth)   strain of the orbital velocity along the short wave
  !   dg/g         = -eps cos(psi)              vertical orbital acceleration
  !   dN/N         =  eps cos(psi)              convergence of the orbital velocity
  !   dsigma/sigma = -(eps/2) cos(psi) sin^2(dth)
  !   dC/C = dCg/Cg = -(eps/2) cos(psi) (1 + cos^2(dth))
  !   U.k          =  eps C_L cos(psi) cos(dth) orbital advection along the short wave
  !
  ! For a spectrum of long waves with independent phases these add, and the
  ! variance of each is a cumulative moment of the longer-wave spectrum
  ! weighted by a power of cos(dth): a single long wave of steepness eps has
  ! the variance eps^2/2 in each. The nonlinear one-dimensional solutions of
  ! the paper are evaluated at the equivalent steepness that gives the same
  ! first-order variance as the spectrum. The phase means of the paper's
  ! solutions are second order in eps, and the quantities that are nonlinear
  ! in the modulated state, the wind input and the breaking dissipation,
  ! respond to the first-order modulation at second order.
  !
  ! The short-wave action, concentrated on the long-wave crests where the
  ! orbital velocity is forward, is transported by the longer waves at their
  ! surface Stokes drift, the phase mean of the orbital advection of the
  ! action, as in the generalized Lagrangian mean wave action equation.

  use, intrinsic :: iso_fortran_env, only: real64

  implicit none

  private
  public :: action_ratio, breaking_modulation, breaking_steepness, breaking_variance, &
            gravity_ratio, frequency_ratio, input_steepness, input_variance, longwave_moments, &
            modulated_group_speed, modulation_factor, modulation_fit, modulation_steepness, &
            num_coef, num_moments, num_sums, phase_means, relative_speed_amplitude, &
            separation_bins, sin_mod, sin_mod_bin, slope_variances, speed_ratio, &
            speed_steepness, speed_variance, steepness_ratio, wavenumber_ratio

  integer, parameter :: num_moments = 12 ! moments of the longer-wave spectrum, see longwave_moments
  integer, parameter :: num_sums = 5     ! direction sums of the spectrum, see longwave_moments
  integer, parameter :: num_coef = 4     ! coefficients of the polynomials of the phase means

  real, parameter :: pi = 3.14159265358979, twopi = 2 * pi

contains

  pure elemental real function wavenumber_ratio(eps, cospsi) result(res)
    ! The wavenumber of a short wave at the phase psi of a long wave of
    ! steepness eps over its unmodulated wavenumber, from the crest
    ! conservation with the orbital velocity at the long-wave surface, (3.2):
    ! exp(x e^x), x = eps cos(psi).
    !$acc routine seq
    real, intent(in) :: eps, cospsi
    real :: x
    x = eps * cospsi
    res = exp(x * exp(x))
  end function wavenumber_ratio


  pure elemental real function action_ratio(eps, cospsi) result(res)
    ! The action of a short wave at the long-wave phase over its unmodulated
    ! action, from the action conservation with the convergence of the
    ! orbital velocity at the long-wave surface, (3.15): the same as the
    ! wavenumber ratio.
    !$acc routine seq
    real, intent(in) :: eps, cospsi
    res = wavenumber_ratio(eps, cospsi)
  end function action_ratio


  pure elemental real function gravity_ratio(eps, cospsi) result(res)
    ! The effective gravity of a short wave at the long-wave phase over g:
    ! the Lagrangian gravity at the long-wave surface, (3.4) with (3.8) and
    ! (3.10), 1 - x e^x + x^2 e^{2x}, x = eps cos(psi), reduced on the crests
    ! by the vertical orbital acceleration.
    !$acc routine seq
    real, intent(in) :: eps, cospsi
    real :: x, ex
    x = eps * cospsi
    ex = exp(x)
    res = 1 - x * ex + x * x * ex * ex
  end function gravity_ratio


  pure elemental real function frequency_ratio(eps, cospsi) result(res)
    ! The intrinsic frequency of a short wave at the long-wave phase over
    ! its unmodulated frequency, sqrt(g~ k~ / (g k)): the first-order
    ! wavenumber and gravity modulations cancel, and the frequency is
    ! modulated at second order.
    !$acc routine seq
    real, intent(in) :: eps, cospsi
    res = sqrt(gravity_ratio(eps, cospsi) * wavenumber_ratio(eps, cospsi))
  end function frequency_ratio


  pure elemental real function speed_ratio(eps, cospsi) result(res)
    ! The phase speed, and in deep water also the group speed, of a short
    ! wave at the long-wave phase over its unmodulated value,
    ! sqrt((g~/g) / (k~/k)): reduced by eps cos(psi) on the long-wave crests.
    !$acc routine seq
    real, intent(in) :: eps, cospsi
    res = sqrt(gravity_ratio(eps, cospsi) / wavenumber_ratio(eps, cospsi))
  end function speed_ratio


  pure elemental real function steepness_ratio(eps, cospsi) result(res)
    ! The steepness of a short wave at the long-wave phase over its
    ! unmodulated steepness, (3.18): (k~/k)^(5/4) (N~/N)^(1/2) (g~/g)^(-1/4),
    ! about 1 + 2 eps cos(psi) at first order.
    !$acc routine seq
    real, intent(in) :: eps, cospsi
    res = wavenumber_ratio(eps, cospsi)**1.75 * gravity_ratio(eps, cospsi)**(-0.25)
  end function steepness_ratio


  pure subroutine phase_means(eps, nphase, sds_power, fc, fks, fds)
    ! Means over the long-wave phase of the speed ratio (fc), of the product
    ! of the wavenumber and frequency ratios (fks), which scales the wind
    ! input, and of the breaking enhancement (fds), the mean of
    ! (s^2 / <s^2>)^n of the steepness ratio s with the saturation exponent
    ! n = sds_power, for the modulated steepness squared over its phase mean
    ! as the local saturation of the bin. Midpoint rule over nphase phases,
    ! exact for trigonometric polynomials of degree below nphase in psi.
    real, intent(in) :: eps, sds_power
    integer, intent(in) :: nphase
    real, intent(out) :: fc, fks, fds
    real :: cospsi, s2, s2_mean, s2n_mean
    integer :: j

    fc = 0
    fks = 0
    s2_mean = 0
    s2n_mean = 0
    do j = 1, nphase
      cospsi = cos((j - 0.5) * twopi / nphase)
      fc = fc + speed_ratio(eps, cospsi)
      fks = fks + wavenumber_ratio(eps, cospsi) * frequency_ratio(eps, cospsi)
      s2 = steepness_ratio(eps, cospsi)**2
      s2_mean = s2_mean + s2
      s2n_mean = s2n_mean + s2**sds_power
    end do
    fc = fc / nphase
    fks = fks / nphase
    s2_mean = s2_mean / nphase
    s2n_mean = s2n_mean / nphase
    fds = s2n_mean / s2_mean**sds_power

  end subroutine phase_means


  pure subroutine modulation_fit(eps_max, nphase, sds_power, coef_c, coef_ks, coef_ds)
    ! Polynomials of degree num_coef in t = (eps / eps_max)^2 for the phase
    ! means, 1 + sum of coef(i) t^i, fitted by least squares to phase_means
    ! at 65 steepnesses from 0 to eps_max, for modulation_factor: the means
    ! are smooth even functions of the steepness, and the polynomials need
    ! neither a square root nor a table; they fit the means to 1e-4.
    integer, intent(in) :: nphase
    real, intent(in) :: eps_max, sds_power
    real, intent(out) :: coef_c(num_coef), coef_ks(num_coef), coef_ds(num_coef)
    integer, parameter :: nodes = 64
    real :: fc, fks, fds
    real(real64) :: t, a(num_coef,num_coef), rhs(num_coef,3), x(num_coef,3)
    integer :: j, i, l

    a = 0
    rhs = 0
    do j = 0, nodes
      t = (real(j, real64) / nodes)**2
      call phase_means(eps_max * j / nodes, nphase, sds_power, fc, fks, fds)
      do i = 1, num_coef
        do l = 1, num_coef
          a(i,l) = a(i,l) + t**(i + l)
        end do
        rhs(i,1) = rhs(i,1) + t**i * (fc - 1)
        rhs(i,2) = rhs(i,2) + t**i * (fks - 1)
        rhs(i,3) = rhs(i,3) + t**i * (fds - 1)
      end do
    end do
    call solve(num_coef, 3, a, rhs, x)
    coef_c = real(x(:,1))
    coef_ks = real(x(:,2))
    coef_ds = real(x(:,3))

  end subroutine modulation_fit


  pure subroutine solve(n, m, a, b, x)
    ! Solves a x = b for the m right-hand sides, by Gaussian elimination
    ! with partial pivoting.
    integer, intent(in) :: n, m
    real(real64), intent(in) :: a(n,n), b(n,m)
    real(real64), intent(out) :: x(n,m)
    real(real64) :: w(n,n+m), row(n+m), factor
    integer :: i, j, p
    w(:,1:n) = a
    w(:,n+1:) = b
    do j = 1, n
      p = j - 1 + maxloc(abs(w(j:n,j)), 1)
      row = w(p,:)
      w(p,:) = w(j,:)
      w(j,:) = row
      do i = j + 1, n
        factor = w(i,j) / w(j,j)
        w(i,:) = w(i,:) - factor * w(j,:)
      end do
    end do
    do j = n, 1, -1
      x(j,:) = (w(j,n+1:) - matmul(w(j,j+1:n), x(j+1:n,:))) / w(j,j)
    end do
  end subroutine solve


  pure elemental real function modulation_factor(c1, c2, c3, c4, xmax, x) result(res)
    ! A phase mean from its polynomial of modulation_fit with the
    ! coefficients c1 to c4 at the squared equivalent steepness x, limited
    ! to xmax = eps_max^2.
    !$acc routine seq
    real, intent(in) :: c1, c2, c3, c4, xmax, x
    real :: t
    t = min(max(x, 0.), xmax) / xmax
    res = 1 + t * (c1 + t * (c2 + t * (c3 + c4 * t)))
  end function modulation_factor


  pure integer function separation_bins(kratio, dlnf) result(res)
    ! The number of bins by which the longer waves of the moments are at
    ! least below a bin, so that their wavenumber is at most 1/kratio of the
    ! bin's, within a twentieth of a bin, with the deep-water wavenumber
    ! k ~ f^2: at least one bin, the next lower bin, as the mean-squared
    ! slope of the breaking dissipation. The default resolution, 37 bins
    ! from 0.0313 to 2 Hz, has a factor 2 in frequency, 4 in wavenumber,
    ! in 6 bins.
    real, intent(in) :: kratio, dlnf
    res = max(1, ceiling(log(max(kratio, 1.)) / (2 * dlnf) - 0.05))
  end function separation_bins


  pure subroutine longwave_moments(nf, nd, sep, e, k3dk, worb, wstk, cth, sth, cth2, sth2, &
                                   cth3, sth3, cth4, sth4, dth, moments, mss, esum)
    ! Cumulative moments of the spectrum of the longer waves, at least sep
    ! bins below each bin, over all directions th':
    !   1-5  sum of e k^3 dk dth {1, cos 2th', sin 2th', cos 4th', sin 4th'}: the slope,
    !        of which 1 is the cumulative mean-squared slope, and 1-3 the slope projected
    !        on a direction,
    !   6-8  sum of e worb dth {1, cos 2th', sin 2th'}: the surface orbital velocity,
    !        worb = sigma^2 coth^2(kd) k dk,
    !   9-12 sum of e wstk dth {cos th', sin th', cos 3th', sin 3th'}: the surface Stokes
    !        drift, wstk = sigma k^2 cosh(2kd) / sinh^2(kd) dk, of which 9 and 10 are its
    !        x- and y-components.
    ! Also the moments of the mean-squared slope of the breaking dissipation,
    ! from the next lower bin, as mss_moments of umwm_source_functions, and
    ! the direction sums of the spectrum of each bin, esum, with the weights
    ! {1, cos 2th', sin 2th', cos 4th', sin 4th'}, for the means over
    ! direction of modulated_group_speed.
    integer, intent(in) :: nf, nd, sep
    real, intent(in) :: e(nf,nd), k3dk(nf), worb(nf), wstk(nf)
    real, intent(in) :: cth(nd), sth(nd), cth2(nd), sth2(nd), cth3(nd), sth3(nd), cth4(nd), sth4(nd)
    real, intent(in) :: dth
    real, intent(out) :: moments(nf,num_moments), mss(nf,3), esum(nf,num_sums)
    real :: d(nf,num_moments)
    real, dimension(nf) :: d0, dc1, ds1, dc2, ds2, dc3, ds3, dc4, ds4 ! direction sums of the spectrum
    integer :: o, p

    d0 = 0
    dc1 = 0
    ds1 = 0
    dc2 = 0
    ds2 = 0
    dc3 = 0
    ds3 = 0
    dc4 = 0
    ds4 = 0
    do p = 1, nd
      do o = 1, nf
        d0(o) = d0(o) + e(o,p)
        dc1(o) = dc1(o) + e(o,p) * cth(p)
        ds1(o) = ds1(o) + e(o,p) * sth(p)
        dc2(o) = dc2(o) + e(o,p) * cth2(p)
        ds2(o) = ds2(o) + e(o,p) * sth2(p)
        dc3(o) = dc3(o) + e(o,p) * cth3(p)
        ds3(o) = ds3(o) + e(o,p) * sth3(p)
        dc4(o) = dc4(o) + e(o,p) * cth4(p)
        ds4(o) = ds4(o) + e(o,p) * sth4(p)
      end do
    end do
    do o = 1, nf
      esum(o,1) = d0(o)
      esum(o,2) = dc2(o)
      esum(o,3) = ds2(o)
      esum(o,4) = dc4(o)
      esum(o,5) = ds4(o)
      d(o,1) = d0(o) * k3dk(o) * dth
      d(o,2) = dc2(o) * k3dk(o) * dth
      d(o,3) = ds2(o) * k3dk(o) * dth
      d(o,4) = dc4(o) * k3dk(o) * dth
      d(o,5) = ds4(o) * k3dk(o) * dth
      d(o,6) = d0(o) * worb(o) * dth
      d(o,7) = dc2(o) * worb(o) * dth
      d(o,8) = ds2(o) * worb(o) * dth
      d(o,9) = dc1(o) * wstk(o) * dth
      d(o,10) = ds1(o) * wstk(o) * dth
      d(o,11) = dc3(o) * wstk(o) * dth
      d(o,12) = ds3(o) * wstk(o) * dth
    end do

    moments(1:min(sep, nf),:) = 0
    do o = sep + 1, nf
      moments(o,:) = moments(o-1,:) + d(o-sep,:)
    end do

    mss(1,:) = 0
    do o = 2, nf
      mss(o,:) = mss(o-1,:) + 0.5 * d(o-1,1:3)
    end do

  end subroutine longwave_moments


  pure elemental subroutine slope_variances(m1, m2, m3, m4, m5, cth2, sth2, cth4, sth4, &
                                            mu_g, mu_gk, mu_k)
    ! The variances of the first-order modulations of a short wave of
    ! direction th by the longer waves, from the slope moments 1-5 of
    ! longwave_moments: of the gravity and the action (mu_g, the cumulative
    ! mean-squared slope), of the wavenumber (mu_k, weighted by cos^4 of the
    ! angle to the short wave), and minus their covariance (mu_gk, weighted
    ! by cos^2), with cos^2 x = (1 + cos 2x)/2 and
    ! cos^4 x = (3 + 4 cos 2x + cos 4x)/8.
    !$acc routine seq
    real, intent(in) :: m1, m2, m3, m4, m5, cth2, sth2, cth4, sth4
    real, intent(out) :: mu_g, mu_gk, mu_k
    real :: p2, p4
    p2 = cth2 * m2 + sth2 * m3
    p4 = cth4 * m4 + sth4 * m5
    mu_g = m1
    mu_gk = 0.5 * (m1 + p2)
    mu_k = 0.125 * (3 * m1 + 4 * p2 + p4)
  end subroutine slope_variances


  pure elemental real function speed_variance(mu_g, mu_gk, mu_k) result(res)
    ! The squared steepness of the equivalent long wave for the modulation
    ! of the phase and group speeds, dC/C = -(eps/2) cos(psi) (1 + cos^2 dth),
    ! with the same variance, (mu_g + 2 mu_gk + mu_k)/4, as the spectrum.
    !$acc routine seq
    real, intent(in) :: mu_g, mu_gk, mu_k
    res = 0.5 * (mu_g + 2 * mu_gk + mu_k)
  end function speed_variance


  pure elemental real function breaking_variance(mu_g, mu_gk, mu_k) result(res)
    ! The squared steepness of the equivalent long wave for the modulation
    ! of the short-wave steepness, eps cos(psi) (5/4 cos^2 dth + 3/4) at
    ! first order from (3.18), with the same variance as the spectrum.
    !$acc routine seq
    real, intent(in) :: mu_g, mu_gk, mu_k
    res = 0.5 * (1.5625 * mu_k + 1.875 * mu_gk + 0.5625 * mu_g)
  end function breaking_variance


  pure elemental real function input_variance(mu_g, mu_gk, mu_k) result(res)
    ! The squared steepness of the equivalent long wave for the modulation
    ! of k sigma, which scales the wind input, eps cos(psi) (3/2 cos^2 dth - 1/2)
    ! at first order, with the same variance as the spectrum.
    !$acc routine seq
    real, intent(in) :: mu_g, mu_gk, mu_k
    res = max(4.5 * mu_k - 3 * mu_gk + 0.5 * mu_g, 0.)
  end function input_variance


  pure elemental real function speed_steepness(mu_g, mu_gk, mu_k) result(res)
    ! The steepness of the equivalent long wave of the speed modulation.
    !$acc routine seq
    real, intent(in) :: mu_g, mu_gk, mu_k
    res = sqrt(speed_variance(mu_g, mu_gk, mu_k))
  end function speed_steepness


  pure elemental real function breaking_steepness(mu_g, mu_gk, mu_k) result(res)
    ! The steepness of the equivalent long wave of the steepness modulation.
    !$acc routine seq
    real, intent(in) :: mu_g, mu_gk, mu_k
    res = sqrt(breaking_variance(mu_g, mu_gk, mu_k))
  end function breaking_steepness


  pure elemental real function input_steepness(mu_g, mu_gk, mu_k) result(res)
    ! The steepness of the equivalent long wave of the modulation of k sigma.
    !$acc routine seq
    real, intent(in) :: mu_g, mu_gk, mu_k
    res = sqrt(input_variance(mu_g, mu_gk, mu_k))
  end function input_steepness


  pure elemental real function relative_speed_amplitude(mu_g, mu_gk, mu_k, m6, m7, m8, m9, m10, &
                                                        m11, m12, cth, sth, cth2, sth2, cth3, &
                                                        sth3, cp) result(res)
    ! The amplitude b of the modulation -b cos(psi) of the wind-wave relative
    ! speed of a short wave of direction th and phase speed cp by the longer
    ! waves: their orbital velocity along the short wave, eps C_L cos(psi)
    ! cos(dth), less the modulation of its phase speed, from the slope
    ! variances and the orbital and Stokes moments 6-12 of longwave_moments.
    ! For a single aligned long wave, b = eps (C_L - cp).
    !$acc routine seq
    real, intent(in) :: mu_g, mu_gk, mu_k, m6, m7, m8, m9, m10, m11, m12
    real, intent(in) :: cth, sth, cth2, sth2, cth3, sth3, cp
    real :: var_u, var_c, cov, sk1, sk3
    var_u = 0.5 * (m6 + cth2 * m7 + sth2 * m8)       ! orbital velocity along the short wave
    var_c = 0.25 * cp * cp * (mu_g + 2 * mu_gk + mu_k) ! phase speed
    sk1 = cth * m9 + sth * m10                        ! Stokes moments projected on the short wave
    sk3 = cth3 * m11 + sth3 * m12
    cov = -0.5 * cp * (1.75 * sk1 + 0.25 * sk3)       ! 2 x their covariance, cos dth + cos^3 dth
    res = sqrt(max(2 * (var_u + var_c + cov), 0.))
  end function relative_speed_amplitude


  pure elemental real function sin_mod_bin(wind, cosd, cp, uc, vc, cth, sth, shelt, sin_fac, &
                                           fieldscale1, fieldscale2, scale, fkovg, fks, lower, b) &
    result(res)
    ! Wind input of one spectral bin (see sin_d12_bin of umwm_source_functions)
    ! averaged over the phase of the longer waves, which modulate the
    ! wind-wave relative speed by -b cos(psi): the mean of dU |dU| over the
    ! phase in closed form, with the coefficient of following winds where
    ! the relative speed is positive and that of opposing winds where it is
    ! negative, and the input scaled by fks, the phase mean of the modulated
    ! k sigma over its unmodulated value; not less than lower.
    !$acc routine seq
    real, intent(in) :: wind, cosd, cp, uc, vc, cth, sth, shelt
    real, intent(in) :: sin_fac, fieldscale1, fieldscale2, scale, fkovg, fks, lower, b
    real :: a, s, coef_neg
    a = wind * cosd - cp - uc * cth - vc * sth
    coef_neg = sin_fac * fieldscale1 * merge(fieldscale2, 1., cosd > 0)
    if (b <= abs(a)) then
      ! the relative speed keeps the sign of its mean over the phase
      s = a * a + 0.5 * b * b
      s = merge(shelt * s, -coef_neg * s, a > 0)
    else
      s = sign_change_mean(a, b, shelt, coef_neg)
    end if
    res = max(scale * s * fkovg * fks, lower)
  end function sin_mod_bin


  pure subroutine sin_mod(nf, nd, no, oc, wind, cosd, cp, uc, vc, cth, sth, cth2, sth2, cth3, sth3, &
                          cth4, sth4, shelt, sin_fac, fieldscale1, fieldscale2, scale, fkovg, &
                          moments, coef_c, coef_ks, xmax, ssin)
    ! Wind input averaged over the phase of the longer waves, evaluated for
    ! frequencies 1:no, with the phase speed of each bin scaled by the phase
    ! mean of its modulation: sin_mod_bin over the column, with the
    ! modulation of each bin from the moments as slope_variances,
    ! speed_variance, input_variance and relative_speed_amplitude, written
    ! out for speed. The bins where the relative speed keeps its sign over
    ! the phase, nearly all, are evaluated in a vectorizable pass, and the
    ! others in a second pass with the inverse cosine of sign_change_mean.
    integer, intent(in) :: nf, nd, no, oc
    real, intent(in) :: wind(nf), cosd(nd), cp(nf), uc, vc
    real, intent(in) :: cth(nd), sth(nd), cth2(nd), sth2(nd), cth3(nd), sth3(nd), cth4(nd), sth4(nd)
    real, intent(in) :: shelt, sin_fac, fieldscale1, fieldscale2, scale, fkovg(nf)
    real, intent(in) :: moments(nf,num_moments), coef_c(num_coef), coef_ks(num_coef), xmax
    real, intent(inout) :: ssin(nf,nd)
    real :: a(nf,nd), b(nf,nd), w(nf,nd) ! relative speed, its modulation amplitude, and the factor of dU |dU|
    real :: lower(nf), wscale(nf), coef_neg(nd)
    real :: inv_xmax, p2, p4, mu_g, mu_gk, mu_k, tc, tks, fc, fks, sk1, sk3, var, s
    integer :: o, p

    inv_xmax = 1 / xmax
    do o = 1, no
      lower(o) = merge(0., -huge(0.), o > oc)
      wscale(o) = scale * fkovg(o)
    end do
    do p = 1, nd
      coef_neg(p) = sin_fac * fieldscale1 * merge(fieldscale2, 1., cosd(p) > 0)
    end do

    do p = 1, nd
      do o = 1, no
        ! the slope variances and the squared equivalent steepnesses
        p2 = cth2(p) * moments(o,2) + sth2(p) * moments(o,3)
        p4 = cth4(p) * moments(o,4) + sth4(p) * moments(o,5)
        mu_g = moments(o,1)
        mu_gk = 0.5 * (mu_g + p2)
        mu_k = 0.125 * (3 * mu_g + 4 * p2 + p4)
        tc = min(0.5 * (mu_g + 2 * mu_gk + mu_k), xmax) * inv_xmax
        tks = min(max(4.5 * mu_k - 3 * mu_gk + 0.5 * mu_g, 0.), xmax) * inv_xmax
        ! the phase means of the phase speed and of k sigma
        fc = 1 + tc * (coef_c(1) + tc * (coef_c(2) + tc * (coef_c(3) + coef_c(4) * tc)))
        fks = 1 + tks * (coef_ks(1) + tks * (coef_ks(2) + tks * (coef_ks(3) + coef_ks(4) * tks)))
        ! the amplitude of the modulation of the relative speed
        sk1 = cth(p) * moments(o,9) + sth(p) * moments(o,10)
        sk3 = cth3(p) * moments(o,11) + sth3(p) * moments(o,12)
        var = 0.5 * (moments(o,6) + cth2(p) * moments(o,7) + sth2(p) * moments(o,8)) &
            + 0.25 * cp(o) * cp(o) * (mu_g + 2 * mu_gk + mu_k) &
            - 0.5 * cp(o) * (1.75 * sk1 + 0.25 * sk3)
        b(o,p) = sqrt(max(2 * var, 0.))
        a(o,p) = wind(o) * cosd(p) - cp(o) * fc - uc * cth(p) - vc * sth(p)
        w(o,p) = wscale(o) * fks
        s = a(o,p) * a(o,p) + 0.5 * b(o,p) * b(o,p)
        s = merge(shelt * s, -coef_neg(p) * s, a(o,p) > 0)
        ssin(o,p) = max(w(o,p) * s, lower(o))
      end do
    end do

    do p = 1, nd
      do o = 1, no
        if (b(o,p) > abs(a(o,p))) &
          ssin(o,p) = max(w(o,p) * sign_change_mean(a(o,p), b(o,p), shelt, coef_neg(p)), lower(o))
      end do
    end do

  end subroutine sin_mod


  pure elemental real function sign_change_mean(a, b, shelt, coef_neg) result(res)
    ! The mean over the phase of dU |dU| with its coefficients, shelt where
    ! dU = a - b cos(psi) is positive and -coef_neg where it is negative,
    ! for b > |a|: negative for |psi| < psi0 = acos(a/b).
    !$acc routine seq
    real, intent(in) :: a, b, shelt, coef_neg
    real :: psi0, ineg
    psi0 = acos(a / b)
    ineg = (a * a * psi0 - 2 * a * b * sin(psi0) + b * b * (0.5 * psi0 + 0.25 * sin(2 * psi0))) / pi
    res = shelt * (a * a + 0.5 * b * b - ineg) - coef_neg * ineg
  end function sign_change_mean


  pure subroutine breaking_modulation(nf, nd, o1, o2, moments, cth2, sth2, cth4, sth4, &
                                      coef_ds, xmax, fds)
    ! The enhancement of the breaking dissipation of frequencies o1:o2 by
    ! the modulation of the short-wave steepness, the phase mean of the
    ! local saturation of each bin raised to the saturation exponent, at
    ! the squared equivalent steepness of breaking_variance, written out
    ! for speed.
    integer, intent(in) :: nf, nd, o1, o2
    real, intent(in) :: moments(nf,num_moments), cth2(nd), sth2(nd), cth4(nd), sth4(nd)
    real, intent(in) :: coef_ds(num_coef), xmax
    real, intent(inout) :: fds(nf,nd)
    real :: inv_xmax, p2, p4, mu_g, t
    integer :: o, p
    inv_xmax = 1 / xmax
    do p = 1, nd
      do o = o1, o2
        p2 = cth2(p) * moments(o,2) + sth2(p) * moments(o,3)
        p4 = cth4(p) * moments(o,4) + sth4(p) * moments(o,5)
        mu_g = moments(o,1)
        ! 0.5 (1.5625 mu_k + 1.875 mu_gk + 0.5625 mu_g) in the slope variances
        t = min(1.04296875 * mu_g + 0.859375 * p2 + 0.09765625 * p4, xmax) * inv_xmax
        fds(o,p) = 1 + t * (coef_ds(1) + t * (coef_ds(2) + t * (coef_ds(3) + coef_ds(4) * t)))
      end do
    end do
  end subroutine breaking_modulation


  pure subroutine modulated_group_speed(nf, esum, moments, coef_c, xmax, cg0, cgm)
    ! The group speed of each frequency scaled by the phase mean of its
    ! modulation by the longer waves at the squared equivalent steepness of
    ! speed_variance averaged over directions, weighted by the spectrum
    ! through its direction sums esum; the unmodulated group speed where
    ! the spectrum is empty.
    integer, intent(in) :: nf
    real, intent(in) :: esum(nf,num_sums), moments(nf,num_moments), coef_c(num_coef), xmax, cg0(nf)
    real, intent(out) :: cgm(nf)
    real :: p2, p4, mu_g, mu_gk, mu_k, t
    integer :: o
    do o = 1, nf
      if (esum(o,1) > 0) then
        p2 = (esum(o,2) * moments(o,2) + esum(o,3) * moments(o,3)) / esum(o,1)
        p4 = (esum(o,4) * moments(o,4) + esum(o,5) * moments(o,5)) / esum(o,1)
        mu_g = moments(o,1)
        mu_gk = 0.5 * (mu_g + p2)
        mu_k = 0.125 * (3 * mu_g + 4 * p2 + p4)
        t = min(0.5 * (mu_g + 2 * mu_gk + mu_k), xmax) / xmax
        cgm(o) = cg0(o) * (1 + t * (coef_c(1) + t * (coef_c(2) + t * (coef_c(3) + coef_c(4) * t))))
      else
        cgm(o) = cg0(o)
      end if
    end do
  end subroutine modulated_group_speed


  pure real function modulation_steepness(nf, nd, o, e, moments, cth2, sth2, cth4, sth4) result(res)
    ! The steepness of the equivalent long wave of the breaking modulation
    ! of frequency o, the mean over directions weighted by the spectrum.
    integer, intent(in) :: nf, nd, o
    real, intent(in) :: e(nf,nd), moments(nf,num_moments), cth2(nd), sth2(nd), cth4(nd), sth4(nd)
    real :: mu_g, mu_gk, mu_k, num, den
    integer :: p
    num = 0
    den = 0
    do p = 1, nd
      call slope_variances(moments(o,1), moments(o,2), moments(o,3), moments(o,4), &
                           moments(o,5), cth2(p), sth2(p), cth4(p), sth4(p), mu_g, mu_gk, mu_k)
      num = num + e(o,p) * breaking_steepness(mu_g, mu_gk, mu_k)
      den = den + e(o,p)
    end do
    res = 0
    if (den > 0) res = num / den
  end function modulation_steepness

end module umwm_modulation
