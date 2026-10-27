module umwm_source_functions
  ! Wave source functions of Donelan et al. (2012), evaluated for one
  ! spectral column (all frequencies and directions of one grid cell).
  ! The prognostic range of frequencies is 1:oc; higher bins are diagnostic.
  ! The elemental functions with the suffix _bin evaluate one spectral bin,
  ! for the column routines and for the GPU kernels.

  use, intrinsic :: iso_fortran_env, only: real64

  implicit none

  private
  public :: mss_moments, power, power_bin, power_threshold, sin_d12, sin_d12_bin, &
            sds_coefficient, sds_coefficient_bin, sds_d12, snl_d12, snl_d12_bin, &
            s_ice, s_ice_bin

  ! sea ice attenuation, Kohout et al. (2014)
  real, parameter :: H_th = 3.0       ! [m]
  real, parameter :: C1   = -5.35e-6  ! [m-1]
  real, parameter :: C2   = C1 * H_th ! []

contains

  pure subroutine sin_d12(nf, nd, no, oc, wind, cosd, cp, uc, vc, cth, sth, &
                          shelt, sin_fac, fieldscale1, fieldscale2, scale, fkovg, ssin)
    ! Wind input function based on Jeffreys's sheltering hypothesis,
    ! evaluated for frequencies 1:no.
    integer, intent(in) :: nf, nd, no, oc
    real, intent(in) :: wind(nf)          ! wind speed at half-wavelength height [m/s]
    real, intent(in) :: cosd(nd)          ! cosine of the wind-wave angle
    real, intent(in) :: cp(nf)            ! intrinsic phase speed [m/s]
    real, intent(in) :: uc, vc            ! ocean current [m/s]
    real, intent(in) :: cth(nd), sth(nd)  ! wave direction cosines
    real, intent(in) :: shelt             ! sheltering coefficient
    real, intent(in) :: sin_fac, fieldscale1, fieldscale2
    real, intent(in) :: scale             ! (1 - fice) * 2 pi * rhoa / rhow
    real, intent(in) :: fkovg(nf)         ! f * k / g
    real, intent(inout) :: ssin(nf,nd)
    real :: lower(no) ! no negative input for the diagnostic tail
    integer :: o, p

    do o = 1, no
      lower(o) = merge(0., -huge(0.), o > oc)
    end do

    do p = 1, nd
      do o = 1, no
        ssin(o,p) = sin_d12_bin(wind(o), cosd(p), cp(o), uc, vc, cth(p), sth(p), shelt, &
                                sin_fac, fieldscale1, fieldscale2, scale, fkovg(o), lower(o))
      end do
    end do

  end subroutine sin_d12


  pure elemental real function sin_d12_bin(wind, cosd, cp, uc, vc, cth, sth, shelt, sin_fac, &
                                           fieldscale1, fieldscale2, scale, fkovg, lower) &
    result(res)
    ! Wind input of one spectral bin (see sin_d12), not less than lower.
    !$acc routine seq
    real, intent(in) :: wind, cosd, cp, uc, vc, cth, sth, shelt
    real, intent(in) :: sin_fac, fieldscale1, fieldscale2, scale, fkovg, lower
    real :: s
    s = wind * cosd - cp - uc * cth - vc * sth
    s = sin_fac * abs(s) * s
    ! variable sheltering coefficient for following winds; damping of
    ! opposing winds, further reduced for swell overrunning the wind
    s = merge(s * shelt / sin_fac, s * fieldscale1 * merge(fieldscale2, 1., cosd > 0), s > 0)
    res = max(scale * s * fkovg, lower)
  end function sin_d12_bin


  pure subroutine mss_moments(nf, nd, e, k3dk, cth2, sth2, dth, moments)
    ! Direction-integrated moments of the mean-squared slope of all longer
    ! waves, such that the slope projected onto direction th(p) is
    ! sum(e(o',p') * k3dk(o') * cos(th(p') - th(p))**2 * dth, o' < o)
    !   = moments(o,1) + cos(2 th(p)) * moments(o,2) + sin(2 th(p)) * moments(o,3),
    ! because cos(x)**2 = (1 + cos(2x)) / 2.
    integer, intent(in) :: nf, nd
    real, intent(in) :: e(nf,nd), k3dk(nf)
    real, intent(in) :: cth2(nd), sth2(nd) ! cos(2 th), sin(2 th)
    real, intent(in) :: dth
    real, intent(out) :: moments(nf,3)
    real :: d0(nf), dc(nf), ds(nf)
    integer :: o, p

    d0 = 0
    dc = 0
    ds = 0
    do p = 1, nd
      do o = 1, nf
        d0(o) = d0(o) + e(o,p)
        dc(o) = dc(o) + e(o,p) * cth2(p)
        ds(o) = ds(o) + e(o,p) * sth2(p)
      end do
    end do

    moments(1,:) = 0
    do o = 2, nf
      moments(o,1) = moments(o-1,1) + 0.5 * dth * d0(o-1) * k3dk(o-1)
      moments(o,2) = moments(o-1,2) + 0.5 * dth * dc(o-1) * k3dk(o-1)
      moments(o,3) = moments(o-1,3) + 0.5 * dth * ds(o-1) * k3dk(o-1)
    end do

  end subroutine mss_moments


  pure subroutine sds_coefficient(nf, nd, o1, o2, f, cth2, sth2, mss_fac, twopisds_fac, &
                                  moments, sds_coef)
    ! Wave breaking dissipation coefficient for frequencies o1:o2,
    ! enhanced by the mean-squared slope of longer waves if mss_fac > 0.
    integer, intent(in) :: nf, nd, o1, o2
    real, intent(in) :: f(nf), cth2(nd), sth2(nd), mss_fac, twopisds_fac, moments(nf,3)
    real, intent(inout) :: sds_coef(nf,nd)
    real :: fac
    integer :: o, p
    fac = max(mss_fac, 0.)
    do p = 1, nd
      do o = o1, o2
        sds_coef(o,p) = sds_coefficient_bin(f(o), cth2(p), sth2(p), fac, twopisds_fac, &
                                            moments(o,1), moments(o,2), moments(o,3))
      end do
    end do
  end subroutine sds_coefficient


  pure elemental real function sds_coefficient_bin(f, cth2, sth2, fac, twopisds_fac, &
                                                   moment1, moment2, moment3) result(res)
    ! Wave breaking dissipation coefficient of one spectral bin (see
    ! sds_coefficient), with fac = max(mss_fac, 0).
    !$acc routine seq
    real, intent(in) :: f, cth2, sth2, fac, twopisds_fac, moment1, moment2, moment3
    res = twopisds_fac * f * (1 + fac * (moment1 + cth2 * moment2 + sth2 * moment3))**2
  end function sds_coefficient_bin


  pure subroutine sds_d12(nf, nd, o1, o2, e, k4, f, cth2, sth2, mss_fac, twopisds_fac, &
                          moments, sds_power, sds)
    ! Wave breaking dissipation for frequencies o1:o2, proportional to a
    ! power of the saturation spectrum. The adjustment for plunging breakers,
    ! coth(0.2 kd), is not included and is applied where sds is used.
    integer, intent(in) :: nf, nd, o1, o2
    real, intent(in) :: e(nf,nd), k4(nf), f(nf), cth2(nd), sth2(nd)
    real, intent(in) :: mss_fac, twopisds_fac, moments(nf,3), sds_power
    real, intent(inout) :: sds(nf,nd)
    real :: saturation(o1:o2,nd), sds_coef(nf,nd)
    integer :: p

    if (o2 < o1) return
    do p = 1, nd
      saturation(:,p) = e(o1:o2,p) * k4(o1:o2)
    end do
    call power(size(saturation), saturation, sds_power)
    call sds_coefficient(nf, nd, o1, o2, f, cth2, sth2, mss_fac, twopisds_fac, moments, sds_coef)
    do p = 1, nd
      sds(o1:o2,p) = sds_coef(o1:o2,p) * saturation(:,p)
    end do

  end subroutine sds_d12


  pure subroutine power(n, x, a)
    ! Raises x to the positive power a in place, in double precision for
    ! accuracy and speed. Inputs and results that would be subnormal in
    ! single precision give zero.
    integer, intent(in) :: n
    real, intent(inout) :: x(n)
    real, intent(in) :: a
    real :: xmin
    xmin = power_threshold(a)
    x = merge(real(exp(a * log(real(max(x, xmin), real64)))), 0., x > xmin)
  end subroutine power


  pure real function power_threshold(a) result(xmin)
    ! The input of power below which the input or the result would be
    ! subnormal in single precision.
    real, intent(in) :: a
    xmin = max(tiny(xmin)**(1 / a), tiny(xmin))
  end function power_threshold


  pure elemental real function power_bin(x, a, xmin) result(res)
    ! Raises x to the positive power a in single precision, zero for x not
    ! above xmin = power_threshold(a); for the GPU, where double precision
    ! is slow.
    !$acc routine seq
    real, intent(in) :: x, a, xmin
    res = merge(max(x, xmin)**a, 0., x > xmin)
  end function power_bin


  pure subroutine snl_d12(nf, nd, oc, e, bf1_renorm, bf2_renorm, snl_fac, sds, snl)
    ! Nonlinear downshifting: spread wave energy to the 2 next longer
    ! wavenumbers exponentially decaying as distance from the donating
    ! wavenumber, and remove the energy from the donating wavenumbers.
    integer, intent(in) :: nf, nd, oc
    real, intent(in) :: e(nf,nd), sds(nf,nd)
    real, intent(in) :: bf1_renorm(nf), bf2_renorm(nf)
    real, intent(in) :: snl_fac
    real, intent(out) :: snl(nf,nd)
    integer :: o, p

    do p = 1, nd
      do o = 1, oc
        snl(o,p) = snl_d12_bin(bf1_renorm(o), bf2_renorm(o), snl_fac, sds(o,p), e(o,p), &
                               sds(o+1,p), e(o+1,p), sds(o+2,p), e(o+2,p))
      end do
      snl(oc+1:nf,p) = 0
    end do

  end subroutine snl_d12


  pure elemental real function snl_d12_bin(bf1_renorm, bf2_renorm, snl_fac, sds, e, &
                                           sds1, e1, sds2, e2) result(res)
    ! Nonlinear downshifting of one spectral bin (see snl_d12), from the
    ! dissipation and spectrum of the bin and of the next two frequencies.
    !$acc routine seq
    real, intent(in) :: bf1_renorm, bf2_renorm, snl_fac, sds, e, sds1, e1, sds2, e2
    res = bf1_renorm * sds1 * e1 + bf2_renorm * sds2 * e2 - snl_fac * sds * e
  end function snl_d12_bin


  pure subroutine s_ice(nf, nd, e, kdk, cg, dth, fice, fice_lth, sice)
    ! Wave attenuation by sea ice, following Kohout et al. (2014).
    integer, intent(in) :: nf, nd
    real, intent(in) :: e(nf,nd), kdk(nf), cg(nf)
    real, intent(in) :: dth, fice, fice_lth
    real, intent(out) :: sice(nf)

    real :: swh
    integer :: p

    sice = 0
    if (fice <= fice_lth) return

    swh = 0
    do p = 1, nd
      swh = swh + sum(e(:,p) * kdk)
    end do
    swh = 4 * sqrt(swh * dth)

    sice = s_ice_bin(cg, swh)

  end subroutine s_ice


  pure elemental real function s_ice_bin(cg, swh) result(res)
    ! Wave attenuation by sea ice of one frequency with group speed cg,
    ! given the significant wave height swh, in its two regimes.
    !$acc routine seq
    real, intent(in) :: cg, swh
    if (swh < H_th) then
      res = 2 * cg * (C1 * swh)
    else
      res = 2 * cg * C2
    end if
  end function s_ice_bin

end module umwm_source_functions
