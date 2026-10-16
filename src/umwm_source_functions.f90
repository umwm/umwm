module umwm_source_functions
  ! Wave source functions of Donelan et al. (2012), evaluated for one
  ! spectral column (all frequencies and directions of one grid cell).
  ! The prognostic range of frequencies is 1:oc; higher bins are diagnostic.

  use, intrinsic :: iso_fortran_env, only: real64

  implicit none

  private
  public :: mss_moments, power, sin_d12, sds_coefficient, sds_d12, snl_d12, s_ice

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
    real :: s, overrun
    integer :: o, p

    do o = 1, no
      lower(o) = merge(0., -huge(0.), o > oc)
    end do

    do p = 1, nd
      overrun = merge(fieldscale2, 1., cosd(p) > 0)
      do o = 1, no
        s = wind(o) * cosd(p) - cp(o) - uc * cth(p) - vc * sth(p)
        s = sin_fac * abs(s) * s
        ! variable sheltering coefficient for following winds; damping of
        ! opposing winds, further reduced for swell overrunning the wind
        s = merge(s * shelt / sin_fac, s * fieldscale1 * overrun, s > 0)
        ssin(o,p) = max(scale * s * fkovg(o), lower(o))
      end do
    end do

  end subroutine sin_d12


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
        sds_coef(o,p) = twopisds_fac * f(o) * (1 + fac &
          * (moments(o,1) + cth2(p) * moments(o,2) + sth2(p) * moments(o,3)))**2
      end do
    end do
  end subroutine sds_coefficient


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
    xmin = max(tiny(x)**(1 / a), tiny(x))
    x = merge(real(exp(a * log(real(max(x, xmin), real64)))), 0., x > xmin)
  end subroutine power


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
        snl(o,p) = bf1_renorm(o) * sds(o+1,p) * e(o+1,p) &
                 + bf2_renorm(o) * sds(o+2,p) * e(o+2,p) &
                 - snl_fac * sds(o,p) * e(o,p)
      end do
      snl(oc+1:nf,p) = 0
    end do

  end subroutine snl_d12


  pure subroutine s_ice(nf, nd, e, kdk, cg, dth, fice, fice_lth, sice)
    ! Wave attenuation by sea ice, following Kohout et al. (2014).
    integer, intent(in) :: nf, nd
    real, intent(in) :: e(nf,nd), kdk(nf), cg(nf)
    real, intent(in) :: dth, fice, fice_lth
    real, intent(out) :: sice(nf)

    real, parameter :: H_th = 3.0       ! [m]
    real, parameter :: C1   = -5.35e-6  ! [m-1]
    real, parameter :: C2   = C1 * H_th ! []

    real :: swh
    integer :: p

    sice = 0
    if (fice <= fice_lth) return

    swh = 0
    do p = 1, nd
      swh = swh + sum(e(:,p) * kdk)
    end do
    swh = 4 * sqrt(swh * dth)

    ! wave attenuation from sea ice in the two SWH regimes
    if (swh < H_th) then
      sice = 2 * cg * (C1 * swh)
    else
      sice = 2 * cg * C2
    end if

  end subroutine s_ice

end module umwm_source_functions
