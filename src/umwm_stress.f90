module umwm_stress

  !! Functions and subroutines to evaluate stresses.

  use umwm_constants, only: rk

  implicit none

  private
  public :: directional_integral, drag_coefficient, drag_coefficient_skin, &
            stress_tail, wind_relative

contains

  real(rk) pure elemental function drag_coefficient(taux, tauy, rhoa, wspd) result(cd)
    !! Computes drag coefficient from stress vector (taux, tauy, N/m^2),
    !! air density (rhoa, kg/m^3), and wind speed (wspd, m/s).
    real(rk), intent(in) :: taux, tauy, rhoa, wspd
    cd = sqrt(taux**2 + tauy**2) / (rhoa * wspd**2)
  end function drag_coefficient

  real(rk) pure elemental function drag_coefficient_skin(cd_form, wspd, wspdrel, z, &
                                                         air_viscosity, von_karman) result(cd)
    !! Computes the skin drag coefficient attenuated by from drag,
    !! given input absolute (wspd) and relative (wspdrel) wind speed (m/s),
    !! height z (m), air viscosity (m^2/s), and Von Karman constant.
    real(rk), intent(in) :: cd_form, wspd, wspdrel, z, air_viscosity, von_karman
    cd = friction_velocity_skin(wspdrel, z, air_viscosity, von_karman)**2 / wspd**2
    cd = cd * (1 + 2 * cd / (cd + cd_form + tiny(cd))) / 3
    if (cd > 1e-2) cd = 1e-2
  end function drag_coefficient_skin

  real(rk) pure elemental function friction_velocity_skin(wspd, z, air_viscosity, von_karman) result(ustar)
    !! Computes the skin friction velocity (flow over smooth surface) in m/s
    !! given input wind speed (m/s), height z (m), air viscosity (m^2/s),
    !! and Von Karman constant. Input wind speed should relative to the
    !! water surface.
    real(rk), intent(in) :: wspd, z, air_viscosity, von_karman
    real(rk) :: z0
    integer :: n
    z0 = 1e-3_rk
    do n = 1, 6
      ustar = von_karman * wspd / log(z / z0)
      z0 = 0.132 * air_viscosity / ustar
    end do
  end function friction_velocity_skin

  pure elemental real(rk) function stress_tail(wspd, kmax) result(tail)
    !! Computes the wind speed-dependent stress tail
    !! from kmax (highest wavenumber bin) to capillary.
    real(rk), intent(in) :: wspd, kmax
    real(rk) :: kmax_pow_tail
    real(rk), parameter :: kcap = 1e3 ! capillary wavenumber
    real(rk), parameter :: a = 0.000112, b = -0.01451, c = -1.0186
    tail = a * wspd**2 + b * wspd + c
    kmax_pow_tail = kmax**tail
    tail = (kcap**(tail + 1) - kmax * kmax_pow_tail)&
         / (kmax_pow_tail * (tail + 1))
  end function stress_tail

  pure subroutine directional_integral(nf, nd, a, b, cth, sth, x, y)
    !! Integrates a(o,p) * b(o,p) over directions, weighted by the x- and
    !! y-components of the propagation direction, for each frequency.
    integer, intent(in) :: nf, nd
    real(rk), intent(in) :: a(nf,nd), b(nf,nd), cth(nd), sth(nd)
    real(rk), intent(out) :: x(nf), y(nf)
    integer :: p
    x = 0
    y = 0
    do p = 1, nd
      x = x + a(:,p) * b(:,p) * cth(p)
      y = y + a(:,p) * b(:,p) * sth(p)
    end do
  end subroutine directional_integral

  pure elemental subroutine wind_relative(wspd, wdir, us, vs, wspdrel, wdirrel)
    !! Computes the wind speed and direction relative to the water surface,
    !! given input absolute wind speed and direction, and u- and v- components
    !! of the surface velocity.
    real(rk), intent(in) :: wspd, wdir, us, vs
    real(rk), intent(out) :: wspdrel, wdirrel
    real(rk) :: urel, vrel
    urel = wspd * cos(wdir) - us
    vrel = wspd * sin(wdir) - vs
    wspdrel = sqrt(urel**2 + vrel**2)
    wdirrel = atan2(vrel, urel)
  end subroutine wind_relative

end module umwm_stress
