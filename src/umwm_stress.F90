module umwm_stress

  !! Module with functions and subroutines to evaluate stresses.

  use umwm_module
  use umwm_advection, only: zerocurrents
  use umwm_constants, only: rk
  use umwm_stokes, only: util

  implicit none

  private
  public :: stress

contains

  subroutine stress(option)

    character(3), intent(in) :: option

    integer :: i, o, p

    real(rk) :: taux_util(om, istart:iend), tauy_util(om, istart:iend)
    real(rk) :: tail(istart:iend)

    ! x- and y-components of surface stokes drift velocities
    real(rk) :: usurf(istart:iend), vsurf(istart:iend)

    ! wind speed and direction relative to surface velocity
    real(rk) :: wspdrel(istart:iend), wdirrel(istart:iend)

    real(rk) :: cd_form(istart:iend), cd_skin(istart:iend)
    real(rk) :: spectrum_integral

    ! evaluate wind speed dependent tail
    tail = stress_tail(wspd(istart:iend), k(om,istart:iend))

    if (option == 'ocn') then

      ! compute stress into ocean top (positive into ocean)
      taux_util = 0
      tauy_util = 0
      do i = istart, iend
        do p = 1, pm

          ! dissipation into currents
          do o = 1, om
            taux_util(o,i) = taux_util(o,i) + e(o,p,i)         &
                           * (sds(o,p,i) + sdt(o,i) + sdv(o,i))&
                           * cth(p) * invcp0(o,i)
            tauy_util(o,i) = tauy_util(o,i) + e(o,p,i)         &
                           * (sds(o,p,i) + sdt(o,i) + sdv(o,i))&
                           * sth(p) * invcp0(o,i)
          end do

          ! correction for the Snl in diagnostic part
          do o = 3, om

            taux_util(o,i) = taux_util(o,i) - snl(o,p,i)     &
                           * (bf1 * cp0(o,i) * invcp0(o-1,i) &
                           +  bf2 * cp0(o,i) * invcp0(o-2,i))&
                           * cth(p) * invcp0(o,i)

            tauy_util(o,i) = tauy_util(o,i) - snl(o,p,i)     &
                           * (bf1 * cp0(o,i) * invcp0(o-1,i) &
                           +  bf2 * cp0(o,i) * invcp0(o-2,i))&
                           * sth(p) * invcp0(o,i)

          end do
        end do
      end do

      do concurrent(i = istart:iend)

        ! compute the tail
        tailocnx(i) = taux_util(om,i) * k(om,i) * tail(i) * rhow(i) * dthg
        tailocny(i) = tauy_util(om,i) * k(om,i) * tail(i) * rhow(i) * dthg

        ! integrate over frequencies
        taux_ocntop(i) = sum(taux_util(:,i) * kdk(:,i), dim=1) * rhow(i) * dthg&
                       + tailocnx(i) + taux_skin(i)
        tauy_ocntop(i) = sum(tauy_util(:,i) * kdk(:,i), dim=1) * rhow(i) * dthg&
                       + tailocny(i) + tauy_skin(i)

      end do

      ! compute stress into ocean bottom (positive into ocean)
      taux_util = 0
      tauy_util = 0
      do i = istart, iend
        do p = 1, pm
          do o = 1, om
            taux_util(o,i) = taux_util(o,i) + e(o,p,i) * sbf(o,i) * cth(p) * invcp0(o,i)
            tauy_util(o,i) = tauy_util(o,i) + e(o,p,i) * sbf(o,i) * sth(p) * invcp0(o,i)
          end do
        end do
      end do

      do concurrent(i = istart:iend)
        taux_ocnbot(i) = sum(taux_util(:,i) * kdk(:,i), dim=1) * rhow(i) * dthg
        tauy_ocnbot(i) = sum(tauy_util(:,i) * kdk(:,i), dim=1) * rhow(i) * dthg
      end do

      ! Snl conserves energy, but not momentum. the following term is
      ! the momentum loss due to non-linear downshifting of energy.
      taux_util = 0
      tauy_util = 0
      do i = istart, iend
        do p = 1, pm
          do o = 3, oc(i)
            taux_util(o,i) = taux_util(o,i) + snl(o,p,i)           &
                           * (bf1 * (1 - cp0(o,i) * invcp0(o-1,i)) &
                           +  bf2 * (1 - cp0(o,i) * invcp0(o-2,i)))&
                           * cth(p) * invcp0(o,i)
            tauy_util(o,i) = tauy_util(o,i) + snl(o,p,i)           &
                           * (bf1 * (1 - cp0(o,i) * invcp0(o-1,i)) &
                           +  bf2 * (1 - cp0(o,i) * invcp0(o-2,i)))&
                           * sth(p) * invcp0(o,i)
          end do
        end do
      end do

      do concurrent(i = istart:iend)
        taux_snl(i) = sum(taux_util(:,i) * kdk(:,i), dim=1) * rhow(i) * dthg
        tauy_snl(i) = sum(tauy_util(:,i) * kdk(:,i), dim=1) * rhow(i) * dthg
      end do

      ! compute energy flux into ocean:
      taux_util = 0
      tauy_util = 0
      do i = istart, iend
        do p = 1, pm

          do o = 1, om
            taux_util(o,i) = taux_util(o,i) + e(o,p,i) * sds(o,p,i) * cth(p)
            tauy_util(o,i) = tauy_util(o,i) + e(o,p,i) * sds(o,p,i) * sth(p)
          end do

          do o = 3, om
            taux_util(o,i) = taux_util(o,i) - snl(o,p,i)     &
                           * (bf1 * cp0(o,i) * invcp0(o-1,i) &
                           +  bf2 * cp0(o,i) * invcp0(o-2,i))&
                           * cth(p)
            tauy_util(o,i) = tauy_util(o,i) - snl(o,p,i)     &
                           * (bf1 * cp0(o,i) * invcp0(o-1,i) &
                           +  bf2 * cp0(o,i) * invcp0(o-2,i))&
                           * sth(p)
          end do

        end do
      end do

      do concurrent(i = istart:iend)
        epsx_ocn(i) = sum(taux_util(:,i) * kdk(:,i), dim=1) * rhow(i) * dthg
        epsy_ocn(i) = sum(tauy_util(:,i) * kdk(:,i), dim=1) * rhow(i) * dthg
      end do

      ! compute energy flux from air:
      taux_util = 0
      tauy_util = 0
      do i = istart, iend
        do p = 1, pm
          do o = 1, om
            taux_util(o,i) = taux_util(o,i) + e(o,p,i) * ssin(o,p,i) * cth(p)
            tauy_util(o,i) = tauy_util(o,i) + e(o,p,i) * ssin(o,p,i) * sth(p)
          end do
        end do
      end do

      do concurrent(i = istart:iend)
        epsx_atm(i) = sum(taux_util(:,i) * kdk(:,i), dim=1) * rhow(i) * dthg
        epsy_atm(i) = sum(tauy_util(:,i) * kdk(:,i), dim=1) * rhow(i) * dthg
      end do

      ! this part calculates the components of form drag.
      ! it has nothing to do with momentum fluxes into ocean,
      ! but we do it here because it is needed only for output.
      taux1 = 0
      tauy1 = 0
      taux2 = 0
      tauy2 = 0
      taux3 = 0
      tauy3 = 0

      do i = istart, iend
        do p = 1, pm
          do o = 1, om

            dummy(o,p,i) = e(o,p,i) * ssin(o,p,i) * invcp0(o,i) * kdk(o,i)

            if(ssin(o,p,i) > 0)then ! positive stress
              ! wind pushing waves
              taux1(i) = taux1(i) + dummy(o,p,i) * cth(p)
              tauy1(i) = tauy1(i) + dummy(o,p,i) * sth(p)
            else ! negative stress, two cases
              if(cos(wdir(i)-th(p)) < 0)then
                ! waves against wind
                taux2(i) = taux2(i) + dummy(o,p,i) * cth(p)
                tauy2(i) = tauy2(i) + dummy(o,p,i) * sth(p)
              else
                ! waves overrunning wind
                taux3(i) = taux3(i) + dummy(o,p,i) * cth(p)
                tauy3(i) = tauy3(i) + dummy(o,p,i) * sth(p)
              end if
            end if

          end do
        end do
      end do

      taux1 = taux1 * rhow(istart:iend) * dthg
      tauy1 = tauy1 * rhow(istart:iend) * dthg
      taux2 = taux2 * rhow(istart:iend) * dthg
      tauy2 = tauy2 * rhow(istart:iend) * dthg
      taux3 = taux3 * rhow(istart:iend) * dthg
      tauy3 = tauy3 * rhow(istart:iend) * dthg

    end if ! if(option=='ocn')

    if (option == 'atm') then

      ! compute form drag from atmosphere: (positive into waves)
      call stress_vector_integral(e(:,:,istart:iend), ssin, cp0, taux_util, tauy_util)

      do concurrent(i = istart:iend)

        ! compute the tail
        tailatmx(i) = taux_util(om,i) * k(om,i) * tail(i) * rhow(i) * dthg
        tailatmy(i) = tauy_util(om,i) * k(om,i) * tail(i) * rhow(i) * dthg

        taux_form(i) = sum(taux_util(:,i) * kdk(:,i), dim=1)&
                     * rhow(i) * dthg + tailatmx(i)
        tauy_form(i) = sum(tauy_util(:,i) * kdk(:,i), dim=1)&
                     * rhow(i) * dthg + tailatmy(i)

        taux_diag(i) = sum(taux_util(oc(i):om,i) * kdk(oc(i):om,i), dim=1)&
                     * rhow(i) * dthg
        tauy_diag(i) = sum(tauy_util(oc(i):om,i) * kdk(oc(i):om,i), dim=1)&
                     * rhow(i) * dthg

      end do

      ! compute surface stokes drift velocities:
      do concurrent(i = istart:iend, p = 1:pm)
        spectrum_integral = sum(util(:,i,1) * e(:,p,i))
        usurf(i) = spectrum_integral * cth(p)
        vsurf(i) = spectrum_integral * sth(p)
      end do

      ! wind speed and direction relative to surface velocity
      call wind_relative(wspd(istart:iend), wdir(istart:iend),&
                         uc(istart:iend) + usurf, vc(istart:iend) + vsurf,&
                         wspdrel, wdirrel)

      ! form-induced drag coefficient
      cd_form = drag_coefficient(taux_form, tauy_form,&
                                 rhoa(istart:iend), wspd(istart:iend))

      ! skin-induced drag coefficient
      cd_skin = drag_coefficient_skin(cd_form, wspd(istart:iend), wspdrel, z, nu_air, kappa)

      taux_skin = rhoa(istart:iend) * cd_skin * wspdrel**2 * cos(wdirrel)
      tauy_skin = rhoa(istart:iend) * cd_skin * wspdrel**2 * sin(wdirrel)

      taux = taux_form + taux_skin
      tauy = tauy_form + tauy_skin

      ! total (form + skin) drag coefficient
      cd = drag_coefficient(taux, tauy, rhoa(istart:iend), wspd(istart:iend))
      where (cd > 5e-3) cd = 5e-3

      ! update friction velocity (form + skin)
      ustar = sqrt(cd) * wspd(istart:iend)

    end if ! if(option=='atm')

  end subroutine stress

  real(rk) pure elemental function drag_coefficient(taux, tauy, rhoa, wspd) result(cd)
    !! Computes drag coefficient from stress vector (taux, tauy, N/m^2),
    !! air density (rhoa, kg/m^3), and wind speed (wspd, m/s).
    real(rk), intent(in) :: taux, tauy, rhoa, wspd
    cd = sqrt(taux**2 + tauy**2) / (rhoa * wspd**2)
  end function drag_coefficient

  real(rk) pure elemental function drag_coefficient_skin(cd_form, wspd, wspdrel, z, air_viscosity, von_karman) result(cd)
    !! Computes the skin drag coefficient attenuated by from drag,
    !! given input absolute (wspd) and relative (wspdrel) wind speed (m/s),
    !! height z (m), air viscosity (m^2/s), and Von Karman constant.
    real(rk), intent(in) :: cd_form, wspd, wspdrel, z, air_viscosity, von_karman
    cd = friction_velocity_skin(wspdrel, z, air_viscosity, von_karman)**2 / wspd**2
    cd = cd * (1 + 2 * cd / (cd + cd_form + tiny(cd))) / 3
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

  pure subroutine stress_vector_integral(e, src, cp, taux, tauy)
    !! Integrates the input source function src over the spectrum with
    !! wave variance e over the directions.
    real(rk), intent(in) :: e(:,:,:), src(:,:,:), cp(:,:)
    real(rk), intent(out) :: taux(:,:), tauy(:,:)
    integer :: i, p, o
    integer :: dim(3)
    dim = shape(e)
    taux = 0
    tauy = 0
    do i = 1, dim(3)
      do p = 1, dim(2)
        do o = 1, dim(1)
          taux(o,i) = taux(o,i) + e(o,p,i) * src(o,p,i) * cth(p) / cp(o,i)
          tauy(o,i) = tauy(o,i) + e(o,p,i) * src(o,p,i) * sth(p) / cp(o,i)
        end do
      end do
    end do
  end subroutine stress_vector_integral

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
