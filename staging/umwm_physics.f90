module umwm_physics
  !use umwm_io, only: currents,seaice
  use umwm_state, only: state_type
  use umwm_forcing, only: forcing_type
  use umwm_spectrum, only: spectrum_type

  implicit none

  private

  public :: &
    wind_input_donelan2012, &
    wave_dissipation_donelan2012, &
    wave_transfer_donelan2012, &
    wave_attenuation_seaice_kohout2014

contains

  subroutine wind_input_donelan2012(state, forcing, spectrum)
    ! Wind input function based on Jeffreys's sheltering hypothesis
    ! and described by Donelan et al. (2012).
    type(state_type), intent(inout) :: state
    type(forcing_type), intent(in) :: forcing
    type(spectrum_type), intent(in) :: spectrum

    integer :: i, o, p

    ! protection against low wind speed values
    wspd = max(wspd, 1e-2)

    ! cut-off frequency (4*pierson-moskowitz peak frequency)
    fcutoff(istart:iend) = 0.53 * g / wspd(istart:iend)

    where (fcutoff > fprog) fcutoff = fprog

    ! search for the cut-off frequency bin:
    do i = istart, iend
      do o = om-2, 1, -1
        oc(i) = o
        !oc(i) = om-2
        if (fcutoff(i) > f(o)) exit
      end do
    end do

    ! compute wind input at height of half wavelenght:
    do concurrent (o = 1:om,p = 1:pm, i = istart:iend)
      ssin(o,p,i) = (wspd(i) + 2.5 * ustar(i) &
                  * (logl2overz(o,i) + psim(i) - psiml2(o,i))) &
                  * cos(wdir(i) - th(p)) - cp0(o,i) &
                  - uc(i) * cth(p) - vc(i) * sth(p)
    end do

    ssin = sin_fac * abs(ssin) * ssin

    ! compute variable sheltering coefficient
    shelt = sheltering_coare35(wspd(istart:iend))

    ! apply variable sheltering coefficient
    do concurrent(o = 1:om, p = 1:pm, i = istart:iend, ssin(o,p,i) > 0)
      ssin(o,p,i) = ssin(o,p,i) * shelt(i) / sin_fac
    end do

    ! adjust input for opposing winds
    do concurrent (o = 1:om, p = 1:pm, i = istart:iend, ssin(o,p,i) < 0)
      ssin(o,p,i) = ssin(o,p,i) * fieldscale1
    end do

    ! further reduce for swell that overruns the wind
    do concurrent (o = 1:om, p = 1:pm, i = istart:iend, ssin(o,p,i) < 0 .and. cos(wdir(i)-th(p)) > 0)
      ssin(o,p,i) = ssin(o,p,i) * fieldscale2
    end do

    do concurrent (o = 1:om, p = 1:pm, i = istart:iend)
      ssin(o,p,i) = (1 - fice(i)) * twopi * rhorat(i) * ssin(o,p,i) * fkovg(o,i)
    end do

    ! prevent negative sin for diagnostic tail
    do i = istart, iend
      do p = 1, pm
        do o = oc(i)+1, om
          ssin(o,p,i) = max(ssin(o,p,i), 0._rk)
        end do
      end do
    end do

  end subroutine wind_input_donelan2012


  subroutine wave_dissipation_donelan2012()
    ! Wave dissipation function described by Donelan et al. (2012).
    integer :: i, o, p

    dummy = 0

    if (mss_fac > 0) then
      do concurrent(p = 1:pm, i = istart:iend)
        do o = 2, om
          dummy(o,p,i) = dummy(o-1,p,i) + sum(e(o-1,:,i) * cth2pp(:,p)) * k3dk(o-1,i)
        end do
      end do
    end if

    dummy = (1 + mss_fac * dummy)**2

    do concurrent(o = 1:om, p = 1:pm, i = istart:iend)
      sds(o,p,i) = twopisds_fac * f(o) * dummy(o,p,i) * (e(o,p,i) * k4(o,i))**sds_power
    end do

  end subroutine wave_dissipation_donelan2012


  subroutine wave_attenuation_seaice_kohout2014()
    ! Wave attenuation by sea ice, following kohout et al. (2014).

    integer :: i, o, p

    ! parameters from Kohout et al. 2014
    real, parameter :: H_th = 3.0       ! [m]
    real, parameter :: C1   = -5.35e-6  ! [m-1]
    real, parameter :: C2   = C1 * H_th ! []

    real, dimension(om,pm) :: spectrumbin

    real :: ht_
  
    sice = 0.0

    do i = istart, iend

      if (fice(i) > fice_lth) then
 
        ht_ = 0.0
 
        do p = 1, pm
          do o = 1, om
            spectrumbin(o,p) = e(o,p,i) * kdk(o,i)
 
            ht_ = ht_ + spectrumbin(o,p)
          end do
        end do
 
        ht_ = 4 * sqrt(ht_ * dth) ! significant wave height
 
        ! wave attenuation from sea ice in the two SWH regimes
        if (ht_ < H_th) then
          sice(:,i) = C1 * ht_
        else
          sice(:,i) = C2
        end if
        
        sice(:,i) = 2 * cg0(:,i) * sice(:,i)
         
      end if
 
    end do

  end subroutine wave_attenuation_seaice_kohout2014

  
  subroutine wave_transfer_donelan2012()

    integer :: o, p, i

    snl = 0

    ! spread wave energy to 2 next longer wavenumbers exponentially decaying
    ! as distance from donating wavenumber, and remove the energy from
    ! donating wavenumbers:
    do concurrent(i = istart:iend)
      do concurrent(o = 1:oc(i), p = 1:pm)
        snl(o,p,i) = bf1_renorm(o,i) * sds(o+1,p,i) * e(o+1,p,i)&
                   + bf2_renorm(o,i) * sds(o+2,p,i) * e(o+2,p,i)&
                   - snl_fac * sds(o,p,i) * e(o,p,i)
      end do
    end do

    ! account for plunging breakers
    do concurrent(o = 1:om, p = 1:pm, i = istart:iend)
      sds(o,p,i) = sds(o,p,i) * cothkd(o,i)
    end do

    ! compute dissipation due to turbulence
    do concurrent(o = 1:om, i = istart:iend)
      sdt(o,i) = sdt_fac * sqrt(rhorat(i)) * ustar(i) * k(o,i)
    end do

  end subroutine wave_transfer_donelan2012


  pure real function polyval(p, x)
    ! Evaluates a polynomial of order p over values of x:
    !   y = p(1)*x**n + p(2)*x**(n-1) + ... + p(n)
    real, intent(in) :: p(:), x
    integer :: n
    polyval = 0
    do n = 1, size(p)
      polyval = polyval + p(n) * x**(size(p) - n)
    end do
  end function polyval


  real elemental function sheltering_coare35(wind_speed) result(res)
    ! Returns the sheltering coefficient given input wind speed
    ! to approximately match the momentum flux of COARE 3.5 algorithm.
    real, intent(in) :: wind_speed ! wind speed [m/s]
    real, parameter :: x1 = 15, x2 = 33, x3 = 60
    real, parameter :: y1 = 0.10, y2 = 0.09, y3 = 0.06
    real, parameter :: intercept = 0.04, curvature = 0.65, decay = 1.6

    real :: a, b, c, s1, slope1, slope2

    ! low range linear growth
    slope1 = (y1 - intercept) / x1
    slope2 = (y3 - y2) / (x3 - x2)

    ! medium range fitting (quadratic)
    c = curvature * (slope2 - slope1) / (x2 - x1)
    b = slope1 - 2 * c * x1
    a = y1 - b * x1 - c * x1**2

    ! high end decay
    s1 = a + b * x2 + c * x2**2

    if (wind_speed <= x1) then
      res = intercept + slope1 * wind_speed
    else if (wind_speed > x1 .and. wind_speed <= x2) then
      res = a + b * wind_speed + c * wind_speed**2
    else
      res = s1 * exp(- (wind_speed - x2) / (decay * wind_speed))
    end if

  end function sheltering_coare35

end module umwm_physics
