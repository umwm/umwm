module umwm_physics

  use umwm_forcing, only: forcing_type
  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type
  use umwm_state, only: state_type

  implicit none

  private

  public :: &
    wind_input_donelan2012, &
    wave_dissipation_donelan2012

  real, parameter :: pi = 4 * atan(1.)

contains

  pure subroutine wind_input_donelan2012(state, forcing, grid, spectrum)
    ! Wind input function based on Jeffreys's sheltering hypothesis
    ! and described by Donelan et al. (2012).
    type(state_type), intent(inout) :: state
    type(forcing_type), intent(in) :: forcing
    type(grid_type), intent(in) :: grid
    type(spectrum_type), intent(in) :: spectrum

    integer :: i, j, nf, nd
    real :: log_lambda_half_over_z

    associate( &
      s_in => state % wind_input, &
      wspd => sqrt(forcing % u_atmosphere**2 + forcing % v_atmosphere**2), &
      wdir => atan2(forcing % v_atmosphere, forcing % u_atmosphere), &
      uc => forcing % u_ocean, &
      vc => forcing % v_ocean, &
      density_ratio => forcing % density_atmosphere / forcing % density_ocean, &
      cp0 => state % phase_speed, &
      theta => spectrum % direction, &
      cos_theta => cos(spectrum % direction), &
      sin_theta => sin(spectrum % direction), &
      omega => 2 * pi * spectrum % frequency, &
      wavelength => 2 * pi / state % wavenumber &
    )

      do j = 1, grid % size_y
        do i = 1, grid % size_x

          ! TODO: prognostic/diagnostic cutoff frequency
          ! based on 4*PM peak frequency
          do nd = 1, spectrum % num_directions
            do nf = 1, spectrum % num_frequencies

              !TODO z is hardcoded to 10 for now
              log_lambda_half_over_z = log(min(wavelength(nf,i,j) / 2, 10.) / state % depth(i,j))

              ! Compute wind input at height of half wavelength:
              ! TODO: atmospheric stability
              s_in(nf,nd,i,j) = ( &
                wspd(i,j) + &
                2.5 * state % ustar(i,j) * log_lambda_half_over_z * cos(wdir(i,j) - theta(nd)) &
                - cp0(nf,i,j) &
                - uc(i,j) * cos_theta(nd) &
                - vc(i,j) * sin_theta(nd) &
              ) / cp0(nf,i,j)**2

              ! TODO ice fraction
              s_in(nf,nd,i,j) = density_ratio(i,j) * abs(s_in(nf,nd,i,j)) * s_in(nf,nd,i,j) * omega(nf)

              ! TODO variable sheltering coefficient; hardcoded for now
              if (s_in(nf,nd,i,j) > 0) then
                ! Growth rate
                s_in(nf,nd,i,j) = 0.11 * s_in(nf,nd,i,j)
              else
                ! Decay rate
                ! TODO differentiate between opposing and outrunning waves
                s_in(nf,nd,i,j) = 0.01 * s_in(nf,nd,i,j)
              end if

            end do
          end do
        end do
      end do

    end associate

  end subroutine wind_input_donelan2012


  pure subroutine wave_dissipation_donelan2012(state, grid, spectrum)
    ! Wave dissipation function described by Donelan et al. (2012).
    type(state_type), intent(inout) :: state
    type(grid_type), intent(in) :: grid
    type(spectrum_type), intent(in) :: spectrum

    integer :: i, j, nf, nd

    real :: cum_mss(spectrum % num_frequencies, spectrum % num_directions)

    real, parameter :: mss_fac = 360 ! FIXME: hardcoded for now
    real, parameter :: sds_fac = 42 ! FIXME: hardcoded for now
    real, parameter :: sds_power = 2.5 ! FIXME: hardcoded for now

    associate( &
      omega => 2 * pi * spectrum % frequency, &
      direction => spectrum % direction, &
      k => state % wavenumber, &
      dk => state % dk, &
      variance => state % variance, &
      dissipation => state % dissipation &
    )

      do j = 1, grid % size_y
        do i = 1, grid % size_x

          cum_mss = 0

          do nd = 1, spectrum % num_directions
            do nf = 2, spectrum % num_frequencies
              cum_mss(nf,nd) = cum_mss(nf-1,nd) + &
                sum(state % variance(nf-1,:,i,j) * cos(direction(:))) * &
                k(nf-1,i,j)**3 * dk(nf-1,i,j)
            end do
          end do

          cum_mss = (1 + mss_fac * cum_mss)**2

          do nd = 1, spectrum % num_directions
            do nf = 1, spectrum % num_frequencies
              dissipation(nf,nd,i,j) = sds_fac * omega(nf) * cum_mss(nf,nd) * &
                (variance(nf,nd,i,j) * k(nf,i,j)**4 * dk(nf,i,j))**sds_power
            end do
          end do

        end do
      end do

    end associate

  end subroutine wave_dissipation_donelan2012
  !
  !
  !subroutine wave_attenuation_seaice_kohout2014()
  !  ! Wave attenuation by sea ice, following kohout et al. (2014).
  !
  !  integer :: i, o, p
  !
  !  ! parameters from Kohout et al. 2014
  !  real, parameter :: H_th = 3.0       ! [m]
  !  real, parameter :: C1   = -5.35e-6  ! [m-1]
  !  real, parameter :: C2   = C1 * H_th ! []
  !
  !  real, dimension(om,pm) :: spectrumbin
  !
  !  real :: ht_
  !  
  !  sice = 0.0
  !
  !  do i = istart, iend
  !
  !    if (fice(i) > fice_lth) then
  ! 
  !      ht_ = 0.0
  ! 
  !      do p = 1, pm
  !        do o = 1, om
  !          spectrumbin(o,p) = e(o,p,i) * kdk(o,i)
  ! 
  !          ht_ = ht_ + spectrumbin(o,p)
  !        end do
  !      end do
  ! 
  !      ht_ = 4 * sqrt(ht_ * dth) ! significant wave height
  ! 
  !      ! wave attenuation from sea ice in the two SWH regimes
  !      if (ht_ < H_th) then
  !        sice(:,i) = C1 * ht_
  !      else
  !        sice(:,i) = C2
  !      end if
  !      
  !      sice(:,i) = 2 * cg0(:,i) * sice(:,i)
  !       
  !    end if
  ! 
  !  end do
  !
  !end subroutine wave_attenuation_seaice_kohout2014
  !
  !
  !subroutine wave_transfer_donelan2012()
  !
  !  integer :: o, p, i
  !
  !  snl = 0
  !
  !  ! spread wave energy to 2 next longer wavenumbers exponentially decaying
  !  ! as distance from donating wavenumber, and remove the energy from
  !  ! donating wavenumbers:
  !  do concurrent(i = istart:iend)
  !    do concurrent(o = 1:oc(i), p = 1:pm)
  !      snl(o,p,i) = bf1_renorm(o,i) * sds(o+1,p,i) * e(o+1,p,i)&
  !                 + bf2_renorm(o,i) * sds(o+2,p,i) * e(o+2,p,i)&
  !                 - snl_fac * sds(o,p,i) * e(o,p,i)
  !    end do
  !  end do
  !
  !  ! account for plunging breakers
  !  do concurrent(o = 1:om, p = 1:pm, i = istart:iend)
  !    sds(o,p,i) = sds(o,p,i) * cothkd(o,i)
  !  end do
  !
  !  ! compute dissipation due to turbulence
  !  do concurrent(o = 1:om, i = istart:iend)
  !    sdt(o,i) = sdt_fac * sqrt(rhorat(i)) * ustar(i) * k(o,i)
  !  end do
  !
  !end subroutine wave_transfer_donelan2012


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
