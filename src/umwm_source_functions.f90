module umwm_source_functions
  ! Module that provides wave source functions.
  use umwm_config, only: config_type
  use umwm_forcing, only: forcing_type, min_wind_speed
  use umwm_grid, only: grid_type
  use umwm_module, only: bf1_renorm, bf2_renorm, cg0, cothkd, cp0, cth, &
                         cth2pp, dth, dummy, e, f, fcutoff, &
                         fieldscale1, fieldscale2, fkovg, &
                         k, k3dk, k4, kdk, logl2overz, &
                         oc, psim, psiml2, sds, &
                         sdt, shelt, sice, &
                         snl, ssin, sth, th, twopi, twopisds_fac, &
                         ustar
  use umwm_constants, only: rk
  use umwm_sheltering, only: sheltering_coare35, sheltering_reynolds
  use umwm_spectrum, only: spectrum_type

  implicit none

  private

  public :: sin_d12, sds_d12, snl_d12, s_ice

contains

  subroutine sin_d12(config, spectrum, grid, forcing)
    ! Wind input function based on Jeffreys's sheltering hypothesis
    ! and described by Donelan et al. (2012).
    type(config_type), intent(in) :: config
    type(spectrum_type), intent(in) :: spectrum
    type(grid_type), intent(in) :: grid
    type(forcing_type), intent(in) :: forcing
    integer :: i, o, p
    real :: wspd_clamped(grid % istart:grid % iend)

    associate(istart => grid % istart, iend => grid % iend)

    ! protection against low wind speed values
    wspd_clamped = max(forcing % wspd(istart:iend), min_wind_speed)

    ! cut-off frequency (4*pierson-moskowitz peak frequency)
    fcutoff(istart:iend) = 0.53 * config % g / wspd_clamped

    where (fcutoff > config % fprog) fcutoff = config % fprog

    ! search for the cut-off frequency bin:
    do i = istart, iend
      do o = spectrum % num_frequencies - 2, 1, -1
        oc(i) = o
        if (fcutoff(i) > f(o)) exit
      end do
    end do

    ! compute wind input at height of half wavelenght:
    do concurrent (o = 1:spectrum % num_frequencies,p = 1:spectrum % num_directions, i = istart:iend)
#ifdef ESMF
      ssin(o,p,i) = (wspd_clamped(i) + 2.5 * ustar(i) * (logl2overz(o,i) + psim(i) - psiml2(o,i)))&
#else
      ssin(o,p,i) = (wspd_clamped(i) + 2.5 * ustar(i) * logl2overz(o,i))&
#endif
                  * cos(forcing % wdir(i) - th(p)) - cp0(o,i)&
                  - forcing % uc(i) * cth(p) - forcing % vc(i) * sth(p)
    end do

    ssin = config % sin_fac * abs(ssin) * ssin

    ! compute variable sheltering coefficient
    shelt = sheltering_coare35(wspd_clamped)

    ! apply variable sheltering coefficient
    do concurrent(o = 1:spectrum % num_frequencies, p = 1:spectrum % num_directions, i = istart:iend, ssin(o,p,i) > 0)
      ssin(o,p,i) = ssin(o,p,i) * shelt(i) / config % sin_fac
    end do

    ! adjust input for opposing winds
    do concurrent (o = 1:spectrum % num_frequencies, p = 1:spectrum % num_directions, i = istart:iend, ssin(o,p,i) < 0)
      ssin(o,p,i) = ssin(o,p,i) * fieldscale1
    end do

    ! further reduce for swell that overruns the wind
    do concurrent (o = 1:spectrum % num_frequencies, p = 1:spectrum % num_directions, &
                   i = istart:iend, ssin(o,p,i) < 0 .and. cos(forcing % wdir(i)-th(p)) > 0)
      ssin(o,p,i) = ssin(o,p,i) * fieldscale2
    end do

    do concurrent (o = 1:spectrum % num_frequencies, p = 1:spectrum % num_directions, i = istart:iend)
      ssin(o,p,i) = (1 - forcing % fice(i)) * twopi * forcing % rhorat(i) * ssin(o,p,i) * fkovg(o,i)
    end do

    ! prevent negative sin for diagnostic tail
    do i = istart, iend
      do p = 1, spectrum % num_directions
        do o = oc(i)+1, spectrum % num_frequencies
          ssin(o,p,i) = max(ssin(o,p,i), 0._rk)
        end do
      end do
    end do

    end associate

  end subroutine sin_d12


  subroutine sds_d12(config, spectrum, grid)
    ! Wave dissipation function described by Donelan et al. (2012).
    type(config_type), intent(in) :: config
    type(spectrum_type), intent(in) :: spectrum
    type(grid_type), intent(in) :: grid
    integer :: i, o, p

    associate(istart => grid % istart, iend => grid % iend)

    dummy = 0

    if (config % mss_fac > 0) then
      do concurrent(p = 1:spectrum % num_directions, i = istart:iend)
        do o = 2, spectrum % num_frequencies
          dummy(o,p,i) = dummy(o-1,p,i) + sum(e(o-1,:,i) * cth2pp(:,p)) * k3dk(o-1,i)
        end do
      end do
    end if

    dummy = (1 + config % mss_fac * dummy)**2

    do concurrent(o = 1:spectrum % num_frequencies, p = 1:spectrum % num_directions, i = istart:iend)
      sds(o,p,i) = twopisds_fac * f(o) * dummy(o,p,i) * (e(o,p,i) * k4(o,i))**config % sds_power
    end do

    end associate

  end subroutine sds_d12


  subroutine s_ice(config, spectrum, grid, forcing)
    ! Wave attenuation by sea ice, following Kohoun et al. (2014).

    type(config_type), intent(in) :: config
    type(spectrum_type), intent(in) :: spectrum
    type(grid_type), intent(in) :: grid
    type(forcing_type), intent(in) :: forcing
    integer :: i, o, p

    ! parameters from Kohout et al. 2014
    real, parameter :: H_th = 3.0       ! [m]
    real, parameter :: C1   = -5.35e-6  ! [m-1]
    real, parameter :: C2   = C1 * H_th ! []

    real, dimension(spectrum % num_frequencies, spectrum % num_directions) :: spectrumbin

    real :: ht_

    associate(istart => grid % istart, iend => grid % iend)

    sice = 0.0

    do i = istart, iend

      if (forcing % fice(i) > config % fice_lth) then

        ht_ = 0.0

        do p = 1, spectrum % num_directions
          do o = 1, spectrum % num_frequencies
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

    end associate

  end subroutine s_ice


  subroutine snl_d12(config, spectrum, grid, forcing)

    type(config_type), intent(in) :: config
    type(spectrum_type), intent(in) :: spectrum
    type(grid_type), intent(in) :: grid
    type(forcing_type), intent(in) :: forcing
    integer :: o, p, i

    associate(istart => grid % istart, iend => grid % iend)

    snl = 0

    ! spread wave energy to 2 next longer wavenumbers exponentially decaying
    ! as distance from donating wavenumber, and remove the energy from
    ! donating wavenumbers:
    do concurrent(i = istart:iend)
      do concurrent(o = 1:oc(i), p = 1:spectrum % num_directions)
        snl(o,p,i) = bf1_renorm(o,i) * sds(o+1,p,i) * e(o+1,p,i)&
                   + bf2_renorm(o,i) * sds(o+2,p,i) * e(o+2,p,i)&
                   - config % snl_fac * sds(o,p,i) * e(o,p,i)
        !snl(o,p,i) = snl_fac * (kdk(o+1,i) * e(o+1,p,i) - kdk(o,i) * e(o,p,i)) / dwn(o,i) ! WIP dk-invariant Snl
      end do
    end do

    ! account for plunging breakers
    do concurrent(o = 1:spectrum % num_frequencies, p = 1:spectrum % num_directions, i = istart:iend)
      sds(o,p,i) = sds(o,p,i) * cothkd(o,i)
    end do

    ! compute dissipation due to turbulence
    do concurrent(o = 1:spectrum % num_frequencies, i = istart:iend)
      sdt(o,i) = config % sdt_fac * sqrt(forcing % rhorat(i)) * ustar(i) * k(o,i)
    end do

    end associate

  end subroutine snl_d12

end module umwm_source_functions
