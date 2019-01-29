module umwm_source_functions
  ! Module that provides wave source functions.
  use umwm_io, only: currents
  use umwm_module
  use umwm_constants, only: rk
  use umwm_sheltering, only: sheltering_coare35, sheltering_reynolds

  implicit none

  private

  public :: sin_d12, sds_d12, snl_d12

contains

  subroutine sin_d12()
    ! Wind input function based on Jeffreys's sheltering hypothesis
    ! and described by Donelan et al. (2012).
    integer :: i, m, n, o, p

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
#ifdef ESMF
      ssin(o,p,i) = (wspd(i) + 2.5 * ustar(i) * (logl2overz(o,i) + psim(i) - psiml2(o,i)))&
#else
      ssin(o,p,i) = (wspd(i) + 2.5 * ustar(i) * logl2overz(o,i))&
#endif
                  * cos(wdir(i) - th(p)) - cp0(o,i)&
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
      ssin(o,p,i) = twopi * rhorat(i) * ssin(o,p,i) * fkovg(o,i)
    end do

    ! prevent negative sin for diagnostic tail
    do i = istart, iend
      do p = 1, pm
        do o = oc(i)+1, om
          ssin(o,p,i) = max(ssin(o,p,i), 0._rk)
        end do
      end do
    end do

  end subroutine sin_d12


  subroutine sds_d12
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

  end subroutine sds_d12

  
  subroutine snl_d12

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

  end subroutine snl_d12

end module umwm_source_functions
