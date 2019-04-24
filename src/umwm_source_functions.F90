module umwm_source_functions
  ! Module that provides wave source functions.
  use umwm_io, only: currents,seaice
!  use umwm_module
  use umwm_constants, only: rk
  use umwm_sheltering, only: sheltering_coare35, sheltering_reynolds

  implicit none

  private

  public :: sin_d12, sds_d12, snl_d12, s_ice

contains

   subroutine sin_d12( wspd, istart, iend, iistart, iiend, fprog, om, oc, pm, f, g, ustar, logl2overz, psim, psiml2, wdir, th, cp0, uc, cth, vc, sth, sin_fac, fieldscale1, fieldscale2, twopi, fkovg, fice, rhorat, fcutoff, shelt, ssin)


    implicit NONE


 ! input parameters


    integer, intent(in)                  :: istart, iend, om, pm, iistart, iiend

    real,   intent(inout)                   :: wspd(istart : iend)

    real,   intent(in)                   :: twopi, g, fprog, fieldscale1, fieldscale2, sin_fac

    real,   intent(in)                   :: fkovg(om,istart:iend), f(om)

    real,   intent(in)                   :: ustar(istart:iend), psim(istart:iend), wdir(istart:iend), fice(istart:iend), rhorat(istart:iend), uc(istart:iend), vc(istart:iend)

    real,   intent(in)                   :: logl2overz(om,istart:iend), psiml2(om,istart:iend)

    real,   intent(in)                   :: cp0(om,iistart-1:iiend)

    real,   intent(in)                   :: th(pm), cth(pm), sth(pm)



!  local variables

    integer                              :: i, m, n, o, p


!  output variables


   integer, dimension(istart:iend), intent(out)      :: oc

    real, intent(inout)                  :: fcutoff(istart:iend)

    real, intent(out)                    :: shelt(istart:iend), ssin(om,pm,istart:iend)

!   Wind input function based on Jeffreys's sheltering hypothesis
    ! and described by Donelan et al. (2012).


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

  end subroutine sin_d12

  subroutine sds_d12 (e, mss_fac, om, pm, iistart,iiend, istart, iend, cth2pp, k3dk, twopisds_fac, f, k4, sds_power, sds, dummy )



    implicit NONE


    ! Input variables

    integer, intent(in)             :: om, pm, istart,iend,iistart,iiend

    real , intent(in)               :: mss_fac, twopisds_fac, sds_power

    real, intent(in)                :: cth2pp(pm, pm), k3dk(om,istart:iend), k4(om,istart:iend), f(om)


    real, intent(in)                :: e(om,pm,iistart-1:iiend)



    ! local variables

    integer                         :: i, o, p


    !output variables

    real, intent(out)             :: sds(om,pm,istart:iend)

    real             :: dummy(om,pm,istart:iend)
    ! Wave dissipation function described by Donelan et al. (2012).

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
  




  subroutine s_ice(istart, iend, iistart, iiend, fice, fice_lth, fice_uth, om, pm , e, kdk, dth, cg0, sice)

    ! Wave attenuation by sea ice, following Kohoun et al. (2014).

     implicit NONE



   ! Input variables

     integer , intent(in)                :: istart, iend, om, pm, iistart, iiend


     real,     intent(in)                :: fice(istart:iend), e(om,pm,iistart-1:iiend)
     real,     intent(in)                :: fice_lth, fice_uth, dth
     real,     intent(in)                :: kdk(om,istart:iend), cg0(om,iistart-1:iiend)


    ! local variables

     integer :: i, o, p

     ! parameters from Kohout et al. 2014
     real, parameter :: H_th = 3.0       ! [m]
     real, parameter :: C1   = -5.35e-6  ! [m-1]
     real, parameter :: C2   = C1 * H_th ! []

     real, dimension(om,pm) :: spectrumbin

     real :: ht_


   ! output variables


     real , intent(inout)                :: sice(om,istart:iend)
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

  end subroutine s_ice

   subroutine snl_d12(iistart,iiend,istart, iend, oc, om, pm, bf1_renorm, bf2_renorm, e, snl_fac, cothkd, sdt_fac, rhorat, ustar, k, sds, sdt, snl  )

    ! input variables


    integer, intent(in)     :: istart, iend, pm, om, iistart,iiend

    real,    intent(in)     :: bf1_renorm(om,istart:iend), bf2_renorm(om,istart:iend), cothkd(om,istart:iend), k(om,istart:iend)

    real,    intent(in)     :: e(om, pm, iistart-1:iiend)



    real, intent(in)        :: snl_fac, sdt_fac

    real, intent(in)        :: rhorat(istart:iend)

    integer, dimension(istart:iend), intent(in)      :: oc

    real, intent(in)        :: ustar(istart:iend)


    ! local variables

    integer :: o, p, i

    ! output variables


    real, intent(inout)   :: sds(om,pm,istart:iend), sdt(om,istart:iend), snl(om,pm,istart:iend)

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
