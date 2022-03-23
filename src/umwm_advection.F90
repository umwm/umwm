module umwm_advection

  ! Subroutines for calculation of wave energy advection in geographical (x, y)
  ! and directional (theta) space. 1st order upstream finite differences.

#if defined(MPI)
  use mpi
#endif
  use umwm_module
#ifndef ESMF
  use umwm_io, only: currents
#endif

  implicit none

  logical :: zerocurrents

contains

  subroutine propagation

    ! 1st order upstream finite difference advection in geographical space.

    integer :: o, p, i
    real :: cge, cgw, cgn, cgs
    real :: feup, fedn, fwup, fwdn, fnup, fndn, fsup, fsdn

    real :: flux(om,pm,istart:iend)

    flux = 0

    do concurrent(i = istart:iend)
      do concurrent(o = 1:oc(i), p = 1:pm)

        !  double group velocity at east, west, north, south cell edges
        cge = cg0(o,i) * cth_curv(p,i) + cg0(o,ie(i)) * cth_curv(p,i)
        cgw = cg0(o,i) * cth_curv(p,i) + cg0(o,iw(i)) * cth_curv(p,i)
        cgs = cg0(o,i) * sth_curv(p,i) + cg0(o,is(i)) * sth_curv(p,i)
        cgn = cg0(o,i) * sth_curv(p,i) + cg0(o,in(i)) * sth_curv(p,i)

        ! advective energy flux in and out of all 4 directions
        flux(o,p,i) = (cge + abs(cge)) * dye(i) * e(o,p,i)     & ! east, out
                    + (cge - abs(cge)) * dye(i) * e(o,p,ie(i)) & ! east, in
                    - (cgw + abs(cgw)) * dyw(i) * e(o,p,iw(i)) & ! west, in
                    - (cgw - abs(cgw)) * dyw(i) * e(o,p,i)     & ! west, out
                    + (cgn + abs(cgn)) * dxn(i) * e(o,p,i)     & ! north, out
                    + (cgn - abs(cgn)) * dxn(i) * e(o,p,in(i)) & ! north, in
                    - (cgs + abs(cgs)) * dxs(i) * e(o,p,is(i)) & ! south, in
                    - (cgs - abs(cgs)) * dxs(i) * e(o,p,i)       ! south, out

      end do
    end do

    ! check if currents are non-zero:
    if (.not. isglobal) then
      zerocurrents = .not. (any(uc(iistart:iiend) /= 0)&
                       .or. any(vc(iistart:iiend) /= 0))
    else
      zerocurrents = .not. (any(uc /= 0) .or. any(vc /= 0))
    end if

    if (.not. zerocurrents) then ! advect wave energy by currents

      do concurrent(i = istart:iend)

        ! x-direction
        feup = (uc(i) + uc(iie(i)) + abs(uc(i) + uc(iie(i)))) * dye(i)
        fedn = (uc(i) + uc(iie(i)) - abs(uc(i) + uc(iie(i)))) * dye(i)
        fwup = (uc(i) + uc(iiw(i)) + abs(uc(i) + uc(iiw(i)))) * dyw(i)
        fwdn = (uc(i) + uc(iiw(i)) - abs(uc(i) + uc(iiw(i)))) * dyw(i)

        ! y-direction
        fnup = (vc(i) + vc(iin(i)) + abs(vc(i) + vc(iin(i)))) * dxn(i)
        fndn = (vc(i) + vc(iin(i)) - abs(vc(i) + vc(iin(i)))) * dxn(i)
        fsup = (vc(i) + vc(iis(i)) + abs(vc(i) + vc(iis(i)))) * dxs(i)
        fsdn = (vc(i) + vc(iis(i)) - abs(vc(i) + vc(iis(i)))) * dxs(i)

        do concurrent(o = 1:oc(i), p = 1:pm)
          flux(o,p,i) = flux(o,p,i)                                &
                      + (feup * e(o,p,i)     + fedn * e(o,p,ie(i)) &
                      -  fwup * e(o,p,iw(i)) - fwdn * e(o,p,i)     &
                      +  fnup * e(o,p,i)     + fndn * e(o,p,in(i)) &
                      -  fsup * e(o,p,is(i)) - fsdn * e(o,p,i))
        end do

      end do

    end if ! .not. zerocurrents

    ! integrate in time
    do concurrent(i = istart:iend)
      do concurrent(o = 1:oc(i), p = 1:pm)
          ef(o,p,i) = ef(o,p,i) - 0.25 * dta * flux(o,p,i) * oneovar(i)

          if (fice(i) > fice_uth) then
            ef(o,p,i) = 0.0
          end if  
      end do
    end do

  end subroutine propagation


  subroutine refraction

    ! 1st order upstream finite difference advection in
    ! directional space -- bottom- and current-induced refraction.

    integer :: i, o, p
    logical :: compute_rotation_tendency
    real :: sendbuffer
    real, save :: dtr_temp

    real :: flux(om,pm,istart:iend)

#ifdef ESMF
    ! always compute in coupled mode:
    compute_rotation_tendency = .true.
#else
    ! compute if varrying currents or first step:
    compute_rotation_tendency = currents .or. first
#endif

    if (compute_rotation_tendency) then

      ! compute rotation
      flux = 0
      do concurrent(i = istart:iend)
        do concurrent(o = 1:oc(i), p = 1:pm)
          flux(o,p,i) = 0.5 * (((cp0(o,ie(i)) - cp0(o,iw(i))) * sth(p) &
                               + vc(iie(i)) - vc(iiw(i))) * oneovdx(i) &
                             - ((cp0(o,in(i)) - cp0(o,is(i))) * cth(p) &
                               + uc(iin(i)) - uc(iis(i))) * oneovdy(i))
        end do
      end do

      ! evaluate rotation at cell edges:
      rotl = 0
      rotr = 0
      do concurrent(i = istart:iend)
        do concurrent(o = 1:oc(i), p = 1:pm)
          rotl(o,p,i) = 0.5 * (flux(o,p,i) + flux(o,pl(p),i))
          rotr(o,p,i) = 0.5 * (flux(o,p,i) + flux(o,pr(p),i))
        end do
      end do

      if (maxval(flux) <= tiny(sendbuffer)) then
        dtr_temp = dts
      else
#ifdef MPI
        sendbuffer = dth / maxval(flux)
        call mpi_allreduce(sendbuffer, dtr_temp, 1, MPI_REAL, mpi_min, MPI_COMM_WORLD, ierr)
#else
        dtr_temp = dth / maxval(flux)
#endif
      end if

    end if

    dtr = min(dtr_temp, dts)

    flux = 0
    do concurrent(i = istart:iend)
      do concurrent(o = 1:oc(i), p = 1:pm)

        ! compute tendencies
        flux(o,p,i) = 0.5 * ((rotl(o,p,i) + abs(rotl(o,p,i))) * e(o,p,i)    &
                           + (rotl(o,p,i) - abs(rotl(o,p,i))) * e(o,pl(p),i)&
                           - (rotr(o,p,i) + abs(rotr(o,p,i))) * e(o,pr(p),i)&
                           - (rotr(o,p,i) - abs(rotr(o,p,i))) * e(o,p,i))   &
                     * oneovdth
        ! integrate
        ef(o,p,i) = ef(o,p,i) - dtr * flux(o,p,i)
        
        if (fice(i) > fice_uth) then
          ef(o,p,i) = 0.0
        end if

      end do
    end do

  end subroutine refraction

end module umwm_advection
