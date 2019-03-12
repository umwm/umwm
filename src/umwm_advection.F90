module umwm_advection

  ! Subroutines for calculation of wave energy advection in geographical (x, y)
  ! and directional (theta) space. 1st order upstream finite differences.

#if defined(MPI)
  use mpi
#endif
  use umwm_module
  use umwm_io, only: currents

  implicit none

  logical :: zerocurrents

contains

  subroutine propagation

    ! 1st order upstream finite difference advection in geographical space.

    integer :: o, p, i
    real :: cge, cgw, cgn, cgs
    real :: feup, fedn, fwup, fwdn, fnup, fndn, fsup, fsdn
    real :: east_in_tx,west_in_tx,south_in_ty, north_in_ty
    real :: flux(om,pm,istart:iend)

    flux = 0

    do concurrent(i = istart:iend)

      do concurrent(o = 1:oc(i), p = 1:pm)


        !  double group velocity at east, west, north, south cell edges
        cge = (cg0(o,i) * cth_curv(p,i) + cg0(o,ie(i)) * cth_curv(p,i))
        cgw = (cg0(o,i) * cth_curv(p,i) + cg0(o,iw(i)) * cth_curv(p,i))
        cgs = (cg0(o,i) * sth_curv(p,i) + cg0(o,is(i)) * sth_curv(p,i))
        cgn = (cg0(o,i) * sth_curv(p,i) + cg0(o,in(i)) * sth_curv(p,i))

        
        east_in_tx    = 2*tx(i)/(1+tx(i))
        west_in_tx    = 0.5*(1+tx(i))
        south_in_ty   = 0.5*(1+ty(i))
        north_in_ty   = 2*ty(i)/(1+ty(i)) 
        
        if(tx(iw(i)) > 0.02 .and. tx(i) > 0.02) then
!
           west_in_tx  =  tx(iw(i)) * (1 + tx(i))/ (1+tx(iw(i)))
 !           west_in_tx = (tx(iw(i)) + tx(i))/ 2

        endif


        if(tx(ie(i)) > 0.02 .and. tx(i) > 0.02) then

           east_in_tx  =  tx(i) * (1 + tx(ie(i)))/ (1+tx(i))
           ! east_in_tx = (tx(ie(i)) + tx(i))/ 2
        endif


        if(ty(is(i)) > 0.02 .and. ty(i) > 0.02) then

           south_in_ty  =  ty(is(i)) * (1 + ty(i))/ (1+ty(is(i)))
           ! south_in_ty = (ty(is(i)) + ty(i))/ 2
        endif



       if(ty(in(i)) > 0.02 .and. ty(i) > 0.02) then
        
           north_in_ty  =  ty(i) * (1 + ty(in(i)))/ (1+ty(i))
           ! north_in_ty = (ty(in(i)) + ty(i))/ 2
       endif
        

        ! advective energy flux in and out of all 4 directions
        flux(o,p,i) = (cge + abs(cge)) * 1.0 * dye(i) * e(o,p,i)  & 
                      ! east, out
          + (cge - abs(cge)) * east_in_tx * dye(i) * e(o,p,ie(i)) & 
                      ! east, in
          - (cgw + abs(cgw)) * west_in_tx * dyw(i) * e(o,p,iw(i)) & 
                      ! west, in
          - (cgw - abs(cgw)) * 1.0 * dyw(i) * e(o,p,i)     & 
                      ! west, out
          + (cgn + abs(cgn)) * 1.0 * dxn(i) * e(o,p,i)     & 
                      ! north, out
          + (cgn - abs(cgn)) * north_in_ty * dxn(i) * e(o,p,in(i)) & 
                      ! north, in
          - (cgs + abs(cgs)) * south_in_ty * dxs(i) * e(o,p,is(i)) & 
                      ! south, in
          - (cgs - abs(cgs)) * 1.0 * dxs(i) * e(o,p,i)       
                      ! south, out
          

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
        feup = (tx(i) * uc(i) + tx(iie(i)) * uc(iie(i)) + abs(tx(i) * uc(i) + tx(iie(i)) * uc(iie(i)))) * dye(i)
        fedn = (tx(i) * uc(i) + tx(iie(i))* uc(iie(i)) - abs(tx(i) * uc(i) + tx(iie(i)) * uc(iie(i)))) * dye(i)
        fwup = (tx(i) * uc(i) + tx(iiw(i)) * uc(iiw(i)) + abs(tx(i) * uc(i) + tx(iiw(i)) * uc(iiw(i)))) * dyw(i)
        fwdn = (tx(i) * uc(i) + tx(iiw(i)) * uc(iiw(i)) - abs(tx(i) * uc(i) + tx(iiw(i)) * uc(iiw(i)))) * dyw(i)

        ! y-direction
        fnup = (ty(i) * vc(i) + ty(iin(i)) * vc(iin(i)) + abs(ty(i) * vc(i) + ty(iin(i)) * vc(iin(i)))) * dxn(i)
        fndn = (ty(i) * vc(i) + ty(iin(i)) * vc(iin(i)) - abs(ty(i) * vc(i) + ty(iin(i)) * vc(iin(i)))) * dxn(i)
        fsup = (ty(i) * vc(i) + ty(iis(i)) * vc(iis(i)) + abs(ty(i) * vc(i) + ty(iis(i)) * vc(iis(i)))) * dxs(i)
        fsdn = (ty(i) * vc(i) + ty(iis(i)) * vc(iis(i)) - abs(ty(i) * vc(i) + ty(iis(i)) * vc(iis(i)))) * dxs(i)

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


          if((tx(i) <= 0.02) .AND. (ty(i) <=0.02)) then
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
       
        if((tx(i) <= 0.02) .AND. (ty(i) <=0.02)) then
          ef(o,p,i) = 0.0
        end if




        if (fice(i) > fice_uth) then
          ef(o,p,i) = 0.0
        end if

      end do
    end do

  end subroutine refraction

end module umwm_advection
