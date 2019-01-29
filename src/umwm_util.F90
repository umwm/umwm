module umwm_util
!======================================================================!
!                                                                      !
! description: a module with utility functions.                        !
!                                                                      !
! contains: remap_i2mn     - remaps an (i) indexed array to (m,n)      !
!           remap_mn2i     - remaps an (m,n) indexed array to (i)      !
!           sigwaveheight  - calculates sig. wave height               !
!           meanwaveperiod - calculates mean wave period               !
!           raiseexception - raises an exception and prints a message  !
!           dealloc        - array deallocation routine                !
!                                                                      !
!======================================================================!
implicit none
!======================================================================!
contains



pure function remap_i2mn(field_i) result(field_mn)
!======================================================================+
!                                                                      !
! remaps an (i) indexed array to (m,n)                                 !
!                                                                      !
!======================================================================+
use umwm_module,only:imm,mm,nm,ii

! arguments
real,dimension(imm),intent(in) :: field_i
real,dimension(mm,nm)          :: field_mn

integer :: m,n
!=======================================================================

do n=1,nm
  do m=1,mm
    field_mn(m,n) = field_i(ii(m,n))
  end do
end do

end function remap_i2mn
!======================================================================>



pure function remap_mn2i(field_mn) result(field_i)
!======================================================================+
!                                                                      !
! remaps an (m,n) indexed array to (i)                                 !
!                                                                      !
!======================================================================+
use umwm_module,only:imm,mm,nm,ii

real,dimension(mm,nm),intent(in) :: field_mn
real,dimension(imm)              :: field_i
integer :: m,n

do n=1,nm
  do m=1,mm
    field_i(ii(m,n)) = field_mn(m,n)
  end do
end do

end function remap_mn2i
!======================================================================>



pure function sigwaveheight(i) result(swh)
!======================================================================+
!                                                                      !
! given a spatial grid index i, returns significant wave height        !
! at that location.                                                    !
!                                                                      !
!======================================================================+
use umwm_module,only:e,kdk,dth,om,pm

! arguments:
integer,intent(in) :: i

integer :: o,p
real    :: swh

!=======================================================================

swh = 0
do p=1,pm
  do o=1,om
    swh = swh+e(o,p,i)*kdk(o,i)
  end do
end do
swh = 4.*sqrt(swh*dth)

end function sigwaveheight
!======================================================================+



pure function meanwaveperiod(i) result(mwp)
!======================================================================+
!                                                                      !
! given a spatial grid index i, returns mean wave period at that       !
! location.                                                            !
!                                                                      !
!======================================================================+
use umwm_module,only:e,f,kdk,om,pm

! arguments:
integer,intent(in) :: i

integer :: o,p
real    :: m0,m2,mwp

m0 = 0
m2 = 0
do p=1,pm
  do o=1,om
    m0 = m0+e(o,p,i)*kdk(o,i)
    m2 = m2+f(o)**2*e(o,p,i)*kdk(o,i)
  end do
end do
mwp = sqrt(m0/m2)

end function meanwaveperiod
!======================================================================+



subroutine raiseexception(exceptiontype,routinename,message,flag)
!======================================================================+
!                                                                      !
! raises an exception, prints a message to stdout and sets the flag    !
! to .false. if present in argument list                               !
!                                                                      !
!======================================================================+

! arguments:
character(len=*),intent(in)    :: exceptiontype
character(len=*),intent(in)    :: routinename
character(len=*),intent(in)    :: message
logical,intent(inout),optional :: flag

!======================================================================+

if(present(flag))then
  flag = .false.
end if

write(unit=*,fmt='(a)')'umwm: '//trim(routinename)//': '  &
                               //trim(exceptiontype)//': '&
                               //trim(message)

endsubroutine raiseexception
!======================================================================+



subroutine dealloc
!======================================================================+
!                                                                      !
! deallocates umwm arrays                                              !
!                                                                      !
!======================================================================>
use umwm_module
!======================================================================>

deallocate(ar_2d,d_2d,dlon,dlat,dx_2d,dy_2d)
deallocate(curv)
deallocate(gustu,gustv)
deallocate(lat,lon)
deallocate(x,y)
deallocate(rhoa_2d,rhow_2d)
deallocate(wspd_2d,wdir_2d)
deallocate(uwb,vwb,uw,vw,uwf,vwf,ucb,uc_2d,ucf,vcb,vc_2d,vcf)
deallocate(dom,f,cth,cth2,sth,th)
deallocate(cth_curv,sth_curv)
deallocate(ar,cd,d,dx,dy,dwd,dwl,dwp,fcutoff)
deallocate(dxn,dxs,dye,dyw)
deallocate(dcp0,dcg0,dcp,dcg)
deallocate(ht,mss,mwd,mwl,mwp)
deallocate(momx,momy)
deallocate(cgmxx,cgmxy,cgmyy)
deallocate(oneovar,oneovdx,oneovdy)
deallocate(psim,psiml2)
deallocate(rhoab,rhoa,rhoaf,rhowb,rhow,rhowf,rhorat)
deallocate(taux,tauy,taux_form,tauy_form,taux_skin,tauy_skin)
deallocate(taux_ocntop,tauy_ocntop,taux_ocnbot,tauy_ocnbot)
deallocate(taux_diag,tauy_diag)
deallocate(taux_snl,tauy_snl)
deallocate(taux1,tauy1,taux2,tauy2,taux3,tauy3)
deallocate(tailatmx,tailatmy)
deallocate(tailocnx,tailocny)
deallocate(epsx_atm, epsy_atm, epsx_ocn, epsy_ocn)
deallocate(uc,vc,ustar)
deallocate(wspd,wdir)
deallocate(shelt)
deallocate(bf1_renorm,bf2_renorm)
deallocate(cg0,cp0,cothkd)
deallocate(dwn,invcp0)
deallocate(fkovg)
deallocate(k,k4,kdk,k3dk,l2,logl2overz,oneoverk4)
deallocate(sbf,sdv,sdt,snl_arg,dummy,e,ef,rotl,rotr,sds,snl,ssin)

endsubroutine dealloc
!======================================================================!
end module umwm_util
