module umwm_util
!======================================================================!
!                                                                      !
! description: a module with utility functions.                        !
!                                                                      !
! contains: sigwaveheight  - calculates sig. wave height               !
!           meanwaveperiod - calculates mean wave period               !
!           raiseexception - raises an exception and prints a message  !
!           dealloc        - array deallocation routine                !
!                                                                      !
!======================================================================!
implicit none
!======================================================================!
contains


pure function sigwaveheight(i, spectrum) result(swh)
!======================================================================+
!                                                                      !
! given a spatial grid index i, returns significant wave height        !
! at that location.                                                    !
!                                                                      !
!======================================================================+
use umwm_module,only:e,kdk,dth
use umwm_spectrum, only: spectrum_type

! arguments:
integer,intent(in) :: i
type(spectrum_type), intent(in) :: spectrum

integer :: o,p
real    :: swh

!=======================================================================

swh = 0
do p=1,spectrum % num_directions
  do o=1,spectrum % num_frequencies
    swh = swh+e(o,p,i)*kdk(o,i)
  end do
end do
swh = 4.*sqrt(swh*dth)

end function sigwaveheight
!======================================================================+



pure function meanwaveperiod(i, spectrum) result(mwp)
!======================================================================+
!                                                                      !
! given a spatial grid index i, returns mean wave period at that       !
! location.                                                            !
!                                                                      !
!======================================================================+
use umwm_module,only:e,f,kdk
use umwm_spectrum, only: spectrum_type

! arguments:
integer,intent(in) :: i
type(spectrum_type), intent(in) :: spectrum

integer :: o,p
real    :: m0,m2,mwp

m0 = 0
m2 = 0
do p=1,spectrum % num_directions
  do o=1,spectrum % num_frequencies
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

deallocate(dom,f,cth,cth2,sth,th)
deallocate(cd,dwd,dwl,dwp,fcutoff)
deallocate(dcp0,dcg0,dcp,dcg)
deallocate(ht,mss,mwd,mwl,mwp)
deallocate(momx,momy)
deallocate(cgmxx,cgmxy,cgmyy)
deallocate(psim,psiml2)
deallocate(taux,tauy,taux_form,tauy_form,taux_skin,tauy_skin)
deallocate(taux_ocntop,tauy_ocntop,taux_ocnbot,tauy_ocnbot)
deallocate(taux_diag,tauy_diag)
deallocate(taux_snl,tauy_snl)
deallocate(taux1,tauy1,taux2,tauy2,taux3,tauy3)
deallocate(tailatmx,tailatmy)
deallocate(tailocnx,tailocny)
deallocate(epsx_atm, epsy_atm, epsx_ocn, epsy_ocn)
deallocate(ustar)
deallocate(shelt)
deallocate(physics_time_step)
deallocate(bf1_renorm,bf2_renorm)
deallocate(cg0,cp0,cothkd)
deallocate(dwn,invcp0)
deallocate(fkovg)
deallocate(k,k4,kdk,k3dk,l2,logl2overz,oneoverk4)
deallocate(sbf,sdv,sdt,snl_arg,dummy,e,ef,rotl,rotr,sds,snl,ssin,sice)

endsubroutine dealloc
!======================================================================!



pure function distance_haversine(lon1, lon2, lat1, lat2) result(distance)
!======================================================================+
!                                                                      !
! calculates shortest distance between two points on the surface       !
! of a sphere using the Haversine formula                              !                                 !
!                                                                      !
!======================================================================+
    
! arguments
real :: distance
real, intent(in) :: lon1, lon2, lat1, lat2

! local
real :: dlon, dlat

dlon = abs(lon2 - lon1)
dlat = abs(lat2 - lat1)

distance = 2 * asin( sqrt( (sin(0.5 * dlat))**2 + &
                     cos(lat1) * cos(lat2) * (sin(0.5 * dlon))**2 ) )

end function distance_haversine
!======================================================================>

end module umwm_util
