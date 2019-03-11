module umwm_physics

  use umwm_module

  implicit none

contains

subroutine source

  ! TODO move to umwm_integration.f90

#ifdef MPI
use mpi
#endif

integer :: i,o,p
real    :: maxarg,send_buff

!real,save :: explim_ramp

!print *, maxval(sice)

! calculate the exponential argument:
ef = 0
do i = istart,iend
  do p = 1,pm
    do o = 1,oc(i)
      ef(o,p,i) = ssin(o,p,i)-sds(o,p,i)*snl_arg(o,i)&
                 -sbf(o,i)-sdt(o,i)-sdv(o,i)+sice(o,p,i)
    end do
  end do
end do

!explim_ramp = explim

! if first global step, apply ramp to explim
!if(firstdtg)then
!  explim_ramp = (0.01+0.99*sumt/dtg)*explim
!end if

! compute the dynamic time step and update in-step model time:
maxarg = maxval(abs(ef))
!dts = min(explim_ramp/maxarg,dtamin,dtg-sumt)
dts = min(explim/maxarg,dtamin,dtg-sumt)

#ifdef MPI
send_buff = dts
call mpi_allreduce(send_buff,dts,1,MPI_REAL,mpi_min,MPI_COMM_WORLD,ierr)
#endif

! if first time step, set time step to zero and
! integrate only the diagnostic part
if (.not. restart .and. first) dts = 0

dta = dts

! increment time:
sumt = sumt + dts

! integrate source terms for the prognostic range (o <= ol)
do i = istart,iend
  do p = 1,pm
    do o = 1,oc(i)
      ef(o,p,i) = e(o,p,i)*exp(dts*(ssin(o,p,i)-sds(o,p,i)      &
                                   -sbf(o,i)-sdt(o,i)-sdv(o,i)  &
                                   +sice(o,p,i)))&
                 +dts*snl(o,p,i)
    end do
  end do
end do

! integrate source terms for the diagnostic range (o > ol)
do i = istart, iend
  do p = 1, pm
    do o = oc(i) + 1, om
      if (ssin(o,p,i) - sdt(o,i) - sdv(o,i) + sice(o,p,i) >= 0) then
        ef(o,p,i) = oneoverk4(o,i) * ((ssin(o,p,i) - sdt(o,i) - sdv(o,i) + sice(o,p,i)) &
                  / (twopisds_fac * f(o) * dummy(o,p,i) * cothkd(o,i)))**inv_sds_power
      end if
    end do
  end do
end do

e(:,:,istart:iend) = 0.5*(e(:,:,istart:iend)+ef(:,:,istart:iend))

endsubroutine source


subroutine diag

! TODO move to umwm_diagnostics.F90

integer              :: o,p,i
integer              :: opeak,ppeak
integer,dimension(2) :: spectrum_peak_loc

real                        :: mag,xcomp,ycomp
real,dimension(istart:iend) :: m0,m2
real,dimension(om,pm)       :: spectrumbin

real :: ekdkovcp

m0 = 0
m2 = 0

do i=istart,iend
  do p=1,pm
    do o=1,om
      m0(i) = m0(i)+e(o,p,i)*kdk(o,i)
      m2(i) = m2(i)+f(o)**2*e(o,p,i)*kdk(o,i)
    end do
  end do
end do

mwp = sqrt(m0 / (m2 + tiny(m2))) ! mean wave period

momx = 0
momy = 0
cgmxx = 0
cgmxy = 0
cgmyy = 0

do i=istart,iend

  ! total wave momentum:
  do p=1,pm
    do o=1,oc(i)

      ekdkovcp = e(o,p,i)*kdk(o,i)*invcp0(o,i)

      momx(i) = momx(i)+ekdkovcp*cth(p)
      momy(i) = momy(i)+ekdkovcp*sth(p)

      cgmxx(i) = cgmxx(i)+cg0(o,i)*ekdkovcp*cth(p)**2
      cgmxy(i) = cgmxy(i)+cg0(o,i)*ekdkovcp*cth(p)*sth(p)
      cgmyy(i) = cgmyy(i)+cg0(o,i)*ekdkovcp*sth(p)**2

    end do
  end do

  momx(i)  = momx(i)*rhow(i)*dthg
  momy(i)  = momy(i)*rhow(i)*dthg
  cgmxx(i) = cgmxx(i)*rhow(i)*dthg
  cgmxy(i) = cgmxy(i)*rhow(i)*dthg
  cgmyy(i) = cgmyy(i)*rhow(i)*dthg

  ! significant wave height:
  ht(i) = 0.
  do p=1,pm
    do o=1,om
      ht(i) = ht(i)+e(o,p,i)*kdk(o,i)
    end do
  end do
  ht(i) = 4*sqrt(ht(i)*dth)

  ! mean wave direction:
  xcomp = 0.
  ycomp = 0.
  do p=1,pm
    mag   = sum(e(:,p,i)*kdk(:,i),dim=1)
    xcomp = xcomp+mag*cth(p)
    ycomp = ycomp+mag*sth(p)
  end do
  mwd(i) = atan2(ycomp,xcomp)

  ! wavenumber spectrum moments:
  m0(i) = 0.
  m2(i) = 0.
  do p=1,pm
    do o=1,om
      m0(i) = m0(i)+e(o,p,i)*kdk(o,i)
      m2(i) = m2(i)+e(o,p,i)*k3dk(o,i)
    end do
  end do

  ! mean-squared slope:
  mss(i) = m2(i) * dth

  mwl(i) = twopi * sqrt(m0(i) / (m2(i) + tiny(m2(i)))) ! mean wavelenght

  do p=1,pm
    do o=1,om
      spectrumbin(o,p) = e(o,p,i)*kdk(o,i)
    end do
  end do

  spectrum_peak_loc = maxloc(spectrumbin)  ! indices of spectrum peak

  opeak = spectrum_peak_loc(1)             ! frequency/wavenumber peak
  ppeak = spectrum_peak_loc(2)             ! direction peak

  dwd(i) = th(ppeak)                       ! dominant wave direction
  dwl(i) = twopi/k(opeak,i)                ! dominant wave length
  dwp(i) = 1./f(opeak)                     ! dominant wave period

  dcp0(i) = cp0(opeak,i)                   ! dominant phase speed, intrinsic
  dcg0(i) = cg0(opeak,i)                   ! dominant group speed, intrinsic

  ! dominant phase and group speed:
  dcp(i) = dcp0(i)+uc(i)*cth(ppeak)+vc(i)*sth(ppeak)
  dcg(i) = dcg0(i)+uc(i)*cth(ppeak)+vc(i)*sth(ppeak)

end do

end subroutine diag

end module umwm_physics
