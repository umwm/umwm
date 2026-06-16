module umwm_init

#ifdef MPI
use mpi
#endif
use umwm_config, only: config_type
use umwm_constants, only: rk
use umwm_dispersion, only: group_speed, wavenumber
use umwm_grid, only: grid_type
use umwm_spectrum, only: spectrum_type
use umwm_module

implicit none

contains


subroutine alloc(grid, spectrum)
  ! Allocates UMWM arrays
  type(grid_type), intent(in) :: grid
  type(spectrum_type), intent(in) :: spectrum

  associate( &
    om => spectrum % num_frequencies, &
    pm => spectrum % num_directions &
  )

  ! 1-d arrays:
  allocate(dom(om),f(om))
  allocate(cth(pm),cth2(pm),pl(pm),pr(pm),sth(pm),th(pm))

  allocate(oc(grid % istart:grid % iend))

  ! 2-d arrays (remapped):
  allocate(cd(grid % istart:grid % iend))        ! air-side drag coefficient
  allocate(fcutoff(grid % istart:grid % iend))   ! cutoff frequency
  allocate(ht(grid % istart:grid % iend))        ! significant wave height
  allocate(mss(grid % istart:grid % iend))       ! mean-squared slope

  allocate(shelt(grid % istart:grid % iend))     ! sheltering coefficient
  shelt = 0

  allocate(physics_time_step(grid % istart:grid % iend))
  physics_time_step = 0

  ! mean spectrum quantities:
  allocate(mwd(grid % istart:grid % iend)) ! direction
  allocate(mwp(grid % istart:grid % iend)) ! period
  allocate(mwl(grid % istart:grid % iend)) ! wavelength

  ! dominant spectrum quantities:
  allocate(dwd(grid % istart:grid % iend))  ! direction
  allocate(dwp(grid % istart:grid % iend))  ! period
  allocate(dwl(grid % istart:grid % iend))  ! wavelength
  allocate(dcp0(grid % istart:grid % iend)) ! intrinsic phase speed
  allocate(dcp(grid % istart:grid % iend))  ! phase speed
  allocate(dcg0(grid % istart:grid % iend)) ! intrinsic group speed
  allocate(dcg(grid % istart:grid % iend))  ! group speed

  dwd  = 0
  dwp  = 0
  dwl  = 0
  dcp0 = 0
  dcg0 = 0
  dcp  = 0
  dcg  = 0

  allocate(momx(grid % istart:grid % iend), momy(grid % istart:grid % iend)) ! total wave momentum
  allocate(cgmxx(grid % istart:grid % iend), cgmxy(grid % istart:grid % iend), &
           cgmyy(grid % istart:grid % iend)) ! cg*m

  ! momentum fluxes:
  allocate(taux(grid % istart:grid % iend), tauy(grid % istart:grid % iend))
  taux = 0; tauy = 0

  allocate(taux_form(grid % istart:grid % iend), tauy_form(grid % istart:grid % iend))
  taux_form = 0; tauy_form = 0

  allocate(taux_skin(grid % istart:grid % iend), tauy_skin(grid % istart:grid % iend))
  taux_skin = 0; tauy_skin = 0

  allocate(taux_diag(grid % istart:grid % iend), tauy_diag(grid % istart:grid % iend))
  taux_diag = 0; tauy_diag = 0

  allocate(taux_ocntop(grid % istart:grid % iend), tauy_ocntop(grid % istart:grid % iend))
  taux_ocntop = 0; tauy_ocntop = 0

  allocate(taux_ocnbot(grid % istart:grid % iend), tauy_ocnbot(grid % istart:grid % iend))
  taux_ocnbot = 0; tauy_ocnbot = 0

  allocate(taux_snl(grid % istart:grid % iend), tauy_snl(grid % istart:grid % iend))
  taux_snl = 0; tauy_snl = 0

  allocate(epsx_atm(grid % istart:grid % iend), epsy_atm(grid % istart:grid % iend))
  epsx_atm = 0; epsy_atm = 0

  allocate(epsx_ocn(grid % istart:grid % iend), epsy_ocn(grid % istart:grid % iend))
  epsx_ocn = 0; epsy_ocn = 0

  allocate(taux1(grid % istart:grid % iend), tauy1(grid % istart:grid % iend))
  allocate(taux2(grid % istart:grid % iend), tauy2(grid % istart:grid % iend))
  allocate(taux3(grid % istart:grid % iend), tauy3(grid % istart:grid % iend))

  taux1 = 0; tauy1 = 0
  taux2 = 0; tauy2 = 0
  taux3 = 0; tauy3 = 0

  allocate(tailatmx(grid % istart:grid % iend), tailatmy(grid % istart:grid % iend))
  allocate(tailocnx(grid % istart:grid % iend), tailocny(grid % istart:grid % iend))

  allocate(ustar(grid % istart:grid % iend)) ! air-side friction velocity

  allocate(psim(grid % imm)) ! integrated stability function for momentum
  psim = 0

  ! 3-d arrays (remapped):
  allocate(bf1_renorm(om,grid % istart:grid % iend))
  allocate(bf2_renorm(om,grid % istart:grid % iend))
  allocate(    cothkd(om,grid % istart:grid % iend))
  allocate(       dwn(om,grid % istart:grid % iend))
  allocate(     fkovg(om,grid % istart:grid % iend))
  allocate(    invcp0(om,grid % istart:grid % iend))
  allocate(        k4(om,grid % istart:grid % iend))
  allocate(       kdk(om,grid % istart:grid % iend))
  allocate(      k3dk(om,grid % istart:grid % iend))
  allocate( oneoverk4(om,grid % istart:grid % iend))
  allocate(       sbf(om,grid % istart:grid % iend))
  allocate(       sdt(om,grid % istart:grid % iend))
  allocate(       sdv(om,grid % istart:grid % iend))
  allocate(   snl_arg(om,grid % istart:grid % iend))
  allocate(        l2(om,grid % istart:grid % iend))
  allocate(logl2overz(om,grid % istart:grid % iend))

  allocate(psiml2(om,grid % istart:grid % iend))
  psiml2 = 0

  ! 4-d arrays (remapped):
  allocate(   ef(om,pm,grid % istart:grid % iend)) ! wave variance spectrum, forward in time
  allocate(dummy(om,pm,grid % istart:grid % iend)) ! dummy array, used in sds, advection and refraction
  allocate( rotl(om,pm,grid % istart:grid % iend)) ! anti-clockwise rotation, used in refraction
  allocate( rotr(om,pm,grid % istart:grid % iend)) ! clockwise rotation, used in refraction
  allocate(  sds(om,pm,grid % istart:grid % iend)) ! wave dissipation sink function
  allocate(  snl(om,pm,grid % istart:grid % iend)) ! wave downshifting source/sink function
  allocate( ssin(om,pm,grid % istart:grid % iend)) ! wind input source/sink function

  allocate( sice(om,grid % istart:grid % iend)) ! wave attenuation by sea ice function

  allocate(e(om,pm,grid % iistart-1:grid % iiend))
  allocate(cp0(om,grid % iistart-1:grid % iiend))
  allocate(cg0(om,grid % iistart-1:grid % iiend))
  allocate(k(om,grid % istart:grid % iend))

  e = tiny(e)

  end associate

end subroutine alloc


subroutine init(config, spectrum, grid, forcing)
! Initialize model variables such as frequencies, direction angles,
! phase speed and group velocity, wave numbers, etc.
use umwm_config, only: config_type
use umwm_forcing, only: forcing_type
use umwm_util,only:raiseexception

type(config_type), intent(in) :: config
type(spectrum_type), intent(in) :: spectrum
type(grid_type), intent(in) :: grid
type(forcing_type), intent(inout) :: forcing
integer :: i, o, p, pp, ind
real :: mindelx

#ifdef MPI
integer :: n
#endif

! initialize legacy spectrum aliases:
om = spectrum % num_frequencies
pm = spectrum % num_directions
dlnf = real(spectrum % dlnf, kind(dlnf))
dth = real(spectrum % dth, kind(dth))

! set frequency bins:
f = real(spectrum % frequency, kind(f))
th = real(spectrum % direction, kind(th))

! define various constants:
dthg          = dth * config % g
oneovdth      = 1./dth
twopisds_fac  = twopi * config % sds_fac
fieldscale1   = config % sin_diss1 / config % sin_fac
fieldscale2   = config % sin_diss2 / config % sin_diss1
inv_sds_power = 1. / config % sds_power

! this limits the Courant number to its theoretical value slightly
! larger than 1/sqrt(2), depending on the number of directional bins;
! this also limits the number of directions to be divisible by 8,
! not 4! (more isotropic in Cartesian projection)

if(mod(pm,8) == 0)then
  cfllim = cos(0.25*pi-0.5*dth)
else
  cfllim = 1./sqrt(2.)
end if

! compute diffusion values in 2 frequenciess:
bf1  = exp(-16*dlnf*dlnf)
bf2  = exp(-64*dlnf*dlnf)
bf1a = bf1/(bf1+bf2)
bf2  = bf2/(bf1+bf2)
bf1  = bf1a

cth = cos(th) ! cosines
sth = sin(th) ! sines

! "left" and "right" directional indices for refraction:
do p=1,pm
  pl(p) = p+1
  pr(p) = p-1
end do
pl(pm) = 1
pr(1)  = pm

dom = twopi*dlnf*f

do p=1,pm
  cth2(p) = cos(dth*(p-1))**2
end do

allocate(cth2pp(pm,pm))

do p=1,pm
  do pp=1,pm
    ind = pp-p+1
    if(ind<=0)ind = pm+ind
    cth2pp(pp,p) = cth2(ind)*dth
  end do
end do

! compute wave numbers, phase speeds, and group velocities:
call dispersion(config, grid)

mindelx = min(minval(grid % dx_2d, grid % mask == 1), minval(grid % dy_2d, grid % mask == 1))
cgmax   = maxval(cg0(:,grid % istart:grid % iend))
dtamin  = 0.98*cfllim*mindelx/cgmax

first    = .true.

if(config % restart)then
  firstdtg = .false.
else
  firstdtg = .true.
end if

! if sea ice from file, update the fice field
if (config % seaice) forcing % fice = grid % remap_mn2i(forcing % ficef)

! if forcing from file, update the wspd field for ustar first guess
if (config % winds) forcing % wspd = grid % remap_mn2i(sqrt(forcing % uwf**2 + forcing % vwf**2))
call forcing % apply_wind_speed_floor()

! initialize drag coefficient (Large and Pond, 1981):
cd = 1.2e-3
do concurrent (i=grid % istart:grid % iend, forcing % wspd(i) > 11)
  cd(i) = (0.49 + 0.065 * forcing % wspd(i)) * 1e-3
end do

! initialize friction velocity:
do i=grid % istart,grid % iend
  ustar(i) = sqrt(cd(i))*forcing % wspd(i)
end do

#ifdef MPI

! figure out which process will print to screen:
iip = grid % ii(config % xpl, config % ypl)
if(grid % mask(config % xpl, config % ypl)==0)then
  if(nproc==0)then
    write(0,*)config % xpl, config % ypl, grid % mask(config % xpl, config % ypl)

    call raiseexception('warning','init',&
                        'land-point chosen for stdout, may go out of bounds')
    call raiseexception('warning','init',&
                        'check xpl and ypl in the output namelist in namelists/main.nml')
    stop
  end if
end if

if(nproc==0)then
  do n=0,mpisize-1
    if(iip>=grid % istart_all(n).and.iip<=grid % iend_all(n))nproc_plot = n
  end do
end if

call mpi_bcast(nproc_plot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(iip,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

#else

iip = grid % ii(config % xpl, config % ypl)

#endif

#ifndef MPI

write(*,fmt=102)
write(*,*)'initialization summary:'
write(*,fmt=103)
write(*,*)'bin  f[hz]    t[s]    min(k)   max(k)   min(c)   max(c)   min(cg)   max(cg)'
write(*,fmt=103)
do o=1,om
  write(*,fmt=101)o,f(o),1./f(o),                                     &
                  minval(k(o,:),dim=1),maxval(k(o,:),dim=1),          &
                  minval(cp0(o,1:grid % im),dim=1),maxval(cp0(o,1:grid % im),dim=1),&
                  minval(cg0(o,1:grid % im),dim=1),maxval(cg0(o,1:grid % im),dim=1)
end do
101 format(1x,i2,8(2x,f7.4))
write(*,fmt=102)

102 format('!',77('='),'!')
103 format('!',77('-'),'!')

#endif

end subroutine init


subroutine dispersion(config, grid)
! Iteratively solve the dispersion relation by iteration,
! and compute phase and group velocities in absolute reference frame (w/ currents)
use umwm_config, only: config_type

#ifdef MPI
integer :: sendcount, recvcount
integer :: sendtag, recvtag
integer :: src, dest
integer :: status(MPI_STATUS_SIZE)
#endif

type(config_type), intent(in) :: config
type(grid_type), intent(in) :: grid
integer :: i,o

real,dimension(om,grid % istart:grid % iend) :: kd

if(nproc==0)write(*,'(a)')'umwm: dispersion: solving for dispersion relationship;'

do concurrent (o=1:om, i=grid % istart:grid % iend)
  k(o,i) = real(wavenumber(real(f(o), rk), real(grid % d(i), rk), real(config % rhow0, rk), &
                            real(config % g, rk), real(config % sfct, rk)), kind(k))
end do

if(nproc==0)write(*,'(a)')'umwm: dispersion: dispersion relationship done;'

do concurrent (o=1:om, i=grid % istart:grid % iend)
  kd(o,i) = k(o,i) * grid % d(i)
end do

! limit kd to avoid floating overflow in transcendental functions:
where(kd>20.)kd = 20.

! phase speed and group velocity:
cp0 = tiny(cp0)
cg0 = tiny(cg0)

do concurrent (o=1:om, i=grid % istart:grid % iend)

  cp0(o,i) = twopi*f(o)/k(o,i)
  cg0(o,i) = real(group_speed(real(k(o,i), rk), real(grid % d(i), rk), &
                              real(config % rhow0, rk), real(config % g, rk), &
                              real(config % sfct, rk)), kind(cg0))
end do

! compute some frequently used arrays:
do concurrent (o=1:om, i=grid % istart:grid % iend)

  dwn(o,i)       = dom(o)/abs(cg0(o,i))                         ! dk
  l2(o,i)        = 0.5*abs(cp0(o,i))/f(o)                       ! lambda/2 (half wavelength)
  k4(o,i)        = k(o,i)**4.                                   ! k^4
  oneoverk4(o,i) = 1./k4(o,i)                                   ! k^-4
  kdk(o,i)       = k(o,i)*dwn(o,i)                              ! k*dk
  k3dk(o,i)      = k(o,i)**3.*dwn(o,i)                          ! k*k*k*dk
  fkovg(o,i)     = f(o)*k(o,i)/config % g                       ! f*k/g
  cothkd(o,i)    = cosh(0.2*kd(o,i))/sinh(0.2*kd(o,i))          ! coth(0.2*kd)
  invcp0(o,i)    = 1./cp0(o,i)                                  ! 1/cp
  sbf(o,i)       = config % sbf_fac*k(o,i)/(sinh(2.*kd(o,i)))&  ! bottom friction
                  +config % sbp_fac*k(o,i)/(cosh(kd(o,i))*cosh(kd(o,i))) ! bottom percolation
  sdv(o,i)       = 4.*config % nu_water*k(o,i)**2.              ! viscosity

end do

! compute renormalization factors for snl:
bf1_renorm = 0.
bf2_renorm = 0.
snl_arg    = 0.

do i=grid % istart,grid % iend
  do o=1,om-2
    bf1_renorm(o,i) = config % snl_fac*bf1*kdk(o+1,i)/kdk(o,i)
    bf2_renorm(o,i) = config % snl_fac*bf2*kdk(o+2,i)/kdk(o,i)
    snl_arg(o,i)    = 1.-(bf1_renorm(o,i)+bf2_renorm(o,i))
  end do
end do

! half-wavelength over z:
logl2overz = log(l2/config % z)

! limit wind input to be at 10 m for l/2 > 10 m:
where(l2>20.)logl2overz = log(20./config % z)

#ifdef MPI
if(nproc<mpisize-1)then ! communicate with process above:

  sendcount = om*(grid % iend-grid % iistart_all(nproc+1)+1) ; dest = nproc+1 ; sendtag = nproc
  recvcount = om*(grid % iiend-grid % iend)                  ; src  = nproc+1 ; recvtag = src

  call mpi_sendrecv(cp0(:,grid % iistart_all(nproc+1):grid % iend),sendcount,&
                    MPI_REAL,dest,sendtag,                  &
                    cp0(:,grid % iend+1:grid % iiend),recvcount,          &
                    MPI_REAL,src,recvtag,                   &
                    MPI_COMM_WORLD,status,ierr)

end if

if(nproc>0)then ! communicate with process below:

  sendcount = om*(grid % iiend_all(nproc-1)-grid % istart+1) ; dest = nproc-1 ; sendtag = nproc
  recvcount = om*(grid % istart-grid % iistart)              ; src  = nproc-1 ; recvtag = src

  call mpi_sendrecv(cp0(:,grid % istart:grid % iiend_all(nproc-1)),sendcount,&
                    MPI_REAL,dest,sendtag,                  &
                    cp0(:,grid % iistart:grid % istart-1),recvcount,      &
                    MPI_REAL,src,recvtag,                   &
                    MPI_COMM_WORLD,status,ierr)

end if

call mpi_barrier(MPI_COMM_WORLD,ierr)

if(nproc<mpisize-1)then ! communicate with process above:

  sendcount = om*(grid % iend-grid % iistart_all(nproc+1)+1) ; dest = nproc+1 ; sendtag = nproc
  recvcount = om*(grid % iiend-grid % iend)                  ; src  = nproc+1 ; recvtag = src

  call mpi_sendrecv(cg0(:,grid % iistart_all(nproc+1):grid % iend),sendcount,&
                    MPI_REAL,dest,sendtag,                  &
                    cg0(:,grid % iend+1:grid % iiend),recvcount,          &
                    MPI_REAL,src,recvtag,                   &
                    MPI_COMM_WORLD,status,ierr)

end if

if(nproc>0)then ! communicate with process below:

  sendcount = om*(grid % iiend_all(nproc-1)-grid % istart+1) ; dest = nproc-1 ; sendtag = nproc
  recvcount = om*(grid % istart-grid % iistart)              ; src  = nproc-1 ; recvtag = src

  call mpi_sendrecv(cg0(:,grid % istart:grid % iiend_all(nproc-1)),sendcount,&
                    MPI_REAL,dest,sendtag,                  &
                    cg0(:,grid % iistart:grid % istart-1),recvcount,      &
                    MPI_REAL,src,recvtag,                   &
                    MPI_COMM_WORLD,status,ierr)

end if

call mpi_barrier(MPI_COMM_WORLD,ierr)

! if periodic domain, connect the east and west:
if(config % isglobal)then

  if(nproc==0)then ! communicate with last tile:

    sendcount = om*grid % first_col_len ; dest = mpisize-1 ; sendtag = nproc
    recvcount = om*grid % last_col_len  ; src  = mpisize-1 ; recvtag = src

    call mpi_sendrecv(cp0(:,grid % istart:(grid % istart+grid % first_col_len-1)),sendcount,&
                      MPI_REAL,dest,sendtag,                           &
                      cp0(:,grid % iistart:grid % istart-1),recvcount,  &
                      MPI_REAL,src,recvtag,                            &
                      MPI_COMM_WORLD,status,ierr)

  end if

  if(nproc==mpisize-1)then ! communicate with first tile:

    sendcount = om*grid % last_col_len  ; dest = 0 ; sendtag = nproc
    recvcount = om*grid % first_col_len ; src  = 0 ; recvtag = src

    call mpi_sendrecv(cp0(:,(grid % iend-grid % last_col_len+1):grid % iend),sendcount,&
                      MPI_REAL,dest,sendtag,                      &
                      cp0(:,grid % iend+1:grid % iiend),recvcount,&
                      MPI_REAL,src,recvtag,                       &
                      MPI_COMM_WORLD,status,ierr)

  end if

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if(nproc==0)then ! communicate with last tile:

    sendcount = om*grid % first_col_len ; dest = mpisize-1 ; sendtag = nproc
    recvcount = om*grid % last_col_len  ; src  = mpisize-1 ; recvtag = src

    call mpi_sendrecv(cg0(:,grid % istart:(grid % istart+grid % first_col_len-1)),sendcount,&
                      MPI_REAL,dest,sendtag,                           &
                      cg0(:,grid % iistart:grid % istart-1),recvcount,  &
                      MPI_REAL,src,recvtag,                            &
                      MPI_COMM_WORLD,status,ierr)

  end if

  if(nproc==mpisize-1)then ! communicate with first tile:

    sendcount = om*grid % last_col_len  ; dest = 0 ; sendtag = nproc
    recvcount = om*grid % first_col_len ; src  = 0 ; recvtag = src

    call mpi_sendrecv(cg0(:,(grid % iend-grid % last_col_len+1):grid % iend),sendcount,&
                      MPI_REAL,dest,sendtag,                      &
                      cg0(:,grid % iend+1:grid % iiend),recvcount,&
                      MPI_REAL,src,recvtag,                       &
                      MPI_COMM_WORLD,status,ierr)

  end if

end if
#endif

! handle land points for cp and cg (needed for advection/refraction)
do o=1,om
  cp0(o,grid % iistart-1) = minval(cp0(o,grid % istart:grid % iend))
  cg0(o,grid % iistart-1) = minval(cg0(o,grid % istart:grid % iend))
end do

end subroutine dispersion

end module umwm_init
