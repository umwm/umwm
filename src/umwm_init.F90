module umwm_init

#ifdef MPI
use mpi
#endif
use umwm_module

implicit none

character :: remap_dir

contains



subroutine environment(option)

#ifdef MPI 
use umwm_mpi, only: istart_, iend_, ilen_, iistart_, iiend_
#endif

character(4), intent(in) :: option

if(option == 'init')then
#ifdef MPI
#ifndef ESMF
  call mpi_init(ierr)                           ! initialize mpi
#endif
  call mpi_comm_rank(MPI_COMM_WORLD,nproc,ierr) ! who am i?
  call mpi_comm_size(MPI_COMM_WORLD,mpisize,ierr)  ! how many processes?
#else
  nproc=0
  mpisize=1
#endif
end if

if(option == 'stop')then
#if defined(MPI) && !defined(ESMF)
  call mpi_barrier(MPI_COMM_WORLD,ierr) ! wait for all
  call mpi_finalize(ierr)               ! finalize mpi
#endif

#ifdef MPI
  deallocate(istart_, iend_, ilen_)
  deallocate(iistart_, iiend_)
#endif
end if

end subroutine environment


subroutine greeting
  use netcdf
  if (nproc == 0) then
    print '(a)'
    print '(a)', ' University of Miami Wave Model v' // version
    print '(a)', '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print '(a)'
    print '(a)', ' compiled using NetCDF library version ' // trim(nf90_inq_libvers())
    print '(a)'
  end if
end subroutine greeting


#ifndef GEOS
subroutine nmlread
! Opens namelist file namelists/main.nml and reads runtime input parameters.
use umwm_module,only:starttimestr_nml,stoptimestr_nml
use umwm_io,  only:winds,currents,air_density,water_density,seaice
use umwm_util,only:raiseexception

logical :: namelistok

! local variables for reading from namelist
character(19) :: starttimestr,stoptimestr

namelist /domain/ isglobal,mm,nm,om,pm,fmin,fmax,fprog,starttimestr,&
stoptimestr,dtg,restart

namelist /physics/ g,nu_air,nu_water,sfct,kappa,z,gustiness,dmin,    &
explim,sin_fac,sin_diss1,sin_diss2,sds_fac,sds_power,mss_fac,snl_fac,&
sdt_fac,sbf_fac,sbp_fac

namelist /grid/ gridfromfile,delx,dely,topofromfile,dpt,fillestuaries,&
filllakes

namelist /forcing/ winds,currents,air_density,water_density,seaice

namelist /forcing_constant/ wspd0,wdir0,uc0,vc0,rhoa0,rhow0,fice0,fice_lth,fice_uth

namelist /output/ outgrid,outspec,outrst,xpl,ypl,stokes

! read simulation parameters from the main namelist:
open(unit=21,file='namelists/main.nml',status='old',&
     form='formatted',access='sequential',action='read')
  read(unit=21,nml=domain)
  read(unit=21,nml=physics)
  read(unit=21,nml=grid)
  read(unit=21,nml=forcing)
  read(unit=21,nml=forcing_constant)
  read(unit=21,nml=output)
close(unit=21)

! check namelist values:
namelistok = .true.
if(nproc == 0)then

  if(mm < 3 .or. nm < 3)&
  call raiseexception('error','nmlread',                             &
                      'bad value in main.nml: mm and nm must be > 2',&
                      namelistok)

  if(mod(pm,4) /= 0)then
    call raiseexception('error','nmlread',                                 &
                        'bad value in main.nml: pm must be divisible by 4',&
                        namelistok)
  else
    if(mod(pm,8) /= 0)&
    call raiseexception('warning','nmlread',&
                        'pm should be divisible by 8 for optimal propagation properties')

  end if

  if(fmin <= 0 .or. fmax <= 0 .or. fprog <= 0)&
  call raiseexception('error','nmlread',                                        &
                      'bad value in main.nml: fmin, fmax and fprog must be > 0',&
                      namelistok)

  if(fmin >= fmax)&
  call raiseexception('error','nmlread',                           &
                      'bad value in main.nml: fmin must be < fmax',&
                      namelistok)

  if(fprog > fmax)&
  call raiseexception('warning','nmlread',&
                      'bad value in main.nml: fprog must be <= fmax; using highest allowed value')

  if(dtg <= 0)&
  call raiseexception('error','nmlread',                       &
                      'bad value in main.nml: dtg must be > 0',&
                      namelistok)

  if( z <= 0)&
  call raiseexception('error','nmlread',                     &
                      'bad value in main.nml: z must be > 0',&
                      namelistok)

  if(gustiness < 0)then
    call raiseexception('error','nmlread',                                  &
                        'bad value in main.nml: gustiness must be positive',&
                        namelistok)
  elseif(gustiness > 0.2)then
    call raiseexception('warning','nmlread',&
                        '(bad) value in main.nml: gustiness > 0.2; proceed with caution')
  end if

  if(dmin <= 0)&
  call raiseexception('error','nmlread',                        &
                      'bad value in main.nml: dmin must be > 0',&
                      namelistok)

  if(.not.gridfromfile)then
    if(delx <= 0 .or. dely <= 0)&
    call raiseexception('error','nmlread',                                 &
                        'bad value in main.nml: delx and dely must be > 0',&
                        namelistok)
  end if

  if(.not.topofromfile)then
    if(dpt <= 0)&
    call raiseexception('error','nmlread',                       &
                        'bad value in main.nml: dpt must be > 0',&
                        namelistok)
  end if

  if(.not.any(outgrid == allowedoutputtimes))&
  call raiseexception('error','nmlread',                                                     &
                      'bad value in main.nml: outgrid must be 0, 1, 2, 3, 4, 6, 8, 12 or 24',&
                      namelistok)

  if(.not.any(outspec == allowedoutputtimes))&
  call raiseexception('error','nmlread',                                                     &
                      'bad value in main.nml: outspec must be 0, 1, 2, 3, 4, 6, 8, 12 or 24',&
                      namelistok)

  if(.not.any(outrst == allowedoutputtimes))&
  call raiseexception('error','nmlread',                                                    &
                      'bad value in main.nml: outrst must be 0, 1, 2, 3, 4, 6, 8, 12 or 24',&
                      namelistok)

  if(.not.namelistok)then
    call raiseexception('abort','nmlread',&
                        'errors encountered in namelist read (see above)')
    stop
  end if

end if

! copy values from namelist into global variables
starttimestr_nml = starttimestr
stoptimestr_nml = stoptimestr

end subroutine nmlread
#endif


subroutine alloc(option)
  ! Allocates UMWM arrays

integer,intent(in) :: option

! allocate 2-d native arrays:
if(option==1)then
#ifdef GEOS
  allocate(ar_2d(mm,nm))
  allocate(curv(mm,nm))
  allocate(d_2d(mm,nm),dx_2d(mm,nm),dy_2d(mm,nm))
  allocate(dlon(mm,nm),dlat(mm,nm))
  allocate(ii(mm,nm))
  allocate(lat(mm,nm),lon(mm,nm))
  allocate(x(mm,nm),y(mm,nm))
  allocate(mask(mm,nm))
  allocate(nproc_out(mm,nm))
#else
  allocate(ar_2d(mm,nm))
  allocate(curv(mm,nm))
  allocate(d_2d(mm,nm),dx_2d(mm,nm),dy_2d(mm,nm))
  allocate(dlon(mm,nm),dlat(mm,nm))
  allocate(gustu(mm,nm),gustv(mm,nm))
  allocate(ii(mm,nm))
  allocate(lat(mm,nm),lon(mm,nm))
  allocate(x(mm,nm),y(mm,nm))
  allocate(mask(mm,nm))
  allocate(nproc_out(mm,nm))
  allocate(rhoa_2d(mm,nm),rhow_2d(mm,nm))
  allocate(wspd_2d(mm,nm))
  allocate(uc_2d(mm,nm),ucb(mm,nm),ucf(mm,nm))
  allocate(uw(mm,nm),uwb(mm,nm),uwf(mm,nm))
  allocate(vc_2d(mm,nm),vcb(mm,nm),vcf(mm,nm))
  allocate(vw(mm,nm),vwb(mm,nm),vwf(mm,nm))
  allocate(wdir_2d(mm,nm))
  allocate(fice_2d(mm,nm),ficeb(mm,nm),ficef(mm,nm))
#endif

! allocate remapped arrays:
elseif(option==2)then
#ifdef GEOS
  allocate(nu_water_(istart:iend))
#endif

  ! 1-d arrays:
  allocate(dom(om),f(om))
  allocate(cth(pm),cth2(pm),pl(pm),pr(pm),sth(pm),th(pm))

  allocate(cth_curv(pm,istart:iend))
  allocate(sth_curv(pm,istart:iend))

  allocate(iw(imm),ie(imm),is(imm),in(imm))
  allocate(iiw(imm),iie(imm),iis(imm),iin(imm))

  allocate(mi(imm),ni(imm))

  allocate(oc(istart:iend))

  ! 2-d arrays (remapped):
  allocate(ar(istart:iend))        ! grid cell surface area
  allocate(cd(istart:iend))        ! air-side drag coefficient
  allocate(d(imm),dx(imm),dy(imm)) ! depth and grid cell increments
  allocate(dxs(istart:iend),dxn(istart:iend)) ! cell edges
  allocate(dyw(istart:iend),dye(istart:iend)) ! cell edges
  allocate(fcutoff(istart:iend))   ! cutoff frequency
  allocate(ht(istart:iend))        ! significant wave height
#ifdef GEOS
  allocate(hts(istart:iend))       ! significant wave height of total swell
  allocate(htw(istart:iend))       ! significant wave height of wind waves
#endif
  allocate(mss(istart:iend))       ! mean-squared slope

  allocate(shelt(istart:iend))     ! sheltering coefficient
  shelt = 0

  allocate(physics_time_step(istart:iend))
  physics_time_step = 0 

  ! mean spectrum quantities:
  allocate(mwd(istart:iend)) ! direction
  allocate(mwp(istart:iend)) ! period
  allocate(mwl(istart:iend)) ! wavelength

  ! dominant spectrum quantities:
  allocate(dwd(istart:iend))  ! direction
  allocate(dwp(istart:iend))  ! period
  allocate(dwl(istart:iend))  ! wavelength
  allocate(dcp0(istart:iend)) ! intrinsic phase speed
  allocate(dcp(istart:iend))  ! phase speed
  allocate(dcg0(istart:iend)) ! intrinsic group speed
  allocate(dcg(istart:iend))  ! group speed

  dwd  = 0
  dwp  = 0
  dwl  = 0
  dcp0 = 0
  dcg0 = 0
  dcp  = 0
  dcg  = 0

  ! inverse area and grid cell increments:
  allocate(oneovar(istart:iend),oneovdx(istart:iend),oneovdy(istart:iend))

  allocate(momx(istart:iend),momy(istart:iend)) ! total wave momentum
  allocate(cgmxx(istart:iend),cgmxy(istart:iend),cgmyy(istart:iend)) ! cg*m

  ! momentum fluxes:
  allocate(taux(istart:iend),tauy(istart:iend))
  taux = 0; tauy = 0

  allocate(taux_form(istart:iend),tauy_form(istart:iend))
  taux_form = 0; tauy_form = 0

  allocate(taux_skin(istart:iend),tauy_skin(istart:iend))
  taux_skin = 0; tauy_skin = 0

  allocate(taux_diag(istart:iend),tauy_diag(istart:iend))
  taux_diag = 0; tauy_diag = 0

  allocate(taux_ocntop(istart:iend),tauy_ocntop(istart:iend))
  taux_ocntop = 0; tauy_ocntop = 0

  allocate(taux_ocnbot(istart:iend),tauy_ocnbot(istart:iend))
  taux_ocnbot = 0; tauy_ocnbot = 0

  allocate(taux_snl(istart:iend),tauy_snl(istart:iend))
  taux_snl = 0; tauy_snl = 0

  allocate(epsx_atm(istart:iend), epsy_atm(istart:iend))
  epsx_atm = 0; epsy_atm = 0

  allocate(epsx_ocn(istart:iend), epsy_ocn(istart:iend))
  epsx_ocn = 0; epsy_ocn = 0

  allocate(taux1(istart:iend),tauy1(istart:iend))
  allocate(taux2(istart:iend),tauy2(istart:iend))
  allocate(taux3(istart:iend),tauy3(istart:iend))

  taux1 = 0; tauy1 = 0
  taux2 = 0; tauy2 = 0
  taux3 = 0; tauy3 = 0

  allocate(tailatmx(istart:iend),tailatmy(istart:iend))
  allocate(tailocnx(istart:iend),tailocny(istart:iend))

  allocate(ustar(istart:iend)) ! air-side friction velocity

  allocate(wspd(imm)) ! wind speed
  wspd = 0

  allocate(fice(imm))
  fice = 0
  
  allocate(wdir(imm))                       ! wind direction
  wdir = 0

  allocate(uc(imm),vc(imm)) ! ocean currents
  uc = 0; vc = 0

#ifdef GEOS
  allocate(rhoa(imm))                       ! air density
  allocate(rhow(imm))                       ! water density
  allocate(rhorat(imm))                     ! air/water density ratio
#else
  allocate(rhoa(imm),rhoab(imm),rhoaf(imm)) ! air density
  allocate(rhow(imm),rhowb(imm),rhowf(imm)) ! water density
  allocate(rhorat(imm))                     ! air/water density ratio
#endif

  allocate(psim(imm)) ! integrated stability function for momentum
  psim = 0

  ! 3-d arrays (remapped):
  allocate(bf1_renorm(om,istart:iend))
  allocate(bf2_renorm(om,istart:iend))
  allocate(    cothkd(om,istart:iend))
  allocate(       dwn(om,istart:iend))
  allocate(     fkovg(om,istart:iend))
  allocate(    invcp0(om,istart:iend))
  allocate(        k4(om,istart:iend))
  allocate(       kdk(om,istart:iend))
  allocate(      k3dk(om,istart:iend))
  allocate( oneoverk4(om,istart:iend))
  allocate(       sbf(om,istart:iend))
  allocate(       sdt(om,istart:iend))
  allocate(       sdv(om,istart:iend))
  allocate(   snl_arg(om,istart:iend))
  allocate(        l2(om,istart:iend))
  allocate(logl2overz(om,istart:iend))

  allocate(psiml2(om,istart:iend))
  psiml2 = 0

  ! 4-d arrays (remapped):
  allocate(   ef(om,pm,istart:iend)) ! wave variance spectrum, forward in time
  allocate(dummy(om,pm,istart:iend)) ! dummy array, used in sds, advection and refraction
  allocate( rotl(om,pm,istart:iend)) ! anti-clockwise rotation, used in refraction
  allocate( rotr(om,pm,istart:iend)) ! clockwise rotation, used in refraction
  allocate(  sds(om,pm,istart:iend)) ! wave dissipation sink function
  allocate(  snl(om,pm,istart:iend)) ! wave downshifting source/sink function
  allocate( ssin(om,pm,istart:iend)) ! wind input source/sink function

  allocate( sice(om,istart:iend)) ! wave attenuation by sea ice function

end if ! if(option)

end subroutine alloc


subroutine grid
! Defines grid spacing and grid cell areas
#ifndef GEOS
use netcdf
use umwm_io,  only:nc_check
use umwm_util,only:raiseexception, distance_haversine
#else
use umwm_util,only:distance_haversine
#endif

logical :: loniscontinuous = .true.

integer :: m,n
integer :: ncid,varid,stat

real,parameter :: r_earth = 6.371009e6

real,dimension(:,:),allocatable :: abscoslat,lon_tmp,rotx,roty
real,dimension(:,:),allocatable :: rlon,rlat

if(gridfromfile)then
#ifndef GEOS
  stat = nf90_open('input/umwm.gridtopo',nf90_nowrite,ncid)

  if(stat/=0)then

    if(nproc==0)then
      call raiseexception('warning','grid',&
                          'input/umwm.gridtopo not found; trying input/umwm.grid')
    end if

    stat = nf90_open('input/umwm.grid',nf90_nowrite,ncid)

    if(stat/=0)then

      if(nproc==0)then
        call raiseexception('abort','grid',&
                            'input/umwm.grid not found. make sure input files are in place')
      end if

      stop

    end if

  end if

  call nc_check(nf90_inq_varid(ncid,'lon',varid))
  call nc_check(nf90_get_var(ncid,varid,lon))
  call nc_check(nf90_inq_varid(ncid,'lat',varid))
  call nc_check(nf90_get_var(ncid,varid,lat))
  call nc_check(nf90_close(ncid))
#endif
  allocate(abscoslat(mm,nm),lon_tmp(mm,nm),rotx(mm,nm),roty(mm,nm))
  allocate(rlon(mm,nm),rlat(mm,nm))

  ! figure out if longitude field is continuous:
  if(minval(lon) <- 175 .and. maxval(lon) >175)loniscontinuous = .false.

  lon_tmp = lon

  ! store original lon array before modifying:
  if(.not.loniscontinuous)where(lon < 0)lon = lon+360

  do n=1,nm
    do m=2,mm-1
      dlon(m,n) = 0.5*(lon(m+1,n)-lon(m-1,n))
    end do
  end do

  dlon(1,:)  = 2*dlon(2,:)-dlon(3,:)
  dlon(mm,:) = 2*dlon(mm-1,:)-dlon(mm-2,:)

  do n=2,nm-1
    do m=1,mm
      dlat(m,n) = 0.5*(lat(m,n+1)-lat(m,n-1))
    end do
  end do
  dlat(:,1)  = 2*dlat(:,2)-dlat(:,3)
  dlat(:,nm) = 2*dlat(:,nm-1)-dlat(:,nm-2)

  abscoslat = abs(cos(dr*lat))

  ! revert to original lon array:
  lon = lon_tmp

  rlon = lon*twopi/360.
  rlat = lat*twopi/360.

  do n=1,nm
    do m=2,mm-1
      dx_2d(m,n) = r_earth * distance_haversine(0.5*(rlon(m-1,n) + rlon(m  ,n)), &
                                                0.5*(rlon(m  ,n) + rlon(m+1,n)), &
                                                0.5*(rlat(m-1,n) + rlat(m  ,n)), &
                                                0.5*(rlat(m  ,n) + rlat(m+1,n)))
    end do
  end do

  dx_2d(1,:)  = 2*dx_2d(2,:)-dx_2d(3,:)
  dx_2d(mm,:) = 2*dx_2d(mm-1,:)-dx_2d(mm-2,:)

  do n=2,nm-1
    do m=1,mm
      dy_2d(m,n) = r_earth * distance_haversine(0.5*(rlon(m,n-1) + rlon(m,n  )), &
                                                0.5*(rlon(m,n  ) + rlon(m,n+1)), &
                                                0.5*(rlat(m,n-1) + rlat(m,n  )), &
                                                0.5*(rlat(m,n  ) + rlat(m,n+1)))
    end do
  end do

  dy_2d(:,1)  = 2*dy_2d(:,2)-dy_2d(:,3)
  dy_2d(:,nm) = 2*dy_2d(:,nm-1)-dy_2d(:,nm-2)

  ! compute grid rotation for great circle propagation
  curv = 0
  do n = 1,nm
    do m = 2,mm-1
      curv(m,n) = atan2(sin(rlon(m+1,n)-rlon(m-1,n))*cos(rlat(m+1,n)),&
                        cos(rlat(m-1,n))*sin(rlat(m+1,n))&
                       -sin(rlat(m-1,n))*cos(rlat(m+1,n))*cos(rlon(m+1,n)-rlon(m-1,n)))
    end do
  end do

  if(isglobal)then
    m = 1
    curv(m,:) = atan2(sin(rlon(m+1,:)-rlon(mm,:))*cos(rlat(m+1,:)),&
                      cos(rlat(mm,:))*sin(rlat(m+1,:))&
                     -sin(rlat(mm,:))*cos(rlat(m+1,:))*cos(rlon(m+1,:)-rlon(mm,:)))
    m = mm
    curv(m,:) = atan2(sin(rlon(1,:)-rlon(m-1,:))*cos(rlat(1,:)),&
                      cos(rlat(m-1,:))*sin(rlat(1,:))&
                     -sin(rlat(m-1,:))*cos(rlat(1,:))*cos(rlon(1,:)-rlon(m-1,:)))
  else
    curv(1,:) = curv(2,:)
    curv(mm,:) = curv(mm-1,:)
  end if

  curv = curv-0.5*pi

  deallocate(abscoslat,lon_tmp,rotx,roty,rlon,rlat)

else ! use constant value from namelist

  dx_2d = delx
  dy_2d = dely

  curv = 0

  lon = 0
  lat = 0
  dlon = 0
  dlat = 0

  x(1,:) = 0
  do m = 2,mm
    x(m,:) = x(m-1,:)+0.5*(dx_2d(m-1,:)+dx_2d(m,:))
  end do

  y(:,1) = 0
  do n = 2,nm
    y(:,n) = y(:,n-1)+0.5*(dy_2d(:,n-1)+dy_2d(:,n))
  end do

end if

if(topofromfile)then ! read depth field from file
#ifndef GEOS
  call nc_check(nf90_open('input/umwm.gridtopo',nf90_nowrite,ncid))
  call nc_check(nf90_inq_varid(ncid,'z',varid))
  call nc_check(nf90_get_var(ncid,varid,d_2d))
  call nc_check(nf90_close(ncid))
#endif
else ! use constant value from namelist

  d_2d = dpt

end if

! compute cell areas and reciprocals:
ar_2d = dx_2d*dy_2d

if(nproc == 0)then
  write(*,fmt=101)'umwm: grid: dx min/max/mean [m]:     ',&
                  minval(dx_2d),maxval(dx_2d),sum(dx_2d)/(mm*nm)
  write(*,fmt=101)'umwm: grid: dy min/max/mean [m]:     ',&
                  minval(dy_2d),maxval(dy_2d),sum(dy_2d)/(mm*nm)
  write(*,fmt=101)'umwm: grid: area min/max/mean [m^2]: ',&
                  minval(ar_2d),maxval(ar_2d),sum(ar_2d)/(mm*nm)
end if

101 format(a,3(f15.2,1x))

end subroutine grid


subroutine masks
! Defines landmasks and optionally removes one-cell wide estuaries and lakes

logical :: iterate

integer :: m, n
integer :: exm, exn
integer :: cnt, fillcount

! set masks:
if(topofromfile)then

  ! set initial seamask everywhere:
  mask = 1

  ! set up boundary points:
  mask(:,1)  = 0
  mask(:,nm) = 0

  ! close e and w edges if limited area:
  if(.not.isglobal)then
    mask(1,:)  = 0
    mask(mm,:) = 0
  end if

  ! set land mask where depth is non-negative, and then
  ! set depth to dmin:
  where(d_2d>=0)
    mask = 0
    d_2d = dmin
  endwhere

  ! make depths positive and limit to dmin:
  d_2d = abs(d_2d)
  where(d_2d<dmin)d_2d = dmin

  ! fill estuaries and isolated sea points:
  if(fillestuaries)then
    fillcount = 0
    iterate = .true.
    do while(iterate) ! iterate as long as there are points to be modified
      iterate = .false.
      do n=2,nm-1
        do m=2,mm-1

          if(mask(m,n)==1)then

            cnt = 0
            if(mask(m-1,n)==1)cnt = cnt+1
            if(mask(m+1,n)==1)cnt = cnt+1
            if(mask(m,n-1)==1)cnt = cnt+1
            if(mask(m,n+1)==1)cnt = cnt+1

            if(cnt<=1)then
              mask(m,n) = 0
              fillcount = fillcount+1
              iterate   = .true.
            end if

          end if ! if(mask(m,n)==1)

        end do
      end do
    end do ! while loop
    if(nproc==0)write(*,fmt=100)'umwm: masks: filled cells with 3-land neighbours,',&
                                fillcount,' cells total.'
  end if

  ! discard lakes or unwanted closed basin from the domain:
  if(filllakes)then
    open(unit=24,file='namelists/exclude.nml')
    do
      read(unit=24,fmt=*,end=107)exm,exn
      fillcount = 0
      call fill(exm,exn,fillcount)
      if(nproc==0)write(*,fmt=101)'umwm: masks: filled closed sea at i,j:',&
                                  exm,exn,', ',fillcount,' cells total.'
    end do
  end if ! if(filllakes)

107 close(unit=24)

100 format(a,i8,a)
101 format(a,2(i5,1x),a,i8,a)

else

  ! set initial mask everywhere:
  mask = 1

  ! set up boundary points
  mask(:, 1) = 0
  mask(:,nm) = 0

  ! close e and w edges if limited area:
  if(.not.isglobal)then
    mask( 1,:)  = 0
    mask(mm,:) = 0
  end if

  where(mask==0)d_2d = dmin

end if

! calculate the upper index for 1-d arrays:
im = count(mask==1)
imm = mm*nm

! print out summary:
if(nproc==0)then
  write(unit=*,fmt=302)imm
  write(unit=*,fmt=303)imm-im,float(imm-im)/float(imm)*100.
  write(unit=*,fmt=304)im,float(im)/float(imm)*100.
  write(unit=*,fmt=305)minval(d_2d,mask==1)
  write(unit=*,fmt=306)maxval(d_2d,mask==1)
end if

302 format('umwm: masks: total number of grid points: ',i9)
303 format('umwm: masks: number of land points:       ',i9,', ',f5.1,'%')
304 format('umwm: masks: number of sea points:        ',i9,', ',f5.1,'%')
305 format('umwm: masks: shallowest point:            ',f9.3,' meters')
306 format('umwm: masks: deepest point:               ',f9.3,' meters')

end subroutine masks



recursive subroutine fill(m, n, fillcount)
! Mask out enclosed seas chosen by user.

integer, intent(in) :: m, n
integer, intent(in out) :: fillcount

! return if reached domain edge:
if (any(m == [1, mm]) .or. any(n == [1, nm])) return

! fill:
mask(m,n) = 0

fillcount = fillcount+1

! recurse in all directions:
if (mask(m-1,n) == 1) call fill(m-1, n, fillcount)
if (mask(m+1,n) == 1) call fill(m+1, n, fillcount)
if (mask(m,n-1) == 1) call fill(m, n-1, fillcount)
if (mask(m,n+1) == 1) call fill(m, n+1, fillcount)

end subroutine fill


subroutine partition
! Partitions the domain for parallel computation.

#ifdef MPI
use umwm_mpi
#endif

integer :: nn

#ifdef ESMF
integer :: i, m, n
#endif

#ifndef MPI
istart  = 1
iend    = im
iistart = istart
iiend   = iend
#else

allocate(istart_(0:mpisize-1),iend_(0:mpisize-1),ilen_(0:mpisize-1))

#if (1)
im_mod = mod(im,mpisize)
if(im_mod==0)then
  ilen   = im/mpisize
  istart = nproc*ilen+1
  iend   = istart+ilen-1
else
  ilen   = floor(float(im)/float(mpisize))
  istart = nproc*ilen+1
  iend   = istart+ilen-1

  if (nproc == mpisize-1) then
    ilen = ilen + im_mod
    iend = iend + im_mod
  end if 
end if

#else
! find out what is the length of my part:
im_mod = mod(im,mpisize)
if(im_mod==0)then
  ilen   = im/mpisize
  istart = nproc*ilen+1
  iend   = istart+ilen-1
else
  ilen   = floor(float(im)/float(mpisize))
  istart = nproc*ilen+1
  iend   = istart+ilen-1
  do nn=1,im_mod
    if(nproc==nn-1)ilen = ilen+1
  end do
end if

! adjust start/end boundaries:
if(nproc<im_mod)then
  if(nproc==0)then
    iend = iend+1
  else
    istart = istart+nproc
    iend   = iend+nproc+1
  end if
else
  istart = istart+im_mod
  iend   = istart+ilen-1
end if
#endif

! Which direction for remapping? (matters only in parallel or global mode)
if (mm >= nm) then
  remap_dir = 'v'
else
  remap_dir = 'h'
end if

if (isglobal) remap_dir = 'v'
if (isglobal) remap_dir = 'v' !!! just testing:: 'h' -> out of bounds; 'v' -> negative counts


! The code below adjusts the start and end indices of each tile
! because currently ESMF DEBlockList accepts only regular rectangular
! domains.

#if defined(ESMF) && !defined(GEOS)
if (remap_dir == 'h') then ! row-major remapping

  ! adjust ends first:
  i = 0
  outer1: do n = 1, nm
    inner1: do m = 1, mm
      if (mask(m,n) == 1) i = i + 1
      if (i == iend .and. m /= mm) then
        iend = iend + count(mask(m+1:mm,n) == 1)
        exit outer1
      end if
    end do inner1
  end do outer1

  ! now adjust beginnings (all but proc 0):
  if (nproc /= 0) then
    i = 0
    outer2: do n = 1, nm
      inner2: do m = 1, mm
        if (mask(m,n) == 1) i = i + 1
        if (i == istart .and. m /= 1 .and. count(mask(1:m-1,n) == 1) > 0) then
          istart = istart + count(mask(m:mm,n) == 1)
          exit outer2
        end if
      end do inner2
    end do outer2
  end if

else if (remap_dir == 'v') then ! column-major remapping

  ! adjust ends first:
  i = 0
  outer3: do m = 1, mm
    inner3: do n = 1, nm
      if (mask(m,n) == 1) i = i + 1
      if (i == iend .and. n /= nm) then
        iend = iend + count(mask(m,n+1:nm) == 1)
        exit outer3
      end if
    end do inner3
  end do outer3

  ! now adjust beginnings (all but proc 0):
  if (nproc /= 0) then
    i = 0
    outer4: do m = 1, mm
      inner4: do n = 1, nm
        if (mask(m,n) == 1) i = i + 1
        if (i == istart .and. n /= 1 .and. count(mask(m,1:n-1) == 1) > 0) then
          istart = istart + count(mask(m,n:nm) == 1)
          exit outer4
        end if
      end do inner4
    end do outer4
  end if

end if

! adjust tile length:
ilen = iend - istart + 1

#endif

! gather tile mpisize information to root process:
call mpi_gather(istart,1,MPI_INTEGER,istart_,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call mpi_gather(iend,1,MPI_INTEGER,iend_,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call mpi_gather(ilen,1,MPI_INTEGER,ilen_,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#if (1)
if(nproc==0)then

  write(*,fmt='(a)')'umwm: partition: tiling summary:'
  write(*,fmt='(a)')'+---------+------------+----------------+----------------+'
  write(*,fmt='(a)')'|  nproc  |    ilen    |     istart     |      iend      |'
  write(*,fmt='(a)')'+---------+------------+----------------+----------------+'
!                      1234567   1234567890   12345678901234   12345678901234
  do nn=0,mpisize-1
    write(*,fmt=307)nn,ilen_(nn),istart_(nn),iend_(nn)
  end do

  write(*,fmt='(a)')'+---------+------------+----------------+----------------+'

end if

307 format('| ',i7,' | ',i10,' | ',2(i14,' | '))
#endif
#endif

end subroutine partition


subroutine remap
! Remaps two-dimensional arrays into one-dimensional
! arrays while leaving out land points. it also assigns
! and links neighboring points that are needed for
! spatial differencing, and builds the communication
! interface between processes.
#ifdef MPI
use umwm_mpi
#endif

integer :: i,m,n,nn
integer :: counter, itemp

integer, dimension(:), allocatable :: n_exchange_indices

! remapping section:
! unroll 2-d array into a contiguous 1-d array of sea-only points.
! first part contains only sea-points; the second part contains only
! land-points, which are necessary for grid mpisize information in advection
! routine.

! in serial mode this does not matter, but we must pick one:
#ifndef MPI
  remap_dir = 'v'
#endif

! first construct (m,n)->(i) transformation:
ii = 0
i  = 0
if(remap_dir == 'h')then ! column-major

  ! sea points:
  do n=1,nm
    do m=1,mm
      if(mask(m,n) == 1)then
        i = i+1
        ii(m,n) = i
        dx(i) = dx_2d(m,n)
        dy(i) = dy_2d(m,n)
        d(i)  = d_2d(m,n)
      end if
    end do
  end do

  ! land points:
  do n=1,nm
    do m=1,mm
      if(mask(m,n) == 0)then
        i = i+1
        ii(m,n) = i
        dx(i) = dx_2d(m,n)
        dy(i) = dy_2d(m,n)
        d(i)  = d_2d(m,n)
      end if
    end do
  end do

elseif(remap_dir == 'v')then ! row-major

  ! sea points:
  do m=1,mm
    do n=1,nm
      if(mask(m,n) == 1)then
        i = i+1
        ii(m,n) = i
        dx(i) = dx_2d(m,n)
        dy(i) = dy_2d(m,n)
        d(i)  = d_2d(m,n)
      end if
    end do
  end do

  ! land points:
  do m=1,mm
    do n=1,nm
      if(mask(m,n) == 0)then
        i = i+1
        ii(m,n) = i
        dx(i) = dx_2d(m,n)
        dy(i) = dy_2d(m,n)
        d(i)  = d_2d(m,n)
      end if
    end do
  end do

else

  stop 'umwm: remap: error - remap_dir must be ''h'' or ''v'''

end if

! now construct (i)->(m,n) transformation:
mi = 0
ni = 0
do n=1,nm
  do m=1,mm
    i = ii(m,n)
    mi(i) = m
    ni(i) = n
  end do
end do

! neighboring point indices for advection and refraction:
is = 0
in = 0
iw = 0
ie = 0
do n=2,nm-1
  do m=2,mm-1
    i = ii(m,n)
    is(i) = ii(m,n-1)
    in(i) = ii(m,n+1)
    iw(i) = ii(m-1,n)
    ie(i) = ii(m+1,n)
  end do
end do

! adjust periodic boundary:
if(isglobal)then
  do n=2,nm-1

    is(ii(1,n)) = ii(1,n-1)
    in(ii(1,n)) = ii(1,n+1)
    iw(ii(1,n)) = ii(mm,n)
    ie(ii(1,n)) = ii(2,n)

    is(ii(mm,n)) = ii(mm,n-1)
    in(ii(mm,n)) = ii(mm,n+1)
    iw(ii(mm,n)) = ii(mm-1,n)
    ie(ii(mm,n)) = ii(1,n)

  end do
end if

! true indices:
! (is,in,iw,ie will be aliased for land points)
iis = is
iin = in
iiw = iw
iie = ie

! cell edges in x and y:
do i = istart,iend
  dxn(i) = 0.5*(dx(i)+dx(iin(i)))
  dxs(i) = 0.5*(dx(i)+dx(iis(i)))
  dye(i) = 0.5*(dy(i)+dy(iie(i)))
  dyw(i) = 0.5*(dy(i)+dy(iiw(i)))
end do

#ifdef MPI
! allocate c, cg, and e, and initialize:
if(nproc==0)then ! root process

  iistart = istart ! start point
  itemp = iend     ! end point (first guess)

  do
    if(remap_dir == 'h')iiend = in(itemp)
    if(remap_dir == 'v')iiend = ie(itemp)
    ! if land point, go one cell back and cycle:
    if(iiend>im)then
      itemp = itemp-1
      cycle
   ! if sea point, exit loop:
    else
      exit
    end if
  end do

elseif(nproc==mpisize-1)then ! last process

  itemp = istart ! start point (first guess)
  iiend = iend   ! end point

  do
    if(remap_dir == 'h')iistart = is(itemp)
    if(remap_dir == 'v')iistart = iw(itemp)
    ! if land point, go one cell forward and cycle:
    if(iistart>im)then
      itemp = itemp+1
      cycle
    ! if sea point, exit loop:
    else
      exit
    end if
  end do

else ! everybody else

  itemp = istart ! start point (first guess)

  do
    if(remap_dir == 'h')iistart = is(itemp)
    if(remap_dir == 'v')iistart = iw(itemp)
    ! if land point, go one cell forward and cycle:
    if(iistart>im)then
      itemp = itemp+1
      cycle
      ! if sea point, exit loop:
    else
      exit
    end if
  end do

  itemp = iend ! end point (first guess)

  do
    if(remap_dir == 'h')iiend = in(itemp)
    if(remap_dir == 'v')iiend = ie(itemp)
    ! if land point, go one cell back and cycle:
    if(iiend>im)then
      itemp = itemp-1
      cycle
    ! if sea point, exit loop:
    else
      exit
    end if
  end do

end if ! if(nproc==0)

if(isglobal)then

  allocate(n_exchange_indices(0))

  first_col_len = 0
  do n = 2, nm-1
    if (mask(1,n) == 1 .and. mask(mm,n) == 1) then
      n_exchange_indices = [n_exchange_indices, n]
      first_col_len = first_col_len + 1
    end if
  end do
  last_col_len = first_col_len

  allocate(i_exchange_indices(first_col_len))

  ! find west neighbor indices on first processor
  if (nproc == 0) then

    ! find unrolled halo exchange indices
    do n = 1, first_col_len
      i_exchange_indices(n) = ii(1,n_exchange_indices(n))
    end do

    ! adjust the start halo index for the number
    ! of water points in easternmost column
    iistart = iistart - last_col_len

    counter = 0
    m = 1
    do n = 2, nm-1
      if (mask(m,n) == 1) then
        i = ii(m,n)
        if (mask(mi(iw(i)),ni(iw(i))) == 0) then
          ! west neighbor is land
          iw(i) = iistart - 1
        else
          ! west neighbor is water
          iw(i) = iistart + counter
          counter = counter + 1
        end if
      end if
    end do

  end if

  ! find east neighbor indices on last processor
  if (nproc == mpisize-1) then

    ! find unrolled halo exchange indices
    do n = 1, first_col_len
      i_exchange_indices(n) = ii(mm,n_exchange_indices(n))
    end do

    ! adjust the end halo index for the number
    ! of water points in easternmost column
    iiend = iiend + first_col_len

    counter = 1
    m = mm
    do n = 2, nm-1
      if (mask(m,n) == 1) then
        i = ii(m,n)
        if (mask(mi(ie(i)),ni(ie(i))) == 0) then
          ! east neighbor is land
          ie(i) = iistart - 1
        else
          ! east neighbor is water
          ie(i) = iend + counter
          counter = counter + 1
        end if
      end if
    end do

  end if

end if

#endif

! allocate and add a ghost land point at the beggining:
allocate(e(om,pm,iistart-1:iiend))
allocate(cp0(om,iistart-1:iiend))
allocate(cg0(om,iistart-1:iiend))
allocate(k(om,istart:iend))

! initialize:
e = tiny(e)

#ifdef MPI
! distribute iistart and iiend to everyone:
allocate(iistart_(0:mpisize-1),iiend_(0:mpisize-1))
call mpi_allgather(iistart,1,MPI_INTEGER,iistart_,1,&
                   MPI_INTEGER,MPI_COMM_WORLD,ierr)
call mpi_allgather(iiend,1,MPI_INTEGER,iiend_,1,&
                   MPI_INTEGER,MPI_COMM_WORLD,ierr)
#endif

! treat land points for cp0, cg0, and e:
! (needed for advection terms)
do i = istart, iend
  if (is(i) < iistart .or. is(i) > iiend) is(i) = iistart - 1
  if (in(i) < iistart .or. in(i) > iiend) in(i) = iistart - 1
  if (ie(i) < iistart .or. ie(i) > iiend) ie(i) = iistart - 1
  if (iw(i) < iistart .or. iw(i) > iiend) iw(i) = iistart - 1
end do

#ifdef MPI
! create the domain partiotioning field for output:
if (nproc == 0) then
nproc_out = -1
  do n = 1, nm
    do m = 1, mm
      do nn = 0, mpisize-1
        if (ii(m,n) >= istart_(nn) .and. ii(m,n) <= iend_(nn)) then
          nproc_out(m,n) = nn
          exit
        end if
      end do
    end do
  end do
end if
#if (1)
if(nproc==0)then

  write(*,fmt='(a)')'umwm: remap: tiling with halo summary:'
  write(*,fmt='(a)')'+---------+------------+----------------+----------------+'
  write(*,fmt='(a)')'|  nproc  |   iilen    |    iistart     |     iiend      |'
  write(*,fmt='(a)')'+---------+------------+----------------+----------------+'
!                      1234567   1234567890   12345678901234   12345678901234
  do nn=0,mpisize-1
    write(*,fmt=310)nn,iiend_(nn)-iistart_(nn)+1,iistart_(nn),iiend_(nn)
  end do

  write(*,fmt='(a)')'+---------+------------+----------------+----------------+'

end if

310 format('| ',i7,' | ',i10,' | ',2(i14,' | '))
#endif
#endif

end subroutine remap


subroutine init
! Initialize model variables such as frequencies, direction angles,
! phase speed and group velocity, wave numbers, etc.
#ifdef MPI
use umwm_mpi, only:istart_,iend_
#endif
use umwm_util,only:raiseexception
#ifndef GEOS
use umwm_io, only: winds,seaice
use umwm_util, only: remap_mn2i
#endif

integer :: i, n, o, p, pp, ind

! set frequency increment:
dlnf = (log(fmax)-log(fmin))/float(om-1)

! set frequency bins:
do o=1,om
  f(o) = exp(log(fmin)+(o-1)*dlnf)
end do

! define various constants:
dth           = twopi/float(pm)
dthg          = dth*g
oneovdth      = 1./dth
log10overz    = log(10./z)
twopisds_fac  = twopi*sds_fac
twonu         = 2.*nu_water
fieldscale1   = sin_diss1/sin_fac
fieldscale2   = sin_diss2/sin_diss1
inv_sds_power = 1./sds_power

! this limits the Courant number to its theoretical value slightly
! larger than 1/sqrt(2), depending on the number of directional bins;
! this also limits the number of directions to be divisible by 8,
! not 4! (more isotropic in Cartesian projection)

if(mod(pm,8) == 0)then
  cfllim = cos(0.25*pi-0.5*dth)
else
  cfllim = 1./sqrt(2.)
end if

! compute grid cell areas and reciprocals:
ar      = dx(istart:iend)*dy(istart:iend)
oneovdx = 1./dx(istart:iend)
oneovdy = 1./dy(istart:iend)
oneovar = 1./ar

! compute diffusion values in 2 frequenciess:
bf1  = exp(-16*dlnf*dlnf)
bf2  = exp(-64*dlnf*dlnf)
bf1a = bf1/(bf1+bf2)
bf2  = bf2/(bf1+bf2)
bf1  = bf1a

do p=1,pm
  th(p) = (p-0.5*(pm+1))*dth ! angles
end do

cth = cos(th) ! cosines
sth = sin(th) ! sines

! calculate wave ray directions adjusted for grid curvature:
do i=istart,iend
  do p=1,pm
    cth_curv(p,i) = cos(th(p)+curv(mi(i),ni(i)))
    sth_curv(p,i) = sin(th(p)+curv(mi(i),ni(i)))
  end do
end do

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
call dispersion(1e-2)

mindelx = min(minval(dx_2d,mask==1),minval(dy_2d,mask==1))
cgmax   = maxval(cg0(:,istart:iend))
dtamin  = 0.98*cfllim*mindelx/cgmax

first    = .true.

if(restart)then
  firstdtg = .false.
else
  firstdtg = .true.
end if

#ifndef GEOS
! if sea ice from file, update the fice field
if (seaice) fice = remap_mn2i(ficef)

! if forcing from file, update the wspd field for ustar first guess
if (winds) wspd = remap_mn2i(sqrt(uwf**2 + vwf**2))

! initialize drag coefficient (Large and Pond, 1981):
cd = 1.2e-3
do concurrent (i=istart:iend, wspd(i) > 11)
  cd(i) = (0.49 + 0.065 * wspd(i)) * 1e-3
end do

! initialize friction velocity:
do i=istart,iend
  ustar(i) = sqrt(cd(i))*wspd(i)
end do
#endif

#ifndef GEOS
#ifdef MPI

! figure out which process will print to screen:
iip = ii(xpl,ypl)
if(mask(xpl,ypl)==0)then
  if(nproc==0)then
    write(0,*)xpl,ypl,mask(xpl,ypl)

    call raiseexception('warning','init',&
                        'land-point chosen for stdout, may go out of bounds')
    call raiseexception('warning','init',&
                        'check xpl and ypl in the output namelist in namelists/main.nml')
    stop
  end if
end if

if(nproc==0)then
  do n=0,mpisize-1
    if(iip>=istart_(n).and.iip<=iend_(n))nproc_plot = n
  end do
end if

call mpi_bcast(nproc_plot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(iip,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

#else

iip = ii(xpl,ypl)

#endif
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
                  minval(cp0(o,1:im),dim=1),maxval(cp0(o,1:im),dim=1),&
                  minval(cg0(o,1:im),dim=1),maxval(cg0(o,1:im),dim=1)
end do
101 format(1x,i2,8(2x,f7.4))
write(*,fmt=102)

102 format('!',77('='),'!')
103 format('!',77('-'),'!')

#endif

end subroutine init


subroutine dispersion(tol)
! Iteratively solve the dispersion relation by iteration,
! and compute phase and group velocities in absolute reference frame (w/ currents)
#ifdef MPI
use umwm_mpi
#endif

integer :: sendcount, recvcount
integer :: sendtag, recvtag
integer :: src, dest

integer :: counter,i,o
real,intent(in) :: tol
real :: dk

real,dimension(istart:iend) :: b
real,dimension(om,istart:iend) :: f_nd,kd,t

if(nproc==0)write(*,'(a)')'umwm: dispersion: solving for dispersion relationship;'

! non-dimesionalize frequencies, and use deep water limit
! as initial guess:
do concurrent (o=1:om, i=istart:iend)
  cp0(o,i)  = twopi*sqrt(d(i)/g)
  f_nd(o,i) = cp0(o,i)*f(o)
  k(o,i)    = f_nd(o,i)*f_nd(o,i)
end do

! non-dimesionalize surface tension:
b = sfct/(rhow0*g*d(istart:iend)**2)

do i=istart,iend
  do o=1,om

    counter = 1
    dk = 2.*tol

    do while(abs(dk) > tol) ! newton-raphson iteration loop

      t(o,i) = tanh(k(o,i))
      dk = -(f_nd(o,i)*f_nd(o,i)-k(o,i)*t(o,i)  &
           *(1.+b(i)*k(o,i)*k(o,i)))            &
           /(3.*b(i)*k(o,i)*k(o,i)*t(o,i)+t(o,i)&
           +k(o,i)*(1.+b(i)*k(o,i)*k(o,i))*(1.-t(o,i)*t(o,i)))
      k(o,i) = k(o,i)-dk

      if(counter == 1000)exit ! escape if stuck
      counter = counter+1

    end do

    k(o,i) = abs(k(o,i))/d(i) ! f(k)=f(-k), so k>0 == k<0 roots

  end do
end do

if(nproc==0)write(*,'(a)')'umwm: dispersion: dispersion relationship done;'

do concurrent (o=1:om, i=istart:iend)
  kd(o,i) = k(o,i) * d(i)
end do

! limit kd to avoid floating overflow in transcendental functions:
where(kd>20.)kd = 20.

! phase speed and group velocity:
cp0 = tiny(cp0)
cg0 = tiny(cg0)

do concurrent (o=1:om, i=istart:iend)

  cp0(o,i) = twopi*f(o)/k(o,i)
  cg0(o,i) = cp0(o,i)*(0.5+k(o,i)*d(i)/sinh(2.*kd(o,i))&
                      +sfct*k(o,i)*k(o,i)/(rhow0*g+sfct*k(o,i)*k(o,i)))
end do

! compute some frequently used arrays:
do concurrent (o=1:om, i=istart:iend)

  dwn(o,i)       = dom(o)/abs(cg0(o,i))                         ! dk
  l2(o,i)        = 0.5*abs(cp0(o,i))/f(o)                       ! lambda/2 (half wavelength)
  k4(o,i)        = k(o,i)**4.                                   ! k^4
  oneoverk4(o,i) = 1./k4(o,i)                                   ! k^-4
  kdk(o,i)       = k(o,i)*dwn(o,i)                              ! k*dk
  k3dk(o,i)      = k(o,i)**3.*dwn(o,i)                          ! k*k*k*dk
  fkovg(o,i)     = f(o)*k(o,i)/g                                ! f*k/g
  cothkd(o,i)    = cosh(0.2*kd(o,i))/sinh(0.2*kd(o,i))          ! coth(0.2*kd)
  invcp0(o,i)    = 1./cp0(o,i)                                  ! 1/cp
  sbf(o,i)       = sbf_fac*k(o,i)/(sinh(2.*kd(o,i)))&           ! bottom friction
                  +sbp_fac*k(o,i)/(cosh(kd(o,i))*cosh(kd(o,i))) ! bottom percolation
#ifdef GEOS
  sdv(o,i)       = 4.*nu_water_(i)*k(o,i)**2.                   ! Prognostic water viscosity
#else
  sdv(o,i)       = 4.*nu_water*k(o,i)**2.                       ! viscosity
#endif
end do

! compute renormalization factors for snl:
bf1_renorm = 0.
bf2_renorm = 0.
snl_arg    = 0.

do i=istart,iend
  do o=1,om-2
    bf1_renorm(o,i) = snl_fac*bf1*kdk(o+1,i)/kdk(o,i)
    bf2_renorm(o,i) = snl_fac*bf2*kdk(o+2,i)/kdk(o,i)
    snl_arg(o,i)    = 1.-(bf1_renorm(o,i)+bf2_renorm(o,i))
  end do
end do

! half-wavelength over z:
logl2overz = log(l2/z)

! limit wind input to be at 10 m for l/2 > 10 m:
where(l2>20.)logl2overz = log(20./z)

#ifdef MPI
if(nproc<mpisize-1)then ! communicate with process above:

  sendcount = om*(iend-iistart_(nproc+1)+1) ; dest = nproc+1 ; sendtag = nproc
  recvcount = om*(iiend-iend)               ; src  = nproc+1 ; recvtag = src

  if (sendcount < 1 .or. recvcount < 1) then
      print *, '***** L1589', nproc, sendcount, recvcount, istart, iend, iistart, iiend, iistart_(nproc+1)
!            ----------------------------------------------------------------------------------------------
!                ***** L1589   0     -1632662    6105       1       116   -163     281    44243
!                ***** L1589   382    5994      -1632662    44062   44176  43900   50     44015
  end if
  call mpi_sendrecv(cp0(:,iistart_(nproc+1):iend),sendcount,&
                    MPI_REAL,dest,sendtag,                  &
                    cp0(:,iend+1:iiend),recvcount,          &
                    MPI_REAL,src,recvtag,                   &
                    MPI_COMM_WORLD,status,ierr)

end if

if(nproc>0)then ! communicate with process below:

  sendcount = om*(iiend_(nproc-1)-istart+1) ; dest = nproc-1 ; sendtag = nproc
  recvcount = om*(istart-iistart)           ; src  = nproc-1 ; recvtag = src

  if (sendcount < 1 .or. recvcount < 1) then
  !    print *, '***** L1605', nproc, sendcount, recvcount, istart, iend, iistart, iiend, iistart_(nproc+1) 
  end if
  call mpi_sendrecv(cp0(:,istart:iiend_(nproc-1)),sendcount,&
                    MPI_REAL,dest,sendtag,                  &
                    cp0(:,iistart:istart-1),recvcount,      &
                    MPI_REAL,src,recvtag,                   &
                    MPI_COMM_WORLD,status,ierr)

end if

call mpi_barrier(MPI_COMM_WORLD,ierr)

if(nproc<mpisize-1)then ! communicate with process above:

  sendcount = om*(iend-iistart_(nproc+1)+1) ; dest = nproc+1 ; sendtag = nproc
  recvcount = om*(iiend-iend)               ; src  = nproc+1 ; recvtag = src

  if (sendcount < 1 .or. recvcount < 1) then
  !    print *, '***** L1623', nproc, sendcount, recvcount, istart, iend, iistart, iiend, iistart_(nproc+1) 
  end if
  call mpi_sendrecv(cg0(:,iistart_(nproc+1):iend),sendcount,&
                    MPI_REAL,dest,sendtag,                  &
                    cg0(:,iend+1:iiend),recvcount,          &
                    MPI_REAL,src,recvtag,                   &
                    MPI_COMM_WORLD,status,ierr)

end if

if(nproc>0)then ! communicate with process below:

  sendcount = om*(iiend_(nproc-1)-istart+1) ; dest = nproc-1 ; sendtag = nproc
  recvcount = om*(istart-iistart)           ; src  = nproc-1 ; recvtag = src

  if (sendcount < 1 .or. recvcount < 1) then
  !    print *, '***** L1639', nproc, sendcount, recvcount, istart, iend, iistart, iiend, iistart_(nproc+1) 
  end if  
  call mpi_sendrecv(cg0(:,istart:iiend_(nproc-1)),sendcount,&
                    MPI_REAL,dest,sendtag,                  &
                    cg0(:,iistart:istart-1),recvcount,      &
                    MPI_REAL,src,recvtag,                   &
                    MPI_COMM_WORLD,status,ierr)

end if

call mpi_barrier(MPI_COMM_WORLD,ierr)

! if periodic domain, connect the east and west:
if(isglobal)then

  if(nproc==0)then ! communicate with last tile:

    sendcount = om*first_col_len ; dest = mpisize-1 ; sendtag = nproc
    recvcount = om*last_col_len  ; src  = mpisize-1 ; recvtag = src

    if (sendcount < 1 .or. recvcount < 1) then
    !  print *, '***** L1660', nproc, sendcount, recvcount, istart, iend, iistart, iiend, iistart_(nproc+1) 
    end if
    call mpi_sendrecv(cp0(:,istart:(istart+first_col_len-1)),sendcount,&
                      MPI_REAL,dest,sendtag,                           &
                      cp0(:,iistart:istart-1),recvcount,               &
                      MPI_REAL,src,recvtag,                            &
                      MPI_COMM_WORLD,status,ierr)

  end if

  if(nproc==mpisize-1)then ! communicate with first tile:

    sendcount = om*last_col_len  ; dest = 0 ; sendtag = nproc
    recvcount = om*first_col_len ; src  = 0 ; recvtag = src

    if (sendcount < 1 .or. recvcount < 1) then
    !  print *, '***** L1676', nproc, sendcount, recvcount, istart, iend, iistart, iiend, iistart_(nproc+1) 
    end if
    call mpi_sendrecv(cp0(:,(iend-last_col_len+1):iend),sendcount,&
                      MPI_REAL,dest,sendtag,                      &
                      cp0(:,iend+1:iiend),recvcount,              &
                      MPI_REAL,src,recvtag,                       &
                      MPI_COMM_WORLD,status,ierr)

  end if

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  if(nproc==0)then ! communicate with last tile:

    sendcount = om*first_col_len ; dest = mpisize-1 ; sendtag = nproc
    recvcount = om*last_col_len  ; src  = mpisize-1 ; recvtag = src

    if (sendcount < 1 .or. recvcount < 1) then
    !  print *, '***** L1694', nproc, sendcount, recvcount, istart, iend, iistart, iiend, iistart_(nproc+1) 
    end if
    call mpi_sendrecv(cg0(:,istart:(istart+first_col_len-1)),sendcount,&
                      MPI_REAL,dest,sendtag,                           &
                      cg0(:,iistart:istart-1),recvcount,               &
                      MPI_REAL,src,recvtag,                            &
                      MPI_COMM_WORLD,status,ierr)

  end if

  if(nproc==mpisize-1)then ! communicate with first tile:

    sendcount = om*last_col_len  ; dest = 0 ; sendtag = nproc
    recvcount = om*first_col_len ; src  = 0 ; recvtag = src

    if (sendcount < 1 .or. recvcount < 1) then
    !  print *, '***** L1710', nproc, sendcount, recvcount, istart, iend, iistart, iiend, iistart_(nproc+1) 
    end if
    call mpi_sendrecv(cg0(:,(iend-last_col_len+1):iend),sendcount,&
                      MPI_REAL,dest,sendtag,                      &
                      cg0(:,iend+1:iiend),recvcount,              &
                      MPI_REAL,src,recvtag,                       &
                      MPI_COMM_WORLD,status,ierr)

  end if

end if
#endif

! handle land points for cp and cg (needed for advection/refraction)
do o=1,om
  cp0(o,iistart-1) = minval(cp0(o,istart:iend))
  cg0(o,iistart-1) = minval(cg0(o,istart:iend))
end do

end subroutine dispersion

end module umwm_init
