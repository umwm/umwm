module umwm_io
! Provides input/output routines for the wave model
use umwm_module
use netcdf

logical :: readfile
logical :: winds,currents,air_density,water_density,seaice

contains



subroutine input_nc(timestr)
! Reads input data files for atmospheric and oceanic fields.
use umwm_util, only: remap_mn2i

character(19), intent(in) :: timestr
character(999) :: nc_infile
character(19) :: readstr
integer :: ncid,varid

readstr = timestr
readstr(11:11) = '_'

! forcing fields are input from this file:
nc_infile = 'input/umwmin_' // readstr // '.nc'

! set the logical switch to .true. only if from file is requested
! for any of the fields:
readfile = any([winds, currents, air_density, water_density, seaice])

if (readfile) call nc_check(nf90_open(trim(nc_infile), nf90_nowrite, ncid))

if(winds)then
  call nc_check(nf90_inq_varid(ncid,'uw',varid))
  call nc_check(nf90_get_var(ncid,varid,uwf))
  call nc_check(nf90_inq_varid(ncid,'vw',varid))
  call nc_check(nf90_get_var(ncid,varid,vwf))
else
  wspd = wspd0
  wdir = wdir0
end if

if(currents)then
  call nc_check(nf90_inq_varid(ncid,'uc',varid))
  call nc_check(nf90_get_var(ncid,varid,ucf))
  call nc_check(nf90_inq_varid(ncid,'vc',varid))
  call nc_check(nf90_get_var(ncid,varid,vcf))
else
  ucf = uc0
  vcf = vc0
  uc  = uc0
  vc  = vc0
end if

if(seaice)then
  call nc_check(nf90_inq_varid(ncid,'fice',varid))
  call nc_check(nf90_get_var(ncid,varid,ficef))
else
  fice_2d = fice0
  fice    = fice0
end if

where(mask == 0)
  ucf = 0
  vcf = 0
endwhere

if(air_density)then
  call nc_check(nf90_inq_varid(ncid,'rhoa',varid))
  call nc_check(nf90_get_var(ncid,varid,rhoa_2d))
else
  rhoa_2d = rhoa0
  rhoa    = rhoa0
end if

if(water_density)then
  call nc_check(nf90_inq_varid(ncid,'rhow',varid))
  call nc_check(nf90_get_var(ncid,varid,rhow_2d))
else
  rhow_2d = rhow0
  rhow    = rhow0
end if

if(readfile)then
  call nc_check(nf90_close(ncid))
end if

! remap to 1-d arrays:
rhoaf = remap_mn2i(rhoa_2d)
rhowf = remap_mn2i(rhow_2d)

end subroutine input_nc


subroutine output_grid
! Outputs grid related fields into a netcdf file.

integer :: ncid
integer :: xdimid,ydimid
integer :: lonid,latid,dlonid,dlatid,dxid,dyid,arid,maskid,did,nprocid
integer :: xid,yid,curvid

if(nproc == 0)then

  call nc_check(nf90_create('output/umwmout.grid',nf90_clobber,ncid))
  call nc_check(nf90_def_dim(ncid,'x',mm,xdimid))
  call nc_check(nf90_def_dim(ncid,'y',nm,ydimid))
  call nc_check(nf90_def_var(ncid,'lon',nf90_float,[xdimid,ydimid],lonid))
  call nc_check(nf90_def_var(ncid,'lat',nf90_float,[xdimid,ydimid],latid))
  call nc_check(nf90_def_var(ncid,'xx',nf90_float,[xdimid,ydimid],xid))
  call nc_check(nf90_def_var(ncid,'yy',nf90_float,[xdimid,ydimid],yid))
  call nc_check(nf90_def_var(ncid,'dlon',nf90_float,[xdimid,ydimid],dlonid))
  call nc_check(nf90_def_var(ncid,'dlat',nf90_float,[xdimid,ydimid],dlatid))
  call nc_check(nf90_def_var(ncid,'dx',nf90_float,[xdimid,ydimid],dxid))
  call nc_check(nf90_def_var(ncid,'dy',nf90_float,[xdimid,ydimid],dyid))
  call nc_check(nf90_def_var(ncid,'curvature',nf90_float,[xdimid,ydimid],curvid))
  call nc_check(nf90_def_var(ncid,'area',nf90_float,[xdimid,ydimid],arid))
  call nc_check(nf90_def_var(ncid,'depth',nf90_float,[xdimid,ydimid],did))
  call nc_check(nf90_def_var(ncid,'seamask',nf90_int,[xdimid,ydimid],maskid))
  call nc_check(nf90_def_var(ncid,'nproc',nf90_int,[xdimid,ydimid],nprocid))
  call nc_check(nf90_enddef(ncid))
  call nc_check(nf90_put_var(ncid,lonid,lon))
  call nc_check(nf90_put_var(ncid,latid,lat))
  call nc_check(nf90_put_var(ncid,xid,x))
  call nc_check(nf90_put_var(ncid,yid,y))
  call nc_check(nf90_put_var(ncid,dlonid,dlon))
  call nc_check(nf90_put_var(ncid,dlatid,dlat))
  call nc_check(nf90_put_var(ncid,dxid,dx_2d))
  call nc_check(nf90_put_var(ncid,dyid,dy_2d))
  call nc_check(nf90_put_var(ncid,curvid,curv))
  call nc_check(nf90_put_var(ncid,arid,ar_2d))
  call nc_check(nf90_put_var(ncid,maskid,mask))
  call nc_check(nf90_put_var(ncid,did,d_2d))
  call nc_check(nf90_put_var(ncid,nprocid,nproc_out))
  call nc_check(nf90_close(ncid))

  write(*,'(a)')'umwm: output_grid: grid definition written in output/umwmout.grid'

end if

end subroutine output_grid


subroutine output_spectrum_nc(timestr)
! Writes out model spectrum output in a netcdf format

! arguments:
character(19),intent(in) :: timestr

character(19),save :: savetimestr

character(9999) :: spectrumoutputfile

character(2) :: coord_input

integer,dimension(2) :: xy_coords

integer :: stat
integer :: ncid
integer :: fdimid,thdimid,tdimid,scalarid
integer :: lon_scalarid,lat_scalarid,wspdid,wdirid
integer :: specid,sinid,sdsid,snlid,freqid,thetaid,wlid
integer :: sdtid,sdvid,sbfid

integer                               :: nn
integer,save                          :: npts = 0
integer,          dimension(999),save :: mspec,nspec,ispec
character(3),dimension(999),save :: cnn
character(40),dimension(999),save :: spectrumid

real :: latspec,lonspec
real :: wspdtmp,wdirtmp

integer,save :: counter = 1
logical,save :: firstrun = .true.

if(firstrun)then

  savetimestr = timestr
  savetimestr(11:11) = '_'

  ! open spectrum points list file:
  open(unit=21,file='namelists/spectrum.nml',status='old',&
       form='formatted',access='sequential',err=100)

  ! are points given in lat/lon or grid indices?
  read(unit=21,fmt='(2a)')coord_input

  ! loop over points in the list:
  do
    npts = npts + 1

    if (any(coord_input == ['xy', 'XY'])) then

      read(21, *, end=100) mspec(npts), nspec(npts), spectrumid(npts)

      if (mspec(npts) < 2 .or. mspec(npts) > mm - 1 .or. &
          nspec(npts) < 2 .or. nspec(npts) > nm - 1) &
          stop 'umwm: output_nc: error: a requested point ' &
            // 'in spectrum.nml is out of bounds'

    else if (any(coord_input == ['ll', 'LL'])) then

      read(21, *, end=100) lonspec, latspec, spectrumid(npts)
      xy_coords = minloc((lonspec - lon)**2 + (latspec - lat)**2)
      mspec(npts) = xy_coords(1)
      nspec(npts) = xy_coords(2)

    else

      stop 'umwm: output_nc: error: first line in namelists/spectrum.nml must' &
        // 'contain "xy" (or "XY") or "ll" (or "LL").'

    end if

    ispec(npts) = ii(mspec(npts),nspec(npts))
    write(unit=cnn(npts),fmt='(i3)')npts
    if(npts<100)cnn(npts) = '0'//adjustl(cnn(npts))
    if(npts< 10)cnn(npts) = '0'//adjustl(cnn(npts))

  end do

  100 npts=npts-1
  close(unit=21)

end if

do nn=1,npts
  if(ispec(nn) >= istart .and. ispec(nn) <= iend)then

    if(firstrun)then

      ! create output file:
      spectrumoutputfile = 'output/umwmspc_'//trim(spectrumid(nn))//'_'//savetimestr//'.nc'

      stat = nf90_create(trim(spectrumoutputfile),nf90_clobber,ncid)

      ! define dimensions:
      stat = nf90_def_dim(ncid,'scalar',1,scalarid)
      stat = nf90_def_dim(ncid,'f',om,fdimid)
      stat = nf90_def_dim(ncid,'th',pm,thdimid)
      stat = nf90_def_dim(ncid,'time',nf90_unlimited,tdimid)

      ! define variables:
      stat = nf90_def_var(ncid,'frequency',nf90_float,[fdimid],freqid)
      stat = nf90_def_var(ncid,'wavenumber',nf90_float,[fdimid],wlid)
      stat = nf90_def_var(ncid,'direction',nf90_float,[thdimid],thetaid)
      stat = nf90_def_var(ncid,'longitude',nf90_float,[scalarid],lon_scalarid)
      stat = nf90_def_var(ncid,'latitude',nf90_float,[scalarid],lat_scalarid)
      stat = nf90_def_var(ncid,'wspd',nf90_float,[scalarid,tdimid],wspdid)
      stat = nf90_def_var(ncid,'wdir',nf90_float,[scalarid,tdimid],wdirid)
      stat = nf90_def_var(ncid,'f',nf90_float,[fdimid,thdimid,tdimid],specid)
      stat = nf90_def_var(ncid,'sin',nf90_float,[fdimid,thdimid,tdimid],sinid)
      stat = nf90_def_var(ncid,'sds',nf90_float,[fdimid,thdimid,tdimid],sdsid)
      stat = nf90_def_var(ncid,'sdt',nf90_float,[fdimid,tdimid],sdtid)
      stat = nf90_def_var(ncid,'sdv',nf90_float,[fdimid],sdvid)
      stat = nf90_def_var(ncid,'sbf',nf90_float,[fdimid],sbfid)
      stat = nf90_def_var(ncid,'snl',nf90_float,[fdimid,thdimid,tdimid],snlid)

      ! end of definition mode:
      stat = nf90_enddef(ncid)

      ! fill in static fields:
      stat = nf90_put_var(ncid,freqid,f,start=[1],count=[om])
      stat = nf90_put_var(ncid,wlid,k(:,ispec(nn)),start=[1],count=[om])
      stat = nf90_put_var(ncid,thetaid,th,start=[1],count=[pm])
      stat = nf90_put_var(ncid,lon_scalarid,lon(mspec(nn),nspec(nn)))
      stat = nf90_put_var(ncid,lat_scalarid,lat(mspec(nn),nspec(nn)))
      stat = nf90_put_var(ncid,sdvid,sdv(:,ispec(nn)),start=[1],count=[om])
      stat = nf90_put_var(ncid,sbfid,sbf(:,ispec(nn)),start=[1],count=[om])

    else

      ! open file for writing:
      spectrumoutputfile = 'output/umwmspc_'//trim(spectrumid(nn))//'_'//savetimestr//'.nc'

      stat = nf90_open(trim(spectrumoutputfile),nf90_write,ncid)

      stat = nf90_inq_varid(ncid,'f',specid)
      stat = nf90_inq_varid(ncid,'sin',sinid)
      stat = nf90_inq_varid(ncid,'sds',sdsid)
      stat = nf90_inq_varid(ncid,'sdt',sdtid)
      stat = nf90_inq_varid(ncid,'snl',snlid)
      stat = nf90_inq_varid(ncid,'wspd',wspdid)
      stat = nf90_inq_varid(ncid,'wdir',wdirid)

    end if

    ! fill in variables:
    stat = nf90_put_var(ncid,specid,e(:,:,ispec(nn)),start=[1,1,counter],count=[om,pm,1])
    stat = nf90_put_var(ncid,sinid,ssin(:,:,ispec(nn)),start=[1,1,counter],count=[om,pm,1])
    stat = nf90_put_var(ncid,sdsid,sds(:,:,ispec(nn)),start=[1,1,counter],count=[om,pm,1])
    stat = nf90_put_var(ncid,sdtid,sdt(:,ispec(nn)),start=[1,counter],count=[om,1])
    stat = nf90_put_var(ncid,snlid,snl(:,:,ispec(nn)),start=[1,1,counter],count=[om,pm,1])

    wspdtmp = wspd(ispec(nn))
    wdirtmp = wdir(ispec(nn))

    stat = nf90_put_var(ncid,wspdid,wspdtmp,start=[1,counter])
    stat = nf90_put_var(ncid,wdirid,wdirtmp,start=[1,counter])

    ! close file:
    stat = nf90_close(ncid)

    write(*,'(a)')'umwm: output_spectrum_nc: spectrum written to '//trim(spectrumoutputfile)

  end if
end do

counter  = counter+1
firstrun = .false.

end subroutine output_spectrum_nc


subroutine output_grid_nc(timestr)
! Writes out model gridded output in a netcdf format
use umwm_stokes,only:depth,lm,us,vs,ds

character(19),intent(in) :: timestr

character(19) :: timestrnew

integer :: stat
integer :: ncid
integer :: xdimid,ydimid,zdimid,fdimid,thdimid,tdimid
integer :: lonid,latid,maskid,depthid,swhid,mwpid
integer :: freqid,thetaid
integer :: wspdid,wdirid
integer :: rhoaid,rhowid
integer :: ficeid
integer :: psimid
integer :: zid,usid,vsid,dsid
integer :: momxid,momyid
integer :: cgmxxid,cgmxyid,cgmyyid
integer :: tfdx1id,tfdy1id
integer :: tfdx2id,tfdy2id
integer :: tfdx3id,tfdy3id
integer :: epsx_atmid, epsy_atmid
integer :: epsx_ocnid, epsy_ocnid
integer :: taux_formid,tauy_formid
integer :: taux_skinid,tauy_skinid
integer :: taux_diagid,tauy_diagid
integer :: taux_ocnid,tauy_ocnid
integer :: taux_botid,tauy_botid
integer :: taux_snlid,tauy_snlid
integer :: tailatmxid,tailatmyid
integer :: tailocnxid,tailocnyid
integer :: dwdid,dwlid,dwpid,mwdid,mwlid,ucid,vcid
integer :: cdid,mssid,ustid,sheltid
integer :: dcpid,dcp0id,dcgid,dcg0id
integer :: physics_time_stepid

integer :: l

real :: output_field(mm,nm)

timestrnew = timestr
timestrnew(11:11) = '_'

! super boring, boiler-plate code follows.

if(nproc == 0)then

  stat = nf90_create('output/umwmout_'//timestrnew//'.nc',nf90_clobber,ncid)

  stat = nf90_def_dim(ncid,'x',mm,xdimid)
  stat = nf90_def_dim(ncid,'y',nm,ydimid)
  stat = nf90_def_dim(ncid,'f',om,fdimid)
  stat = nf90_def_dim(ncid,'th',pm,thdimid)
  stat = nf90_def_dim(ncid,'time',nf90_unlimited,tdimid)

  if(stokes)then

    stat = nf90_def_dim(ncid,'z',lm,zdimid)

    stat = nf90_def_var(ncid,'z',nf90_float,[zdimid],zid)
    stat = nf90_put_att(ncid,zid,name='description',values='depth')
    stat = nf90_put_att(ncid,zid,name='units',values='m')

    stat = nf90_def_var(ncid,'u_stokes',nf90_float,[xdimid,ydimid,zdimid,tdimid],usid)
    stat = nf90_put_att(ncid,usid,name='description',values='stokes drift x-component')
    stat = nf90_put_att(ncid,usid,name='units',values='m/s')

    stat = nf90_def_var(ncid,'v_stokes',nf90_float,[xdimid,ydimid,zdimid,tdimid],vsid)
    stat = nf90_put_att(ncid,vsid,name='description',values='stokes drift y-component')
    stat = nf90_put_att(ncid,vsid,name='units',values='m/s')

    stat = nf90_def_var(ncid,'d_stokes',nf90_float,[xdimid,ydimid,tdimid],dsid)
    stat = nf90_put_att(ncid,dsid,name='description',values='stokes drift e-folding depth')
    stat = nf90_put_att(ncid,dsid,name='units',values='m')

  end if

  stat = nf90_def_var(ncid,'frequency',nf90_float,[fdimid],freqid)
  stat = nf90_put_att(ncid,freqid,name='description',values='frequency')
  stat = nf90_put_att(ncid,freqid,name='units',values='hz')

  stat = nf90_def_var(ncid,'theta',nf90_float,[thdimid],thetaid)
  stat = nf90_put_att(ncid,thetaid,name='description',values='directions')
  stat = nf90_put_att(ncid,thetaid,name='units',values='rad')

  stat = nf90_def_var(ncid,'lon',nf90_float,[xdimid,ydimid,tdimid],lonid)
  stat = nf90_put_att(ncid,lonid,name='description',values='longitude')
  stat = nf90_put_att(ncid,lonid,name='units',values='degrees east')

  stat = nf90_def_var(ncid,'lat',nf90_float,[xdimid,ydimid,tdimid],latid)
  stat = nf90_put_att(ncid,latid,name='description',values='latitude')
  stat = nf90_put_att(ncid,latid,name='units',values='degrees north')

  stat = nf90_def_var(ncid,'seamask',nf90_int,[xdimid,ydimid,tdimid],maskid)
  stat = nf90_put_att(ncid,maskid,name='description',values='seamask')
  stat = nf90_put_att(ncid,maskid,name='units',values='non-dimensional')

  stat = nf90_def_var(ncid,'depth',nf90_float,[xdimid,ydimid,tdimid],depthid)
  stat = nf90_put_att(ncid,depthid,name='description',values='ocean depth')
  stat = nf90_put_att(ncid,depthid,name='units',values='m')

  stat = nf90_def_var(ncid,'wspd',nf90_float,[xdimid,ydimid,tdimid],wspdid)
  stat = nf90_put_att(ncid,wspdid,name='description',values='wind speed')
  stat = nf90_put_att(ncid,wspdid,name='units',values='m/s')

  stat = nf90_def_var(ncid,'wdir',nf90_float,[xdimid,ydimid,tdimid],wdirid)
  stat = nf90_put_att(ncid,wdirid,name='description',values='wind direction')
  stat = nf90_put_att(ncid,wdirid,name='units',values='rad')

  stat = nf90_def_var(ncid,'uc',nf90_float,[xdimid,ydimid,tdimid],ucid)
  stat = nf90_put_att(ncid,ucid,name='description',values='ocean current, x-component')
  stat = nf90_put_att(ncid,ucid,name='units',values='m/s')

  stat = nf90_def_var(ncid,'vc',nf90_float,[xdimid,ydimid,tdimid],vcid)
  stat = nf90_put_att(ncid,vcid,name='description',values='ocean current, y-component')
  stat = nf90_put_att(ncid,vcid,name='units',values='m/s')

  stat = nf90_def_var(ncid,'rhoa',nf90_float,[xdimid,ydimid,tdimid],rhoaid)
  stat = nf90_put_att(ncid,rhoaid,name='description',values='air density')
  stat = nf90_put_att(ncid,rhoaid,name='units',values='kg/m^3')

  stat = nf90_def_var(ncid,'rhow',nf90_float,[xdimid,ydimid,tdimid],rhowid)
  stat = nf90_put_att(ncid,rhowid,name='description',values='water density')
  stat = nf90_put_att(ncid,rhowid,name='units',values='kg/m^3')

  stat = nf90_def_var(ncid,'fice',nf90_float,[xdimid,ydimid,tdimid],ficeid)
  stat = nf90_put_att(ncid,ficeid,name='description',values='seaice fraction')
  stat = nf90_put_att(ncid,ficeid,name='units',values='non-dimensional')

  stat = nf90_def_var(ncid,'psim',nf90_float,[xdimid,ydimid,tdimid],psimid)
  stat = nf90_put_att(ncid,psimid,name='description',values='universal stability function for momentum')
  stat = nf90_put_att(ncid,psimid,name='units',values='non-dimensional')

  stat = nf90_def_var(ncid,'momx',nf90_float,[xdimid,ydimid,tdimid],momxid)
  stat = nf90_put_att(ncid,momxid,name='description',values='momentum, x-component')
  stat = nf90_put_att(ncid,momxid,name='units',values='kgm/s')

  stat = nf90_def_var(ncid,'momy',nf90_float,[xdimid,ydimid,tdimid],momyid)
  stat = nf90_put_att(ncid,momyid,name='description',values='momentum, y-component')
  stat = nf90_put_att(ncid,momyid,name='units',values='kgm/s')

  stat = nf90_def_var(ncid,'cgmxx',nf90_float,[xdimid,ydimid,tdimid],cgmxxid)
  stat = nf90_put_att(ncid,cgmxxid,name='description',values='cg*momentum, xx-component')
  stat = nf90_put_att(ncid,cgmxxid,name='units',values='kgm^2/s^2')

  stat = nf90_def_var(ncid,'cgmxy',nf90_float,[xdimid,ydimid,tdimid],cgmxyid)
  stat = nf90_put_att(ncid,cgmxyid,name='description',values='cg*momentum, xy-component')
  stat = nf90_put_att(ncid,cgmxyid,name='units',values='kgm^2/s^2')

  stat = nf90_def_var(ncid,'cgmyy',nf90_float,[xdimid,ydimid,tdimid],cgmyyid)
  stat = nf90_put_att(ncid,cgmyyid,name='description',values='cg*momentum, yy-component')
  stat = nf90_put_att(ncid,cgmyyid,name='units',values='kgm^2/s^2')

  stat = nf90_def_var(ncid,'shelt',nf90_float,[xdimid,ydimid,tdimid],sheltid)
  stat = nf90_put_att(ncid,sheltid,name='description',values='sheltering coefficient')
  stat = nf90_put_att(ncid,sheltid,name='units',values='non-dimensional')

  stat = nf90_def_var(ncid,'epsx_atm',nf90_float,[xdimid,ydimid,tdimid],epsx_atmid)
  stat = nf90_put_att(ncid,epsx_atmid,name='description',values='wave energy growth flux, x-component')
  stat = nf90_put_att(ncid,epsx_atmid,name='units',values='kg/s^3')

  stat = nf90_def_var(ncid,'epsy_atm',nf90_float,[xdimid,ydimid,tdimid],epsy_atmid)
  stat = nf90_put_att(ncid,epsy_atmid,name='description',values='wave energy growth flux, y-component')
  stat = nf90_put_att(ncid,epsy_atmid,name='units',values='kg/s^3')

  stat = nf90_def_var(ncid,'epsx_ocn',nf90_float,[xdimid,ydimid,tdimid],epsx_ocnid)
  stat = nf90_put_att(ncid,epsx_ocnid,name='description',values='wave energy dissipation flux, x-component')
  stat = nf90_put_att(ncid,epsx_ocnid,name='units',values='kg/s^3')

  stat = nf90_def_var(ncid,'epsy_ocn',nf90_float,[xdimid,ydimid,tdimid],epsy_ocnid)
  stat = nf90_put_att(ncid,epsy_ocnid,name='description',values='wave energy dissipation flux, y-component')
  stat = nf90_put_att(ncid,epsy_ocnid,name='units',values='kg/s^3')

  stat = nf90_def_var(ncid,'taux_form',nf90_float,[xdimid,ydimid,tdimid],taux_formid)
  stat = nf90_put_att(ncid,taux_formid,name='description',values='form drag, x-component')
  stat = nf90_put_att(ncid,taux_formid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_form',nf90_float,[xdimid,ydimid,tdimid],tauy_formid)
  stat = nf90_put_att(ncid,tauy_formid,name='description',values='form drag, y-component')
  stat = nf90_put_att(ncid,tauy_formid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_form_1',nf90_float,[xdimid,ydimid,tdimid],tfdx1id)
  stat = nf90_put_att(ncid,tfdx1id,name='description',values='form drag, part 1, x-component')
  stat = nf90_put_att(ncid,tfdx1id,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_form_1',nf90_float,[xdimid,ydimid,tdimid],tfdy1id)
  stat = nf90_put_att(ncid,tfdy1id,name='description',values='form drag, part 1, y-component')
  stat = nf90_put_att(ncid,tfdy1id,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_form_2',nf90_float,[xdimid,ydimid,tdimid],tfdx2id)
  stat = nf90_put_att(ncid,tfdx2id,name='description',values='form drag, part 2, x-component')
  stat = nf90_put_att(ncid,tfdx2id,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_form_2',nf90_float,[xdimid,ydimid,tdimid],tfdy2id)
  stat = nf90_put_att(ncid,tfdy2id,name='description',values='form drag, part 2, y-component')
  stat = nf90_put_att(ncid,tfdy2id,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_form_3',nf90_float,[xdimid,ydimid,tdimid],tfdx3id)
  stat = nf90_put_att(ncid,tfdx3id,name='description',values='form drag, part 3, x-component')
  stat = nf90_put_att(ncid,tfdx3id,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_form_3',nf90_float,[xdimid,ydimid,tdimid],tfdy3id)
  stat = nf90_put_att(ncid,tfdy3id,name='description',values='form drag, part 3, y-component')
  stat = nf90_put_att(ncid,tfdy3id,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_skin',nf90_float,[xdimid,ydimid,tdimid],taux_skinid)
  stat = nf90_put_att(ncid,taux_skinid,name='description',values='skin drag, x-component')
  stat = nf90_put_att(ncid,taux_skinid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_skin',nf90_float,[xdimid,ydimid,tdimid],tauy_skinid)
  stat = nf90_put_att(ncid,tauy_skinid,name='description',values='skin drag, y-component')
  stat = nf90_put_att(ncid,tauy_skinid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_diag',nf90_float,[xdimid,ydimid,tdimid],taux_diagid)
  stat = nf90_put_att(ncid,taux_diagid,name='description',values='diagnostic form drag, x-component')
  stat = nf90_put_att(ncid,taux_diagid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_diag',nf90_float,[xdimid,ydimid,tdimid],tauy_diagid)
  stat = nf90_put_att(ncid,tauy_diagid,name='description',values='diagnostic form drag, y-component')
  stat = nf90_put_att(ncid,tauy_diagid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_ocn',nf90_float,[xdimid,ydimid,tdimid],taux_ocnid)
  stat = nf90_put_att(ncid,taux_ocnid,name='description',&
                      values='momentum flux from breaking waves to ocean top, x-component')
  stat = nf90_put_att(ncid,taux_ocnid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_ocn',nf90_float,[xdimid,ydimid,tdimid],tauy_ocnid)
  stat = nf90_put_att(ncid,tauy_ocnid,name='description',&
                      values='momentum flux from breaking waves to ocean top, y-component')
  stat = nf90_put_att(ncid,tauy_ocnid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_bot',nf90_float,[xdimid,ydimid,tdimid],taux_botid)
  stat = nf90_put_att(ncid,taux_botid,name='description',&
                      values='momentum flux from waves to ocean bottom, x-component')
  stat = nf90_put_att(ncid,taux_botid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_bot',nf90_float,[xdimid,ydimid,tdimid],tauy_botid)
  stat = nf90_put_att(ncid,tauy_botid,name='description',&
                      values='momentum flux from waves to ocean bottom, y-component')
  stat = nf90_put_att(ncid,tauy_botid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_snl',nf90_float,[xdimid,ydimid,tdimid],taux_snlid)
  stat = nf90_put_att(ncid,taux_snlid,name='description',&
                      values='momentum flux due to snl, x-component')
  stat = nf90_put_att(ncid,taux_snlid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_snl',nf90_float,[xdimid,ydimid,tdimid],tauy_snlid)
  stat = nf90_put_att(ncid,tauy_snlid,name='description',&
                      values='momentum flux due to snl, y-component')
  stat = nf90_put_att(ncid,tauy_snlid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tailatmx',nf90_float,[xdimid,ydimid,tdimid],tailatmxid)
  stat = nf90_put_att(ncid,tailatmxid,name='description',&
                      values='atmosphere tail stress part, x-component')
  stat = nf90_put_att(ncid,tailatmxid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tailatmy',nf90_float,[xdimid,ydimid,tdimid],tailatmyid)
  stat = nf90_put_att(ncid,tailatmyid,name='description',&
                      values='atmosphere tail stress part, y-component')
  stat = nf90_put_att(ncid,tailatmyid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tailocnx',nf90_float,[xdimid,ydimid,tdimid],tailocnxid)
  stat = nf90_put_att(ncid,tailocnxid,name='description',&
                      values='ocean tail stress part, x-component')
  stat = nf90_put_att(ncid,tailocnxid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tailocny',nf90_float,[xdimid,ydimid,tdimid],tailocnyid)
  stat = nf90_put_att(ncid,tailocnyid,name='description',&
                      values='ocean tail stress part, y-component')
  stat = nf90_put_att(ncid,tailocnyid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'cd',nf90_float,[xdimid,ydimid,tdimid],cdid)
  stat = nf90_put_att(ncid,cdid,name='description',values='drag coefficient of air')
  stat = nf90_put_att(ncid,cdid,name='units',values='non-dimensional')

  stat = nf90_def_var(ncid,'ust',nf90_float,[xdimid,ydimid,tdimid],ustid)
  stat = nf90_put_att(ncid,ustid,name='description',values='friction velocity of air')
  stat = nf90_put_att(ncid,ustid,name='units',values='m/s')

  stat = nf90_def_var(ncid,'swh',nf90_float,[xdimid,ydimid,tdimid],swhid)
  stat = nf90_put_att(ncid,swhid,name='description',values='significant wave height')
  stat = nf90_put_att(ncid,swhid,name='units',values='m')

  stat = nf90_def_var(ncid,'mss',nf90_float,[xdimid,ydimid,tdimid],mssid)
  stat = nf90_put_att(ncid,mssid,name='description',values='mean-squared slope')
  stat = nf90_put_att(ncid,mssid,name='units',values='non-dimensional')

  stat = nf90_def_var(ncid,'mwp',nf90_float,[xdimid,ydimid,tdimid],mwpid)
  stat = nf90_put_att(ncid,mwpid,name='description',values='mean wave period')
  stat = nf90_put_att(ncid,mwpid,name='units',values='s')

  stat = nf90_def_var(ncid,'mwl',nf90_float,[xdimid,ydimid,tdimid],mwlid)
  stat = nf90_put_att(ncid,mwlid,name='description',values='mean wavelength')
  stat = nf90_put_att(ncid,mwlid,name='units',values='m')

  stat = nf90_def_var(ncid,'mwd',nf90_float,[xdimid,ydimid,tdimid],mwdid)
  stat = nf90_put_att(ncid,mwdid,name='description',values='mean wave direction')
  stat = nf90_put_att(ncid,mwdid,name='units',values='rad')

  stat = nf90_def_var(ncid,'dwp',nf90_float,[xdimid,ydimid,tdimid],dwpid)
  stat = nf90_put_att(ncid,dwpid,name='description',values='dominant wave period')
  stat = nf90_put_att(ncid,dwpid,name='units',values='s')

  stat = nf90_def_var(ncid,'dwl',nf90_float,[xdimid,ydimid,tdimid],dwlid)
  stat = nf90_put_att(ncid,dwlid,name='description',values='dominant wavelength')
  stat = nf90_put_att(ncid,dwlid,name='units',values='m')

  stat = nf90_def_var(ncid,'dwd',nf90_float,[xdimid,ydimid,tdimid],dwdid)
  stat = nf90_put_att(ncid,dwdid,name='description',values='dominant wave direction')
  stat = nf90_put_att(ncid,dwdid,name='units',values='rad')

  stat = nf90_def_var(ncid,'dcp0',nf90_float,[xdimid,ydimid,tdimid],dcp0id)
  stat = nf90_put_att(ncid,dcp0id,name='description',values='dominant phase speed, intrinsic')
  stat = nf90_put_att(ncid,dcp0id,name='units',values='m/s')

  stat = nf90_def_var(ncid,'dcg0',nf90_float,[xdimid,ydimid,tdimid],dcg0id)
  stat = nf90_put_att(ncid,dcg0id,name='description',values='dominant group speed, intrinsic')
  stat = nf90_put_att(ncid,dcg0id,name='units',values='m/s')

  stat = nf90_def_var(ncid,'dcp',nf90_float,[xdimid,ydimid,tdimid],dcpid)
  stat = nf90_put_att(ncid,dcpid,name='description',values='dominant phase speed')
  stat = nf90_put_att(ncid,dcpid,name='units',values='m/s')

  stat = nf90_def_var(ncid,'dcg',nf90_float,[xdimid,ydimid,tdimid],dcgid)
  stat = nf90_put_att(ncid,dcgid,name='description',values='dominant group speed')
  stat = nf90_put_att(ncid,dcgid,name='units',values='m/s')

  stat = nf90_def_var(ncid,'physics_time_step',nf90_float,[xdimid,ydimid,tdimid],physics_time_stepid)
  stat = nf90_put_att(ncid,physics_time_stepid,name='description',values='Physics time step')
  stat = nf90_put_att(ncid,physics_time_stepid,name='units',values='s')

  stat = nf90_enddef(ncid)

end if

if(nproc == 0)then

  stat = nf90_put_var(ncid,freqid,f,start=[1],count=[om])
  stat = nf90_put_var(ncid,thetaid,th,start=[1],count=[pm])
  stat = nf90_put_var(ncid,lonid,lon,start=[1,1,1],count=[mm,nm,1])
  stat = nf90_put_var(ncid,latid,lat,start=[1,1,1],count=[mm,nm,1])
  stat = nf90_put_var(ncid,maskid,mask,start=[1,1,1],count=[mm,nm,1])
  stat = nf90_put_var(ncid,depthid,d_2d,start=[1,1,1],count=[mm,nm,1])

end if

call gatherfield(wspd(istart:iend),output_field)
if(nproc == 0)stat = nf90_put_var(ncid,wspdid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(wdir(istart:iend),output_field)
if(nproc == 0)stat = nf90_put_var(ncid,wdirid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(uc(istart:iend),output_field)
if(nproc == 0)stat = nf90_put_var(ncid,ucid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(vc(istart:iend),output_field)
if(nproc == 0)stat = nf90_put_var(ncid,vcid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(rhoa(istart:iend),output_field)
if(nproc == 0)stat = nf90_put_var(ncid,rhoaid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(rhow(istart:iend),output_field)
if(nproc == 0)stat = nf90_put_var(ncid,rhowid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(fice(istart:iend),output_field)
if(nproc == 0)stat = nf90_put_var(ncid,ficeid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(psim(istart:iend),output_field)
if(nproc == 0)stat = nf90_put_var(ncid,psimid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(momx,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,momxid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(momy,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,momyid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(cgmxx,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,cgmxxid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(cgmxy,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,cgmxyid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(cgmyy,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,cgmyyid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(shelt,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,sheltid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(epsx_atm,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,epsx_atmid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(epsy_atm,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,epsy_atmid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(epsx_ocn,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,epsx_ocnid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(epsy_ocn,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,epsy_ocnid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(physics_time_step, output_field)
if(nproc == 0)stat = nf90_put_var(ncid,physics_time_stepid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(taux_form,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,taux_formid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(tauy_form,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tauy_formid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(taux1,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tfdx1id,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(tauy1,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tfdy1id,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(taux2,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tfdx2id,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(tauy2,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tfdy2id,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(taux3,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tfdx3id,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(tauy3,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tfdy3id,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(taux_skin,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,taux_skinid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(tauy_skin,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tauy_skinid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(taux_diag,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,taux_diagid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(tauy_diag,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tauy_diagid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(taux_ocntop,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,taux_ocnid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(tauy_ocntop,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tauy_ocnid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(taux_ocnbot,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,taux_botid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(tauy_ocnbot,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tauy_botid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(taux_snl,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,taux_snlid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(tauy_snl,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tauy_snlid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(tailatmx,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tailatmxid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(tailatmy,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tailatmyid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(tailocnx,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tailocnxid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(tailocny,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,tailocnyid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(cd,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,cdid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(ustar,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,ustid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(ht,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,swhid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(mss,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,mssid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(mwp,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,mwpid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(mwl,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,mwlid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(mwd,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,mwdid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(dwp,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,dwpid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(dwl,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,dwlid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(dwd,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,dwdid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(dcp0,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,dcp0id,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(dcg0,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,dcg0id,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(dcp,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,dcpid,output_field,start=[1,1,1],count=[mm,nm,1])

call gatherfield(dcg,output_field)
if(nproc == 0)stat = nf90_put_var(ncid,dcgid,output_field,start=[1,1,1],count=[mm,nm,1])

if(stokes)then

  if(nproc == 0)stat = nf90_put_var(ncid,zid,depth,start=[1],count=[lm])

  do l=1,lm

    call gatherfield(us(istart:iend,l),output_field)
    if(nproc == 0)stat = nf90_put_var(ncid,usid,output_field,start=[1,1,l,1],count=[mm,nm,1,1])

    call gatherfield(vs(istart:iend,l),output_field)
    if(nproc == 0)stat = nf90_put_var(ncid,vsid,output_field,start=[1,1,l,1],count=[mm,nm,1,1])

  end do

  call gatherfield(ds(istart:iend),output_field)
  if(nproc == 0)stat = nf90_put_var(ncid,dsid,output_field,start=[1,1,1],count=[mm,nm,1])

end if

if(nproc == 0)then
  stat = nf90_close(ncid)
  write(unit=*,fmt='(a)')'umwm: output_nc: output written to output/umwmout_'//timestrnew//'.nc'
end if

end subroutine output_grid_nc


subroutine gatherfield(field, field_mn)
! This subroutine gathers a field on root processor 
! and remaps it on a 2-d array.
#ifdef MPI
use mpi
use umwm_mpi
#endif
use umwm_util,only:remap_i2mn
use, intrinsic :: ieee_arithmetic

real, intent(in) :: field(istart:iend)
real, intent(out) :: field_mn(mm,nm)

real :: field_ii(imm)
integer :: nn
real :: nan

field_ii = ieee_value(nan, ieee_quiet_nan)

#ifdef MPI
if(mpiisblocking)then

  if(nproc/=0)then
    call mpi_send(field(istart:iend),ilen,MPI_REAL,&
                  0,nproc,MPI_COMM_WORLD,ierr)
  else
    do nn=1,mpisize-1
      call mpi_recv(field_ii(istart_(nn):iend_(nn)),ilen_(nn),&
                    MPI_REAL,nn,nn,MPI_COMM_WORLD,status,ierr)
    end do
  end if

else

 ! non-blocking gather:
  call gather_array(field, field_ii(1:im))

end if
#endif

if(nproc == 0)then
  field_ii(istart:iend) = field(istart:iend)
  field_mn = remap_i2mn(field_ii)
end if

end subroutine gatherfield


subroutine nc_check(stat)
! Checks for netcdf errors and if any, print and abort.

integer,intent(in) :: stat

if(stat /= nf90_noerr)then
  write(*,*)'error in netcdf i/o'
  write(*,*)trim(nf90_strerror(stat))
  stop
end if

end subroutine nc_check

end module umwm_io
