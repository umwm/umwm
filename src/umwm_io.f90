module umwm_io
! Provides input/output routines for the wave model
use, intrinsic :: iso_fortran_env, only: real64
use umwm_forcing, only: forcing_type
use umwm_grid, only: grid_type
use umwm_module
use umwm_spectrum, only: spectrum_type
use umwm_config, only: config_type
use datetime_module, only: date2num, datetime, strptime
use netcdf

contains


subroutine output_grid(grid)
! Outputs grid related fields into a netcdf file.
type(grid_type), intent(in) :: grid
integer :: ncid
integer :: xdimid,ydimid
integer :: lonid,latid,dlonid,dlatid,dxid,dyid,arid,maskid,did,nprocid
integer :: xid,yid,curvid

if(nproc == 0)then

  call nc_check(nf90_create('output/umwmout.grid',nf90_clobber,ncid))
  call nc_check(nf90_def_dim(ncid,'x',grid % mm,xdimid))
  call nc_check(nf90_def_dim(ncid,'y',grid % nm,ydimid))
  call nc_check(nf90_def_var(ncid,'lon',NF90_FLOAT,[xdimid,ydimid],lonid))
  call nc_check(nf90_def_var(ncid,'lat',NF90_FLOAT,[xdimid,ydimid],latid))
  call nc_check(nf90_def_var(ncid,'xx',NF90_FLOAT,[xdimid,ydimid],xid))
  call nc_check(nf90_def_var(ncid,'yy',NF90_FLOAT,[xdimid,ydimid],yid))
  call nc_check(nf90_def_var(ncid,'dlon',NF90_FLOAT,[xdimid,ydimid],dlonid))
  call nc_check(nf90_def_var(ncid,'dlat',NF90_FLOAT,[xdimid,ydimid],dlatid))
  call nc_check(nf90_def_var(ncid,'dx',NF90_FLOAT,[xdimid,ydimid],dxid))
  call nc_check(nf90_def_var(ncid,'dy',NF90_FLOAT,[xdimid,ydimid],dyid))
  call nc_check(nf90_def_var(ncid,'curvature',NF90_FLOAT,[xdimid,ydimid],curvid))
  call nc_check(nf90_def_var(ncid,'area',NF90_FLOAT,[xdimid,ydimid],arid))
  call nc_check(nf90_def_var(ncid,'depth',NF90_FLOAT,[xdimid,ydimid],did))
  call nc_check(nf90_def_var(ncid,'seamask',nf90_int,[xdimid,ydimid],maskid))
  call nc_check(nf90_def_var(ncid,'nproc',nf90_int,[xdimid,ydimid],nprocid))
  call nc_check(nf90_enddef(ncid))
  call nc_check(nf90_put_var(ncid,lonid,grid % lon))
  call nc_check(nf90_put_var(ncid,latid,grid % lat))
  call nc_check(nf90_put_var(ncid,xid,grid % x))
  call nc_check(nf90_put_var(ncid,yid,grid % y))
  call nc_check(nf90_put_var(ncid,dlonid,grid % dlon))
  call nc_check(nf90_put_var(ncid,dlatid,grid % dlat))
  call nc_check(nf90_put_var(ncid,dxid,grid % dx_2d))
  call nc_check(nf90_put_var(ncid,dyid,grid % dy_2d))
  call nc_check(nf90_put_var(ncid,curvid,grid % curv))
  call nc_check(nf90_put_var(ncid,arid,grid % ar_2d))
  call nc_check(nf90_put_var(ncid,maskid,grid % mask))
  call nc_check(nf90_put_var(ncid,did,grid % d_2d))
  call nc_check(nf90_put_var(ncid,nprocid,grid % nproc_out))
  call nc_check(nf90_close(ncid))

  write(*,'(a)')'umwm: output_grid: grid definition written in output/umwmout.grid'

end if

end subroutine output_grid


subroutine output_spectrum_nc(config, timestr, spectrum, grid, forcing)
! Writes out model spectrum output in a netcdf format
type(config_type), intent(in) :: config
character(19),intent(in) :: timestr
type(spectrum_type), intent(in) :: spectrum
type(grid_type), intent(in) :: grid
type(forcing_type), intent(in) :: forcing

character(19),save :: savetimestr

character(9999) :: spectrumoutputfile

character(2) :: coord_input

integer,dimension(2) :: xy_coords

integer :: stat
integer :: ncid
integer :: fdimid,thdimid,tdimid,scalarid
integer :: lon_scalarid,lat_scalarid,timeid,wspdid,wdirid
integer :: specid,sinid,sdsid,snlid,freqid,thetaid,wlid
integer :: sdtid,sdvid,sbfid

integer                               :: nn
integer,save                          :: npts = 0
integer,          dimension(999),save :: mspec,nspec,ispec
character(3),dimension(999),save :: cnn
character(40),dimension(999),save :: spectrumid

real :: latspec,lonspec
real :: wspdtmp,wdirtmp
real(real64) :: time_value(1)

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

      if (mspec(npts) < 2 .or. mspec(npts) > config % mm - 1 .or. &
          nspec(npts) < 2 .or. nspec(npts) > config % nm - 1) &
          stop 'umwm: output_nc: error: a requested point ' &
            // 'in spectrum.nml is out of bounds'

    else if (any(coord_input == ['ll', 'LL'])) then

      read(21, *, end=100) lonspec, latspec, spectrumid(npts)
      xy_coords = minloc((lonspec - grid % lon)**2 + (latspec - grid % lat)**2)
      mspec(npts) = xy_coords(1)
      nspec(npts) = xy_coords(2)

    else

      stop 'umwm: output_nc: error: first line in namelists/spectrum.nml must' &
        // 'contain "xy" (or "XY") or "ll" (or "LL").'

    end if

    ispec(npts) = grid % ii(mspec(npts),nspec(npts))
    write(unit=cnn(npts),fmt='(i3)')npts
    if(npts<100)cnn(npts) = '0'//adjustl(cnn(npts))
    if(npts< 10)cnn(npts) = '0'//adjustl(cnn(npts))

  end do

  100 npts=npts-1
  close(unit=21)

end if

do nn=1,npts
  if(ispec(nn) >= grid % istart .and. ispec(nn) <= grid % iend)then

    if(firstrun)then

      ! create output file:
      spectrumoutputfile = 'output/umwmspc_'//trim(spectrumid(nn))//'_'//savetimestr//'.nc'

      stat = nf90_create(trim(spectrumoutputfile),nf90_clobber,ncid)

      ! define dimensions:
      stat = nf90_def_dim(ncid, 'scalar',1, scalarid)
      stat = nf90_def_dim(ncid, 'frequency', spectrum % num_frequencies, fdimid)
      stat = nf90_def_dim(ncid, 'direction', spectrum % num_directions, thdimid)
      stat = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, tdimid)

      ! define variables:
      stat = nf90_def_var(ncid, 'time', NF90_DOUBLE, [tdimid], timeid)
      stat = nf90_def_var(ncid, 'frequency', NF90_FLOAT, [fdimid], freqid)
      stat = nf90_def_var(ncid, 'wavenumber', NF90_FLOAT, [fdimid], wlid)
      stat = nf90_def_var(ncid, 'direction', NF90_FLOAT, [thdimid], thetaid)
      stat = nf90_def_var(ncid, 'longitude', NF90_FLOAT, [scalarid], lon_scalarid)
      stat = nf90_def_var(ncid, 'latitude', NF90_FLOAT, [scalarid], lat_scalarid)
      stat = nf90_def_var(ncid, 'wspd', NF90_FLOAT, [scalarid,tdimid], wspdid)
      stat = nf90_def_var(ncid, 'wdir', NF90_FLOAT, [scalarid,tdimid], wdirid)
      stat = nf90_def_var(ncid, 'F', NF90_FLOAT, [fdimid,thdimid,tdimid], specid)
      stat = nf90_def_var(ncid, 'Sin', NF90_FLOAT, [fdimid,thdimid,tdimid], sinid)
      stat = nf90_def_var(ncid, 'Sds', NF90_FLOAT, [fdimid,thdimid,tdimid], sdsid)
      stat = nf90_def_var(ncid, 'Sdt', NF90_FLOAT, [fdimid,tdimid], sdtid)
      stat = nf90_def_var(ncid, 'Sdv', NF90_FLOAT, [fdimid], sdvid)
      stat = nf90_def_var(ncid, 'Sbf', NF90_FLOAT, [fdimid], sbfid)
      stat = nf90_def_var(ncid, 'Snl', NF90_FLOAT, [fdimid,thdimid,tdimid], snlid)

      call put_time_metadata(ncid, timeid, config % reftimestr)
      call put_description(ncid, freqid, 'Frequency')
      call put_description(ncid, wlid, 'Wavenumber')
      call put_description(ncid, thetaid, 'Direction')
      call put_description(ncid, lon_scalarid, 'Longitude')
      call put_description(ncid, lat_scalarid, 'Latitude')
      call put_description(ncid, wspdid, 'Wind speed')
      call put_description(ncid, wdirid, 'Wind direction')
      call put_description(ncid, specid, 'Surface elevation variance spectrum')
      call put_description(ncid, sinid, 'Wind input/export rate')
      call put_description(ncid, sdsid, 'Wave breaking decay rate')
      call put_description(ncid, sdtid, 'Turbulent dissipation rate')
      call put_description(ncid, sdvid, 'Viscous dissipation rate')
      call put_description(ncid, sbfid, 'Combined bottom friction and bottom percolation rate')
      call put_description(ncid, snlid, 'Nonlinear wave-wave interaction absolute source/sink tendency')

      ! end of definition mode:
      stat = nf90_enddef(ncid)

      ! fill in static fields:
      stat = nf90_put_var(ncid,freqid,spectrum % frequency,start=[1],count=[om])
      stat = nf90_put_var(ncid,wlid,k(:,ispec(nn)),start=[1],count=[om])
      stat = nf90_put_var(ncid,thetaid,spectrum % direction,start=[1],count=[pm])
      stat = nf90_put_var(ncid,lon_scalarid,grid % lon(mspec(nn),nspec(nn)))
      stat = nf90_put_var(ncid,lat_scalarid,grid % lat(mspec(nn),nspec(nn)))
      stat = nf90_put_var(ncid,sdvid,sdv(:,ispec(nn)),start=[1],count=[om])
      stat = nf90_put_var(ncid,sbfid,sbf(:,ispec(nn)),start=[1],count=[om])

    else

      ! open file for writing:
      spectrumoutputfile = 'output/umwmspc_'//trim(spectrumid(nn))//'_'//savetimestr//'.nc'

      stat = nf90_open(trim(spectrumoutputfile),nf90_write,ncid)

      stat = nf90_inq_varid(ncid, 'F', specid)
      stat = nf90_inq_varid(ncid, 'Sin', sinid)
      stat = nf90_inq_varid(ncid, 'Sds', sdsid)
      stat = nf90_inq_varid(ncid, 'Sdt', sdtid)
      stat = nf90_inq_varid(ncid, 'Snl', snlid)
      stat = nf90_inq_varid(ncid, 'time', timeid)
      stat = nf90_inq_varid(ncid, 'wspd', wspdid)
      stat = nf90_inq_varid(ncid, 'wdir', wdirid)

    end if

    ! fill in variables:
    time_value(1) = seconds_since_reference(timestr, config % reftimestr)
    stat = nf90_put_var(ncid,timeid,time_value,start=[counter],count=[1])
    stat = nf90_put_var(ncid,specid,e(:,:,ispec(nn)),start=[1,1,counter],count=[om,pm,1])
    stat = nf90_put_var(ncid,sinid,ssin(:,:,ispec(nn)),start=[1,1,counter],count=[om,pm,1])
    stat = nf90_put_var(ncid,sdsid,sds(:,:,ispec(nn)),start=[1,1,counter],count=[om,pm,1])
    stat = nf90_put_var(ncid,sdtid,sdt(:,ispec(nn)),start=[1,counter],count=[om,1])
    stat = nf90_put_var(ncid,snlid,snl(:,:,ispec(nn)),start=[1,1,counter],count=[om,pm,1])

    wspdtmp = forcing % wspd(ispec(nn))
    wdirtmp = forcing % wdir(ispec(nn))

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


subroutine output_grid_nc(config, timestr, spectrum, grid, forcing)
! Writes out model gridded output in a netcdf format
use umwm_stokes,only:depth,lm,us,vs,ds

type(config_type), intent(in) :: config
character(19),intent(in) :: timestr
type(spectrum_type), intent(in) :: spectrum
type(grid_type), intent(in) :: grid
type(forcing_type), intent(in) :: forcing

character(19) :: timestrnew

integer :: stat
integer :: ncid
integer :: xdimid,ydimid,zdimid,fdimid,thdimid,tdimid
integer :: lonid,latid,timeid,maskid,depthid,swhid,mwpid
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

real :: output_field(grid % mm,grid % nm)
real(real64) :: time_value(1)

timestrnew = timestr
timestrnew(11:11) = '_'

associate(istart => grid % istart, iend => grid % iend)

! super boring, boiler-plate code follows.

if(nproc == 0)then

  stat = nf90_create('output/umwmout_'//timestrnew//'.nc',nf90_clobber,ncid)

  stat = nf90_def_dim(ncid,'x',grid % mm,xdimid)
  stat = nf90_def_dim(ncid,'y',grid % nm,ydimid)
  stat = nf90_def_dim(ncid,'f',spectrum % num_frequencies,fdimid)
  stat = nf90_def_dim(ncid,'th',spectrum % num_directions,thdimid)
  stat = nf90_def_dim(ncid,'time',NF90_UNLIMITED,tdimid)

  stat = nf90_def_var(ncid,'time',NF90_DOUBLE,[tdimid],timeid)
  call put_time_metadata(ncid, timeid, config % reftimestr)

  if(config % stokes)then

    stat = nf90_def_dim(ncid,'z',lm,zdimid)

    stat = nf90_def_var(ncid,'z',NF90_FLOAT,[zdimid],zid)
    call put_description(ncid,zid,'depth')
    stat = nf90_put_att(ncid,zid,name='units',values='m')

    stat = nf90_def_var(ncid,'u_stokes',NF90_FLOAT,[xdimid,ydimid,zdimid,tdimid],usid)
    call put_description(ncid,usid,'stokes drift x-component')
    stat = nf90_put_att(ncid,usid,name='units',values='m/s')

    stat = nf90_def_var(ncid,'v_stokes',NF90_FLOAT,[xdimid,ydimid,zdimid,tdimid],vsid)
    call put_description(ncid,vsid,'stokes drift y-component')
    stat = nf90_put_att(ncid,vsid,name='units',values='m/s')

    stat = nf90_def_var(ncid,'d_stokes',NF90_FLOAT,[xdimid,ydimid,tdimid],dsid)
    call put_description(ncid,dsid,'stokes drift e-folding depth')
    stat = nf90_put_att(ncid,dsid,name='units',values='m')

  end if

  stat = nf90_def_var(ncid,'frequency',NF90_FLOAT,[fdimid],freqid)
  call put_description(ncid,freqid,'frequency')
  stat = nf90_put_att(ncid,freqid,name='units',values='hz')

  stat = nf90_def_var(ncid,'theta',NF90_FLOAT,[thdimid],thetaid)
  call put_description(ncid,thetaid,'directions')
  stat = nf90_put_att(ncid,thetaid,name='units',values='rad')

  stat = nf90_def_var(ncid,'lon',NF90_FLOAT,[xdimid,ydimid,tdimid],lonid)
  call put_description(ncid,lonid,'longitude')
  stat = nf90_put_att(ncid,lonid,name='units',values='degrees east')

  stat = nf90_def_var(ncid,'lat',NF90_FLOAT,[xdimid,ydimid,tdimid],latid)
  call put_description(ncid,latid,'latitude')
  stat = nf90_put_att(ncid,latid,name='units',values='degrees north')

  stat = nf90_def_var(ncid,'seamask',nf90_int,[xdimid,ydimid,tdimid],maskid)
  call put_description(ncid,maskid,'seamask')
  stat = nf90_put_att(ncid,maskid,name='units',values='non-dimensional')

  stat = nf90_def_var(ncid,'depth',NF90_FLOAT,[xdimid,ydimid,tdimid],depthid)
  call put_description(ncid,depthid,'ocean depth')
  stat = nf90_put_att(ncid,depthid,name='units',values='m')

  stat = nf90_def_var(ncid,'wspd',NF90_FLOAT,[xdimid,ydimid,tdimid],wspdid)
  call put_description(ncid,wspdid,'wind speed')
  stat = nf90_put_att(ncid,wspdid,name='units',values='m/s')

  stat = nf90_def_var(ncid,'wdir',NF90_FLOAT,[xdimid,ydimid,tdimid],wdirid)
  call put_description(ncid,wdirid,'wind direction')
  stat = nf90_put_att(ncid,wdirid,name='units',values='rad')

  stat = nf90_def_var(ncid,'uc',NF90_FLOAT,[xdimid,ydimid,tdimid],ucid)
  call put_description(ncid,ucid,'ocean current, x-component')
  stat = nf90_put_att(ncid,ucid,name='units',values='m/s')

  stat = nf90_def_var(ncid,'vc',NF90_FLOAT,[xdimid,ydimid,tdimid],vcid)
  call put_description(ncid,vcid,'ocean current, y-component')
  stat = nf90_put_att(ncid,vcid,name='units',values='m/s')

  stat = nf90_def_var(ncid,'rhoa',NF90_FLOAT,[xdimid,ydimid,tdimid],rhoaid)
  call put_description(ncid,rhoaid,'air density')
  stat = nf90_put_att(ncid,rhoaid,name='units',values='kg/m^3')

  stat = nf90_def_var(ncid,'rhow',NF90_FLOAT,[xdimid,ydimid,tdimid],rhowid)
  call put_description(ncid,rhowid,'water density')
  stat = nf90_put_att(ncid,rhowid,name='units',values='kg/m^3')

  stat = nf90_def_var(ncid,'fice',NF90_FLOAT,[xdimid,ydimid,tdimid],ficeid)
  call put_description(ncid,ficeid,'seaice fraction')
  stat = nf90_put_att(ncid,ficeid,name='units',values='non-dimensional')

  stat = nf90_def_var(ncid,'psim',NF90_FLOAT,[xdimid,ydimid,tdimid],psimid)
  call put_description(ncid,psimid,'universal stability function for momentum')
  stat = nf90_put_att(ncid,psimid,name='units',values='non-dimensional')

  stat = nf90_def_var(ncid,'momx',NF90_FLOAT,[xdimid,ydimid,tdimid],momxid)
  call put_description(ncid,momxid,'momentum, x-component')
  stat = nf90_put_att(ncid,momxid,name='units',values='kgm/s')

  stat = nf90_def_var(ncid,'momy',NF90_FLOAT,[xdimid,ydimid,tdimid],momyid)
  call put_description(ncid,momyid,'momentum, y-component')
  stat = nf90_put_att(ncid,momyid,name='units',values='kgm/s')

  stat = nf90_def_var(ncid,'cgmxx',NF90_FLOAT,[xdimid,ydimid,tdimid],cgmxxid)
  call put_description(ncid,cgmxxid,'cg*momentum, xx-component')
  stat = nf90_put_att(ncid,cgmxxid,name='units',values='kgm^2/s^2')

  stat = nf90_def_var(ncid,'cgmxy',NF90_FLOAT,[xdimid,ydimid,tdimid],cgmxyid)
  call put_description(ncid,cgmxyid,'cg*momentum, xy-component')
  stat = nf90_put_att(ncid,cgmxyid,name='units',values='kgm^2/s^2')

  stat = nf90_def_var(ncid,'cgmyy',NF90_FLOAT,[xdimid,ydimid,tdimid],cgmyyid)
  call put_description(ncid,cgmyyid,'cg*momentum, yy-component')
  stat = nf90_put_att(ncid,cgmyyid,name='units',values='kgm^2/s^2')

  stat = nf90_def_var(ncid,'shelt',NF90_FLOAT,[xdimid,ydimid,tdimid],sheltid)
  call put_description(ncid,sheltid,'sheltering coefficient')
  stat = nf90_put_att(ncid,sheltid,name='units',values='non-dimensional')

  stat = nf90_def_var(ncid,'epsx_atm',NF90_FLOAT,[xdimid,ydimid,tdimid],epsx_atmid)
  call put_description(ncid,epsx_atmid,'wave energy growth flux, x-component')
  stat = nf90_put_att(ncid,epsx_atmid,name='units',values='kg/s^3')

  stat = nf90_def_var(ncid,'epsy_atm',NF90_FLOAT,[xdimid,ydimid,tdimid],epsy_atmid)
  call put_description(ncid,epsy_atmid,'wave energy growth flux, y-component')
  stat = nf90_put_att(ncid,epsy_atmid,name='units',values='kg/s^3')

  stat = nf90_def_var(ncid,'epsx_ocn',NF90_FLOAT,[xdimid,ydimid,tdimid],epsx_ocnid)
  call put_description(ncid,epsx_ocnid,'wave energy dissipation flux, x-component')
  stat = nf90_put_att(ncid,epsx_ocnid,name='units',values='kg/s^3')

  stat = nf90_def_var(ncid,'epsy_ocn',NF90_FLOAT,[xdimid,ydimid,tdimid],epsy_ocnid)
  call put_description(ncid,epsy_ocnid,'wave energy dissipation flux, y-component')
  stat = nf90_put_att(ncid,epsy_ocnid,name='units',values='kg/s^3')

  stat = nf90_def_var(ncid,'taux_form',NF90_FLOAT,[xdimid,ydimid,tdimid],taux_formid)
  call put_description(ncid,taux_formid,'form drag, x-component')
  stat = nf90_put_att(ncid,taux_formid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_form',NF90_FLOAT,[xdimid,ydimid,tdimid],tauy_formid)
  call put_description(ncid,tauy_formid,'form drag, y-component')
  stat = nf90_put_att(ncid,tauy_formid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_form_1',NF90_FLOAT,[xdimid,ydimid,tdimid],tfdx1id)
  call put_description(ncid,tfdx1id,'form drag, part 1, x-component')
  stat = nf90_put_att(ncid,tfdx1id,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_form_1',NF90_FLOAT,[xdimid,ydimid,tdimid],tfdy1id)
  call put_description(ncid,tfdy1id,'form drag, part 1, y-component')
  stat = nf90_put_att(ncid,tfdy1id,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_form_2',NF90_FLOAT,[xdimid,ydimid,tdimid],tfdx2id)
  call put_description(ncid,tfdx2id,'form drag, part 2, x-component')
  stat = nf90_put_att(ncid,tfdx2id,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_form_2',NF90_FLOAT,[xdimid,ydimid,tdimid],tfdy2id)
  call put_description(ncid,tfdy2id,'form drag, part 2, y-component')
  stat = nf90_put_att(ncid,tfdy2id,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_form_3',NF90_FLOAT,[xdimid,ydimid,tdimid],tfdx3id)
  call put_description(ncid,tfdx3id,'form drag, part 3, x-component')
  stat = nf90_put_att(ncid,tfdx3id,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_form_3',NF90_FLOAT,[xdimid,ydimid,tdimid],tfdy3id)
  call put_description(ncid,tfdy3id,'form drag, part 3, y-component')
  stat = nf90_put_att(ncid,tfdy3id,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_skin',NF90_FLOAT,[xdimid,ydimid,tdimid],taux_skinid)
  call put_description(ncid,taux_skinid,'skin drag, x-component')
  stat = nf90_put_att(ncid,taux_skinid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_skin',NF90_FLOAT,[xdimid,ydimid,tdimid],tauy_skinid)
  call put_description(ncid,tauy_skinid,'skin drag, y-component')
  stat = nf90_put_att(ncid,tauy_skinid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_diag',NF90_FLOAT,[xdimid,ydimid,tdimid],taux_diagid)
  call put_description(ncid,taux_diagid,'diagnostic form drag, x-component')
  stat = nf90_put_att(ncid,taux_diagid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_diag',NF90_FLOAT,[xdimid,ydimid,tdimid],tauy_diagid)
  call put_description(ncid,tauy_diagid,'diagnostic form drag, y-component')
  stat = nf90_put_att(ncid,tauy_diagid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_ocn',NF90_FLOAT,[xdimid,ydimid,tdimid],taux_ocnid)
  call put_description(ncid,taux_ocnid,'momentum flux from breaking waves to ocean top, x-component')
  stat = nf90_put_att(ncid,taux_ocnid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_ocn',NF90_FLOAT,[xdimid,ydimid,tdimid],tauy_ocnid)
  call put_description(ncid,tauy_ocnid,'momentum flux from breaking waves to ocean top, y-component')
  stat = nf90_put_att(ncid,tauy_ocnid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_bot',NF90_FLOAT,[xdimid,ydimid,tdimid],taux_botid)
  call put_description(ncid,taux_botid,'momentum flux from waves to ocean bottom, x-component')
  stat = nf90_put_att(ncid,taux_botid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_bot',NF90_FLOAT,[xdimid,ydimid,tdimid],tauy_botid)
  call put_description(ncid,tauy_botid,'momentum flux from waves to ocean bottom, y-component')
  stat = nf90_put_att(ncid,tauy_botid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'taux_snl',NF90_FLOAT,[xdimid,ydimid,tdimid],taux_snlid)
  call put_description(ncid,taux_snlid,'momentum flux due to snl, x-component')
  stat = nf90_put_att(ncid,taux_snlid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tauy_snl',NF90_FLOAT,[xdimid,ydimid,tdimid],tauy_snlid)
  call put_description(ncid,tauy_snlid,'momentum flux due to snl, y-component')
  stat = nf90_put_att(ncid,tauy_snlid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tailatmx',NF90_FLOAT,[xdimid,ydimid,tdimid],tailatmxid)
  call put_description(ncid,tailatmxid,'atmosphere tail stress part, x-component')
  stat = nf90_put_att(ncid,tailatmxid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tailatmy',NF90_FLOAT,[xdimid,ydimid,tdimid],tailatmyid)
  call put_description(ncid,tailatmyid,'atmosphere tail stress part, y-component')
  stat = nf90_put_att(ncid,tailatmyid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tailocnx',NF90_FLOAT,[xdimid,ydimid,tdimid],tailocnxid)
  call put_description(ncid,tailocnxid,'ocean tail stress part, x-component')
  stat = nf90_put_att(ncid,tailocnxid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'tailocny',NF90_FLOAT,[xdimid,ydimid,tdimid],tailocnyid)
  call put_description(ncid,tailocnyid,'ocean tail stress part, y-component')
  stat = nf90_put_att(ncid,tailocnyid,name='units',values='n/m^2')

  stat = nf90_def_var(ncid,'cd',NF90_FLOAT,[xdimid,ydimid,tdimid],cdid)
  call put_description(ncid,cdid,'drag coefficient of air')
  stat = nf90_put_att(ncid,cdid,name='units',values='non-dimensional')

  stat = nf90_def_var(ncid,'ust',NF90_FLOAT,[xdimid,ydimid,tdimid],ustid)
  call put_description(ncid,ustid,'friction velocity of air')
  stat = nf90_put_att(ncid,ustid,name='units',values='m/s')

  stat = nf90_def_var(ncid,'swh',NF90_FLOAT,[xdimid,ydimid,tdimid],swhid)
  call put_description(ncid,swhid,'significant wave height')
  stat = nf90_put_att(ncid,swhid,name='units',values='m')

  stat = nf90_def_var(ncid,'mss',NF90_FLOAT,[xdimid,ydimid,tdimid],mssid)
  call put_description(ncid,mssid,'mean-squared slope')
  stat = nf90_put_att(ncid,mssid,name='units',values='non-dimensional')

  stat = nf90_def_var(ncid,'mwp',NF90_FLOAT,[xdimid,ydimid,tdimid],mwpid)
  call put_description(ncid,mwpid,'mean wave period')
  stat = nf90_put_att(ncid,mwpid,name='units',values='s')

  stat = nf90_def_var(ncid,'mwl',NF90_FLOAT,[xdimid,ydimid,tdimid],mwlid)
  call put_description(ncid,mwlid,'mean wavelength')
  stat = nf90_put_att(ncid,mwlid,name='units',values='m')

  stat = nf90_def_var(ncid,'mwd',NF90_FLOAT,[xdimid,ydimid,tdimid],mwdid)
  call put_description(ncid,mwdid,'mean wave direction')
  stat = nf90_put_att(ncid,mwdid,name='units',values='rad')

  stat = nf90_def_var(ncid,'dwp',NF90_FLOAT,[xdimid,ydimid,tdimid],dwpid)
  call put_description(ncid,dwpid,'dominant wave period')
  stat = nf90_put_att(ncid,dwpid,name='units',values='s')

  stat = nf90_def_var(ncid,'dwl',NF90_FLOAT,[xdimid,ydimid,tdimid],dwlid)
  call put_description(ncid,dwlid,'dominant wavelength')
  stat = nf90_put_att(ncid,dwlid,name='units',values='m')

  stat = nf90_def_var(ncid,'dwd',NF90_FLOAT,[xdimid,ydimid,tdimid],dwdid)
  call put_description(ncid,dwdid,'dominant wave direction')
  stat = nf90_put_att(ncid,dwdid,name='units',values='rad')

  stat = nf90_def_var(ncid,'dcp0',NF90_FLOAT,[xdimid,ydimid,tdimid],dcp0id)
  call put_description(ncid,dcp0id,'dominant phase speed, intrinsic')
  stat = nf90_put_att(ncid,dcp0id,name='units',values='m/s')

  stat = nf90_def_var(ncid,'dcg0',NF90_FLOAT,[xdimid,ydimid,tdimid],dcg0id)
  call put_description(ncid,dcg0id,'dominant group speed, intrinsic')
  stat = nf90_put_att(ncid,dcg0id,name='units',values='m/s')

  stat = nf90_def_var(ncid,'dcp',NF90_FLOAT,[xdimid,ydimid,tdimid],dcpid)
  call put_description(ncid,dcpid,'dominant phase speed')
  stat = nf90_put_att(ncid,dcpid,name='units',values='m/s')

  stat = nf90_def_var(ncid,'dcg',NF90_FLOAT,[xdimid,ydimid,tdimid],dcgid)
  call put_description(ncid,dcgid,'dominant group speed')
  stat = nf90_put_att(ncid,dcgid,name='units',values='m/s')

  stat = nf90_def_var(ncid,'physics_time_step',NF90_FLOAT,[xdimid,ydimid,tdimid],physics_time_stepid)
  call put_description(ncid,physics_time_stepid,'Physics time step')
  stat = nf90_put_att(ncid,physics_time_stepid,name='units',values='s')

  stat = nf90_enddef(ncid)

end if

if(nproc == 0)then

  time_value(1) = seconds_since_reference(timestr, config % reftimestr)
  stat = nf90_put_var(ncid,timeid,time_value,start=[1],count=[1])
  stat = nf90_put_var(ncid,freqid,spectrum % frequency,start=[1],count=[om])
  stat = nf90_put_var(ncid,thetaid,spectrum % direction,start=[1],count=[pm])
  stat = nf90_put_var(ncid,lonid,grid % lon,start=[1,1,1],count=[grid % mm,grid % nm,1])
  stat = nf90_put_var(ncid,latid,grid % lat,start=[1,1,1],count=[grid % mm,grid % nm,1])
  stat = nf90_put_var(ncid,maskid,grid % mask,start=[1,1,1],count=[grid % mm,grid % nm,1])
  stat = nf90_put_var(ncid,depthid,grid % d_2d,start=[1,1,1],count=[grid % mm,grid % nm,1])

end if

call gatherfield(forcing % wspd(istart:iend),output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,wspdid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(forcing % wdir(istart:iend),output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,wdirid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(forcing % uc(istart:iend),output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,ucid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(forcing % vc(istart:iend),output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,vcid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(forcing % rhoa(istart:iend),output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,rhoaid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(forcing % rhow(istart:iend),output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,rhowid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(forcing % fice(istart:iend),output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,ficeid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(psim(istart:iend),output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,psimid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(momx,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,momxid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(momy,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,momyid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(cgmxx,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,cgmxxid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(cgmxy,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,cgmxyid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(cgmyy,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,cgmyyid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(shelt,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,sheltid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(epsx_atm,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,epsx_atmid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(epsy_atm,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,epsy_atmid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(epsx_ocn,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,epsx_ocnid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(epsy_ocn,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,epsy_ocnid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(physics_time_step, output_field, grid)
if(nproc == 0)stat = nf90_put_var(ncid,physics_time_stepid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(taux_form,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,taux_formid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(tauy_form,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tauy_formid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(taux1,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tfdx1id,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(tauy1,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tfdy1id,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(taux2,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tfdx2id,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(tauy2,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tfdy2id,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(taux3,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tfdx3id,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(tauy3,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tfdy3id,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(taux_skin,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,taux_skinid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(tauy_skin,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tauy_skinid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(taux_diag,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,taux_diagid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(tauy_diag,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tauy_diagid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(taux_ocntop,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,taux_ocnid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(tauy_ocntop,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tauy_ocnid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(taux_ocnbot,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,taux_botid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(tauy_ocnbot,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tauy_botid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(taux_snl,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,taux_snlid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(tauy_snl,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tauy_snlid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(tailatmx,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tailatmxid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(tailatmy,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tailatmyid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(tailocnx,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tailocnxid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(tailocny,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,tailocnyid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(cd,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,cdid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(ustar,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,ustid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(ht,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,swhid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(mss,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,mssid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(mwp,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,mwpid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(mwl,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,mwlid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(mwd,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,mwdid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(dwp,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,dwpid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(dwl,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,dwlid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(dwd,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,dwdid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(dcp0,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,dcp0id,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(dcg0,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,dcg0id,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(dcp,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,dcpid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

call gatherfield(dcg,output_field,grid)
if(nproc == 0)stat = nf90_put_var(ncid,dcgid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

if(config % stokes)then

  if(nproc == 0)stat = nf90_put_var(ncid,zid,depth,start=[1],count=[lm])

  do l=1,lm

    call gatherfield(us(istart:iend,l),output_field,grid)
    if(nproc == 0)stat = nf90_put_var(ncid,usid,output_field,start=[1,1,l,1],count=[grid % mm,grid % nm,1,1])

    call gatherfield(vs(istart:iend,l),output_field,grid)
    if(nproc == 0)stat = nf90_put_var(ncid,vsid,output_field,start=[1,1,l,1],count=[grid % mm,grid % nm,1,1])

  end do

  call gatherfield(ds(istart:iend),output_field,grid)
  if(nproc == 0)stat = nf90_put_var(ncid,dsid,output_field,start=[1,1,1],count=[grid % mm,grid % nm,1])

end if

if(nproc == 0)then
  stat = nf90_close(ncid)
  write(unit=*,fmt='(a)')'umwm: output_nc: output written to output/umwmout_'//timestrnew//'.nc'
end if

end associate

end subroutine output_grid_nc


subroutine put_description(ncid, varid, description)
  integer, intent(in) :: ncid
  integer, intent(in) :: varid
  character(*), intent(in) :: description
  integer :: stat

  stat = nf90_put_att(ncid, varid, name='description', values=description)
  stat = nf90_put_att(ncid, varid, name='long_name', values=description)

end subroutine put_description


subroutine put_time_metadata(ncid, timeid, reftimestr)
  integer, intent(in) :: ncid
  integer, intent(in) :: timeid
  character(*), intent(in) :: reftimestr
  integer :: stat

  call put_description(ncid, timeid, 'time')
  stat = nf90_put_att(ncid, timeid, name='cartesian_axis', values='T')
  stat = nf90_put_att(ncid, timeid, name='units', values='seconds since '//trim(reftimestr))
  stat = nf90_put_att(ncid, timeid, name='time_origin', values=trim(reftimestr))

end subroutine put_time_metadata


real(real64) function seconds_since_reference(timestr, reftimestr) result(seconds)
  character(*), intent(in) :: timestr
  character(*), intent(in) :: reftimestr
  type(datetime) :: output_time
  type(datetime) :: reference_time

  output_time = strptime(normalized_timestamp(timestr), '%Y-%m-%d %H:%M:%S', tz=0._real64)
  reference_time = strptime(normalized_timestamp(reftimestr), '%Y-%m-%d %H:%M:%S', tz=0._real64)
  seconds = anint((date2num(output_time) - date2num(reference_time)) * 86400._real64)

end function seconds_since_reference


character(19) function normalized_timestamp(timestr) result(normalized)
  character(*), intent(in) :: timestr

  normalized = timestr
  if (normalized(11:11) == '_') normalized(11:11) = ' '

end function normalized_timestamp


subroutine gatherfield(field, field_mn, grid)
! This subroutine gathers a field on root processor
! and remaps it on a 2-d array.
#ifdef MPI
use mpi
use umwm_mpi, only: gather_array
#endif
use, intrinsic :: ieee_arithmetic

type(grid_type), intent(in) :: grid
real, intent(in) :: field(grid % istart:grid % iend)
real, intent(out) :: field_mn(grid % mm,grid % nm)

real :: field_ii(grid % imm)
real :: nan

#ifdef MPI
integer :: nn
integer :: status(MPI_STATUS_SIZE)
#endif

field_ii = ieee_value(nan, ieee_quiet_nan)

#ifdef MPI
if(mpiisblocking)then

  if(nproc/=0)then
    call mpi_send(field(grid % istart:grid % iend),grid % ilen,MPI_REAL,&
                  0,nproc,MPI_COMM_WORLD,ierr)
  else
    do nn=1,mpisize-1
      call mpi_recv(field_ii(grid % istart_all(nn):grid % iend_all(nn)),grid % ilen_all(nn),&
                    MPI_REAL,nn,nn,MPI_COMM_WORLD,status,ierr)
    end do
  end if

else

 ! non-blocking gather:
  call gather_array(field, field_ii(1:grid % im), grid)

end if
#endif

if(nproc == 0)then
  field_ii(grid % istart:grid % iend) = field(grid % istart:grid % iend)
  field_mn = grid % remap_i2mn(field_ii)
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
