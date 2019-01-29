!
! University of Miami Wave Model (UMWM)
! wrf2umwmin
! Copyright (C) 2012  Milan Curcic
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!======================================================================!
PROGRAM wrf2umwmin
!=======================================================================
!
! DESCRIPTION: Takes a WRF file name as a command line argument, and 
! outputs a UMWM input file based on data in the WRF file. Supported
! WRF files are input (wrfinput*), output (wrfout*) or restart 
! (wrfrst*). Conversion of fields goes as:
!
!     WRF            UMWM        ACTION
!-----------------------------------------------------------------------
!    XLONG        -> lon     - Longitude, copied
!    XLAT         -> lat     - Latitude,  copied
!    LANDMASK     -> seamask - Seamask
!    U10          -> uw      - 10-m wind, x-component, copied
!    V10          -> vw      - 10-m wind, y-component, copied 
!    PSFC, T2, Q2 -> rhoa    - Moist air density, calculated from
!                              surface pressure, temperature and 
!                              specific humidity
! 
! Note: Water density (rhow) and surface current fields (uc,vc)
! are not available from WRF files.
!
! REVISION HISTORY:
!              
! July 2012, 1.0.1, First public release
! 
!=======================================================================
USE netcdf
USE umwm_utils,ONLY:nc_check
!=======================================================================
IMPLICIT NONE

INTEGER :: i,j
INTEGER :: idm,jdm

INTEGER :: arglen,stat

INTEGER :: ncid,xdimid,ydimid
INTEGER :: varid,lonid,latid,maskid,uwid,vwid,rhoaid

CHARACTER(LEN=9999)          :: wrfoutFile
CHARACTER(LEN=29)            :: umwminFile
CHARACTER(LEN=19)            :: time
CHARACTER(LEN=NF90_MAX_NAME) :: dimname

REAL,DIMENSION(:,:),ALLOCATABLE :: lon,lat,u10,v10,t2,q2,psfc,rhoa

INTEGER,DIMENSION(:,:),ALLOCATABLE :: mask

REAL,PARAMETER :: Rd = 286.9
REAL,PARAMETER :: Rw = 461.5

!=======================================================================

! 1) Get command line argument:
CALL GET_COMMAND_ARGUMENT(1,VALUE=wrfoutFile,STATUS=stat)

! Error checking:
IF(stat/=0)THEN
  WRITE(*,*)'wrf2umwmin: Error: Cannot retrieve argument'
  STOP
ENDIF

time = wrfoutFile(LEN_TRIM(wrfoutFile)-18:)

umwminFile = 'umwmin_'//time//'.nc'

!=======================================================================

! 2) Read input file:

CALL nc_check(nf90_open(TRIM(wrfoutFile),nf90_nowrite,ncid))
CALL nc_check(nf90_inq_dimid(ncid,'west_east',xdimid))
CALL nc_check(nf90_inq_dimid(ncid,'south_north',ydimid))
CALL nc_check(nf90_inquire_dimension(ncid,xdimid,dimname,idm))
CALL nc_check(nf90_inquire_dimension(ncid,ydimid,dimname,jdm))

ALLOCATE(lon(idm,jdm),lat(idm,jdm),mask(idm,jdm),u10(idm,jdm),v10(idm,jdm))
ALLOCATE(psfc(idm,jdm),t2(idm,jdm),q2(idm,jdm),rhoa(idm,jdm))

CALL nc_check(nf90_inq_varid(ncid,'XLONG',varid))
CALL nc_check(nf90_get_var(ncid,varid,lon,start=[1,1,1],count=[idm,jdm,1]))

CALL nc_check(nf90_inq_varid(ncid,'XLAT',varid))
CALL nc_check(nf90_get_var(ncid,varid,lat,start=[1,1,1],count=[idm,jdm,1]))

CALL nc_check(nf90_inq_varid(ncid,'LANDMASK',varid))
CALL nc_check(nf90_get_var(ncid,varid,mask,start=[1,1,1],count=[idm,jdm,1]))

CALL nc_check(nf90_inq_varid(ncid,'U10',varid))
CALL nc_check(nf90_get_var(ncid,varid,u10,start=[1,1,1],count=[idm,jdm,1]))

CALL nc_check(nf90_inq_varid(ncid,'V10',varid))
CALL nc_check(nf90_get_var(ncid,varid,v10,start=[1,1,1],count=[idm,jdm,1]))

CALL nc_check(nf90_inq_varid(ncid,'PSFC',varid))
CALL nc_check(nf90_get_var(ncid,varid,psfc,start=[1,1,1],count=[idm,jdm,1]))

CALL nc_check(nf90_inq_varid(ncid,'T2',varid))
CALL nc_check(nf90_get_var(ncid,varid,t2,start=[1,1,1],count=[idm,jdm,1]))

CALL nc_check(nf90_inq_varid(ncid,'Q2',varid))
CALL nc_check(nf90_get_var(ncid,varid,q2,start=[1,1,1],count=[idm,jdm,1]))

CALL nc_check(nf90_close(ncid))

!=======================================================================

! Compute air density at the surface:
rhoa = psfc*(1+q2)/(t2*(Rd+Rw*q2))

!=======================================================================

! 3) Generate new umwm.grid file:

CALL nc_check(nf90_create(umwminFile,nf90_clobber,ncid))
CALL nc_check(nf90_def_dim(ncid,'x',idm,xdimid))
CALL nc_check(nf90_def_dim(ncid,'y',jdm,ydimid))
CALL nc_check(nf90_def_var(ncid,'lon',nf90_float,[xdimid,ydimid],lonid))
CALL nc_check(nf90_def_var(ncid,'lat',nf90_float,[xdimid,ydimid],latid))
CALL nc_check(nf90_def_var(ncid,'seamask',nf90_int,[xdimid,ydimid],maskid))
CALL nc_check(nf90_def_var(ncid,'uw',nf90_float,[xdimid,ydimid],uwid))
CALL nc_check(nf90_def_var(ncid,'vw',nf90_float,[xdimid,ydimid],vwid))
CALL nc_check(nf90_def_var(ncid,'rhoa',nf90_float,[xdimid,ydimid],rhoaid))
CALL nc_check(nf90_put_att(ncid,NF90_GLOBAL,name='history',values='Generated by wrfout2umwmin'))
CALL nc_check(nf90_enddef(ncid))
CALL nc_check(nf90_put_var(ncid,lonid,lon))
CALL nc_check(nf90_put_var(ncid,latid,lat))
CALL nc_check(nf90_put_var(ncid,maskid,1-mask))
CALL nc_check(nf90_put_var(ncid,uwid,u10))
CALL nc_check(nf90_put_var(ncid,vwid,v10))
CALL nc_check(nf90_put_var(ncid,rhoaid,rhoa))
CALL nc_check(nf90_close(ncid))

!=======================================================================
ENDPROGRAM wrf2umwmin
