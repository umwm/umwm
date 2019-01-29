!
! University of Miami Wave Model (UMWM)
! wrf2umwmgrid
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
PROGRAM wrf2umwmgrid
!=======================================================================
!
! DESCRIPTION: Utility program that reads lon, lat and mask fields
! from a WRF file, and writes them into a umwm.grid file which can be
! read from umwm_topogen for topography generation.
!              
! REVISION HISTORY:
!              
! Oct 2013, 1.1.0, Added the option to refine or coarsen grid by an
!                  integer factor
! Jun 2012, 1.0.1, First public release
! 
!=======================================================================
USE netcdf
USE umwm_utils,ONLY:nc_check
USE grid_utils,ONLY:refineField,coarsenField
!=======================================================================
IMPLICIT NONE

INTEGER :: i,j
INTEGER :: idm_src,jdm_src
INTEGER :: idm_dst,jdm_dst

INTEGER :: stat
INTEGER :: ncid,xdimid,ydimid,lonid,latid,maskid

CHARACTER(LEN=9999) :: wrfInputFile
CHARACTER(LEN=2)    :: refinementFactorStr

CHARACTER(LEN=NF90_MAX_NAME) :: dimname

REAL,DIMENSION(:,:),ALLOCATABLE :: lon_src,lat_src
REAL,DIMENSION(:,:),ALLOCATABLE :: lon_dst,lat_dst

INTEGER,DIMENSION(:,:),ALLOCATABLE :: mask_src
REAL,   DIMENSION(:,:),ALLOCATABLE :: mask_dst

INTEGER :: refinementFactor = 1

!=======================================================================

! 1) Get command line argument:
CALL GET_COMMAND_ARGUMENT(1,VALUE=wrfInputFile,STATUS=stat)

! Error checking:
IF(stat /=0 )THEN
  WRITE(*,*)'wrf2umwmgrid: Error: Cannot retrieve argument 1'
  STOP
ELSE
  WRITE(*,*)'wrf2umwmgrid: Processing file '//TRIM(wrfInputFile)
ENDIF

CALL GET_COMMAND_ARGUMENT(2,VALUE=refinementFactorStr,STATUS=stat)
IF(stat == 0)THEN
  READ(UNIT=refinementFactorStr,FMT=*)refinementFactor
ENDIF

IF(refinementFactor == 0)THEN
  WRITE(*,*)'wrf2umwmgrid: Error: refinementFactor must be a non-zero integer'
  STOP
ENDIF

IF(refinementFactor >= 1)THEN
  WRITE(*,FMT='(A,I2)')'wrf2umwmgrid: Refinement factor: ',refinementFactor
ELSEIF(refinementFactor < 0)THEN
  WRITE(*,FMT='(A,I2)')'wrf2umwmgrid: Refinement factor: 1 /',-refinementFactor
ENDIF

!=======================================================================

! 2) Read input grid:
WRITE(*,*)'wrf2umwmgrid: Reading input (WRF) grid.'

CALL nc_check(nf90_open(TRIM(wrfInputFile),nf90_nowrite,ncid))
CALL nc_check(nf90_inq_dimid(ncid,'west_east',xdimid))
CALL nc_check(nf90_inq_dimid(ncid,'south_north',ydimid))
CALL nc_check(nf90_inquire_dimension(ncid,xdimid,dimname,idm_src))
CALL nc_check(nf90_inquire_dimension(ncid,ydimid,dimname,jdm_src))

WRITE(*,*)'wrf2umwmgrid: WRF grid idm/jdm:',idm_src,jdm_src

ALLOCATE(lon_src(idm_src,jdm_src))
ALLOCATE(lat_src(idm_src,jdm_src))
ALLOCATE(mask_src(idm_src,jdm_src))

CALL nc_check(nf90_inq_varid(ncid,'XLONG',lonid))
CALL nc_check(nf90_inq_varid(ncid,'XLAT',latid))
CALL nc_check(nf90_inq_varid(ncid,'XLAND',maskid))
CALL nc_check(nf90_get_var(ncid,lonid,lon_src))
CALL nc_check(nf90_get_var(ncid,latid,lat_src))
CALL nc_check(nf90_get_var(ncid,maskid,mask_src))
CALL nc_check(nf90_close(ncid))

WHERE(lon_src < 0)lon_src = lon_src+360

WRITE(*,*)'wrf2umwmgrid: WRF grid longitude min/max:',&
          MINVAL(lon_src),MAXVAL(lon_src)
WRITE(*,*)'wrf2umwmgrid: WRF grid latitude min/max:',&
          MINVAL(lat_src),MAXVAL(lat_src)
WRITE(*,*)'wrf2umwmgrid: WRF grid number of sea points:',&
          COUNT(mask_src == 0)
WRITE(*,*)'wrf2umwmgrid: WRF grid number of land points:',&
          COUNT(mask_src == 1)

! Adjust land/sea values:
mask_src = mask_src-1

!=======================================================================

! 3) Refine or coarsen fields if requested:

IF(refinementFactor > 1)THEN

  idm_dst = (idm_src-1)*refinementFactor+1
  jdm_dst = (jdm_src-1)*refinementFactor+1

ELSEIF(refinementFactor < 0)THEN

  idm_dst = (idm_src-1)/FLOAT(ABS(refinementFactor))+1
  jdm_dst = (jdm_src-1)/FLOAT(ABS(refinementFactor))+1

ELSE
  
  idm_dst = idm_src
  jdm_dst = jdm_src

ENDIF

WRITE(*,*)'wrf2umwmgrid: UMWM grid idm/jdm:',idm_dst,jdm_dst

ALLOCATE(lon_dst(idm_dst,jdm_dst))
ALLOCATE(lat_dst(idm_dst,jdm_dst))
ALLOCATE(mask_dst(idm_dst,jdm_dst))

lon_dst  = 0
lat_dst  = 0
mask_dst = 0

IF(refinementFactor > 1)THEN ! refine:

  CALL refineField(lon_src,lon_dst,refinementFactor)
  CALL refineField(lat_src,lat_dst,refinementFactor)
  CALL refineField(FLOAT(mask_src),mask_dst,refinementFactor)

ELSEIF(refinementFactor == 1)THEN ! copy:

  lon_dst  = lon_src
  lat_dst  = lat_src
  mask_dst = mask_src

ELSE ! coarsen:

  CALL coarsenField(lon_src,lon_dst,-refinementFactor)
  CALL coarsenField(lat_src,lat_dst,-refinementFactor)
  CALL coarsenField(FLOAT(mask_src),mask_dst,-refinementFactor)

ENDIF

WRITE(*,*)'wrf2umwmgrid: UMWM grid longitude min/max:',&
          MINVAL(lon_dst),MAXVAL(lon_dst)
WRITE(*,*)'wrf2umwmgrid: UMWM grid latitude min/max:',&
          MINVAL(lat_dst),MAXVAL(lat_dst)
WRITE(*,*)'wrf2umwmgrid: UMWM grid number of sea points:',&
          COUNT(NINT(mask_dst) == 1)
WRITE(*,*)'wrf2umwmgrid: UMWM grid number of land points:',&
          COUNT(NINT(mask_dst) == 0)

WHERE(lon_dst > 180)lon_dst = lon_dst-360

!=======================================================================

! 4) Generate new umwm.grid file:
WRITE(*,*)'wrf2umwmgrid: Writing umwm.grid file.'

CALL nc_check(nf90_create('umwm.grid',nf90_clobber,ncid))
CALL nc_check(nf90_def_dim(ncid,'x',idm_dst,xdimid))
CALL nc_check(nf90_def_dim(ncid,'y',jdm_dst,ydimid))
CALL nc_check(nf90_def_var(ncid,'lon',nf90_float,[xdimid,ydimid],lonid))
CALL nc_check(nf90_def_var(ncid,'lat',nf90_float,[xdimid,ydimid],latid))
CALL nc_check(nf90_def_var(ncid,'seamask',nf90_int,[xdimid,ydimid],maskid))
CALL nc_check(nf90_put_att(ncid,NF90_GLOBAL,name='history',values='Generated by wrf2umwmgrid'))
CALL nc_check(nf90_enddef(ncid))
CALL nc_check(nf90_put_var(ncid,lonid,lon_dst))
CALL nc_check(nf90_put_var(ncid,latid,lat_dst))
CALL nc_check(nf90_put_var(ncid,maskid,NINT(mask_dst)))
CALL nc_check(nf90_close(ncid))

!=======================================================================

! 5) Clean up:
DEALLOCATE(lon_src,lat_src,mask_src)
DEALLOCATE(lon_dst,lat_dst,mask_dst)

WRITE(*,*)'wrf2umwmgrid: Success.'

!=======================================================================
ENDPROGRAM wrf2umwmgrid
