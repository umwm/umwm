!
! University of Miami Wave Model (UMWM)
! umwm_topogen
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
PROGRAM umwm_topogen
!=======================================================================
!
! DESCRIPTION: Topography generation utility for UMWM. Reads UMWM input
! grid file and ETOPO topography field. Outputs grid and topography
! in NetCDF format (umwm.gridtopo) that can be used by the wave model. 
!     
! REVISION HISTORY:
!         
! Jun 2012, 1.0.1, Added error checking for NetCDF function calls
!                  Added Seamask input capability to match in order
!                  to generate mask-matching bathymetry
! Apr 2012, 1.0.0, First public release
! 
!=======================================================================
USE netcdf
USE umwm_utils,ONLY:nc_check
!=======================================================================
IMPLICIT NONE

INTEGER :: i,j
INTEGER :: i0,j0
INTEGER :: ic,jc
INTEGER :: is,ie,js,je

INTEGER,DIMENSION(2) :: coords

LOGICAL :: useSeamask

INTEGER :: STATUS
INTEGER :: ncid,xdimid,ydimid,lonid,latid,zid,maskid

CHARACTER(LEN=NF90_MAX_NAME) :: dimname

CHARACTER(LEN=9999) :: topoInputFile,umwmInputFile

TYPE grid
  INTEGER                         :: idm,jdm
  REAL,DIMENSION(:),  ALLOCATABLE :: lon1d,lat1d
  REAL,DIMENSION(:,:),ALLOCATABLE :: lon,lat,z
ENDTYPE grid 

TYPE(grid) :: topo,umwm

INTEGER,DIMENSION(:,:),ALLOCATABLE :: mask

NAMELIST /topogen/ umwmInputFile,topoInputFile,useSeamask
!=======================================================================

! 1) Read namelist:
OPEN(UNIT=10,FILE='namelists/topogen.nml',STATUS='OLD',&
     FORM='FORMATTED',ACCESS='SEQUENTIAL')
READ(UNIT=10,NML=topogen)
CLOSE(UNIT=10)

!=======================================================================

! 2) Read input grid:
WRITE(*,*)'umwm_topogen: Reading input (UMWM) grid.'

CALL nc_check(nf90_open(TRIM(umwmInputFile),nf90_nowrite,ncid))
CALL nc_check(nf90_inq_dimid(ncid,'x',xdimid))
CALL nc_check(nf90_inq_dimid(ncid,'y',ydimid))
CALL nc_check(nf90_inquire_dimension(ncid,xdimid,dimname,umwm%idm))
CALL nc_check(nf90_inquire_dimension(ncid,ydimid,dimname,umwm%jdm))

ALLOCATE(umwm%lon(umwm%idm,umwm%jdm))
ALLOCATE(umwm%lat(umwm%idm,umwm%jdm))
ALLOCATE(umwm%z(umwm%idm,umwm%jdm))

IF(useSeamask)ALLOCATE(mask(umwm%idm,umwm%jdm))

CALL nc_check(nf90_inq_varid(ncid,'lon',lonid))
CALL nc_check(nf90_inq_varid(ncid,'lat',latid))
CALL nc_check(nf90_get_var(ncid,lonid,umwm%lon))
CALL nc_check(nf90_get_var(ncid,latid,umwm%lat))

IF(useSeamask)THEN
  CALL nc_check(nf90_inq_varid(ncid,'seamask',maskid))
  CALL nc_check(nf90_get_var(ncid,maskid,mask))
ENDIF

CALL nc_check(nf90_close(ncid))

!=======================================================================

! 3) Read ETOPO01 grid and topography:
WRITE(*,*)'umwm_topogen: Reading input topography data.'

CALL nc_check(nf90_open(TRIM(topoInputFile),nf90_nowrite,ncid))
CALL nc_check(nf90_inq_dimid(ncid,'x',xdimid))
CALL nc_check(nf90_inq_dimid(ncid,'y',ydimid))
CALL nc_check(nf90_inquire_dimension(ncid,xdimid,dimname,topo%idm))
CALL nc_check(nf90_inquire_dimension(ncid,ydimid,dimname,topo%jdm))

ALLOCATE(topo%lon1d(topo%idm))
ALLOCATE(topo%lat1d(topo%jdm))
ALLOCATE(topo%z(topo%idm,topo%jdm))

CALL nc_check(nf90_inq_varid(ncid,'x',lonid))
CALL nc_check(nf90_inq_varid(ncid,'y',latid))
CALL nc_check(nf90_inq_varid(ncid,'z',zid))

CALL nc_check(nf90_get_var(ncid,lonid,topo%lon1d))
CALL nc_check(nf90_get_var(ncid,latid,topo%lat1d))
CALL nc_check(nf90_get_var(ncid,zid,topo%z))

CALL nc_check(nf90_close(ncid))

!=======================================================================

! 4) Loop over target (UMWM) grid and sample ETOPO01 values:
WRITE(*,*)'umwm_topogen: Sampling topography on target grid.'

DO j=1,umwm%jdm
  WRITE(*,*)'umwm_topogen: Sampling row',j
  DO i=1,umwm%idm

    i0 = MINLOC(ABS(umwm%lon(i,j)-topo%lon1d),DIM=1)
    j0 = MINLOC(ABS(umwm%lat(i,j)-topo%lat1d),DIM=1)

    umwm%z(i,j) = topo%z(i0,j0)

  ENDDO
ENDDO

! Adjust bathymetry to match the input seamask:
IF(useSeamask)THEN
  WHERE(mask==0.AND.umwm%z <0)umwm%z = 0
  WHERE(mask==1.AND.umwm%z>=0)umwm%z = -1
ENDIF

!=======================================================================

! 5) Generate new umwm.gridtopo file:
WRITE(*,*)'umwm_topogen: Writing umwm.gridtopo file.'

CALL nc_check(nf90_create('umwm.gridtopo',nf90_clobber,ncid))
CALL nc_check(nf90_def_dim(ncid,'x',umwm%idm,xdimid))
CALL nc_check(nf90_def_dim(ncid,'y',umwm%jdm,ydimid))
CALL nc_check(nf90_def_var(ncid,'lon',nf90_float,[xdimid,ydimid],lonid))
CALL nc_check(nf90_def_var(ncid,'lat',nf90_float,[xdimid,ydimid],latid))
IF(useSeamask)CALL nc_check(nf90_def_var(ncid,'seamask',nf90_int,[xdimid,ydimid],maskid))
CALL nc_check(nf90_def_var(ncid,'z',nf90_float,[xdimid,ydimid],zid))
CALL nc_check(nf90_put_att(ncid,NF90_GLOBAL,name='history',values='Generated by umwm_topogen'))
CALL nc_check(nf90_enddef(ncid))
CALL nc_check(nf90_put_var(ncid,lonid,umwm%lon))
CALL nc_check(nf90_put_var(ncid,latid,umwm%lat))
IF(useSeamask)CALL nc_check(nf90_put_var(ncid,maskid,mask))
CALL nc_check(nf90_put_var(ncid,zid,umwm%z))
CALL nc_check(nf90_close(ncid))

WRITE(*,*)'umwm_topogen: Success.'

!=======================================================================
ENDPROGRAM umwm_topogen
