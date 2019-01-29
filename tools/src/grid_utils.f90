! University of Miami Wave Model (UMWM)
! grid_utils
! Copyright (C) 2013  Milan Curcic
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
MODULE grid_utils
!=======================================================================

IMPLICIT NONE

PRIVATE

PUBLIC :: refineField
PUBLIC :: coarsenField

!=======================================================================
CONTAINS


SUBROUTINE refineField(field_src,field_dst,ratio)
!=======================================================================

! Dummy arguments:
REAL,DIMENSION(:,:),INTENT(IN)  :: field_src
REAL,DIMENSION(:,:),INTENT(OUT) :: field_dst
INTEGER,            INTENT(IN)  :: ratio

INTEGER :: i,j
INTEGER :: i2,j2
INTEGER :: isub,jsub
INTEGER :: im,jm
INTEGER :: im_dst,jm_dst
INTEGER :: im_src,jm_src

REAL :: a,b
REAL :: w1,w2,w3,w4

!=======================================================================

IF(ratio == 1)THEN
  field_dst = field_src
  RETURN
ENDIF

! Get size of source field:
im_src = SIZE(field_src,DIM=1)
jm_src = SIZE(field_src,DIM=2)

! Compute the size of target field:
im = (im_src-1)*ratio+1
jm = (jm_src-1)*ratio+1

! Get size of the target field dummy argument:
im_dst = SIZE(field_dst,DIM=1)
jm_dst = SIZE(field_dst,DIM=2)

IF(im_dst >im .OR. jm_dst > jm)THEN
  WRITE(*,*)'RefineField: ERROR: &
             Source field smaller than required by the target field size:'
  WRITE(*,*)'Source field:                im_src,jm_src = ',im_src,jm_src
  WRITE(*,*)'Interpolates to:                    im,jm  = ',im,jm
  WRITE(*,*)'Target field dummy argument: im_dst,jm_dst = ',im_dst,jm_dst
  RETURN
ENDIF

field_dst = 0

! Outer loop over the coarse grid:
DO j=1,jm_src-1
  DO i=1,im_src-1
    i2 = (i-1)*ratio+1
    j2 = (j-1)*ratio+1
    DO jsub=0,ratio-1
      DO isub=0,ratio-1
        a = 1.-isub/FLOAT(ratio)
        b = 1.-jsub/FLOAT(ratio)
        w1 = a*b
        w2 = b*(1.-a)
        w3 = a*(1.-b)
        w4 = (1.-a)*(1.-b)
        field_dst(i2+isub,j2+jsub) = w1*field_src(i  ,j  )&
                                    +w2*field_src(i+1,j  )&
                                    +w3*field_src(i  ,j+1)&
                                    +w4*field_src(i+1,j+1)
      ENDDO
    ENDDO
  ENDDO
ENDDO

! Close last row:
DO i=1,im_src-1
  i2 = (i-1)*ratio+1
  DO isub=0,ratio-1
    w1 = 1.-isub/FLOAT(ratio)
    w2 = 1.-w1
    field_dst(i2+isub,jm) = w1*field_src(i,jm_src)&
                           +w2*field_src(i+1,jm_src)
  ENDDO
ENDDO

! Close last column:
DO j=1,jm_src-1
  j2 = (j-1)*ratio+1
  DO jsub=0,ratio-1
    w1 = 1.-jsub/FLOAT(ratio)
    w2 = 1.-w1
    field_dst(im,j2+jsub) = w1*field_src(im_src,j)&
                           +w2*field_src(im_src,j+1)
  ENDDO
ENDDO

! Close top-right corner:
field_dst(im,jm) = field_src(im_src,jm_src)

ENDSUBROUTINE refineField
!=======================================================================



SUBROUTINE coarsenField(field_src,field_dst,ratio)
!=======================================================================

! Dummy arguments:
REAL,DIMENSION(:,:),INTENT(IN)  :: field_src
REAL,DIMENSION(:,:),INTENT(OUT) :: field_dst
INTEGER,            INTENT(IN)  :: ratio

INTEGER :: i,j
INTEGER :: isub,jsub
INTEGER :: im,jm

INTEGER :: stride

REAL :: invNumPoints

!=======================================================================

invNumPoints = (1./DBLE(ratio))**2.
stride = ratio-1

im = SIZE(field_dst,DIM=1)
jm = SIZE(field_dst,DIM=2)

! Loop over the coarse grid and average values from the fine grid:
DO j=1,jm
  jsub = (j-1)*ratio+1
  DO i=1,im
    isub = (i-1)*ratio+1
    field_dst(i,j) = SUM(field_src(isub:isub+stride, &
                                   jsub:jsub+stride))&
                    *invNumPoints
  ENDDO
ENDDO

ENDSUBROUTINE coarsenField
!=======================================================================
ENDMODULE grid_utils
