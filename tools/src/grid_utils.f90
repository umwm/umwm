module grid_utils

implicit none

private
public :: refineField, coarsenField

contains


subroutine refineField(field_src,field_dst,ratio)
  real, intent(in) :: field_src(:,:)
  real, intent(out) :: field_dst(:,:)
  integer, intent(in) :: ratio

  integer :: i,j
  integer :: i2,j2
  integer :: isub,jsub
  integer :: im,jm
  integer :: im_dst,jm_dst
  integer :: im_src,jm_src

  real :: a,b
  real :: w1,w2,w3,w4

  if (ratio == 1) then
    field_dst = field_src
    return
  end if

  ! Get size of source field:
  im_src = size(field_src, dim=1)
  jm_src = size(field_src, dim=2)

  ! Compute the size of target field:
  im = (im_src-1) *ratio + 1
  jm = (jm_src-1) *ratio + 1

  ! Get size of the target field dummy argument:
  im_dst = size(field_dst, dim=1)
  jm_dst = size(field_dst, dim=2)

  if (im_dst > im .or. jm_dst > jm) then
    write(*,*) 'RefineField: ERROR: &
                Source field smaller than required by the target field size:'
    write(*,*) 'Source field:                im_src,jm_src = ',im_src,jm_src
    write(*,*) 'Interpolates to:                    im,jm  = ',im,jm
    write(*,*) 'Target field dummy argument: im_dst,jm_dst = ',im_dst,jm_dst
    return
  end if

  field_dst = 0

  ! Outer loop over the coarse grid:
  do j=1,jm_src-1
    do i=1,im_src-1
      i2 = (i-1)*ratio+1
      j2 = (j-1)*ratio+1
      do jsub=0,ratio-1
        do isub=0,ratio-1
          a = 1.-isub/real(ratio)
          b = 1.-jsub/real(ratio)
          w1 = a*b
          w2 = b*(1.-a)
          w3 = a*(1.-b)
          w4 = (1.-a)*(1.-b)
          field_dst(i2+isub,j2+jsub) = w1*field_src(i  ,j  )&
                                      +w2*field_src(i+1,j  )&
                                      +w3*field_src(i  ,j+1)&
                                      +w4*field_src(i+1,j+1)
        end do
      end do
    end do
  end do

  ! Close last row:
  do i=1,im_src-1
    i2 = (i-1)*ratio+1
    do isub=0,ratio-1
      w1 = 1.-isub/real(ratio)
      w2 = 1.-w1
      field_dst(i2+isub,jm) = w1*field_src(i,jm_src)&
                             +w2*field_src(i+1,jm_src)
    end do
  end do

  ! Close last column:
  do j=1,jm_src-1
    j2 = (j-1)*ratio+1
    do jsub=0,ratio-1
      w1 = 1.-jsub/real(ratio)
      w2 = 1.-w1
      field_dst(im,j2+jsub) = w1*field_src(im_src,j)&
                             +w2*field_src(im_src,j+1)
    end do
  end do

  ! Close top-right corner:
  field_dst(im,jm) = field_src(im_src, jm_src)

end subroutine refineField


subroutine coarsenField(field_src, field_dst, ratio)
  real, intent(in) :: field_src(:,:)
  real, intent(out) :: field_dst(:,:)
  integer, intent(in) :: ratio

  integer :: i,j
  integer :: isub,jsub
  integer :: im,jm

  integer :: stride

  real :: invNumPoints

  invNumPoints = (1 / dble(ratio))**2
  stride = ratio - 1

  im = size(field_dst, dim=1)
  jm = size(field_dst, dim=2)

  ! Loop over the coarse grid and average values from the fine grid:
  do j = 1, jm
    jsub = (j - 1) * ratio + 1
    do i = 1, im
      isub = (i - 1) * ratio + 1
      field_dst(i,j) = sum(field_src(isub:isub+stride,jsub:jsub+stride)) * invNumPoints
    end do
  end do

end subroutine coarsenField

end module grid_utils
