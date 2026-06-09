module umwm_diagnostics

  use umwm_constants, only: rk, twopi

  implicit none

  private
  public :: significant_wave_height, spectral_moment

contains

  pure function significant_wave_height(Fk, k, dk, th) result(res)
    real(rk), intent(in) :: Fk(:,:,:), k(:,:), dk(:), th(:)
    real(rk), allocatable :: res(:)
    real(rk), allocatable :: m0(:)

    m0 = spectral_moment(Fk, k, dk, th, 0)
    res = 4._rk * sqrt(m0)
  end function significant_wave_height

  pure function spectral_moment(Fk, k, dk, th, order) result(res)
    real(rk), intent(in) :: Fk(:,:,:), k(:,:), dk(:), th(:)
    integer, intent(in) :: order
    real(rk), allocatable :: res(:)
    integer :: i, o, p
    real(rk) :: dth

    allocate(res(size(Fk, 3)))
    res = 0._rk

    if (size(th) == 0) return

    if (size(th) == 1) then
      dth = twopi
    else
      dth = abs(th(2) - th(1))
    end if

    do i = 1, size(Fk, 3)
      do p = 1, size(Fk, 2)
        do o = 1, size(Fk, 1)
          res(i) = res(i) + Fk(o,p,i) * k(o,i)**(order + 1) * dk(o) * dth
        end do
      end do
    end do
  end function spectral_moment

end module umwm_diagnostics
