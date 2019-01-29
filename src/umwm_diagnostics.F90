module umwm_diagnostics

  use umwm_constants, only: rk

  implicit none

  private
  public :: significant_wave_height, spectral_moment

contains

  pure function significant_wave_height(Fk, k, dk, th) result(res)
    real(rk), intent(in) :: Fk(:,:,:), k(:,:), dk(:), th(:)
    real(rk), allocatable :: res(:)
  end function significant_wave_height

  pure function spectral_moment(Fk, k, dk, th, order) result(res)
    real(rk), intent(in) :: Fk(:,:,:), k(:,:), dk(:), th(:)
    integer, intent(in) :: order
    real(rk), allocatable :: res(:)
  end function spectral_moment

end module umwm_diagnostics
