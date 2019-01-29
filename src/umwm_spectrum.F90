module umwm_spectrum

  use umwm_constants, only: rk

  implicit none

  private
  public :: Spectrum_type

  type :: Spectrum_type
    real(rk), allocatable :: f(:), th(:)
    real(rk), allocatable :: df(:), dth(:)
  end type Spectrum_type

contains

end module umwm_spectrum
