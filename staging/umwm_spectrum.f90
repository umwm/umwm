module umwm_spectrum

  implicit none

  private
  public :: spectrum_type

  type :: spectrum_type
    real, allocatable :: f(:), th(:)
    real, allocatable :: df(:), dth(:)
  end type spectrum_type

contains

end module umwm_spectrum
