module umwm_domain

  ! Under construction. The Domain_type will supersede 
  ! the global namespace in umwm_module.

  use umwm_constants, only: rk
  !use umwm_grid, only: Grid_type
  use umwm_spectrum, only: Spectrum_type

  implicit none

  private
  public :: Domain_type

  type :: Domain_type
    !type(Grid_type) :: grid
    type(Spectrum_type) :: spectrum
  end type Domain_type

contains

end module umwm_domain
