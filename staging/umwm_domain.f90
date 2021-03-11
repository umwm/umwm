module umwm_domain

  !! domain_type is the main UMWM derived type.
  !! Its components are an instance of grid_type, which defines its geographical
  !! grid, and an instance of spectrum_type, which defines its spectral space.

  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type

  implicit none

  private
  public :: domain_type

  type :: domain_type
    type(grid_type) :: grid
    type(spectrum_type) :: spectrum
    real, allocatable :: variance(:,:)
  end type domain_type

contains

end module umwm_domain
