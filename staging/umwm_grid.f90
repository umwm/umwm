module umwm_grid

  implicit none

  private
  public :: grid_type

  type :: grid_type
    real, allocatable :: lon(:,:), lat(:,:)
  end type grid_type

contains

end module umwm_grid
