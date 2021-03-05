module umwm_grid

  use umwm_constants, only: rk

  implicit none

  private
  public :: grid_type

  type :: grid_type
    real(rk), allocatable :: lon(:,:), lat(:,:)
  end type grid_type

contains

end module umwm_grid
