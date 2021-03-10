module umwm_grid

  use umwm_config, only: config_type

  implicit none

  private
  public :: grid_type

  type :: grid_type
    integer :: size_x
    integer :: size_y
    real, allocatable :: lon(:,:)
    real, allocatable :: lat(:,:)
  end type grid_type

  interface grid_type
    module procedure :: grid_type_cons
  end interface grid_type

contains

  pure type(grid_type) function grid_type_cons(config) result(res)
    type(config_type), intent(in) :: config
    integer :: stat

    res % size_x = config % grid_size_x
    res % size_y = config % grid_size_y

    allocate(res % lon(res % size_x, res % size_y), stat=stat)
    res % lon = 0

    allocate(res % lat(res % size_x, res % size_y), stat=stat)
    res % lat = 0

  end function grid_type_cons

end module umwm_grid
