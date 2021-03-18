module umwm_grid

  use umwm_config, only: config_type
  use umwm_parallel, only: tile_indices

  implicit none

  private
  public :: grid_type

  type :: grid_type
    integer :: size_x, size_y
    integer :: tile_size_x, tile_size_y
    integer :: is, ie, js, je
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
    integer :: tile_start_end(4) ! [is, ie, js, je]

    res % size_x = config % grid_size_x
    res % size_y = config % grid_size_y

    allocate(res % lon(res % size_x, res % size_y), stat=stat)
    res % lon = 0

    allocate(res % lat(res % size_x, res % size_y), stat=stat)
    res % lat = 0

    tile_start_end = tile_indices([res % size_x, res % size_y])
    res % is = tile_start_end(1)
    res % ie = tile_start_end(2)
    res % js = tile_start_end(3)
    res % je = tile_start_end(4)

  end function grid_type_cons

end module umwm_grid
