module umwm_grid

  use umwm_parallel, only: tile_indices

  implicit none

  private
  public :: grid_type

  type :: grid_type
    integer :: size_x, size_y
    integer :: tile_size_x, tile_size_y
    integer :: is, ie, js, je
    real :: dx, dy
    real, allocatable :: x(:,:)
    real, allocatable :: y(:,:)
    real, allocatable :: lon(:,:)
    real, allocatable :: lat(:,:)
  contains
    procedure :: load_latlon_from_netcdf
  end type grid_type

  interface grid_type
    module procedure :: grid_type_cons
  end interface grid_type

contains

  pure type(grid_type) function grid_type_cons(grid_size_x, grid_size_y, dx, dy) result(res)
    integer, intent(in) :: grid_size_x
    integer, intent(in) :: grid_size_y
    real, intent(in) :: dx
    real, intent(in) :: dy

    integer :: i, j
    integer :: stat
    integer :: tile_start_end(4) ! [is, ie, js, je]

    res % size_x = grid_size_x
    res % size_y = grid_size_y
    res % dx = dx
    res % dy = dy

    allocate(res % x(res % size_x, res % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating grid % x.'
    do i = 1, res % size_x
      res % x(i,:) = (i - 1) * dx
    end do

    allocate(res % y(res % size_x, res % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating grid % y.'
    do j = 1, res % size_y
      res % y(:,j) = (j - 1) * dy
    end do

    tile_start_end = tile_indices([res % size_x, res % size_y])
    res % is = tile_start_end(1)
    res % ie = tile_start_end(2)
    res % js = tile_start_end(3)
    res % je = tile_start_end(4)

  end function grid_type_cons


  subroutine load_latlon_from_netcdf(self, filename)
    class(grid_type), intent(inout) :: self
    character(len=*), intent(in) :: filename

    integer :: stat

    allocate(self % lon(self % size_x, self % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating grid % lon.'
    self % lon = 0

    allocate(self % lat(self % size_x, self % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating grid % lat.'
    self % lat = 0
    ! TODO
  end subroutine load_latlon_from_netcdf

end module umwm_grid
