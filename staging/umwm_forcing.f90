module umwm_forcing

  use umwm_grid, only: grid_type

  implicit none

  private
  public :: forcing_type

  type :: forcing_type
    real, allocatable :: u_atmosphere(:,:)
    real, allocatable :: v_atmosphere(:,:)
    real, allocatable :: density_atmosphere(:,:)
    real, allocatable :: u_ocean(:,:)
    real, allocatable :: v_ocean(:,:)
    real, allocatable :: density_ocean(:,:)
  contains
    procedure :: load_from_netcdf
  end type forcing_type

  interface forcing_type
    module procedure :: forcing_type_cons
  end interface forcing_type

contains

  pure type(forcing_type) function forcing_type_cons(grid) result(res)
    type(grid_type), intent(in) :: grid
    integer :: stat

    allocate(res % u_atmosphere(grid % size_x, grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating forcing % u_atmosphere.'
    res % u_atmosphere = 0

    allocate(res % v_atmosphere(grid % size_x, grid % size_y), stat=stat) 
    if (stat /= 0) error stop 'Error allocating forcing % v_atmosphere.'
    res % v_atmosphere = 0

    allocate(res % density_atmosphere(grid % size_x, grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating forcing % density_atmosphere.'
    res % density_atmosphere = 0

    allocate(res % u_ocean(grid % size_x, grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating forcing % u_ocean.'
    res % u_ocean = 0

    allocate(res % v_ocean(grid % size_x, grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating forcing % v_ocean.'
    res % v_ocean = 0

    allocate(res % density_ocean(grid % size_x, grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating forcing % density_ocean.'
    res % density_ocean = 0

  end function forcing_type_cons

  subroutine load_from_netcdf(self, filename)
    class(forcing_type), intent(inout) :: self
    character(len=*), intent(in) :: filename
    ! TODO
  end subroutine load_from_netcdf

end module umwm_forcing
