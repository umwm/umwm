module umwm_forcing

  use datetime_module, only: datetime
  use umwm_config, only: config_type
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
    real, allocatable :: seaice_fraction(:,:)
  contains
    procedure :: load
    procedure :: set
    procedure :: update
  end type forcing_type

  interface forcing_type
    module procedure :: forcing_type_cons
  end interface forcing_type

contains

  pure type(forcing_type) function forcing_type_cons(grid, config) result(res)
    type(grid_type), intent(in) :: grid
    type(config_type), intent(in) :: config
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

    allocate(res % seaice_fraction(grid % size_x, grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating forcing % seaice_fraction.'
    res % seaice_fraction = 0

    ! Set initial forcing from config
    if (.not. config % forcing_from_file) then
      call res % set( &
        wind_speed = config % wind_speed, &
        current_speed = config % current_speed, &
        air_density = config % air_density, &
        water_density = config % water_density &
      )
    end if

  end function forcing_type_cons


  subroutine load(self, time)
    ! Load from NetCDF file.
    class(forcing_type), intent(inout) :: self
    type(datetime), intent(in) :: time
    ! TODO
  end subroutine load


  pure subroutine set(self, wind_speed, current_speed, air_density, water_density)
    class(forcing_type), intent(inout) :: self
    real, intent(in), optional :: wind_speed, current_speed, air_density, water_density
    if (present(wind_speed)) self % u_atmosphere = wind_speed
    if (present(current_speed)) self % u_ocean = current_speed
    if (present(air_density)) self % density_atmosphere = air_density
    if (present(water_density)) self % density_ocean = water_density
  end subroutine set


  pure subroutine update(self, time)
    class(forcing_type), intent(inout) :: self
    type(datetime), intent(in) :: time
    ! TODO
  end subroutine update

end module umwm_forcing
