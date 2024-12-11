module umwm_state

  use umwm_config, only: config_type
  use umwm_dispersion, only: wavenumber, group_speed
  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type

  implicit none

  private
  public :: state_type

  real, parameter :: twopi = 2 * acos(-1.)

  type :: state_type
    real, allocatable :: depth(:,:) ! Mean water depth; may vary in time.
    real, allocatable :: variance(:,:,:,:)
    real, allocatable :: wavenumber(:,:,:)
    real, allocatable :: phase_speed(:,:,:)
    real, allocatable :: group_speed(:,:,:)
    real, allocatable :: dk(:,:,:) ! Wavenumber spacing
    real, allocatable :: ustar(:,:)
    real, allocatable :: wind_input(:,:,:,:)
    real, allocatable :: dissipation(:,:,:,:)
  end type state_type

  interface state_type
    module procedure :: state_type_cons
  end interface state_type

contains

  elemental type(state_type) function state_type_cons(grid, spectrum, config) result(res)

    type(grid_type), intent(in) :: grid
    type(spectrum_type), intent(in) :: spectrum
    type(config_type), intent(in) :: config
    integer :: i, j
    integer :: stat
   
    allocate(res % depth(grid % size_x, grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % state % action.'
    res % depth = 1000 ! FIXME load from file, CLI, or config file

    allocate(res % variance(spectrum % num_frequencies, &
                            spectrum % num_directions, &
                            grid % size_x, &
                            grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % state % action.'
    res % variance = 0

    allocate(res % wind_input(spectrum % num_frequencies, &
                              spectrum % num_directions, &
                              grid % size_x, &
                              grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % state % wind_input.'
    res % wind_input = 0

    allocate(res % dissipation(spectrum % num_frequencies, &
                               spectrum % num_directions, &
                               grid % size_x, &
                               grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % state % dissipation.'
    res % dissipation = 0

    allocate(res % wavenumber(spectrum % num_frequencies, &
                              grid % size_x, &
                              grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % state % wavenumber.'

    allocate(res % phase_speed(spectrum % num_frequencies, &
                               grid % size_x, &
                               grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % state % phase_speed.'

    allocate(res % group_speed(spectrum % num_frequencies, &
                                grid % size_x, &
                                grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % state % group_speed.'

    allocate(res % dk(spectrum % num_frequencies, &
                       grid % size_x, &
                       grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % state % dk.'

    allocate(res % ustar(grid % size_x, grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % state % ustar.'
    res % ustar = 0

    ! Initialize wave parameters by solving the dispersion relationship.
    do j = 1, grid % size_y
      do i = 1, grid % size_x

        res % wavenumber(:,i,j) = wavenumber( &
          spectrum % frequency, &
          res % depth(i,j), &
          1e3, &
          config % gravity, &
          config % surface_tension &
        )

        res % phase_speed(:,i,j) = twopi * spectrum % frequency / res % wavenumber(:,i,j)

        res % group_speed(:,i,j) = group_speed( &
          res % wavenumber(:,i,j), &
          res % depth(i,j), &
          1e3, &
          config % gravity, &
          config % surface_tension &
        )

        res % dk(:,i,j) = twopi * spectrum % frequency * spectrum % dlnf / res % group_speed(:,i,j)

      end do
    end do

  end function state_type_cons

end module umwm_state
