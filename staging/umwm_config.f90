module umwm_config

  use, intrinsic :: iso_fortran_env, only: stdout => output_unit, stderr => error_unit
  use datetime_module, only: datetime, strptime, timedelta
  use tomlf, only: get_value, toml_parse, toml_table

  implicit none

  private
  public :: config_type

  type :: config_type
    type(datetime) :: start_time
    type(datetime) :: stop_time
    type(timedelta) :: interval
    character(:), allocatable :: name
    integer :: grid_size_x
    integer :: grid_size_y
    real :: dx
    real :: dy
    integer :: num_frequencies
    integer :: num_directions
    real :: frequency_min
    real :: frequency_max
    real :: gravity
    real :: surface_tension
    logical :: forcing_from_file
    real :: wind_speed
    real :: current_speed
    real :: air_density
    real :: water_density
  end type config_type

  interface config_type
    module procedure :: config_type_cons
  end interface config_type

contains

  type(config_type) function config_type_cons(filename) result(res)
    character(*), intent(in), optional :: filename
    character(:), allocatable :: fn
    integer :: unit
    type(toml_table), allocatable :: table
    type(toml_table), pointer :: domain_table, spectrum_table, physics_table, forcing_table
    !type(toml_key), allocatable :: keys(:)
    character(:), allocatable :: start_time_str, stop_time_str
    integer :: output_interval_seconds
    logical :: ok = .true.

    ! If file name is provided we'll use that,
    ! otherwise we default to umwm.toml.
    if (present(filename)) then
      fn = filename
    else
      fn = 'umwm.toml'
    end if

    ! Open a file in read-only mode. It must already exist.
    open(newunit=unit, file=fn, status='old', action='read')
    call toml_parse(table, unit)
    close(unit)

    call get_value(table, 'name', res % name)
    call get_value(table, 'domain', domain_table)
    call get_value(table, 'spectrum', spectrum_table)
    call get_value(table, 'physics', physics_table)
    call get_value(table, 'forcing', forcing_table)

    call get_value(domain_table, 'start_time', start_time_str)
    call get_value(domain_table, 'stop_time', stop_time_str)
    call get_value(domain_table, 'output_interval_seconds', output_interval_seconds)
    call get_value(domain_table, 'grid_size_x', res % grid_size_x)
    call get_value(domain_table, 'grid_size_y', res % grid_size_y)
    call get_value(domain_table, 'dx', res % dx)
    call get_value(domain_table, 'dy', res % dy)

    call get_value(spectrum_table, 'num_frequencies', res % num_frequencies)
    call get_value(spectrum_table, 'num_directions', res % num_directions)
    call get_value(spectrum_table, 'frequency_min', res % frequency_min)
    call get_value(spectrum_table, 'frequency_max', res % frequency_max)

    call get_value(physics_table, 'gravity', res % gravity)
    call get_value(physics_table, 'surface_tension', res % surface_tension)

    call get_value(forcing_table, 'from_file', res % forcing_from_file)
    call get_value(forcing_table, 'wind_speed', res % wind_speed)
    call get_value(forcing_table, 'current_speed', res % current_speed)
    call get_value(forcing_table, 'air_density', res % air_density)
    call get_value(forcing_table, 'water_density', res % water_density)

    res % start_time = strptime(start_time_str, '%Y-%m-%d %H:%M:%S')
    res % stop_time = strptime(stop_time_str, '%Y-%m-%d %H:%M:%S')
    res % interval = timedelta(seconds=output_interval_seconds)

    if (len(res % name) == 0) then
      ok = .false.
      write(stderr, '(a)') 'Error: Lenght of name in ' // fn // ' must be > 0.'
    end if

    if (res % stop_time < res % start_time) then
      ok = .false.
      write(stderr, '(a)') 'Error: stop_time in ' // fn // ' must be >= start_time.'
    end if

    if (output_interval_seconds <= 0) then
      ok = .false.
      write(stderr, '(a)') 'Error: output_interval_seconds in ' // fn // ' must be > 0.'
    end if

    if (res % grid_size_x < 1) then
      ok = .false.
      write(stderr, '(a)') 'Error: grid_size_x in ' // fn // ' must be > 0.'
    end if

    if (res % grid_size_x < 1) then
      ok = .false.
      write(stderr, '(a)') 'Error: grid_size_y in ' // fn // ' must be > 0.'
    end if

    if (res % dx <= 0) then
      ok = .false.
      write(stderr, '(a)') 'Error: dx in ' // fn // ' must be > 0.'
    end if

    if (res % dy <= 0) then
      ok = .false.
      write(stderr, '(a)') 'Error: dy in ' // fn // ' must be > 0.'
    end if

    if (res % num_frequencies < 1) then
      ok = .false.
      write(stderr, '(a)') 'Error: num_frequencies in ' // fn // ' must be > 0.'
    end if

    if (res % num_directions < 1) then
      ok = .false.
      write(stderr, '(a)') 'Error: num_directions in ' // fn // ' must be > 0.'
    end if

    if (res % frequency_min > res % frequency_max) then
      ok = .false.
      write(stderr, '(a)') 'Error: frequency_min in ' // fn // ' must be <= frequency_max.'
    end if

    if (res % gravity <= 0) then
      ok = .false.
      write(stderr, '(a)') 'Error: gravity in ' // fn // ' must be > 0.'
    end if

    if (res % surface_tension <= 0) then
      ok = .false.
      write(stderr, '(a)') 'Error: surface_tension in ' // fn // ' must be > 0.'
    end if

    if (res % air_density <= 0) then
      ok = .false.
      write(stderr, '(a)') 'Error: air_density in ' // fn // ' must be > 0.'
    end if

    if (res % water_density <= 0) then
      ok = .false.
      write(stderr, '(a)') 'Error: water_density in ' // fn // ' must be > 0.'
    end if

    if (.not. ok) error stop 1

  end function config_type_cons

end module umwm_config
