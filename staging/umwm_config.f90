module umwm_config

  use tomlf, only: get_value, toml_parse, toml_table

  implicit none

  private
  public :: config_type

  type :: config_type
    character(:), allocatable :: name
    integer :: grid_size_x
    integer :: grid_size_y
    integer :: num_frequencies
    integer :: num_directions
    real :: frequency_min
    real :: frequency_max
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
    type(toml_table), pointer :: domain_table, spectrum_table
    !type(toml_key), allocatable :: keys(:)

    ! If file name is provided we'll use that,
    ! otherwise we default to umwm.toml.
    if (present(filename)) then
      fn = filename
    else
      fn = 'umwm.toml'
    end if

    ! Open a file in read-only mode. It must already exist.
    open(newunit=unit, file='umwm.toml', status='old', action='read')
    call toml_parse(table, unit)
    close(unit)

    !print *, 'Keys in umwm.toml:'
    !call table % get_keys(keys)
    !do n = 1, size(keys)
    !  print *, n, keys(n) % key
    !end do

    call get_value(table, 'name', res % name)
    ! TODO CHECK len(res % name) > 0
    call get_value(table, 'domain', domain_table)
    call get_value(table, 'spectrum', spectrum_table)
    call get_value(domain_table, 'grid_size_x', res % grid_size_x)
    ! TODO CHECK 1 < res % grid_size_x
    call get_value(domain_table, 'grid_size_y', res % grid_size_y)
    ! TODO CHECK 1 < res % grid_size_y
    call get_value(spectrum_table, 'num_frequencies', res % num_frequencies)
    call get_value(spectrum_table, 'num_directions', res % num_directions)
    call get_value(spectrum_table, 'frequency_min', res % frequency_min)
    call get_value(spectrum_table, 'frequency_max', res % frequency_max)

  end function config_type_cons

end module umwm_config
