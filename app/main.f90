program umwm
  
  use tomlf, only: toml_parse, toml_table, toml_array, get_value, toml_key
  use umwm_domain, only: domain_type
  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type
  
  implicit none

  integer :: n, unit
  integer :: grid_size_x, grid_size_y
  type(toml_table), allocatable :: table
  type(toml_table), pointer :: domain_table
  type(toml_key), allocatable :: keys(:)
  type(domain_type) :: domain
  character(:), allocatable :: name

  open(newunit=unit, file='umwm.toml', status='old')
  call toml_parse(table, unit)
  close(unit)

  print *, 'Keys in umwm.toml:'
  call table % get_keys(keys)
  do n = 1, size(keys)
    print *, n, keys(n) % key
  end do

  call get_value(table, 'name', name)
  call get_value(table, 'domain', domain_table)
  call get_value(domain_table, 'grid_size_x', grid_size_x)
  call get_value(domain_table, 'grid_size_y', grid_size_y)

  print *, name
  print *, grid_size_x
  print *, grid_size_y
  
  domain = domain_type(grid_type(), spectrum_type())

end program umwm
