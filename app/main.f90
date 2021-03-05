program umwm
  
  use tomlf, only: toml_error, toml_parse, toml_table
  use umwm_domain, only: domain_type
  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type
  
  implicit none

  integer :: unit
  integer :: gridsize_x, gridsize_y
  type(toml_error), allocatable :: error
  type(toml_table), allocatable :: table
  type(domain_type) :: domain

  open(newunit=unit, file='umwm.toml', status='old')
  call toml_parse(table, unit, error)

  if (allocated(error)) then
    print *, error % stat, error % message
    stop
  end if

  if (allocated(table)) then
  end if

  close(unit)

  domain = domain_type(grid_type(), spectrum_type())

  print *, 'Hello from umwm.'

end program umwm
