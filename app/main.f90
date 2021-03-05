program umwm
  
  use umwm_domain, only: domain_type
  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type
  
  implicit none

  type(domain_type) :: domain

  domain = domain_type(grid_type(), spectrum_type())

  print *, 'Hello from umwm.'

end program umwm
