program umwm
  
  use umwm_config, only: config_type
  use umwm_domain, only: domain_type
  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type

  implicit none

  type(config_type) :: config
  type(domain_type) :: domain

  config = config_type()
  print *, config % name
  print *, config % grid_size_x
  print *, config % grid_size_y

  domain = domain_type(grid_type(config), spectrum_type(config))

  print *, domain % spectrum % frequency

end program umwm
