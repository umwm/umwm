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

  domain = domain_type(config % start_time, &
                       config % stop_time, &
                       grid_type(config), &
                       spectrum_type(config))

  print *, domain % start_time % strftime('%Y-%m-%d %H:%M:%S')
  print *, domain % stop_time % strftime('%Y-%m-%d %H:%M:%S')
  print *, domain % spectrum % frequency
  print *, domain % spectrum % direction

end program umwm
