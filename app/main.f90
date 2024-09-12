program umwm
  
  use umwm_clock, only: clock_type
  use umwm_config, only: config_type
  use umwm_domain, only: domain_type
  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type

  implicit none

  type(config_type) :: config
  type(clock_type) :: clock
  type(grid_type) :: grid
  type(spectrum_type) :: spectrum
  type(domain_type) :: domain

  config = config_type()

  clock = clock_type(config % start_time, config % stop_time)

  grid = grid_type(config % grid_size_x, config % grid_size_y)

  spectrum = spectrum_type( &
    config % num_frequencies, &
    config % num_directions, &
    config % frequency_min, &
    config % frequency_max &
  )

  domain = domain_type(clock, grid, spectrum)

  call domain % run()

end program umwm
