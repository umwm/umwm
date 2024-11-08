program test_domain
  use, intrinsic :: iso_fortran_env, only: stdout => output_unit, stderr => error_unit
  use datetime_module, only: datetime
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
  logical :: ok = .true.

  ! Set up test configuration
  config = config_type('test/umwm-test.toml')
  clock = clock_type(config % start_time, config % stop_time, config % interval)
  grid = grid_type(config % grid_size_x, config % grid_size_y, config % dx, config % dy)
  spectrum = spectrum_type( &
    config % num_frequencies, &
    config % num_directions, &
    config % frequency_min, &
    config % frequency_max &
  )

  ! Create domain instance
  domain = domain_type(clock, grid, spectrum, config)

  ! Test domain components are properly initialized
  if (.not. allocated(domain % forcing % u_atmosphere)) then
    write(stderr, '(a)') 'domain % forcing % u_atmosphere is allocated.. failed'
    ok = .false.
  end if

  if (.not. allocated(domain % state % variance)) then
    write(stderr, '(a)') 'domain % state % variance is allocated.. failed'
    ok = .false.
  end if

  if (domain % clock % current /= config % start_time) then
    write(stderr, '(a)') 'domain % clock % current matches start_time.. failed'
    ok = .false.
  end if

  if (domain % grid % size_x /= config % grid_size_x) then
    write(stderr, '(a)') 'domain % grid % size_x matches config.. failed'
    ok = .false.
  end if

  if (domain % spectrum % num_frequencies /= config % num_frequencies) then
    write(stderr, '(a)') 'domain % spectrum % num_frequencies matches config.. failed'
    ok = .false.
  end if

  ! Test domain step
  call domain % step()
  if (domain % clock % current /= config % start_time + config % interval) then
    write(stderr, '(a)') 'domain step advances clock.. failed'
    ok = .false.
  end if

  if (ok) then
    write(stdout, '(a)') 'test_domain: All tests passed.'
  else
    write(stderr, '(a)') 'test_domain: One or more tests failed.'
    error stop 1
  end if

end program test_domain
