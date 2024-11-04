program test_config
  use, intrinsic :: iso_fortran_env, only: stdout => output_unit, stderr => error_unit
  use datetime_module, only: datetime
  use umwm_config, only: config_type

  implicit none

  type(config_type) :: config
  logical :: ok = .true.

  config = config_type('test/umwm-test.toml')

  if (config % name /= 'Blackwater Bay') then
    write(stderr, '(a)') 'config % name is expected.. failed'
    ok = .false.
  end if

  if (config % start_time /= datetime(2024, 1, 1)) then
    write(stderr, '(a)') 'config % start_time is expected.. failed'
    ok = .false.
  end if

  if (config % stop_time /= datetime(2024, 1, 2)) then
    write(stderr, '(a)') 'config % stop_time is expected.. failed'
    ok = .false.
  end if

  if (config % grid_size_x /= 80) then
    write(stderr, '(a)') 'config % grid_size_x is expected.. failed'
    ok = .false.
  end if

  if (config % grid_size_y /= 60) then
    write(stderr, '(a)') 'config % grid_size_y is expected.. failed'
    ok = .false.
  end if

  if (config % num_frequencies /= 50) then
    write(stderr, '(a)') 'config % num_frequencies is expected.. failed'
    ok = .false.
  end if

  if (config % num_directions /= 36) then
    write(stderr, '(a)') 'config % num_directions is expected.. failed'
    ok = .false.
  end if

  if (config % frequency_min /= 0.04) then
    write(stderr, '(a)') 'config % frequency_min is expected.. failed'
    ok = .false.
  end if

  if (config % frequency_max /= 2) then
    write(stderr, '(a)') 'config % frequency_max is expected.. failed'
    ok = .false.
  end if

  if (config % gravity /= 9.8) then
    write(stderr, '(a)') 'config % gravity is expected.. failed'
    ok = .false.
  end if

  if (config % surface_tension /= 0.074) then
    write(stderr, '(a)') 'config % surface_tension is expected.. failed'
    ok = .false.
  end if

  if (ok) then
    write(stdout, '(a)') 'test_config: All tests passed.'
  else
    write(stderr, '(a)') 'test_config: One or more tests failed.'
    error stop 1
  end if

end program test_config