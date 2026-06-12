program test_config
  use tuff, only: test, test_result, nearly_equal
  use umwm_config, only: config_type

  implicit none

  type(test_result) :: suite

  suite = test('test_config', [ &
    test(valid_config_read), &
    test(invalid_grid_size), &
    test(invalid_pm), &
    test(invalid_frequency_ordering), &
    test(invalid_output_interval), &
    test(invalid_seaice_thresholds), &
    test(invalid_stokes_depths), &
    test(unreadable_namelist) &
  ])

  if (.not. suite % ok) error stop 1

contains

  function valid_config_read() result(res)
    type(test_result) :: res
    type(config_type) :: config
    logical :: ok

    config = valid_config()
    ok = config % validate(stop_on_error=.false., rank=1)

    res = test('valid_config_read', ok .and. &
      (.not. config % isglobal) .and. &
      config % mm == 21 .and. config % nm == 11 .and. &
      config % om == 37 .and. config % pm == 36 .and. &
      nearly_equal(config % fmin, 0.0313) .and. &
      nearly_equal(config % fmax, 2.0) .and. &
      trim(config % starttimestr) == '2012-01-01 00:00:00' .and. &
      trim(config % stoptimestr) == '2012-01-01 06:00:00' .and. &
      nearly_equal(config % g, 9.80665) .and. &
      nearly_equal(config % delx, 10000.0) .and. &
      (.not. config % winds) .and. &
      nearly_equal(config % wspd0, 10.0) .and. &
      config % outgrid == 1 .and. config % outspec == 0 .and. &
      config % outrst == 6 .and. config % stokes .and. &
      allocated(config % stokes_depths) .and. &
      size(config % stokes_depths) == 42 .and. &
      nearly_equal(config % stokes_depths(1), 0.1) .and. &
      nearly_equal(config % stokes_depths(42), 100.0))
  end function valid_config_read


  function invalid_grid_size() result(res)
    type(test_result) :: res
    type(config_type) :: config

    config = valid_config()
    config % mm = 2
    res = test('invalid_grid_size', .not. config % validate(stop_on_error=.false., rank=1))
  end function invalid_grid_size


  function invalid_pm() result(res)
    type(test_result) :: res
    type(config_type) :: config

    config = valid_config()
    config % pm = 10
    res = test('invalid_pm', .not. config % validate(stop_on_error=.false., rank=1))
  end function invalid_pm


  function invalid_frequency_ordering() result(res)
    type(test_result) :: res
    type(config_type) :: config

    config = valid_config()
    config % fmin = config % fmax
    res = test('invalid_frequency_ordering', &
      .not. config % validate(stop_on_error=.false., rank=1))
  end function invalid_frequency_ordering


  function invalid_output_interval() result(res)
    type(test_result) :: res
    type(config_type) :: config

    config = valid_config()
    config % outgrid = 5
    res = test('invalid_output_interval', &
      .not. config % validate(stop_on_error=.false., rank=1))
  end function invalid_output_interval


  function invalid_seaice_thresholds() result(res)
    type(test_result) :: res
    type(config_type) :: config

    config = valid_config()
    config % fice_lth = 0.8
    config % fice_uth = 0.3
    res = test('invalid_seaice_thresholds', &
      .not. config % validate(stop_on_error=.false., rank=1))
  end function invalid_seaice_thresholds


  function invalid_stokes_depths() result(res)
    type(test_result) :: res
    type(config_type) :: config

    config = valid_config()
    deallocate(config % stokes_depths)
    allocate(config % stokes_depths(0))
    res = test('invalid_stokes_depths', &
      .not. config % validate(stop_on_error=.false., rank=1))
  end function invalid_stokes_depths


  function unreadable_namelist() result(res)
    type(test_result) :: res
    type(config_type) :: config

    config = config_type('/tmp/umwm-test-config-does-not-exist.nml')
    res = test('unreadable_namelist', &
      .not. config % validate(stop_on_error=.false., rank=1))
  end function unreadable_namelist


  function valid_config() result(config)
    type(config_type) :: config

    config = config_type('../namelists/main.nml')
  end function valid_config

end program test_config
