program test_config
  use iso_fortran_env, only: iostat_end
  use tuff, only: test, test_result, nearly_equal
  use umwm_config, only: config_type

  implicit none

  type(test_result) :: suite

  suite = test('test_config', [ &
    test(valid_config_read), &
    test(invalid_grid_size), &
    test(invalid_pm), &
    test(invalid_frequency_ordering), &
    test(reference_time_read), &
    test(invalid_reference_time), &
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
      trim(config % reftimestr) == '1970-01-01 00:00:00' .and. &
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


  function reference_time_read() result(res)
    type(test_result) :: res
    type(config_type) :: config
    character(len=*), parameter :: test_path = '/tmp/umwm-test-config-reftime.nml'
    character(len=*), parameter :: ref_time = '2010-05-06 07:08:09'
    logical :: ok

    ok = write_config_with_ref_time(test_path, ref_time)
    config = config_type(test_path)
    ok = ok .and. config % validate(stop_on_error=.false., rank=1)

    res = test('reference_time_read', ok .and. &
      trim(config % reftimestr) == ref_time)
  end function reference_time_read


  function invalid_reference_time() result(res)
    type(test_result) :: res
    type(config_type) :: config

    config = valid_config()
    config % reftimestr = '2012-99-01 00:00:00'
    res = test('invalid_reference_time', &
      .not. config % validate(stop_on_error=.false., rank=1))
  end function invalid_reference_time


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


  logical function write_config_with_ref_time(path, ref_time) result(ok)
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: ref_time

    character(len=512) :: line
    integer :: input_unit
    integer :: output_unit
    integer :: stat
    logical :: in_domain

    ok = .false.
    open(newunit=input_unit, file='../namelists/main.nml', status='old', &
         form='formatted', access='sequential', action='read', iostat=stat)
    if (stat /= 0) return

    open(newunit=output_unit, file=path, status='replace', form='formatted', &
         access='sequential', action='write', iostat=stat)
    if (stat /= 0) then
      close(input_unit)
      return
    end if

    in_domain = .false.
    do
      read(input_unit, '(a)', iostat=stat) line
      if (stat /= 0) exit

      if (trim(adjustl(line)) == '&DOMAIN') in_domain = .true.
      if (in_domain .and. trim(adjustl(line)) == '/') then
        write(output_unit, '(a)') "  refTimeStr  = '" // ref_time // "'"
        in_domain = .false.
      end if
      write(output_unit, '(a)') trim(line)
    end do

    close(input_unit)
    close(output_unit)
    ok = stat == iostat_end
  end function write_config_with_ref_time

end program test_config
