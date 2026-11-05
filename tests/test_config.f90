program test_config
  use iso_fortran_env, only: iostat_end
  use tuff, only: test, test_result, nearly_equal
  use umwm_config, only: config_type, sibling_namelist
#ifdef MPI
  use umwm_env, only: env_init, env_stop
#endif

  implicit none

  type(test_result) :: suite

#ifdef MPI
  call env_init()
#endif
  suite = test('test_config', [ &
    test(valid_config_read), &
    test(invalid_grid_size), &
    test(invalid_pm), &
    test(invalid_frequency_ordering), &
    test(invalid_frequency_range), &
    test(reference_time_read), &
    test(invalid_reference_time), &
    test(invalid_output_interval), &
    test(invalid_seaice_thresholds), &
    test(invalid_stokes_depths), &
    test(cutoff_factor), &
    test(sibling_paths), &
    test(lonlat_grid_settings), &
    test(waves_and_source), &
    test(unreadable_namelist) &
  ])

#ifdef MPI
  call env_stop()
#endif
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
      trim(config % starttimestr) == '2026-01-01 00:00:00' .and. &
      trim(config % stoptimestr) == '2026-01-01 06:00:00' .and. &
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


  function invalid_frequency_range() result(res)
    ! The frequency range must be wider than the shift of the downshifting.
    type(test_result) :: res
    type(config_type) :: config

    config = valid_config()
    config % fmin = 0.9 * config % fmax
    res = test('invalid_frequency_range', &
      .not. config % validate(stop_on_error=.false., rank=1))
  end function invalid_frequency_range


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


  function cutoff_factor() result(res)
    ! fcut_fac is read, defaults to 0.53 when not set, and must be positive.
    type(test_result) :: res
    type(config_type) :: config, without_cutoff
    logical :: invalid

    config = valid_config()
    without_cutoff = config_type('namelists/forcing_rankine.nml')
    config % fcut_fac = 0
    invalid = .not. config % validate(stop_on_error=.false., rank=1)
    config = valid_config()

    res = test('cutoff_factor', invalid .and. &
      nearly_equal(config % fcut_fac, 0.53) .and. &
      nearly_equal(without_cutoff % fcut_fac, 0.53) .and. &
      without_cutoff % validate(stop_on_error=.false., rank=1))
  end function cutoff_factor


  function sibling_paths() result(res)
    ! Files beside the main namelist, such as spectrum.nml, with the prefix
    ! of a main namelist named *main.nml, such as the copies of an output
    ! directory; in the input directory, input/ by default; and in the
    ! output directory, the current directory by default.
    type(test_result) :: res
    type(config_type) :: config, here, copy

    config = valid_config()
    here % path = 'main.nml'
    here % input_dir = '../examples/seamount/input'
    here % output_dir = 'runs/seamount'
    copy % path = 'runs/seamount/umwmout_main.nml'
    res = test('sibling_paths', &
      config % sibling_path('spectrum.nml') == '../namelists/spectrum.nml' .and. &
      here % sibling_path('exclude.nml') == 'exclude.nml' .and. &
      copy % sibling_path('spectrum.nml') == 'runs/seamount/umwmout_spectrum.nml' .and. &
      sibling_namelist('cases/lake.nml', 'spectrum.nml') == 'cases/spectrum.nml' .and. &
      sibling_namelist('cases/domain.nml', 'spectrum.nml') == 'cases/spectrum.nml' .and. &
      sibling_namelist('_main.nml', 'exclude.nml') == '_exclude.nml' .and. &
      config % input_path('umwm.grid') == 'input/umwm.grid' .and. &
      here % input_path('umwm.gridtopo') == '../examples/seamount/input/umwm.gridtopo' .and. &
      config % output_path('umwmout.grid') == './umwmout.grid' .and. &
      here % output_path('umwmrst_2026-01-01_00:00:00.nc') == 'runs/seamount/umwmrst_2026-01-01_00:00:00.nc' &
      .and. .not. config % steady)
  end function sibling_paths


  function lonlat_grid_settings() result(res)
    ! A lon-lat grid needs both spacings, latitudes within [-90, 90], and,
    ! if global, longitudes that close the circle; the namelist default is
    ! no lon-lat grid.
    type(test_result) :: res
    type(config_type) :: config
    logical :: ok

    config = valid_config()
    ok = .not. config % lonlat_grid() .and. config % dlon == 0 .and. config % dlat == 0
    config % dlon = 0.25
    ok = ok .and. config % lonlat_grid() .and. .not. config % validate(stop_on_error=.false., rank=-1)
    config % dlat = 0.25
    config % lat0 = 88 ! to 90.5
    ok = ok .and. .not. config % validate(stop_on_error=.false., rank=-1)
    config % lat0 = 10
    ok = ok .and. config % validate(stop_on_error=.false., rank=-1)
    config % gridfromfile = .true.
    ok = ok .and. .not. config % lonlat_grid()

    res = test('lonlat_grid_settings', ok)
  end function lonlat_grid_settings


  function waves_and_source() result(res)
    ! The wave partitions are read from file with waves, and the input files
    ! come from a data source unless it is 'none' or empty, the default.
    type(test_result) :: res
    type(config_type) :: config
    logical :: ok

    config = valid_config()
    ok = .not. config % waves .and. .not. config % has_source() .and. &
         .not. config % forcing_from_file()
    config % waves = .true.
    ok = ok .and. config % forcing_from_file()
    config % source = 'era5'
    ok = ok .and. config % has_source()
    config % source = ' None'
    ok = ok .and. .not. config % has_source()
    config % source = ''
    ok = ok .and. .not. config % has_source()

    res = test('waves_and_source', ok)
  end function waves_and_source


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
