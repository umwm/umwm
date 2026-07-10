program test_forcing
  use tuff, only: test, test_result, nearly_equal, all_nearly_equal
  use umwm_config, only: config_type
  use umwm_forcing, only: forcing_type, min_wind_speed
  use umwm_grid, only: grid_type
  use umwm_module, only: sumt

  implicit none

  type(test_result) :: suite

  suite = test('test_forcing', [ &
    test(interpolation_snapshot), &
    test(wind_speed_floor), &
    test(stored_wind_speed_floor), &
    test(file_backed_rankine_forcing) &
  ])

  if (.not. suite % ok) error stop 1

contains

  function interpolation_snapshot() result(res)
    type(test_result) :: res
    type(config_type) :: config
    type(grid_type) :: grid
    type(forcing_type) :: forcing
    logical :: ok
    real, allocatable :: expected_wspd(:,:), expected_wdir(:,:), expected_uc(:,:), expected_vc(:,:)
    real, allocatable :: expected_wspd_i(:), expected_wdir_i(:), expected_uc_i(:), expected_vc_i(:)
    real, allocatable :: expected_rhoa(:), expected_rhow(:), expected_rhorat(:), expected_fice(:,:)
    real, allocatable :: expected_fice_i(:)

    config = valid_config()
    config % winds = .true.
    config % currents = .true.
    config % air_density = .true.
    config % water_density = .true.
    config % seaice = .true.
    config % dtg = 10.0
    config % gustiness = 0.0

    call grid % initialize(config)
    call forcing % initialize(config, grid)

    allocate(expected_wspd(config % mm, config % nm))
    allocate(expected_wdir(config % mm, config % nm))
    allocate(expected_uc(config % mm, config % nm))
    allocate(expected_vc(config % mm, config % nm))
    allocate(expected_wspd_i(grid % imm))
    allocate(expected_wdir_i(grid % imm))
    allocate(expected_uc_i(grid % imm))
    allocate(expected_vc_i(grid % imm))
    allocate(expected_rhoa(grid % imm))
    allocate(expected_rhow(grid % imm))
    allocate(expected_rhorat(grid % imm))
    allocate(expected_fice(config % mm, config % nm))
    allocate(expected_fice_i(grid % imm))

    forcing % uwb = 1.0
    forcing % vwb = 0.0
    forcing % uwf = 2.0 * forcing % uwb
    forcing % vwf = 0.0
    forcing % ucb = 1.0
    forcing % vcb = 2.0
    forcing % ucf = 3.0
    forcing % vcf = 4.0
    forcing % ficeb = 0.1
    forcing % ficef = 0.9
    forcing % rhoab = 1.2
    forcing % rhoaf = 1.4
    forcing % rhowb = 1.0
    forcing % rhowf = 1.1

    sumt = 5.0
    call forcing % interpolate(config, grid)

    expected_wspd = 1.5
    expected_wdir = 0.0
    expected_uc = 2.0
    expected_vc = 3.0
    expected_wspd_i = 1.5
    expected_wdir_i = 0.0
    expected_uc_i = 2.0
    expected_vc_i = 3.0
    expected_rhoa = 1.3
    expected_rhow = 1.05
    expected_rhorat = 1.3 / 1.05
    expected_fice = 0.5
    expected_fice_i = 0.5

    ok = all_nearly_equal(reshape(forcing % wspd_2d, [size(forcing % wspd_2d)]), &
                           reshape(expected_wspd, [size(expected_wspd)])) .and. &
         all_nearly_equal(reshape(forcing % wdir_2d, [size(forcing % wdir_2d)]), &
                           reshape(expected_wdir, [size(expected_wdir)])) .and. &
         all_nearly_equal(reshape(forcing % uc_2d, [size(forcing % uc_2d)]), &
                           reshape(expected_uc, [size(expected_uc)])) .and. &
         all_nearly_equal(reshape(forcing % vc_2d, [size(forcing % vc_2d)]), &
                           reshape(expected_vc, [size(expected_vc)])) .and. &
         all_nearly_equal(forcing % wspd, expected_wspd_i) .and. &
         all_nearly_equal(forcing % wdir, expected_wdir_i) .and. &
         all_nearly_equal(forcing % uc, expected_uc_i) .and. &
         all_nearly_equal(forcing % vc, expected_vc_i) .and. &
         all_nearly_equal(forcing % rhoa, expected_rhoa) .and. &
         all_nearly_equal(forcing % rhow, expected_rhow) .and. &
         all_nearly_equal(forcing % rhorat, expected_rhorat) .and. &
         all_nearly_equal(reshape(forcing % fice_2d, [size(forcing % fice_2d)]), &
                           reshape(expected_fice, [size(expected_fice)])) .and. &
         all_nearly_equal(forcing % fice, expected_fice_i)

    call forcing % finalize()
    call grid % finalize()

    res = test('interpolation_snapshot', ok)
  end function interpolation_snapshot

  function wind_speed_floor() result(res)
    type(test_result) :: res
    type(config_type) :: config
    type(grid_type) :: grid
    type(forcing_type) :: forcing
    logical :: ok

    config = valid_config()
    config % winds = .true.
    config % dtg = 10.0
    config % gustiness = 0.0

    call grid % initialize(config)
    call forcing % initialize(config, grid)

    forcing % uwb = 0.0
    forcing % vwb = 0.0
    forcing % uwf = 0.0
    forcing % vwf = 0.0

    sumt = 5.0
    call forcing % interpolate(config, grid)

    ok = all(abs(forcing % wspd - min_wind_speed) < 1e-6)

    call forcing % finalize()
    call grid % finalize()

    res = test('wind_speed_floor', ok)
  end function wind_speed_floor

  function stored_wind_speed_floor() result(res)
    type(test_result) :: res
    type(config_type) :: config
    type(grid_type) :: grid
    type(forcing_type) :: forcing
    logical :: ok

    config = valid_config()
    config % winds = .false.
    config % wspd0 = 0.0

    call grid % initialize(config)
    call forcing % initialize(config, grid)

    ok = all(abs(forcing % wspd - min_wind_speed) < 1e-6)

    call forcing % load(config, config % starttimestr, grid)
    ok = ok .and. all(abs(forcing % wspd - min_wind_speed) < 1e-6)

    forcing % wspd = 0.0
    call forcing % apply_wind_speed_floor()
    ok = ok .and. all(abs(forcing % wspd - min_wind_speed) < 1e-6)

    call forcing % finalize()
    call grid % finalize()

    res = test('stored_wind_speed_floor', ok)
  end function stored_wind_speed_floor

  function file_backed_rankine_forcing() result(res)
    type(test_result) :: res
    type(config_type) :: config
    type(grid_type) :: grid
    type(forcing_type) :: forcing
    real, allocatable :: uw0(:,:), vw0(:,:), uc0(:,:), vc0(:,:)
    real, allocatable :: rhoa0(:), rhow0(:)
    real, allocatable :: expected_uw(:,:), expected_vw(:,:)
    real, allocatable :: expected_uc_2d(:,:), expected_vc_2d(:,:)
    real, allocatable :: expected_wspd_2d(:,:), expected_wdir_2d(:,:)
    real, allocatable :: expected_wspd(:), expected_wdir(:)
    real, allocatable :: expected_uc(:), expected_vc(:)
    real, allocatable :: expected_rhoa(:), expected_rhow(:), expected_rhorat(:)
    logical :: ok

    config = rankine_config()

    call grid % initialize(config)
    call forcing % initialize(config, grid)

    call forcing % load(config, config % starttimestr, grid)
    uw0 = forcing % uwf
    vw0 = forcing % vwf
    uc0 = forcing % ucf
    vc0 = forcing % vcf
    rhoa0 = forcing % rhoaf
    rhow0 = forcing % rhowf

    call forcing % update(config, '2026-01-01 01:00:00', grid)

    expected_uw = 0.5 * (uw0 + forcing % uwf)
    expected_vw = 0.5 * (vw0 + forcing % vwf)
    expected_uc_2d = 0.5 * (uc0 + forcing % ucf)
    expected_vc_2d = 0.5 * (vc0 + forcing % vcf)
    expected_wspd_2d = sqrt(expected_uw**2 + expected_vw**2)
    expected_wdir_2d = atan2(expected_vw, expected_uw)
    expected_wspd = max(grid % remap_mn2i(expected_wspd_2d), 1e-2)
    expected_wdir = grid % remap_mn2i(expected_wdir_2d)
    expected_uc = grid % remap_mn2i(expected_uc_2d)
    expected_vc = grid % remap_mn2i(expected_vc_2d)
    expected_rhoa = 0.5 * (rhoa0 + forcing % rhoaf)
    expected_rhow = 0.5 * (rhow0 + forcing % rhowf)
    expected_rhorat = expected_rhoa / expected_rhow

    sumt = 0.5 * config % dtg
    call forcing % interpolate(config, grid)

    ok = maxval(abs(uw0 - forcing % uwf)) > 1.0 .and. &
         maxval(abs(vw0 - forcing % vwf)) > 1.0 .and. &
         maxval(expected_wspd_2d) > 20.0 .and. &
         all_nearly_equal(reshape(forcing % uw, [size(forcing % uw)]), &
                          reshape(expected_uw, [size(expected_uw)])) .and. &
         all_nearly_equal(reshape(forcing % vw, [size(forcing % vw)]), &
                          reshape(expected_vw, [size(expected_vw)])) .and. &
         all_nearly_equal(reshape(forcing % wspd_2d, [size(forcing % wspd_2d)]), &
                          reshape(expected_wspd_2d, [size(expected_wspd_2d)])) .and. &
         all_nearly_equal(reshape(forcing % wdir_2d, [size(forcing % wdir_2d)]), &
                          reshape(expected_wdir_2d, [size(expected_wdir_2d)])) .and. &
         all_nearly_equal(reshape(forcing % uc_2d, [size(forcing % uc_2d)]), &
                          reshape(expected_uc_2d, [size(expected_uc_2d)])) .and. &
         all_nearly_equal(reshape(forcing % vc_2d, [size(forcing % vc_2d)]), &
                          reshape(expected_vc_2d, [size(expected_vc_2d)])) .and. &
         all_nearly_equal(forcing % wspd, expected_wspd) .and. &
         all_nearly_equal(forcing % wdir, expected_wdir) .and. &
         all_nearly_equal(forcing % uc, expected_uc) .and. &
         all_nearly_equal(forcing % vc, expected_vc) .and. &
         all_nearly_equal(forcing % rhoa, expected_rhoa) .and. &
         all_nearly_equal(forcing % rhow, expected_rhow) .and. &
         all_nearly_equal(forcing % rhorat, expected_rhorat)

    call forcing % finalize()
    call grid % finalize()

    res = test('file_backed_rankine_forcing', ok)
  end function file_backed_rankine_forcing

  function valid_config() result(config)
    type(config_type) :: config
    config = config_type('../namelists/main.nml')
  end function valid_config

  function rankine_config() result(config)
    type(config_type) :: config
    config = config_type('namelists/forcing_rankine.nml')
  end function rankine_config

end program test_forcing
