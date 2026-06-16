program test_forcing
  use tuff, only: test, test_result, nearly_equal, all_nearly_equal
  use umwm_config, only: config_type
  use umwm_forcing, only: forcing_type
  use umwm_grid, only: grid_type
  use umwm_module, only: sumt

  implicit none

  type(test_result) :: suite

  suite = test('test_forcing', [ &
    test(interpolation_snapshot) &
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

  function valid_config() result(config)
    type(config_type) :: config
    config = config_type('../namelists/main.nml')
  end function valid_config

end program test_forcing
