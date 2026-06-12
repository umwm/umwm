program umwm

  use umwm_config, only: config_type
  use umwm_constants, only: rk
  use umwm_env, only: env_init
  use umwm_grid, only: grid_type
  use umwm_module, only: nproc
  use umwm_spectrum, only: spectrum_type
  use umwm_top, only: umwm_initialize, umwm_run, umwm_finalize

  type(config_type) :: config
  type(grid_type) :: grid
  type(spectrum_type) :: spectrum
  logical :: ok

  call env_init()
  config = config_type('namelists/main.nml')
  ok = config % validate(rank=nproc)
  if (.not. ok) error stop 1
  spectrum = spectrum_type(config % om, config % pm, &
    real(config % fmin, rk), real(config % fmax, rk))

  call grid % initialize(config)
  call grid % print_diagnostics()
  call grid % initialize_direction_projection(spectrum)
  call umwm_initialize(config, spectrum, grid)
  call umwm_run(config, spectrum, grid)
  call umwm_finalize(grid)

end program umwm
