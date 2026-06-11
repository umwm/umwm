program umwm

  use umwm_config, only: config_type
  use umwm_constants, only: rk
  use umwm_env, only: env_init
  use umwm_module, only: nproc
  use umwm_spectrum, only: spectrum_type
  use umwm_top, only: umwm_initialize, umwm_run, umwm_finalize

  type(config_type) :: config
  type(spectrum_type) :: spectrum
  logical :: ok

  call env_init()
  config = config_type('namelists/main.nml', rank=nproc)
  ok = config % validate(rank=nproc)
  if (.not. ok) error stop 1
  spectrum = spectrum_type(config % om, config % pm, &
    real(config % fmin, rk), real(config % fmax, rk))

  call umwm_initialize(config, spectrum)
  call umwm_run(config, spectrum)
  call umwm_finalize()

end program umwm
