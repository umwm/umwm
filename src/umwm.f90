program umwm

  use umwm_domain, only: domain_type
  use umwm_env, only: env_init, env_stop

  implicit none

  type(domain_type) :: domain

  call env_init()
  call domain % initialize('namelists/main.nml')
  call domain % run()
  call domain % finalize()
  call env_stop()

end program umwm
