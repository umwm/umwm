program umwm

  use umwm_cli, only: parse_command_line
  use umwm_domain, only: domain_type
  use umwm_env, only: env_init, env_stop

  implicit none

  type(domain_type) :: domain
  character(:), allocatable :: namelist
  logical :: run

  call env_init()
  call parse_command_line(namelist, run)
  if (run) then
    call domain % initialize(namelist)
    call domain % run()
    call domain % finalize()
  end if
  call env_stop()

end program umwm
