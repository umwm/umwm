program umwm

  use umwm_cli, only: parse_command_line, prepare_inputs, write_namelist_copies
  use umwm_domain, only: domain_type
  use umwm_env, only: env_init, env_stop

  implicit none

  type(domain_type) :: domain
  character(:), allocatable :: namelist, input_dir, output_dir
  logical :: run

  call env_init()
  call parse_command_line(namelist, input_dir, output_dir, run)
  if (run) then
    call prepare_inputs(namelist, input_dir)
    call domain % initialize(namelist, input_dir=input_dir, output_dir=output_dir)
    call write_namelist_copies(namelist, output_dir)
    call domain % run()
    call domain % finalize()
  end if
  call env_stop()

end program umwm
