program umwm

  use umwm_module, only: starttimestr => starttimestr_nml,&
                         stoptimestr  => stoptimestr_nml
  use umwm_top,only: umwm_initialize, umwm_run, umwm_finalize

  call umwm_initialize()
  call umwm_run(starttimestr, stoptimestr)
  call umwm_finalize()

end program umwm
