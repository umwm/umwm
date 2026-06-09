program umwm

  use umwm_module, only: starttimestr => starttimestr_nml,&
                         stoptimestr  => stoptimestr_nml
  use umwm_spectrum, only: spectrum_type
  use umwm_top, only: umwm_initialize, umwm_run, umwm_finalize

  type(spectrum_type) :: spectrum

  call umwm_initialize(spectrum)
  call umwm_run(starttimestr, stoptimestr, spectrum)
  call umwm_finalize()

end program umwm
