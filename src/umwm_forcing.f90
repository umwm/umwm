module umwm_forcing

  use umwm_config, only: config_type
  use umwm_io, only: input_nc, readfile
  use umwm_module, only: fice, fice_2d, ficeb, ficef, &
                         gustu, gustv, rhoa, rhoab, rhoaf, rhorat, &
                         rhow, rhowb, rhowf, sumt, uc, uc_2d, ucb, ucf, &
                         uw, uwb, uwf, vc, vc_2d, vcb, vcf, vw, vwb, &
                         vwf, wdir, wdir_2d, wspd, wspd_2d
  use umwm_util, only: remap_mn2i

  implicit none

contains

  subroutine forcinginput(config, timestr)
    ! Reads atmospheric and oceanic input forcing fields
    type(config_type), intent(in) :: config
    character(len=19), intent(in) :: timestr

    ! save wind at time level n
    if (config % winds) then
      uwb = uwf
      vwb = vwf
    end if

    ! save currents at time level n
    if (config % currents) then
      ucb = ucf
      vcb = vcf
    end if

    ! save air and water density at time level n
    if (config % air_density) rhoab = rhoaf
    if (config % water_density) rhowb = rhowf

    ! save sea ice fraction at time level n
    if (config % seaice) ficeb = ficef

    ! load input fields at time level n+1
    if (readfile) call input_nc(config, timestr)

  end subroutine forcinginput


  subroutine forcinginterpolate(config)
    ! interpolates in time atmospheric and oceanic input forcing fields
    type(config_type), intent(in) :: config

    if (config % winds) then

      uw = uwb * (1 - sumt / config % dtg) + uwf * sumt / config % dtg
      vw = vwb * (1 - sumt / config % dtg) + vwf * sumt / config % dtg

      ! add wind gustiness if requested in the namelist;
      ! uw, vw will get a uniformly distributed gust component
      if (config % gustiness > 0) then

        ! get random numbers [0, 1]
        call random_number(gustu)
        call random_number(gustv)

        gustu = config % gustiness * (2 * gustu - 1)
        gustv = config % gustiness * (2 * gustv - 1)

        ! add gustiness to the wind fields
        uw = uw * (1 + gustu)
        vw = vw * (1 + gustv)

      end if

      wspd_2d = sqrt(uw**2 + vw**2)
      wdir_2d = atan2(vw, uw)

      ! remap to 1-d arrays:
      wspd = remap_mn2i(wspd_2d)
      wdir = remap_mn2i(wdir_2d)

    end if ! winds

    if (config % seaice) then

      fice_2d = ficeb * (1 - sumt / config % dtg) + ficef * sumt/config % dtg

      ! remap to 1-D arrays:
      fice = remap_mn2i(fice_2d)

     end if
    

    if (config % currents) then

      uc_2d = ucb * (1 - sumt / config % dtg) + ucf * sumt / config % dtg
      vc_2d = vcb * (1 - sumt / config % dtg) + vcf * sumt / config % dtg

      ! remap to 1-d arrays:
      uc = remap_mn2i(uc_2d)
      vc = remap_mn2i(vc_2d)

    end if

    if (config % air_density) rhoa = rhoab * (1 - sumt / config % dtg) + rhoaf * sumt / config % dtg
    if (config % water_density) rhow = rhowb * (1 - sumt / config % dtg) + rhowf * sumt / config % dtg

    rhorat = rhoa / rhow

  end subroutine forcinginterpolate

end module umwm_forcing
