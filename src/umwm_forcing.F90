module umwm_forcing

  use umwm_io, only: input_nc, readfile, winds, currents, air_density, water_density, seaice
  use umwm_module
  use umwm_util, only: remap_mn2i

  implicit none

contains

  subroutine forcinginput(timestr)
    ! Reads atmospheric and oceanic input forcing fields
    character(len=19), intent(in) :: timestr

    ! save wind at time level n
    if (winds) then
      uwb = uwf
      vwb = vwf
    end if

    ! save currents at time level n
    if (currents) then
      ucb = ucf
      vcb = vcf
    end if

    ! save air and water density at time level n
    if (air_density) rhoab = rhoaf
    if (water_density) rhowb = rhowf

    ! save sea ice fraction at time level n
    if (seaice) ficeb = ficef

    ! load input fields at time level n+1
    if (readfile) call input_nc(timestr)

  end subroutine forcinginput


  subroutine forcinginterpolate()
    ! interpolates in time atmospheric and oceanic input forcing fields
    if (winds) then

      uw = uwb * (1 - sumt / dtg) + uwf * sumt / dtg
      vw = vwb * (1 - sumt / dtg) + vwf * sumt / dtg

      ! add wind gustiness if requested in the namelist;
      ! uw, vw will get a uniformly distributed gust component
      if (gustiness > 0) then

        ! get random numbers [0, 1]
        call random_number(gustu)
        call random_number(gustv)

        gustu = gustiness * (2 * gustu - 1)
        gustv = gustiness * (2 * gustv - 1)

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

    if (seaice) then

      fice_2d = ficeb * (1 - sumt / dtg) + ficef * sumt/dtg

      ! remap to 1-D arrays:
      fice = remap_mn2i(fice_2d)

     end if
    

    if (currents) then

      uc_2d = ucb * (1 - sumt / dtg) + ucf * sumt / dtg
      vc_2d = vcb * (1 - sumt / dtg) + vcf * sumt / dtg

      ! remap to 1-d arrays:
      uc = remap_mn2i(uc_2d)
      vc = remap_mn2i(vc_2d)

    end if

    if (air_density) rhoa = rhoab * (1 - sumt / dtg) + rhoaf * sumt / dtg
    if (water_density) rhow = rhowb * (1 - sumt / dtg) + rhowf * sumt / dtg

    rhorat = rhoa / rhow

  end subroutine forcinginterpolate

end module umwm_forcing
