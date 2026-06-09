module umwm_dispersion

  use umwm_constants, only: rk, twopi

  implicit none

  private
  public :: angular_frequency, group_speed, wavenumber

contains

  elemental real(rk) function angular_frequency( &
    wavenumber, depth, water_density, gravity, surface_tension &
  )
    ! Return angular frequency using the linear water wave dispersion relation.
    real(rk), intent(in) :: wavenumber
    real(rk), intent(in) :: depth
    real(rk), intent(in) :: water_density
    real(rk), intent(in) :: gravity
    real(rk), intent(in) :: surface_tension

    angular_frequency = sqrt( &
      (gravity * wavenumber + surface_tension * wavenumber**3 / water_density) &
      * tanh(wavenumber * depth) &
    )

  end function angular_frequency


  elemental real(rk) function group_speed( &
    wavenumber, depth, water_density, gravity, surface_tension &
  )
    real(rk), intent(in) :: wavenumber
    real(rk), intent(in) :: depth
    real(rk), intent(in) :: water_density
    real(rk), intent(in) :: gravity
    real(rk), intent(in) :: surface_tension

    real(rk) :: cp, kd, kd_limited, omega

    omega = angular_frequency(wavenumber, depth, water_density, gravity, surface_tension)
    cp = omega / wavenumber
    kd = wavenumber * depth
    kd_limited = min(kd, 20._rk)

    group_speed = cp * ( &
      0.5_rk &
      + kd_limited / sinh(2._rk * kd_limited) &
      + surface_tension * wavenumber**2 &
        / (water_density * gravity + surface_tension * wavenumber**2) &
    )

  end function group_speed


  elemental real(rk) function wavenumber( &
    frequency, depth, water_density, gravity, surface_tension &
  ) result(k)
    ! Solve the linear water wave dispersion relationship using
    ! Newton-Raphson iteration.
    real(rk), intent(in) :: frequency
    real(rk), intent(in) :: depth
    real(rk), intent(in) :: water_density
    real(rk), intent(in) :: gravity
    real(rk), intent(in) :: surface_tension

    real(rk) :: dk, b, frequency_nondim, tanhk
    integer, parameter :: max_iter = 100
    real(rk), parameter :: eps = 1e-8_rk
    integer :: count

    frequency_nondim = twopi * frequency * sqrt(depth / gravity)
    k = frequency_nondim**2
    b = surface_tension / (water_density * gravity * depth**2)

    count = 0
    dk = 2e-3_rk
    do while (count < max_iter)
      tanhk = tanh(k)
      dk = - (frequency_nondim**2 - k * tanhk * (1._rk + b * k**2)) &
        / (3._rk * b * k**2 * tanhk + tanhk &
        + k * (1._rk + b * k**2) * (1._rk - tanhk**2))
      k = k - dk

      if (abs(dk) < eps) exit
      count = count + 1
    end do

    k = abs(k) / depth

  end function wavenumber

end module umwm_dispersion
