module umwm_dispersion

  implicit none

  private
  public :: frequency, wavenumber

  real, parameter :: twopi = 2 * acos(-1.0d0)

contains

  elemental real function wavenumber( &
    frequency, depth, water_density, gravity, surface_tension &
  ) result(k)
    ! Solve the linear water wave dispersion relationship using
    ! Newton-Raphson iteration.
    real, intent(in) :: frequency
    real, intent(in) :: depth
    real, intent(in) :: water_density
    real, intent(in) :: gravity
    real, intent(in) :: surface_tension

    real :: dk, b, frequency_nondim, tanhk
    integer, parameter :: max_iter = 100
    real, parameter :: eps = 1e-8
    integer :: count

    frequency_nondim = twopi * frequency * sqrt(depth / gravity)
    k = frequency_nondim**2
    b = surface_tension / (water_density * gravity * depth**2)

    count = 0
    dk = 2e-3
    do while (count < max_iter)
      tanhk = tanh(k)
      dk = - (frequency_nondim**2 - k * tanhk * (1 + b * k**2)) &
        / (3 * b * k**2 * tanhk + tanhk + k * (1 + b * k**2) * (1 - tanhk**2))
      k = k - dk

      if (abs(dk) < eps) exit
      count = count + 1
    end do

    k = k / depth

  end function wavenumber


  elemental real function frequency( &
    wavenumber, depth, water_density, gravity, surface_tension &
  )
    ! Return the (non-angular) frequency using the linear water wave dispersion
    ! relationship.
    real, intent(in) :: wavenumber
    real, intent(in) :: depth
    real, intent(in) :: water_density
    real, intent(in) :: gravity
    real, intent(in) :: surface_tension

    frequency = sqrt( &
      (gravity * wavenumber + surface_tension * wavenumber**3 / water_density) &
      * tanh(wavenumber * depth) &
    ) / twopi

  end function frequency

end module umwm_dispersion