module umwm_dispersion

  implicit none

  private
  public :: wavenumber

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
    real, parameter :: eps = 1e-5
    integer :: counter

    real, parameter :: twopi = 2 * acos(-1.)

    frequency_nondim = twopi * frequency * sqrt(depth / gravity)
    k = frequency_nondim**2
    b = surface_tension / (water_density * gravity * depth**2)

    counter = 1
    dk = 2e-3
    do
      tanhk = tanh(k)
      dk = - (frequency_nondim**2 - k * tanhk * (1 + b * k**2)) &
        / (3 * b * k**2 * tanhk + tanhk + k * (1 + b * k**2) * (1 - tanhk**2))
      k = k - dk

      if (abs(dk) < eps .or. counter > 100) exit
      counter = counter+1

    end do

    k = k / depth

  end function wavenumber

end module umwm_dispersion