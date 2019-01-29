module umwm_constants

  use iso_fortran_env,  only: real32, real64, real128, output_unit, error_unit

  implicit none

#ifdef REAL128
  integer, parameter :: rk = real128
#elif REAL64
  integer, parameter :: rk = real64
#else
  integer, parameter :: rk = real32
#endif

  real(rk), parameter :: one = 1
  real(rk), parameter :: two = 2
  real(rk), parameter :: three = 3
  real(rk), parameter :: half = 0.5_rk
  real(rk), parameter :: quart = 0.25_rk
  real(rk), parameter :: onethird = one / three
  real(rk), parameter :: pi = 4 * atan(1._rk) ! pi
  real(rk), parameter :: euler = exp(1d0)     ! e
  real(rk), parameter :: eulerinv = 1 / euler ! 1/e
  real(rk), parameter :: invpi  = one / pi ! 1/pi
  real(rk), parameter :: twopi  = two * pi ! 2*pi
  real(rk), parameter :: twopi2 = twopi**2 ! 4*pi^2
  real(rk), parameter :: dr     = pi / 180 ! deg -> rad
  real(rk), parameter :: r_earth = 6.371009e6_rk ! earth radius
  real(rk), parameter :: tiny_real = tiny(one)
  real(rk), parameter :: huge_real = huge(one)

  integer, parameter :: stdout = output_unit
  integer, parameter :: stderr = error_unit

end module umwm_constants
