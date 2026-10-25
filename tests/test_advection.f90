program test_advection
  use tuff, only: test, test_result
  use umwm_advection, only: current_rotation_rate, depth_rotation_rate, refraction
  use umwm_constants, only: rk, stderr, twopi
  use umwm_dispersion, only: group_speed, wavenumber
#ifdef MPI
  use umwm_env, only: env_init, env_stop
#endif

  implicit none

  integer, parameter :: nf = 3, nd = 36
  real, parameter :: dth = real(twopi) / nd

  type(test_result) :: suite

#ifdef MPI
  call env_init()
#endif
  suite = test('test_advection', [ &
    test(current_rotation_follows_ray_equation), &
    test(current_rotation_is_half_vorticity_in_solid_body_eddy), &
    test(depth_rotation_follows_snell_law), &
    test(refraction_substeps_conserve_and_stay_positive), &
    test(refraction_substeps_rotate_by_the_whole_step) &
  ])
#ifdef MPI
  call env_stop()
#endif
  if (.not. suite % ok) error stop 1

contains

  function current_rotation_follows_ray_equation() result(res)
    ! The rotation of the wavenumber direction from dk/dt = -grad(k.U), for
    ! a current with vorticity and strain.
    type(test_result) :: res
    real, parameter :: du_dx = 3e-5, du_dy = -7e-5, dv_dx = 1.1e-4, dv_dy = 5e-5
    real :: th, kx, ky, dkx_dt, dky_dt, expected
    integer :: p
    logical :: ok
    ok = .true.
    do p = 1, nd
      th = (p - 1) * dth
      kx = cos(th)
      ky = sin(th)
      dkx_dt = -(kx * du_dx + ky * dv_dx)
      dky_dt = -(kx * du_dy + ky * dv_dy)
      expected = kx * dky_dt - ky * dkx_dt
      ok = ok .and. abs(current_rotation_rate(du_dx, du_dy, dv_dx, dv_dy, kx, ky) - expected) &
                    <= 1e-5 * abs(dv_dx)
    end do
    res = test('current_rotation_follows_ray_equation', ok)
  end function current_rotation_follows_ray_equation


  function current_rotation_is_half_vorticity_in_solid_body_eddy() result(res)
    ! A current in solid-body rotation, (u, v) = omega * (-y, x), rotates
    ! the waves of all directions at omega, half its vorticity.
    type(test_result) :: res
    real, parameter :: omega = 1e-4
    real :: th
    integer :: p
    logical :: ok
    ok = .true.
    do p = 1, nd
      th = (p - 1) * dth
      ok = ok .and. abs(current_rotation_rate(0., -omega, omega, 0., cos(th), sin(th)) - omega) &
                    <= 1e-5 * omega
    end do
    res = test('current_rotation_is_half_vorticity_in_solid_body_eddy', ok)
  end function current_rotation_is_half_vorticity_in_solid_body_eddy


  function depth_rotation_follows_snell_law() result(res)
    ! Integrates a ray of 10 s waves at 45 degrees to the depth gradient
    ! across straight and parallel depth contours, from 100 m to 5 m depth,
    ! with the rotation rate from the phase speed gradient at fixed
    ! frequency times cg/cp, as in the model. Snell's law, sin(th)/cp =
    ! const, holds along the ray.
    type(test_result) :: res
    integer, parameter :: num_steps = 2000
    real(rk), parameter :: frequency = 0.1_rk, length = 20e3_rk
    real(rk) :: x, dx, th, th_mid, snell_start, snell_end
    integer :: n
    dx = length / num_steps
    x = 0
    th = twopi / 8
    snell_start = sin(th) / phase_speed(x)
    do n = 1, num_steps
      th_mid = th + 0.5_rk * dx * dth_dx(x, th)
      th = th + dx * dth_dx(x + 0.5_rk * dx, th_mid)
      x = x + dx
    end do
    snell_end = sin(th) / phase_speed(x)
    res = test('depth_rotation_follows_snell_law', abs(snell_end / snell_start - 1) < 1e-3_rk)
    if (.not. res % ok) write(stderr, '(a,f10.6)') &
      'test_advection: ratio of sin(th)/cp at the end and start of the ray: ', &
      snell_end / snell_start
  end function depth_rotation_follows_snell_law


  real(rk) function dth_dx(x, th)
    ! Rate of change of the ray direction th with x along the ray, which
    ! travels at the group speed.
    real(rk), intent(in) :: x, th
    real(rk), parameter :: h = 10
    real(rk) :: dcp_dx, rate
    dcp_dx = (phase_speed(x + h) - phase_speed(x - h)) / (2 * h)
    rate = depth_rotation_rate(real(group_speed_at(x) / phase_speed(x) * dcp_dx), 0., &
                               real(cos(th)), real(sin(th)))
    dth_dx = rate / (group_speed_at(x) * cos(th))
  end function dth_dx


  real(rk) function depth(x)
    real(rk), intent(in) :: x
    depth = 100 - 95 * x / 20e3_rk
  end function depth


  real(rk) function phase_speed(x)
    real(rk), intent(in) :: x
    phase_speed = twopi * 0.1_rk / wavenumber(0.1_rk, depth(x), 1025._rk, 9.80665_rk, 0._rk)
  end function phase_speed


  real(rk) function group_speed_at(x)
    real(rk), intent(in) :: x
    group_speed_at = group_speed(wavenumber(0.1_rk, depth(x), 1025._rk, 9.80665_rk, 0._rk), &
                                 depth(x), 1025._rk, 9.80665_rk, 0._rk)
  end function group_speed_at


  function refraction_substeps_conserve_and_stay_positive() result(res)
    ! With a current gradient that needs several substeps in a time step,
    ! refraction keeps the energy of each frequency and a positive spectrum,
    ! and leaves frequencies above the highest refracted bin unchanged.
    type(test_result) :: res
    real :: e(nf,nd), ef(nf,nd), cth(nd), sth(nd)
    real :: zero(nf)
    integer :: o, p
    logical :: ok
    call directions(cth, sth)
    do p = 1, nd
      do o = 1, nf
        e(o,p) = 1 + 0.5 * sin(3. * p + o) + merge(10., 0., p == 7)
      end do
    end do
    zero = 0
    ef = e
    ! strain and shear of 1e-3 1/s, which turn the waves by up to 9 bins in 1000 s
    call refraction(nf, nd, nf - 1, e, zero, zero, 1e-3, 0., 1e-3, -1e-3, cth, sth, &
                    1000., dth, ef)
    ok = all(ef >= 0) .and. all(ef(nf,:) == e(nf,:)) .and. any(ef(1,:) /= e(1,:))
    do o = 1, nf - 1
      ok = ok .and. abs(sum(ef(o,:)) / sum(e(o,:)) - 1) < 1e-5
    end do
    res = test('refraction_substeps_conserve_and_stay_positive', ok)
  end function refraction_substeps_conserve_and_stay_positive


  function refraction_substeps_rotate_by_the_whole_step() result(res)
    ! A solid-body eddy rotates the waves of all directions at the same
    ! rate, anti-clockwise for an anti-clockwise eddy. In a time step of
    ! 2.5 bins, which takes 3 substeps, the mean direction of a narrow
    ! spectrum turns by 2.5 bins.
    type(test_result) :: res
    real, parameter :: omega = 1e-4
    real :: e(nf,nd), ef(nf,nd), cth(nd), sth(nd)
    real :: zero(nf), dt, shift
    integer :: p
    call directions(cth, sth)
    e = 0
    e(:,10:12) = 1
    zero = 0
    ef = e
    dt = 2.5 * dth / omega
    call refraction(nf, nd, nf, e, zero, zero, 0., -omega, omega, 0., cth, sth, dt, dth, ef)
    shift = sum([(p * ef(1,p), p = 1, nd)]) / sum(ef(1,:)) &
          - sum([(p * e(1,p), p = 1, nd)]) / sum(e(1,:))
    res = test('refraction_substeps_rotate_by_the_whole_step', abs(shift - 2.5) < 1e-4)
    if (.not. res % ok) write(stderr, '(a,f10.6)') &
      'test_advection: mean direction turned by bins: ', shift
  end function refraction_substeps_rotate_by_the_whole_step


  subroutine directions(cth, sth)
    real, intent(out) :: cth(nd), sth(nd)
    integer :: p
    do p = 1, nd
      cth(p) = cos((p - 1) * dth)
      sth(p) = sin((p - 1) * dth)
    end do
  end subroutine directions

end program test_advection
