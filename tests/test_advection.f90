program test_advection
  use tuff, only: test, test_result
  use umwm_advection, only: current_rotation_rate, depth_rotation_rate, frequency_shift, &
                            propagation, refraction
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
    test(refraction_substeps_rotate_by_the_whole_step), &
    test(frequency_shift_conserves_action_and_shifts_by_the_whole_step), &
    test(frequency_shift_exchanges_action_with_the_diagnostic_tail), &
    test(propagation_shoals_by_linear_theory), &
    test(action_flux_is_conserved_on_opposing_current) &
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


  function frequency_shift_conserves_action_and_shifts_by_the_whole_step() result(res)
    ! Deep-water waves against a current that strengthens along them, with
    ! du/dx = -1e-4 1/s, shorten at d(ln f)/dt = (cg/cp) 1e-4 = 5e-5 1/s. In a
    ! time step of 2.5 bins, which takes 3 substeps, the action of a narrow
    ! spectrum away from the ends of the frequency range is conserved, and
    ! its mean ln(f) rises by 2.5 bins.
    type(test_result) :: res
    integer, parameter :: nf = 30
    real, parameter :: dlnf = 0.05, du_dx = -1e-4
    real :: e(nf,1), ef(nf,1), kovcg(nf), invcp(nf), cg(nf), zero(nf), a(nf), a_new(nf)
    real :: sigma, k, dt, shift
    integer :: o
    do o = 1, nf
      sigma = twopi * 0.1 * exp((o - 1) * dlnf)
      k = sigma**2 / 9.81
      invcp(o) = k / sigma
      cg(o) = 0.5 * sigma / k
      kovcg(o) = k / cg(o)
    end do
    e = 0
    e(10:12,1) = 1 / kovcg(10:12)
    ef = e
    zero = 0
    dt = 2.5 * dlnf / (0.5 * abs(du_dx))
    call frequency_shift(nf, 1, nf, e, kovcg, invcp, cg, zero, zero, 0., 0., du_dx, 0., 0., 0., &
                         [1.], [0.], dt, dlnf, ef)
    a = kovcg * e(:,1)
    a_new = kovcg * ef(:,1)
    shift = sum([(o * a_new(o), o = 1, nf)]) / sum(a_new) - sum([(o * a(o), o = 1, nf)]) / sum(a)
    res = test('frequency_shift_conserves_action_and_shifts_by_the_whole_step', &
               abs(sum(a_new) / sum(a) - 1) < 1e-5 .and. abs(shift - 2.5) < 1e-3 &
               .and. all(ef >= 0))
    if (.not. res % ok) write(stderr, '(a,2f10.6)') &
      'test_advection: action ratio and shift in bins: ', sum(a_new) / sum(a), shift
  end function frequency_shift_conserves_action_and_shifts_by_the_whole_step


  function frequency_shift_exchanges_action_with_the_diagnostic_tail() result(res)
    ! With 6 prognostic bins of 10, deep-water waves against a current that
    ! strengthens along them, du/dx = -1e-4 1/s, shorten at d(ln f)/dt =
    ! 5e-5 1/s: in a time step of 0.2 bins, the highest prognostic bin loses
    ! 0.2 of its action to the diagnostic tail. With the current reversed,
    ! the waves lengthen, and it gains 0.2 of the action of the lowest bin
    ! of the tail. The bins of the tail are left unchanged.
    type(test_result) :: res
    integer, parameter :: nf = 10, no = 6
    real, parameter :: dlnf = 0.05, rate = 5e-5
    real :: e(nf,1), ef(nf,1), kovcg(nf), invcp(nf), cg(nf), zero(nf), sigma, k, dt
    integer :: o
    logical :: ok
    do o = 1, nf
      sigma = twopi * 0.1 * exp((o - 1) * dlnf)
      k = sigma**2 / 9.81
      invcp(o) = k / sigma
      cg(o) = 0.5 * sigma / k
      kovcg(o) = k / cg(o)
    end do
    zero = 0
    dt = 0.2 * dlnf / rate

    ! shortening: out of bin no into the tail
    e = 0
    e(no,1) = 1
    ef = e
    call frequency_shift(nf, 1, no, e, kovcg, invcp, cg, zero, zero, 0., 0., -2 * rate, 0., &
                         0., 0., [1.], [0.], dt, dlnf, ef)
    ok = abs(ef(no,1) - 0.8) < 1e-5 .and. all(ef(1:no-1,1) == 0) .and. all(ef(no+1:,1) == 0)

    ! lengthening: into bin no from the tail
    e = 0
    e(no+1:,1) = 1
    ef = e
    call frequency_shift(nf, 1, no, e, kovcg, invcp, cg, zero, zero, 0., 0., 2 * rate, 0., &
                         0., 0., [1.], [0.], dt, dlnf, ef)
    ok = ok .and. abs(ef(no,1) * kovcg(no) / kovcg(no+1) - 0.2) < 1e-5 &
         .and. all(ef(1:no-1,1) == 0) .and. all(ef(no+1:,1) == e(no+1:,1))

    res = test('frequency_shift_exchanges_action_with_the_diagnostic_tail', ok)
    if (.not. res % ok) write(stderr, '(a,10es11.3)') 'test_advection: spectrum: ', ef(:,1)
  end function frequency_shift_exchanges_action_with_the_diagnostic_tail


  function propagation_shoals_by_linear_theory() result(res)
    ! Steady propagation of 10 s waves at normal incidence over a slope
    ! from 100 m to 8 m. Without dissipation the energy flux is conserved,
    ! and the wave height follows linear shoaling, sqrt(cg0/cg), only if
    ! the action of each frequency bin, e k/cg, is propagated rather than e.
    type(test_result) :: res
    integer, parameter :: n = 200
    real, parameter :: dx = 500
    real :: e(1,1,0:n+1), ef(1,1,0:n+1), k(1,0:n+1), cg(1,0:n+1), kovcg(1,0:n+1)
    real :: d, dt, height, height_linear, worst
    integer :: i, step
    do i = 0, n + 1
      d = 100 - 92 * real(i) / (n + 1)
      k(1,i) = real(wavenumber(0.1_rk, real(d, rk), 1025._rk, 9.81_rk, 0._rk))
      cg(1,i) = real(group_speed(real(k(1,i), rk), real(d, rk), 1025._rk, 9.81_rk, 0._rk))
    end do
    kovcg = k / cg
    dt = 0.5 * dx / maxval(cg)
    e = 0
    e(1,1,0) = 1
    do step = 1, 3000
      ef = e
      do i = 1, n
        call propagation(1, 1, 1, e(:,:,i), e(:,:,i+1), e(:,:,i-1), e(:,:,i), e(:,:,i), &
                         cg(:,i), cg(:,i+1), cg(:,i-1), cg(:,i), cg(:,i), kovcg(:,i), &
                         kovcg(:,i+1), kovcg(:,i-1), kovcg(:,i), kovcg(:,i), [1.], [0.], &
                         dx, dx, dx, dx, .false., 0., 0., 0., 0., 0.25 * dt / dx**2, ef(:,:,i))
      end do
      e(:,:,1:n) = ef(:,:,1:n)
    end do
    ! the variance per unit frequency is e k/cg
    worst = 0
    do i = 1, n - 1
      height = sqrt(e(1,1,i) * kovcg(1,i) / (e(1,1,0) * kovcg(1,0)))
      height_linear = sqrt(cg(1,0) / cg(1,i))
      worst = max(worst, abs(height / height_linear - 1))
    end do
    res = test('propagation_shoals_by_linear_theory', worst < 0.01)
    if (.not. res % ok) write(stderr, '(a,f10.6)') &
      'test_advection: largest relative error of the shoaling wave height: ', worst
  end function propagation_shoals_by_linear_theory


  function action_flux_is_conserved_on_opposing_current() result(res)
    ! Deep-water waves with a narrow spectrum around 8 s propagate against
    ! a current that strengthens from 0 to 1 m/s over 100 km. In steady
    ! state, along each ray the absolute frequency sigma + k U and the
    ! action flux (cg + U) E/sigma are conserved, so the variance rises by
    ! the mean over the spectrum of (sigma/sigma0) cg0 / (cg + U), from the
    ! shortening of the waves; it would rise by only cg0 / (cg0 + U) with
    ! the energy conserved at a fixed frequency.
    type(test_result) :: res
    integer, parameter :: n = 100, nf = 60
    real, parameter :: dx = 1000, dlnf = 0.02, g = 9.81
    real :: e(nf,1,0:n+1), ef(nf,1,0:n+1), sigma(nf), k(nf), cg(nf), kovcg(nf), invcp(nf)
    real :: u(0:n+1), zero(nf), dt, variance0, variance, ratio, expected, s, weight
    integer :: i, o, step
    do o = 1, nf
      sigma(o) = twopi * 0.08 * exp((o - 1) * dlnf)
      k(o) = sigma(o)**2 / g
      cg(o) = 0.5 * g / sigma(o)
    end do
    kovcg = k / cg
    invcp = k / sigma
    zero = 0
    do i = 0, n + 1
      u(i) = -min(max(real(i - 1) / (n - 1), 0.), 1.)
    end do
    e = 0
    do o = 1, nf
      e(o,1,0) = exp(-0.5 * (log(sigma(o) / (twopi / 8)) / 0.05)**2) / kovcg(o)
    end do
    dt = 0.5 * dx / maxval(cg)
    do step = 1, 4000
      ef = e
      do i = 1, n
        call propagation(nf, 1, nf, e(:,:,i), e(:,:,i+1), e(:,:,i-1), e(:,:,i), e(:,:,i), &
                         cg, cg, cg, cg, cg, kovcg, kovcg, kovcg, kovcg, kovcg, [1.], [0.], &
                         dx, dx, dx, dx, .true., u(i) + u(i+1), u(i) + u(i-1), 0., 0., &
                         0.25 * dt / dx**2, ef(:,:,i))
        call frequency_shift(nf, 1, nf, e(:,:,i), kovcg, invcp, cg, zero, zero, u(i), 0., &
                             (u(i+1) - u(i-1)) / (2 * dx), 0., 0., 0., [1.], [0.], dt, dlnf, &
                             ef(:,:,i))
      end do
      e(:,:,1:n) = ef(:,:,1:n)
    end do
    ! variance per bin: e k dk = e (k/cg) sigma dlnf
    i = n - 5
    variance0 = sum(e(:,1,0) * kovcg * sigma)
    variance = sum(e(:,1,i) * kovcg * sigma)
    ! ray theory: sigma from sigma0 = sigma + (sigma**2/g) U
    expected = 0
    weight = 0
    do o = 1, nf
      s = (1 - sqrt(1 + 4 * u(i) * sigma(o) / g)) / (-2 * u(i) / g)
      expected = expected + e(o,1,0) * kovcg(o) * sigma(o) * (s / sigma(o)) * cg(o) &
                 / (0.5 * g / s + u(i))
      weight = weight + e(o,1,0) * kovcg(o) * sigma(o)
    end do
    expected = expected / weight
    ratio = variance / variance0
    res = test('action_flux_is_conserved_on_opposing_current', abs(ratio / expected - 1) < 0.02)
    if (.not. res % ok) write(stderr, '(a,3f10.4)') &
      'test_advection: variance ratio, expected, at fixed frequency: ', ratio, expected, &
      sum(cg * e(:,1,0) * kovcg * sigma) / sum((cg + u(i)) * e(:,1,0) * kovcg * sigma)
  end function action_flux_is_conserved_on_opposing_current


  subroutine directions(cth, sth)
    real, intent(out) :: cth(nd), sth(nd)
    integer :: p
    do p = 1, nd
      cth(p) = cos((p - 1) * dth)
      sth(p) = sin((p - 1) * dth)
    end do
  end subroutine directions

end program test_advection
