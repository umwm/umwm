program test_advection
  use tuff, only: test, test_result
  use umwm_advection, only: current_rotation_rate, depth_rotation_rate, frequency_shift, &
                            great_circle_rotation_rate, propagation, refraction
  use umwm_config, only: config_type
  use umwm_constants, only: dr, r_earth, rk, stderr, twopi
  use umwm_dispersion, only: group_speed, wavenumber
  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type
#ifdef MPI
  use umwm_env, only: env_init, env_stop
#endif

  implicit none

  integer, parameter :: nf = 3, nd = 36
  real, parameter :: dth = real(twopi) / nd
  ! The direction of the one-dimensional tests, along x. As named arrays,
  ! unlike array constructors, they need no temporary at each call.
  real, parameter :: cth_x(1) = [1.], sth_x(1) = [0.]

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
    test(action_flux_is_conserved_on_opposing_current), &
    test(great_circle_rotation_follows_clairaut), &
    test(uniform_isotropic_swell_stays_so_on_a_sphere), &
    test(swell_follows_great_circles), &
    test(swell_follows_great_circles_on_rotated_grid) &
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
    call refraction(nf, nd, nf - 1, e, zero, zero, zero, 0., 0., 1e-3, 0., 1e-3, -1e-3, &
                    cth, sth, 1000., dth, ef)
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
    call refraction(nf, nd, nf, e, zero, zero, zero, 0., 0., 0., -omega, omega, 0., cth, sth, &
                    dt, dth, ef)
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
                         cth_x, sth_x, dt, dlnf, ef)
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
                         0., 0., cth_x, sth_x, dt, dlnf, ef)
    ok = abs(ef(no,1) - 0.8) < 1e-5 .and. all(ef(1:no-1,1) == 0) .and. all(ef(no+1:,1) == 0)

    ! lengthening: into bin no from the tail
    e = 0
    e(no+1:,1) = 1
    ef = e
    call frequency_shift(nf, 1, no, e, kovcg, invcp, cg, zero, zero, 0., 0., 2 * rate, 0., &
                         0., 0., cth_x, sth_x, dt, dlnf, ef)
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
                         kovcg(:,i+1), kovcg(:,i-1), kovcg(:,i), kovcg(:,i), cth_x, sth_x, &
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
                         cg, cg, cg, cg, cg, kovcg, kovcg, kovcg, kovcg, kovcg, cth_x, sth_x, &
                         dx, dx, dx, dx, .true., u(i) + u(i+1), u(i) + u(i-1), 0., 0., &
                         0.25 * dt / dx**2, ef(:,:,i))
        call frequency_shift(nf, 1, nf, e(:,:,i), kovcg, invcp, cg, zero, zero, u(i), 0., &
                             (u(i+1) - u(i-1)) / (2 * dx), 0., 0., 0., cth_x, sth_x, dt, dlnf, &
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


  function great_circle_rotation_follows_clairaut() result(res)
    ! Integrates a ray along the surface of the sphere from 40N at 30
    ! degrees north of east, turning at the rate of the great circles
    ! relative to the parallels of a lon-lat grid, of geodesic curvature
    ! tan(lat)/R. Clairaut's relation, cos(lat) cos(th) = const, holds along
    ! the ray, which crosses its highest latitude, where it heads east, and
    ! comes back down to 40N heading 30 degrees south of east.
    type(test_result) :: res
    real(rk), parameter :: cg = 10, ds = 1e3_rk
    real(rk) :: lat, th, lat_mid, th_mid, clairaut_start, clairaut_max, lat_max
    integer :: n
    logical :: ok
    lat = 40 * dr
    th = 30 * dr
    clairaut_start = cos(lat) * cos(th)
    clairaut_max = 0
    lat_max = lat
    do n = 1, 20000
      lat_mid = lat + 0.5_rk * ds * sin(th) / r_earth
      th_mid = th + 0.5_rk * ds * ray_turning(cg, lat, th)
      lat = lat + ds * sin(th_mid) / r_earth
      th = th + ds * ray_turning(cg, lat_mid, th_mid)
      clairaut_max = max(clairaut_max, abs(cos(lat) * cos(th) / clairaut_start - 1))
      lat_max = max(lat_max, lat)
      if (lat < 40 * dr) exit
    end do
    ok = clairaut_max < 1e-4_rk .and. abs(th / dr + 30) < 0.1_rk &
         .and. abs(cos(lat_max) / clairaut_start - 1) < 1e-4_rk
    res = test('great_circle_rotation_follows_clairaut', ok)
    if (.not. res % ok) write(stderr, '(a,es10.3,2f10.4)') &
      'test_advection: Clairaut error, final direction, highest latitude [deg]: ', &
      clairaut_max, th / dr, lat_max / dr
  end function great_circle_rotation_follows_clairaut


  real(rk) function ray_turning(cg, lat, th)
    ! The rate of change of the direction th with the distance along a ray
    ! at the latitude lat on the sphere, at the group speed cg, relative to
    ! the parallels, of geodesic curvature tan(lat)/R.
    real(rk), intent(in) :: cg, lat, th
    ray_turning = great_circle_rotation_rate(real(cg), real(tan(lat) / r_earth), 0., &
                                             real(cos(th)), real(sin(th))) / cg
  end function ray_turning


  function uniform_isotropic_swell_stays_so_on_a_sphere() result(res)
    ! On a sphere, a uniform and isotropic spectrum stays so. On a lon-lat
    ! grid, the meridians converge toward the pole, so the propagation on
    ! its own gains the action of the waves toward the pole, and loses that
    ! of the waves toward the equator, at (cg/R) tan(lat) sin(th). The
    ! turning of the great circles relative to the parallels, in the
    ! refraction, takes it back, to within the error of the differences over
    ! direction, (1 - sin(dth)/dth) of it.
    type(test_result) :: res
    real, parameter :: cg = 11.7, dt = 1800
    type(grid_type) :: grid
    real, allocatable :: e(:,:,:), ef(:,:,:), cgv(:,:), kovcg(:,:)
    real :: change_propagation, change_both, expected
    integer :: i, p
    logical :: ok
    call lonlat_grid(grid, -10., 25., 1., 8, 30)
    allocate(e(1,nd,0:grid % imm), ef(1,nd,0:grid % imm), cgv(1,0:grid % imm), &
             kovcg(1,0:grid % imm))
    e = 1
    cgv = cg
    kovcg = 1
    ok = .true.
    change_propagation = 0
    change_both = 0
    do i = grid % istart, grid % iend
      ef(:,:,i) = e(:,:,i)
      call propagate(grid, i, e, cgv, kovcg, dt, ef(:,:,i))
      expected = dt * cg * tan(dr * grid % lat(grid % mi(i), grid % ni(i))) / r_earth
      do p = 1, nd
        ok = ok .and. abs(ef(1,p,i) - 1 - expected * grid % sth_curv(p,i)) < 0.02 * expected
      end do
      change_propagation = max(change_propagation, maxval(abs(ef(:,:,i) - 1)))
      call refraction(1, nd, 1, e(:,:,i), cgv(:,i), [0.], [0.], grid % kappa_x(i), &
                      grid % kappa_y(i), 0., 0., 0., 0., grid % cth_curv(:,i), &
                      grid % sth_curv(:,i), dt, dth, ef(:,:,i))
      change_both = max(change_both, maxval(abs(ef(:,:,i) - 1)))
    end do
    ok = ok .and. change_both < 0.01 * change_propagation
    call grid % finalize()
    res = test('uniform_isotropic_swell_stays_so_on_a_sphere', ok)
    if (.not. res % ok) write(stderr, '(a,2es12.4)') &
      'test_advection: largest change by propagation, and with refraction: ', &
      change_propagation, change_both
  end function uniform_isotropic_swell_stays_so_on_a_sphere


  function swell_follows_great_circles() result(res)
    ! Packets of 15 s swell in deep water, of group speed 11.7 m/s, from
    ! 40N, 40E on a one-degree lon-lat grid, one launched east and one west,
    ! travel 4000 km in 95 hours. Both follow the great circle, which turns
    ! toward the equator, to 31.3N at 83.5E and 3.5W, within the numerical
    ! diffusion of the upstream scheme, and not the parallel, 970 km north.
    type(test_result) :: res
    real, parameter :: cg = 11.7, distance = 4e6, lon0 = 40, lat0 = 40
    type(grid_type) :: grid
    real :: lat_east, lon_east, lat_west, lon_west, lat_gc, lon_gc_east, lon_gc_west
    real :: error_east, error_west
    logical :: ok
    call lonlat_grid(grid, -25., 15., 1., 131, 46)
    call propagate_packet(grid, cg, lon0, lat0, [18, 19], distance, lat_east, lon_east)
    call propagate_packet(grid, cg, lon0, lat0, [1, nd], distance, lat_west, lon_west)
    call grid % finalize()
    call great_circle(lon0, lat0, 0., distance, lon_gc_east, lat_gc)
    call great_circle(lon0, lat0, real(twopi) / 2, distance, lon_gc_west, lat_gc)
    error_east = central_angle(lon_east, lat_east, lon_gc_east, lat_gc) * r_earth
    error_west = central_angle(lon_west, lat_west, lon_gc_west, lat_gc) * r_earth
    ok = error_east < 0.04 * distance .and. error_west < 0.04 * distance &
         .and. abs(lat_east - lat_west) < 0.05 .and. abs(lon_east + lon_west - 2 * lon0) < 0.05
    res = test('swell_follows_great_circles', ok)
    if (.not. res % ok) write(stderr, '(a,8f9.3)') &
      'test_advection: east packet and great circle (lat, lon), west packet and great circle: ', &
      lat_east, lon_east, lat_gc, lon_gc_east, lat_west, lon_west, lat_gc, lon_gc_west
  end function swell_follows_great_circles


  function swell_follows_great_circles_on_rotated_grid() result(res)
    ! The packets of swell_follows_great_circles from 45N, 0E, at the center
    ! of a one-degree rotated-pole grid, whose x axis points east there and
    ! turns by up to 40 degrees across it, travel 3000 km in 71 hours. Both
    ! follow the great circle to 39.1N at 35.8E and 35.8W, within the
    ! numerical diffusion of the upstream scheme, and not the parallel, 660
    ! km north, along which the waves would keep their direction relative to
    ! east as the grid rotates under them.
    type(test_result) :: res
    real, parameter :: cg = 11.7, distance = 3e6, lon0 = 0, lat0 = 45
    type(grid_type) :: grid
    real :: lat_east, lon_east, lat_west, lon_west, lat_gc, lon_gc_east, lon_gc_west
    real :: error_east, error_west
    logical :: ok
    call rotated_pole_grid(grid, 90 - lat0, 71, 41)
    call propagate_packet(grid, cg, lon0, lat0, [18, 19], distance, lat_east, lon_east)
    call propagate_packet(grid, cg, lon0, lat0, [1, nd], distance, lat_west, lon_west)
    call grid % finalize()
    call great_circle(lon0, lat0, 0., distance, lon_gc_east, lat_gc)
    call great_circle(lon0, lat0, real(twopi) / 2, distance, lon_gc_west, lat_gc)
    error_east = central_angle(lon_east, lat_east, lon_gc_east, lat_gc) * r_earth
    error_west = central_angle(lon_west, lat_west, lon_gc_west, lat_gc) * r_earth
    ok = error_east < 0.05 * distance .and. error_west < 0.05 * distance
    res = test('swell_follows_great_circles_on_rotated_grid', ok)
    if (.not. res % ok) write(stderr, '(a,8f9.3)') &
      'test_advection: east packet and great circle (lat, lon), west packet and great circle: ', &
      lat_east, lon_east, lat_gc, lon_gc_east, lat_west, lon_west, lat_gc, lon_gc_west
  end function swell_follows_great_circles_on_rotated_grid


  subroutine rotated_pole_grid(grid, tilt, mm, nm)
    ! A one-degree grid of mm x nm cells, regular in the longitude and
    ! latitude about a pole tilted by tilt [deg] toward 180E, centered at
    ! the rotated equator and prime meridian, at (0E, 90 - tilt N), and
    ! read from file, with the directions of the spectrum projected onto it.
    use netcdf
    type(grid_type), intent(out) :: grid
    real, intent(in) :: tilt
    integer, intent(in) :: mm, nm
    type(config_type) :: config
    real :: lon(mm,nm), lat(mm,nm), x, y, z
    integer :: ncid, xdim, ydim, lonid, latid, m, n
    do n = 1, nm
      do m = 1, mm
        x = cos(dr * (n - (nm + 1) / 2)) * cos(dr * (m - (mm + 1) / 2))
        y = cos(dr * (n - (nm + 1) / 2)) * sin(dr * (m - (mm + 1) / 2))
        z = sin(dr * (n - (nm + 1) / 2))
        lon(m,n) = atan2(y, x * cos(dr * tilt) - z * sin(dr * tilt)) / dr
        lat(m,n) = asin(x * sin(dr * tilt) + z * cos(dr * tilt)) / dr
      end do
    end do
    call execute_command_line('mkdir -p output/test_advection')
    call nc_check(nf90_create('output/test_advection/umwm.grid', nf90_clobber, ncid))
    call nc_check(nf90_def_dim(ncid, 'x', mm, xdim))
    call nc_check(nf90_def_dim(ncid, 'y', nm, ydim))
    call nc_check(nf90_def_var(ncid, 'lon', nf90_float, [xdim, ydim], lonid))
    call nc_check(nf90_def_var(ncid, 'lat', nf90_float, [xdim, ydim], latid))
    call nc_check(nf90_enddef(ncid))
    call nc_check(nf90_put_var(ncid, lonid, lon))
    call nc_check(nf90_put_var(ncid, latid, lat))
    call nc_check(nf90_close(ncid))
    config = config_type('namelists/default_main.nml')
    config % input_dir = 'output/test_advection'
    config % gridfromfile = .true.
    config % mm = mm
    config % nm = nm
    config % pm = nd
    call grid % initialize(config)
    call grid % initialize_direction_projection(spectrum_type(nf, nd, 0.05_rk, 0.1_rk))
  end subroutine rotated_pole_grid


  subroutine nc_check(stat)
    use netcdf, only: nf90_noerr, nf90_strerror
    integer, intent(in) :: stat
    if (stat /= nf90_noerr) then
      write(stderr, '(a)') 'test_advection: netCDF error: ' // trim(nf90_strerror(stat))
      error stop 1
    end if
  end subroutine nc_check


  subroutine lonlat_grid(grid, lon0, lat0, spacing, mm, nm)
    ! A regional lon-lat grid of mm x nm cells from (lon0, lat0) by spacing
    ! [deg], with the directions of the spectrum projected onto it.
    type(grid_type), intent(out) :: grid
    real, intent(in) :: lon0, lat0, spacing
    integer, intent(in) :: mm, nm
    type(config_type) :: config
    config = config_type('namelists/default_main.nml')
    config % mm = mm
    config % nm = nm
    config % pm = nd
    config % lon0 = lon0
    config % lat0 = lat0
    config % dlon = spacing
    config % dlat = spacing
    call grid % initialize(config)
    call grid % initialize_direction_projection(spectrum_type(nf, nd, 0.05_rk, 0.1_rk))
  end subroutine lonlat_grid


  subroutine propagate(grid, i, e, cg, kovcg, dt, ef)
    ! Integrates the propagation of the spectrum e of one frequency at the
    ! cell i of the grid over dt into ef.
    type(grid_type), intent(in) :: grid
    integer, intent(in) :: i
    real, intent(in) :: e(:,:,0:), cg(:,0:), kovcg(:,0:), dt
    real, intent(inout) :: ef(1,nd)
    associate(ie => grid % ie(i), iw => grid % iw(i), in => grid % in(i), is => grid % is(i))
    call propagation(1, nd, 1, e(:,:,i), e(:,:,ie), e(:,:,iw), e(:,:,in), e(:,:,is), &
                     cg(:,i), cg(:,ie), cg(:,iw), cg(:,in), cg(:,is), kovcg(:,i), kovcg(:,ie), &
                     kovcg(:,iw), kovcg(:,in), kovcg(:,is), grid % cth_curv(:,i), &
                     grid % sth_curv(:,i), grid % dxn(i), grid % dxs(i), grid % dye(i), &
                     grid % dyw(i), .false., 0., 0., 0., 0., 0.25 * dt * grid % oneovar(i), ef)
    end associate
  end subroutine propagate


  subroutine propagate_packet(grid, cg, lon0, lat0, bins, distance, lat, lon)
    ! Propagates and refracts a Gaussian packet of e-folding radius 2.5
    ! degrees at (lon0, lat0) [deg], in the direction bins, at the uniform
    ! group speed cg over distance, and returns the center of its energy on
    ! the sphere (lat, lon) [deg].
    type(grid_type), intent(in) :: grid
    real, intent(in) :: cg, lon0, lat0, distance
    integer, intent(in) :: bins(:)
    real, intent(out) :: lat, lon
    real, allocatable :: e(:,:,:), ef(:,:,:), cgv(:,:), kovcg(:,:)
    real :: dt, w, x, y, z, clat
    integer :: i, n, num_steps
    allocate(e(1,nd,0:grid % imm), ef(1,nd,0:grid % imm), cgv(1,0:grid % imm), &
             kovcg(1,0:grid % imm))
    cgv = cg
    kovcg = 1
    e = 0
    do i = 1, grid % im
      e(1,bins,i) = exp(-(central_angle(lon0, lat0, grid % lon(grid % mi(i), grid % ni(i)), &
                                         grid % lat(grid % mi(i), grid % ni(i))) &
                          / (2.5 * dr))**2)
    end do
    ! at 0.9 of the largest stable time step, cg dt (1/dx + 1/dy) = 1
    num_steps = ceiling(distance * (1 / minval(grid % dx(1:grid % im)) &
                                    + 1 / minval(grid % dy(1:grid % im))) / 0.9)
    dt = distance / cg / num_steps
    do n = 1, num_steps
      ef = e
      do i = grid % istart, grid % iend
        call propagate(grid, i, e, cgv, kovcg, dt, ef(:,:,i))
        call refraction(1, nd, 1, e(:,:,i), cgv(:,i), [0.], [0.], grid % kappa_x(i), &
                        grid % kappa_y(i), 0., 0., 0., 0., grid % cth_curv(:,i), &
                        grid % sth_curv(:,i), dt, dth, ef(:,:,i))
      end do
      e(:,:,1:grid % im) = ef(:,:,1:grid % im)
    end do
    x = 0
    y = 0
    z = 0
    do i = 1, grid % im
      w = sum(e(1,:,i)) * grid % ar(i)
      clat = cos(dr * grid % lat(grid % mi(i), grid % ni(i)))
      x = x + w * clat * cos(dr * grid % lon(grid % mi(i), grid % ni(i)))
      y = y + w * clat * sin(dr * grid % lon(grid % mi(i), grid % ni(i)))
      z = z + w * sin(dr * grid % lat(grid % mi(i), grid % ni(i)))
    end do
    lat = atan2(z, sqrt(x**2 + y**2)) / dr
    lon = atan2(y, x) / dr
  end subroutine propagate_packet


  real function central_angle(lon1, lat1, lon2, lat2) result(res)
    ! The angle [rad] between two points (lon, lat) [deg] on the sphere.
    real, intent(in) :: lon1, lat1, lon2, lat2
    res = 2 * asin(sqrt(sin(0.5 * dr * (lat2 - lat1))**2 &
                        + cos(dr * lat1) * cos(dr * lat2) * sin(0.5 * dr * (lon2 - lon1))**2))
  end function central_angle


  subroutine great_circle(lon1, lat1, th, distance, lon2, lat2)
    ! The point (lon2, lat2) [deg] at distance [m] from (lon1, lat1) along
    ! the great circle that leaves it in the direction th [rad], anti-clockwise
    ! from east.
    real, intent(in) :: lon1, lat1, th, distance
    real, intent(out) :: lon2, lat2
    real :: azimuth, angle
    azimuth = real(twopi) / 4 - th
    angle = distance / r_earth
    lat2 = asin(sin(dr * lat1) * cos(angle) + cos(dr * lat1) * sin(angle) * cos(azimuth))
    lon2 = lon1 + atan2(sin(azimuth) * sin(angle) * cos(dr * lat1), &
                        cos(angle) - sin(dr * lat1) * sin(lat2)) / dr
    lat2 = lat2 / dr
  end subroutine great_circle


  subroutine directions(cth, sth)
    real, intent(out) :: cth(nd), sth(nd)
    integer :: p
    do p = 1, nd
      cth(p) = cos((p - 1) * dth)
      sth(p) = sin((p - 1) * dth)
    end do
  end subroutine directions

end program test_advection
