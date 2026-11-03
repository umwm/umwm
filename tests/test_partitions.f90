program test_partitions
  ! Tests the parametric spectra of wave partitions: the discrete spectrum
  ! of a partition has its significant wave height exactly, and its mean
  ! period, mean direction and directional spread, in deep and shallow
  ! water; partitions add up, absent ones add nothing; and the partitions
  ! are read from an input file.
  use netcdf
  use tuff, only: test, test_result
  use umwm_constants, only: pi, rk, twopi
  use umwm_dispersion, only: group_speed, wavenumber
#ifdef MPI
  use umwm_env, only: env_init, env_stop
#endif
  use umwm_partitions, only: partitions_type, partition_spectrum, period_mean, period_peak

  implicit none

  integer, parameter :: om = 37, pm = 36
  real, parameter :: fmin = 0.0313, fmax = 2.0
  real :: f(om), th(pm), dth, dlnf
  type(test_result) :: suite
  integer :: o, p

#ifdef MPI
  call env_init()
#endif

  dlnf = log(fmax / fmin) / (om - 1)
  f = [(exp(log(fmin) + (o - 1) * dlnf), o = 1, om)]
  dth = twopi / pm
  th = [((p - 0.5 * (pm + 1)) * dth, p = 1, pm)]

  suite = test('test_partitions', [ &
    test(deep_water_partition), &
    test(shallow_water_partition), &
    test(narrow_and_broad_spreading), &
    test(peak_period_partition), &
    test(partitions_add_up), &
    test(absent_partitions), &
    test(period_beyond_frequencies), &
    test(read_from_file) &
  ])

#ifdef MPI
  call env_stop()
#endif
  if (.not. suite % ok) error stop 1

contains

  function deep_water_partition() result(res)
    type(test_result) :: res
    real :: e(om,pm), kdk(om)
    kdk = bin_weights(1000.)
    e = 0
    call partition_spectrum(om, pm, f, th, dth, kdk, 2.5, 10., period_mean, 0.7, 0.4, 3.3, e)
    res = test('deep_water_partition', matches(e, kdk, 2.5, 10., 0.7, 0.4))
  end function deep_water_partition


  function shallow_water_partition() result(res)
    ! The frequency spectrum, and so the mean period, does not depend on the
    ! depth, while the spectrum e(k) and the weights k dk do.
    type(test_result) :: res
    real :: e(om,pm), kdk(om)
    kdk = bin_weights(8.)
    e = 0
    call partition_spectrum(om, pm, f, th, dth, kdk, 1.2, 7., period_mean, -2.9, 0.6, 1., e)
    res = test('shallow_water_partition', matches(e, kdk, 1.2, 7., -2.9, 0.6))
  end function shallow_water_partition


  function narrow_and_broad_spreading() result(res)
    ! A narrow swell resolved by a few directions, and a nearly isotropic sea.
    type(test_result) :: res
    real :: narrow(om,pm), broad(om,pm), kdk(om)
    kdk = bin_weights(4000.)
    narrow = 0
    broad = 0
    call partition_spectrum(om, pm, f, th, dth, kdk, 3., 14., period_mean, 1.0, 0.15, 3.3, narrow)
    call partition_spectrum(om, pm, f, th, dth, kdk, 1., 5., period_mean, 0.3, 1.3, 3.3, broad)
    res = test('narrow_and_broad_spreading', &
      matches(narrow, kdk, 3., 14., 1.0, 0.15, spread_tolerance=0.03) &
      .and. matches(broad, kdk, 1., 5., 0.3, 1.3, direction_tolerance=0.05))
  end function narrow_and_broad_spreading


  function peak_period_partition() result(res)
    ! With a peak period, the variance peaks in the bin of 1/Tp or the next,
    ! where f S(f) of the JONSWAP spectrum is largest on bins 12% apart.
    type(test_result) :: res
    real :: e(om,pm), kdk(om), variance(om)
    integer :: peak
    kdk = bin_weights(1000.)
    e = 0
    call partition_spectrum(om, pm, f, th, dth, kdk, 2., 12., period_peak, 0., 0.3, 3.3, e)
    variance = sum(e, dim=2) * kdk * dth
    peak = maxloc(variance, dim=1)
    res = test('peak_period_partition', &
      abs(log(f(peak) * 12.)) <= dlnf .and. abs(hs(e, kdk) - 2.) < 1e-4 * 2.)
  end function peak_period_partition


  function partitions_add_up() result(res)
    ! The variance of several partitions is the sum of theirs.
    type(test_result) :: res
    real :: e(om,pm), kdk(om)
    kdk = bin_weights(50.)
    e = 0
    call partition_spectrum(om, pm, f, th, dth, kdk, 1.5, 4., period_mean, 2.0, 0.5, 3.3, e)
    call partition_spectrum(om, pm, f, th, dth, kdk, 2.0, 12., period_mean, -1.0, 0.3, 3.3, e)
    call partition_spectrum(om, pm, f, th, dth, kdk, 0.5, 16., period_mean, -2.0, 0.2, 3.3, e)
    res = test('partitions_add_up', abs(hs(e, kdk) - sqrt(1.5**2 + 2.0**2 + 0.5**2)) < 1e-4 * 2.5)
  end function partitions_add_up


  function absent_partitions() result(res)
    ! Partitions without energy, with missing data, or without a period add
    ! nothing; a spread or gamma that is not positive takes its default.
    type(test_result) :: res
    real :: e(om,pm), defaults(om,pm), kdk(om)
    kdk = bin_weights(1000.)
    e = 0
    call partition_spectrum(om, pm, f, th, dth, kdk, 0., 10., period_mean, 0., 0.4, 3.3, e)
    call partition_spectrum(om, pm, f, th, dth, kdk, -1., 10., period_mean, 0., 0.4, 3.3, e)
    call partition_spectrum(om, pm, f, th, dth, kdk, 1., 0., period_mean, 0., 0.4, 3.3, e)
    defaults = 0
    call partition_spectrum(om, pm, f, th, dth, kdk, 1., 8., period_mean, 0., 0., 0., defaults)
    res = test('absent_partitions', all(e == 0) &
      .and. matches(defaults, kdk, 1., 8., 0., 0.5))
  end function absent_partitions


  function period_beyond_frequencies() result(res)
    ! A mean period longer than the frequencies of the model can resolve
    ! keeps the energy of the partition at the lowest frequencies.
    type(test_result) :: res
    real :: e(om,pm), kdk(om), variance(om)
    kdk = bin_weights(1000.)
    e = 0
    call partition_spectrum(om, pm, f, th, dth, kdk, 1., 60., period_mean, 0., 0.3, 3.3, e)
    variance = sum(e, dim=2) * kdk * dth
    res = test('period_beyond_frequencies', abs(hs(e, kdk) - 1.) < 1e-4 &
      .and. maxloc(variance, dim=1) <= 2)
  end function period_beyond_frequencies


  function read_from_file() result(res)
    ! Reads two partitions, one with a mean period and default spreading,
    ! one with a peak period and spread, and missing data at one cell.
    type(test_result) :: res
    type(partitions_type) :: partitions
    integer, parameter :: mm = 4, nm = 3
    character(*), parameter :: path = 'output/test_partitions.nc'
    real :: e(om,pm), kdk(om)
    integer :: ncid, dims(3)
    logical :: ok

    call check(nf90_create(path, nf90_clobber, ncid))
    call check(nf90_def_dim(ncid, 'x', mm, dims(1)))
    call check(nf90_def_dim(ncid, 'y', nm, dims(2)))
    call check(nf90_def_dim(ncid, 'time', nf90_unlimited, dims(3)))
    call check(nf90_enddef(ncid))
    call put(ncid, dims, 'phs0', 1.0)
    call put(ncid, dims, 'ptm0', 5.0)
    call put(ncid, dims, 'pdir0', 0.5)
    call put(ncid, dims, 'phs2', 2.0, missing_at=[2, 2])
    call put(ncid, dims, 'ptp2', 13.0)
    call put(ncid, dims, 'pdir2', -1.0)
    call put(ncid, dims, 'pspr2', 0.25)
    call check(nf90_close(ncid))

    call check(nf90_open(path, nf90_nowrite, ncid))
    call partitions % read(ncid, path, mm, nm)
    call check(nf90_close(ncid))

    ok = partitions % num == 2 .and. all(partitions % id(1:2) == [0, 2]) &
      .and. all(partitions % period_kind(1:2) == [period_mean, period_peak]) &
      .and. all(partitions % spread(:,:,1) == 0.5) .and. all(partitions % spread(:,:,2) == 0.25) &
      .and. all(partitions % gamma == 3.3) &
      .and. partitions % hs(2,2,2) < 0 .and. partitions % valid(2,2) .and. partitions % valid(1,1)

    kdk = bin_weights(1000.)
    call partitions % spectrum(3, 2, om, pm, f, th, dth, kdk, e)
    ok = ok .and. abs(hs(e, kdk) - sqrt(5.)) < 1e-4 * sqrt(5.)
    call partitions % spectrum(2, 2, om, pm, f, th, dth, kdk, e)
    ok = ok .and. abs(hs(e, kdk) - 1.) < 1e-4
    call partitions % finalize()

    res = test('read_from_file', ok)
  end function read_from_file


  subroutine put(ncid, dims, name, value, missing_at)
    ! Writes a 4 by 3 field of a value, -1 at a missing cell, to the open
    ! file ncid.
    integer, intent(in) :: ncid, dims(3)
    character(*), intent(in) :: name
    real, intent(in) :: value
    integer, intent(in), optional :: missing_at(2)
    real :: values(4,3,1)
    integer :: varid
    values = value
    if (present(missing_at)) values(missing_at(1),missing_at(2),1) = -1
    call check(nf90_redef(ncid))
    call check(nf90_def_var(ncid, name, nf90_float, dims, varid))
    call check(nf90_enddef(ncid))
    call check(nf90_put_var(ncid, varid, values, start=[1, 1, 1], count=[4, 3, 1]))
  end subroutine put


  subroutine check(stat)
    integer, intent(in) :: stat
    if (stat /= nf90_noerr) then
      write(*, '(a)') 'test_partitions: ' // trim(nf90_strerror(stat))
      error stop 1
    end if
  end subroutine check


  function bin_weights(depth) result(kdk)
    ! k dk of the frequency bins at the depth, as the model computes them.
    real, intent(in) :: depth
    real :: kdk(om), k(om), cg(om)
    k = real(wavenumber(real(f, rk), real(depth, rk), 1025._rk, 9.80665_rk, 0.07_rk))
    cg = real(group_speed(real(k, rk), real(depth, rk), 1025._rk, 9.80665_rk, 0.07_rk))
    kdk = k * twopi * dlnf * f / cg
  end function bin_weights


  real function hs(e, kdk)
    real, intent(in) :: e(om,pm), kdk(om)
    hs = 4 * sqrt(sum(sum(e, dim=2) * kdk) * dth)
  end function hs


  logical function matches(e, kdk, hs_expected, tm10, direction, spread, &
                           direction_tolerance, spread_tolerance)
    ! Whether the spectrum e has the significant wave height (relative 1e-4),
    ! mean period Tm-1,0 (relative 1e-3), mean direction, and directional
    ! spread sqrt(2 (1 - M1)) (relative 1e-2 by default) given.
    real, intent(in) :: e(om,pm), kdk(om), hs_expected, tm10, direction, spread
    real, intent(in), optional :: direction_tolerance, spread_tolerance
    real :: variance(om,pm), m0, m_1, c, s, m1, mean_direction, spread_found, tolerance(2)
    integer :: q
    tolerance = [1e-3, 1e-2]
    if (present(direction_tolerance)) tolerance(1) = direction_tolerance
    if (present(spread_tolerance)) tolerance(2) = spread_tolerance
    do q = 1, pm
      variance(:,q) = e(:,q) * kdk * dth
    end do
    m0 = sum(variance)
    m_1 = sum(sum(variance, dim=2) / f)
    c = sum(sum(variance, dim=1) * cos(th))
    s = sum(sum(variance, dim=1) * sin(th))
    m1 = sqrt(c**2 + s**2) / m0
    mean_direction = atan2(s, c)
    spread_found = sqrt(2 * (1 - m1))
    matches = abs(4 * sqrt(m0) - hs_expected) <= 1e-4 * hs_expected &
      .and. abs(m_1 / m0 - tm10) <= 1e-3 * tm10 &
      .and. abs(modulo(mean_direction - direction + pi, twopi) - pi) <= tolerance(1) &
      .and. abs(spread_found - spread) <= tolerance(2) * spread
    if (.not. matches) write(*, '(a,4(1x,f9.5),a,4(1x,f9.5))') 'test_partitions: found', &
      4 * sqrt(m0), m_1 / m0, mean_direction, spread_found, ', expected', &
      hs_expected, tm10, direction, spread
  end function matches

end program test_partitions
