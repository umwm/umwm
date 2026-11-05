program test_grid
  use tuff, only: test, test_result, nearly_equal
  use umwm_config, only: config_type
#ifdef MPI
  use umwm_env, only: env_init, env_stop
#endif
  use umwm_grid, only: grid_type

  implicit none

  type(test_result) :: suite

#ifdef MPI
  call env_init()
#endif
  suite = test('test_grid', [ &
    test(constant_limited_grid), &
    test(remap_round_trip), &
    test(limited_neighbor_aliases), &
    test(global_periodic_neighbors), &
    test(lonlat_grid_across_dateline), &
    test(global_lonlat_grid), &
    test(input_files_for_another_grid) &
  ])

#ifdef MPI
  call env_stop()
#endif
  if (.not. suite % ok) error stop 1

contains

  function constant_limited_grid() result(res)
    type(test_result) :: res
    type(config_type) :: config
    type(grid_type) :: grid
    logical :: ok

    config = valid_config()
    call grid % initialize(config)

    ok = grid % mm == config % mm .and. &
         grid % nm == config % nm .and. &
         grid % imm == config % mm * config % nm .and. &
         grid % im == (config % mm - 2) * (config % nm - 2) .and. &
         grid % istart >= 1 .and. grid % iend <= grid % im .and. &
         grid % istart <= grid % iend .and. &
         grid % iistart <= grid % istart .and. grid % iiend >= grid % iend .and. &
         all(grid % mask(2:config % mm-1,2:config % nm-1) == 1) .and. &
         all(grid % mask(:,1) == 0) .and. all(grid % mask(:,config % nm) == 0) .and. &
         all(grid % mask(1,:) == 0) .and. all(grid % mask(config % mm,:) == 0) .and. &
         all(nearly_equal(grid % dx_2d, config % delx)) .and. &
         all(nearly_equal(grid % dy_2d, config % dely)) .and. &
         all(nearly_equal(grid % ar, config % delx * config % dely)) .and. &
         all(nearly_equal(grid % oneovdx, 1.0 / config % delx)) .and. &
         all(nearly_equal(grid % oneovdy, 1.0 / config % dely)) .and. &
         all(nearly_equal(grid % oneovar, 1.0 / (config % delx * config % dely))) .and. &
         all(nearly_equal(grid % d_2d(2:config % mm-1,2:config % nm-1), config % dpt)) .and. &
         all(nearly_equal(grid % d_2d(:,1), config % dmin)) .and. &
         grid % istart == 1 .and. grid % iend == grid % im .and. &
         grid % iistart == grid % istart .and. grid % iiend == grid % iend

    call grid % finalize()

    res = test('constant_limited_grid', ok)
  end function constant_limited_grid


  function remap_round_trip() result(res)
    type(test_result) :: res
    type(config_type) :: config
    type(grid_type) :: grid
    real, allocatable :: field_mn(:,:), field_back(:,:), field_i(:)
    integer :: m, n
    logical :: ok

    config = valid_config()
    call grid % initialize(config)

    allocate(field_mn(grid % mm,grid % nm))
    do n = 1, grid % nm
      do m = 1, grid % mm
        field_mn(m,n) = real(m + 100 * n)
      end do
    end do

    field_i = grid % remap_mn2i(field_mn)
    field_back = grid % remap_i2mn(field_i)

    ok = size(field_i) == grid % imm .and. all(nearly_equal(field_back, field_mn))

    call grid % finalize()

    res = test('remap_round_trip', ok)
  end function remap_round_trip


  function limited_neighbor_aliases() result(res)
    type(test_result) :: res
    type(config_type) :: config
    type(grid_type) :: grid
    integer :: i
    logical :: ok

    config = valid_config()
    call grid % initialize(config)

    i = grid % ii(2,2)
    if (i >= grid % istart .and. i <= grid % iend) then
      ok = grid % iw(i) == grid % iistart - 1 .and. &
           grid % is(i) == grid % iistart - 1 .and. &
           grid % iiw(i) == grid % ii(1,2) .and. &
           grid % iis(i) == grid % ii(2,1) .and. &
           grid % ie(i) == grid % ii(3,2) .and. &
           grid % in(i) == grid % ii(2,3) .and. &
           nearly_equal(grid % dxs(i), config % delx) .and. &
           nearly_equal(grid % dxn(i), config % delx) .and. &
           nearly_equal(grid % dyw(i), config % dely) .and. &
           nearly_equal(grid % dye(i), config % dely)
    else
      ok = .false.
    end if

    call grid % finalize()

    res = test('limited_neighbor_aliases', ok)
  end function limited_neighbor_aliases


  function global_periodic_neighbors() result(res)
    type(test_result) :: res
    type(config_type) :: config
    type(grid_type) :: grid
    integer :: west_edge, east_edge
    logical :: ok

    config = valid_config()
    config % isglobal = .true.
    call grid % initialize(config)

    west_edge = grid % ii(1,2)
    east_edge = grid % ii(grid % mm,2)

    ok = grid % im == config % mm * (config % nm - 2) .and. &
         all(grid % mask(1,2:config % nm-1) == 1) .and. &
         all(grid % mask(config % mm,2:config % nm-1) == 1)
#ifdef MPI
    ! MPI remapping replaces off-tile neighbor indices with halo locations;
    ! the mask checks above verify that both periodic edge columns are active.
#else
    ok = ok .and. &
         grid % iw(west_edge) == east_edge .and. &
         grid % ie(east_edge) == west_edge
#endif

    call grid % finalize()

    res = test('global_periodic_neighbors', ok)
  end function global_periodic_neighbors


  function lonlat_grid_across_dateline() result(res)
    ! A regular lon-lat grid made by the model from lon0, lat0, dlon and
    ! dlat, across the dateline: longitudes in [-180, 180), and the spacing
    ! of a sphere, with the grid rotation of great circles.
    use umwm_constants, only: dr, r_earth
    type(test_result) :: res
    type(config_type) :: config
    type(grid_type) :: grid
    integer :: n
    logical :: ok

    config = valid_config()
    config % lon0 = 170
    config % lat0 = 10
    config % dlon = 0.5
    config % dlat = 0.25
    ok = config % lonlat_grid() .and. config % validate(stop_on_error=.false., rank=-1)
    call grid % initialize(config)

    ok = ok .and. nearly_equal(grid % lon(1,1), 170.) .and. nearly_equal(grid % lon(20,1), 179.5) &
         .and. nearly_equal(grid % lon(21,1), -180.) .and. nearly_equal(grid % lat(1,11), 12.5) &
         .and. all(nearly_equal(grid % lon(:,1), grid % lon(:,11)))
    do n = 2, config % nm - 1
      ok = ok .and. all(abs(grid % dx_2d(2:20,n) / (r_earth * dr * 0.5 * cos(dr * grid % lat(2,n))) - 1) < 1e-3) &
              .and. all(abs(grid % dy_2d(2:20,n) / (r_earth * dr * 0.25) - 1) < 1e-3) &
              .and. all(grid % curv(2:20,n) < 0 .and. grid % curv(2:20,n) > -0.01)
    end do

    call grid % finalize()

    res = test('lonlat_grid_across_dateline', ok)
  end function lonlat_grid_across_dateline


  function global_lonlat_grid() result(res)
    ! A global lon-lat grid is periodic in longitude, also in its spacing.
    type(test_result) :: res
    type(config_type) :: config
    type(grid_type) :: grid
    logical :: ok

    config = valid_config()
    config % isglobal = .true.
    config % mm = 72
    config % nm = 31
    config % xpl = 10
    config % lon0 = 0
    config % lat0 = -75
    config % dlon = 5
    config % dlat = 5
    ok = config % validate(stop_on_error=.false., rank=-1)
    call grid % initialize(config)

    ok = ok .and. grid % im == config % mm * (config % nm - 2) &
         .and. all(abs(grid % dx_2d(1,:) / grid % dx_2d(2,:) - 1) < 1e-4) &
         .and. all(abs(grid % dx_2d(72,:) / grid % dx_2d(2,:) - 1) < 1e-4) &
         .and. all(abs(grid % curv(1,:) - grid % curv(2,:)) < 1e-5) &
         .and. nearly_equal(grid % lon(37,1), -180.)

    call grid % finalize()

    ! the longitudes of a global grid must close the circle
    config % dlon = 4
    ok = ok .and. .not. config % validate(stop_on_error=.false., rank=-1)

    res = test('global_lonlat_grid', ok)
  end function global_lonlat_grid


  function input_files_for_another_grid() result(res)
    ! An input file is for another grid if its fields are of another size,
    ! or, if it has lon and lat, its first and last cells are not those of
    ! the lon-lat grid of the configuration, modulo 360 in longitude.
    use umwm_grid, only: input_grid_problem
    type(test_result) :: res
    type(config_type) :: config
    character(*), parameter :: with_coordinates = 'output/test_grid_input.nc'
    character(*), parameter :: without_coordinates = 'output/test_grid_input_nocoords.nc'
    logical :: ok

    call write_input(with_coordinates, .true.)
    call write_input(without_coordinates, .false.)
    config = valid_config()
    config % gridfromfile = .false.
    config % mm = 4
    config % nm = 3
    config % lon0 = -98.5
    config % lat0 = 17.5
    config % dlon = 0.25
    config % dlat = 0.25
    ok = len(input_grid_problem(config, with_coordinates)) == 0 &
         .and. len(input_grid_problem(config, without_coordinates)) == 0

    ! the same grid, 360 degrees of longitude apart
    config % lon0 = 261.5
    ok = ok .and. len(input_grid_problem(config, with_coordinates)) == 0

    ! another region of the same size, which only the coordinates tell
    config % lon0 = -90
    ok = ok .and. index(input_grid_problem(config, with_coordinates), &
                        'its first and last cells are at (-98.500, 17.500) and (-97.750, 18.000)') > 0 &
            .and. len(input_grid_problem(config, without_coordinates)) == 0

    ! another size
    config % lon0 = -98.5
    config % mm = 5
    ok = ok .and. index(input_grid_problem(config, with_coordinates), &
                        'lon is 4 x 3, but the grid is 5 x 3 (mm x nm)') > 0 &
            .and. index(input_grid_problem(config, without_coordinates), &
                        'uw is 4 x 3, but the grid is 5 x 3 (mm x nm)') > 0

    ! a Cartesian grid has no coordinates to compare
    config % mm = 4
    config % lon0 = -90
    config % dlon = 0
    config % dlat = 0
    ok = ok .and. len(input_grid_problem(config, with_coordinates)) == 0

    res = test('input_files_for_another_grid', ok)
  end function input_files_for_another_grid


  subroutine write_input(path, coordinates)
    ! An input file of the wind on a 4 x 3 lon-lat grid from (-98.5, 17.5)
    ! by 0.25 degrees, with its lon and lat if coordinates.
    use netcdf
    character(*), intent(in) :: path
    logical, intent(in) :: coordinates
    integer :: ncid, xdim, ydim, tdim, lonid, latid, uwid, m, n
    real :: lon(4,3), lat(4,3)

    do n = 1, 3
      do m = 1, 4
        lon(m,n) = -98.5 + (m - 1) * 0.25
        lat(m,n) = 17.5 + (n - 1) * 0.25
      end do
    end do
    call nc_check(nf90_create(path, nf90_clobber, ncid))
    call nc_check(nf90_def_dim(ncid, 'x', 4, xdim))
    call nc_check(nf90_def_dim(ncid, 'y', 3, ydim))
    call nc_check(nf90_def_dim(ncid, 'time', nf90_unlimited, tdim))
    if (coordinates) then
      call nc_check(nf90_def_var(ncid, 'lon', nf90_float, [xdim, ydim], lonid))
      call nc_check(nf90_def_var(ncid, 'lat', nf90_float, [xdim, ydim], latid))
    end if
    call nc_check(nf90_def_var(ncid, 'uw', nf90_float, [xdim, ydim, tdim], uwid))
    call nc_check(nf90_enddef(ncid))
    if (coordinates) then
      call nc_check(nf90_put_var(ncid, lonid, lon))
      call nc_check(nf90_put_var(ncid, latid, lat))
    end if
    call nc_check(nf90_put_var(ncid, uwid, reshape(0 * lon, [4, 3, 1])))
    call nc_check(nf90_close(ncid))
  end subroutine write_input


  subroutine nc_check(stat)
    use netcdf, only: nf90_noerr, nf90_strerror
    integer, intent(in) :: stat
    if (stat /= nf90_noerr) then
      write(*, '(a)') 'test_grid: netCDF error: ' // trim(nf90_strerror(stat))
      error stop 1
    end if
  end subroutine nc_check


  function valid_config() result(config)
    type(config_type) :: config

    config = config_type('../namelists/main.nml')
  end function valid_config

end program test_grid
