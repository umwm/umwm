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
    test(global_periodic_neighbors) &
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


  function valid_config() result(config)
    type(config_type) :: config

    config = config_type('../namelists/main.nml')
  end function valid_config

end program test_grid
