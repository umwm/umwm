program test_forcing
  use, intrinsic :: iso_fortran_env, only: stdout => output_unit, stderr => error_unit
  use umwm_config, only: config_type
  use umwm_grid, only: grid_type
  use umwm_forcing, only: forcing_type

  implicit none

  type(config_type) :: config
  type(grid_type) :: grid
  type(forcing_type) :: forcing
  integer :: grid_size_x, grid_size_y
  real :: dx, dy
  logical :: ok = .true.

  config = config_type('test/umwm-test.toml')

  grid_size_x = config % grid_size_x
  grid_size_y = config % grid_size_y
  dx = config % dx
  dy = config % dy

  grid = grid_type(grid_size_x, grid_size_y, dx, dy)
  forcing = forcing_type(grid, config)

  if (.not. allocated(forcing % u_atmosphere)) then
    write(stderr, '(a)') 'forcing % u_atmosphere is allocated.. failed'
    ok = .false.
  end if

  if (.not. allocated(forcing % v_atmosphere)) then
    write(stderr, '(a)') 'forcing % v_atmosphere is allocated.. failed'
    ok = .false.
  end if

  if (.not. allocated(forcing % density_atmosphere)) then
    write(stderr, '(a)') 'forcing % density_atmosphere is allocated.. failed'
    ok = .false.
  end if

  if (.not. allocated(forcing % u_ocean)) then
    write(stderr, '(a)') 'forcing % u_ocean is allocated.. failed'
    ok = .false.
  end if

  if (.not. allocated(forcing % v_ocean)) then
    write(stderr, '(a)') 'forcing % v_ocean is allocated.. failed'
    ok = .false.
  end if

  if (.not. allocated(forcing % density_ocean)) then
    write(stderr, '(a)') 'forcing % density_ocean is allocated.. failed'
    ok = .false.
  end if

  if (.not. all(shape(forcing % u_atmosphere) == [grid_size_x, grid_size_y])) then
    write(stderr, '(a)') 'forcing % u_atmosphere is expected shape.. failed'
    ok = .false.
  end if

  if (.not. all(shape(forcing % v_atmosphere) == [grid_size_x, grid_size_y])) then
    write(stderr, '(a)') 'forcing % v_atmosphere is expected shape.. failed'
    ok = .false.
  end if

  if (.not. all(shape(forcing % density_atmosphere) == [grid_size_x, grid_size_y])) then
    write(stderr, '(a)') 'forcing % density_atmosphere is expected shape.. failed'
    ok = .false.
  end if

  if (.not. all(shape(forcing % u_ocean) == [grid_size_x, grid_size_y])) then
    write(stderr, '(a)') 'forcing % u_ocean is expected shape.. failed'
    ok = .false.
  end if

  if (.not. all(shape(forcing % v_ocean) == [grid_size_x, grid_size_y])) then
    write(stderr, '(a)') 'forcing % v_ocean is expected shape.. failed'
    ok = .false.
  end if

  if (.not. all(shape(forcing % density_ocean) == [grid_size_x, grid_size_y])) then
    write(stderr, '(a)') 'forcing % density_ocean is expected shape.. failed'
    ok = .false.
  end if

  if (ok) then
    write(stdout, '(a)') 'test_forcing: All tests passed.'
  else
    write(stderr, '(a)') 'test_forcing: One or more tests failed.'
    error stop 1
  end if

end program test_forcing
