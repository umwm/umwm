program test_grid
  use, intrinsic :: iso_fortran_env, only: stdout => output_unit, stderr => error_unit
  use umwm_grid, only: grid_type

  implicit none

  type(grid_type) :: grid
  integer :: grid_size_x, grid_size_y
  real :: dx, dy
  logical :: ok = .true.

  grid_size_x = 80
  grid_size_y = 60
  dx = 1000
  dy = 1000

  grid = grid_type(grid_size_x, grid_size_y, dx, dy)

  if (grid % size_x /= grid_size_x) then
    write(stderr, '(a)') 'grid % size_x is expected.. failed'
    ok = .false.
  end if

  if (grid % size_y /= grid_size_y) then
    write(stderr, '(a)') 'grid % size_y is expected.. failed'
    ok = .false.
  end if

  if (.not. allocated(grid % x)) then
    write(stderr, '(a)') 'grid % x is allocated.. failed'
    ok = .false.
  end if

  if (.not. allocated(grid % y)) then
    write(stderr, '(a)') 'grid % y is allocated.. failed'
    ok = .false.
  end if

  if (.not. all(shape(grid % x) == [grid_size_x, grid_size_y])) then
    write(stderr, '(a)') 'grid % x is expected shape.. failed'
    ok = .false.
  end if

  if (.not. all(shape(grid % y) == [grid_size_x, grid_size_y])) then
    write(stderr, '(a)') 'grid % y is expected shape.. failed'
    ok = .false.
  end if

  ! 1-d case
  grid = grid_type(grid_size_x, 1, dx, dy)

  if (grid % size_x /= grid_size_x) then
    write(stderr, '(a)') 'grid % size_x is expected.. failed'
    ok = .false.
  end if

  if (grid % size_y /= 1) then
    write(stderr, '(a)') 'grid % size_y is expected.. failed'
    ok = .false.
  end if

  grid = grid_type(1, grid_size_y, dx, dy)

  if (grid % size_x /= 1) then
    write(stderr, '(a)') 'grid % size_x is expected.. failed'
    ok = .false.
  end if

  if (grid % size_y /= grid_size_y) then
    write(stderr, '(a)') 'grid % size_y is expected.. failed'
    ok = .false.
  end if

  if (ok) then
    write(stdout, '(a)') 'test_grid: All tests passed.'
  else
    write(stderr, '(a)') 'test_grid: One or more tests failed.'
    error stop 1
  end if

end program test_grid