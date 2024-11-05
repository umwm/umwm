program test_state
  use, intrinsic :: iso_fortran_env, only: stdout => output_unit, stderr => error_unit
  use umwm_state, only: state_type
  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type

  implicit none

  type(grid_type) :: grid
  type(spectrum_type) :: spectrum
  type(state_type) :: state

  integer :: num_frequencies, num_directions
  real :: frequency_min, frequency_max
  integer :: grid_size_x, grid_size_y
  logical :: ok = .true.

  grid_size_x = 80
  grid_size_y = 60
  num_frequencies = 50
  num_directions = 36
  frequency_min = 0.04
  frequency_max = 2

  grid = grid_type(grid_size_x, grid_size_y)

  spectrum = spectrum_type( &
    num_frequencies, num_directions, frequency_min, frequency_max &
  )

  state = state_type(grid, spectrum)
  
  if (.not. allocated(state % depth)) then
    write(stderr, '(a)') 'state % depth is allocated.. failed'
    ok = .false.
  end if

  if (.not. allocated(state % variance)) then
    write(stderr, '(a)') 'state % variance is allocated.. failed'
    ok = .false.
  end if

  if (.not. allocated(state % wavenumber)) then
    write(stderr, '(a)') 'state % wavenumber is allocated.. failed'
    ok = .false.
  end if

  if (.not. all(shape(state % variance) == [num_frequencies, num_directions, grid % size_x, grid % size_y])) then
    write(stderr, '(a)') 'state % variance has expected shape.. failed'
    ok = .false.
  end if

  if (.not. all(shape(state % wavenumber) == [num_frequencies, grid % size_x, grid % size_y])) then
    write(stderr, '(a)') 'state % wavenumber has expected shape.. failed'
    ok = .false.
  end if

  if (.not. all(state % wavenumber == state % wavenumber)) then
    write(stderr, '(a)') 'All state % wavenumber values are finite.. failed'
    ok = .false.
  end if

  if (.not. all(state % wavenumber > 0)) then
    write(stderr, '(a)') 'All state % wavenumber values are positive.. failed'
    ok = .false.
  end if

  if (.not. all(state % phase_speed > 0)) then
    write(stderr, '(a)') 'All state % phase_speed values are positive.. failed'
    ok = .false.
  end if

  if (.not. all(state % group_speed > 0)) then
    write(stderr, '(a)') 'All state % group_speed values are positive.. failed'
    ok = .false.
  end if

  if (ok) then
    write(stdout, '(a)') 'test_state: All tests passed.'
  else
    write(stderr, '(a)') 'test_state: One or more tests failed.'
    error stop 1
  end if

end program test_state