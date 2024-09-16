module umwm_state

  use umwm_dispersion, only: wavenumber
  use umwm_grid, only: grid_type
  use umwm_spectrum, only: spectrum_type

  implicit none

  private
  public :: state_type

  type :: state_type
    real, allocatable :: depth(:,:) ! Mean water depth; may vary in time.
    real, allocatable :: variance(:,:,:,:)
    real, allocatable :: wavenumber(:,:,:)
  end type state_type

  interface state_type
    module procedure :: state_type_cons
  end interface state_type

contains

  elemental type(state_type) function state_type_cons(grid, spectrum) result(res)

    type(grid_type), intent(in) :: grid
    type(spectrum_type), intent(in) :: spectrum
    integer :: i, j
    integer :: stat
   
    allocate(res % depth(grid % size_x, grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % state % action.'
    res % depth = 100 ! FIXME load from file, CLI, or config file

    allocate(res % variance(spectrum % num_frequencies, &
                            spectrum % num_directions, &
                            grid % size_x, &
                            grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % state % action.'
    res % variance = 0

    allocate(res % wavenumber(spectrum % num_frequencies, &
                              grid % size_x, &
                              grid % size_y), stat=stat)
    if (stat /= 0) error stop 'Error allocating domain % state % wavenumber.'

    ! Initialize by solving the dispersion relationship.
    do j = 1, grid % size_y
      do i = 1, grid % size_x
        res % wavenumber(:,i,j) = wavenumber(spectrum % frequency, res % depth(i,j), 1e3, 9.8, 0.074)
      end do
    end do

  end function state_type_cons

end module umwm_state
