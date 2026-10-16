module umwm_mpi
  ! Data exchange between processes in parallel mode.

#ifdef MPI

  use mpi
  use umwm_env, only: mpisize, nproc
  use umwm_grid, only: grid_type

  implicit none

  private
  public :: exchange_halo, gather_array

contains

  subroutine exchange_halo(field, n, grid, isglobal)
    ! Exchange the 1-cell wide halo of field(n, iistart-1:iiend) with the
    ! neighboring tiles, and across the east-west seam if the domain is global.
    integer, intent(in) :: n
    type(grid_type), intent(in) :: grid
    real, intent(inout) :: field(n, grid % iistart-1:grid % iiend)
    logical, intent(in) :: isglobal

    real, allocatable :: seam(:,:)
    integer :: requests(6), nreq, ierr, last, seam_tag

    last = mpisize - 1
    seam_tag = mpisize ! seam messages use tags mpisize + rank
    nreq = 0

    if (nproc < last) then ! tile above
      associate(lb => grid % iistart_all(nproc + 1))
        call irecv(field(:, grid % iend+1:grid % iiend), nproc + 1, nproc + 1)
        call isend(field(:, lb:grid % iend), nproc + 1, nproc)
      end associate
    end if

    if (nproc > 0) then ! tile below
      associate(ub => grid % iiend_all(nproc - 1))
        call irecv(field(:, grid % iistart:grid % istart-1), nproc - 1, nproc - 1)
        call isend(field(:, grid % istart:ub), nproc - 1, nproc)
      end associate
    end if

    if (isglobal .and. (nproc == 0 .or. nproc == last)) then ! periodic seam
      seam = field(:, grid % i_exchange_indices)
      if (nproc == 0) then
        call irecv(field(:, grid % iistart:grid % istart-1), last, seam_tag + last)
        call isend(seam, last, seam_tag)
      end if
      if (nproc == last) then
        call irecv(field(:, grid % iend+1:grid % iiend), 0, seam_tag)
        call isend(seam, 0, seam_tag + last)
      end if
    end if

    call mpi_waitall(nreq, requests, MPI_STATUSES_IGNORE, ierr)

  contains

    subroutine irecv(buf, src, tag)
      real, intent(inout), contiguous :: buf(:,:)
      integer, intent(in) :: src, tag
      nreq = nreq + 1
      call mpi_irecv(buf, size(buf), MPI_REAL, src, tag, MPI_COMM_WORLD, requests(nreq), ierr)
    end subroutine irecv

    subroutine isend(buf, dest, tag)
      real, intent(in), contiguous :: buf(:,:)
      integer, intent(in) :: dest, tag
      nreq = nreq + 1
      call mpi_isend(buf, size(buf), MPI_REAL, dest, tag, MPI_COMM_WORLD, requests(nreq), ierr)
    end subroutine isend

  end subroutine exchange_halo


  subroutine gather_array(field, field_all, grid)
    ! Gathers a tile field into field_all(1:im) on the root process.
    type(grid_type), intent(in) :: grid
    real, intent(in) :: field(grid % istart:grid % iend)
    real, intent(inout) :: field_all(grid % im)
    integer :: ierr

    call mpi_gatherv(field, grid % ilen, MPI_REAL, &
                     field_all, grid % ilen_all, grid % istart_all - 1, MPI_REAL, &
                     0, MPI_COMM_WORLD, ierr)

  end subroutine gather_array

#endif

end module umwm_mpi
