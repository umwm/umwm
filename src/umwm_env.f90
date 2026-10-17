module umwm_env
  ! Parallel execution environment: process rank and number of processes.

#ifdef MPI
  use mpi
#endif

  implicit none

  private
  public :: env_init, env_stop, nproc, mpisize

  integer, protected :: nproc = 0 ! process rank
  integer, protected :: mpisize = 1 ! number of processes

contains

  subroutine env_init()
    ! Initializes MPI, unless the caller (for example mpi4py) already did,
    ! and queries the process rank and the number of processes.
#ifdef MPI
    integer :: ierr
#ifndef ESMF
    logical :: initialized
    call mpi_initialized(initialized, ierr)
    if (.not. initialized) call mpi_init(ierr)         ! initialize mpi
#endif
    call mpi_comm_rank(MPI_COMM_WORLD, nproc, ierr)    ! who am i?
    call mpi_comm_size(MPI_COMM_WORLD, mpisize, ierr)  ! how many processes?
#endif
  end subroutine env_init


  subroutine env_stop()
    ! Finalizes MPI, unless it is already finalized.
#if defined(MPI) && !defined(ESMF)
    integer :: ierr
    logical :: finalized
    call mpi_finalized(finalized, ierr)
    if (finalized) return
    call mpi_barrier(MPI_COMM_WORLD, ierr) ! wait for all
    call mpi_finalize(ierr)                ! finalize mpi
#endif
  end subroutine env_stop

end module umwm_env
