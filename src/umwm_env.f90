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
#ifdef MPI
    integer :: ierr
#ifndef ESMF
    call mpi_init(ierr)                                ! initialize mpi
#endif
    call mpi_comm_rank(MPI_COMM_WORLD, nproc, ierr)    ! who am i?
    call mpi_comm_size(MPI_COMM_WORLD, mpisize, ierr)  ! how many processes?
#endif
  end subroutine env_init


  subroutine env_stop()
#if defined(MPI) && !defined(ESMF)
    integer :: ierr
    call mpi_barrier(MPI_COMM_WORLD, ierr) ! wait for all
    call mpi_finalize(ierr)                ! finalize mpi
#endif
  end subroutine env_stop

end module umwm_env
