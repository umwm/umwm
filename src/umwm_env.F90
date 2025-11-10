module umwm_env

#ifdef MPI
  use mpi
#endif

  use umwm_module, only: nproc, mpisize, ierr

  implicit none

  private
  public :: env_init, env_stop

contains

  subroutine env_init()
#ifdef MPI
#ifndef ESMF
      call mpi_init(ierr)                           ! initialize mpi
#endif
      call mpi_comm_rank(MPI_COMM_WORLD, nproc, ierr) ! who am i?
      call mpi_comm_size(MPI_COMM_WORLD, mpisize, ierr)  ! how many processes?
#else
      nproc = 0
      mpisize = 1
#endif
  end subroutine env_init


  subroutine env_stop()
#if defined(MPI) && !defined(ESMF)
    call mpi_barrier(MPI_COMM_WORLD, ierr) ! wait for all
    call mpi_finalize(ierr)                ! finalize mpi
#endif
  end subroutine env_stop

end module umwm_env