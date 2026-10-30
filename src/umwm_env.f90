module umwm_env
  ! Parallel execution environment: process rank and number of processes,
  ! and in the GPU build the GPU of the process.

#ifdef MPI
  use mpi
#ifdef _OPENACC
  use openacc, only: acc_device_nvidia, acc_get_num_devices, acc_set_device_num
#endif
#endif

  implicit none

  private
  public :: env_init, env_stop, nproc, mpisize

  integer, protected :: nproc = 0 ! process rank
  integer, protected :: mpisize = 1 ! number of processes

contains

  subroutine env_init()
    ! Initializes MPI, unless the caller (for example mpi4py) already did,
    ! and queries the process rank and the number of processes; in the GPU
    ! build, selects the GPU of the process.
#ifdef MPI
    integer :: ierr
#ifndef ESMF
    logical :: initialized
    call mpi_initialized(initialized, ierr)
    if (.not. initialized) call mpi_init(ierr)         ! initialize mpi
#endif
    call mpi_comm_rank(MPI_COMM_WORLD, nproc, ierr)    ! who am i?
    call mpi_comm_size(MPI_COMM_WORLD, mpisize, ierr)  ! how many processes?
#ifdef _OPENACC
    call select_device()
#endif
#endif
  end subroutine env_init


#if defined(MPI) && defined(_OPENACC)
  subroutine select_device()
    ! Distributes the processes of a node over its visible GPUs, in blocks of
    ! consecutive ranks when there are more processes than GPUs, so that
    ! neighboring tiles share a GPU. Mapping the processes to the NUMA
    ! domains in rank order (for example mpirun --map-by ppr:2:numa) then
    ! keeps each process next to its GPU on nodes where the GPUs are
    ! numbered by NUMA domain.
    integer :: node_comm, node_rank, node_size, num_devices, ierr
    call mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &
                             node_comm, ierr)
    call mpi_comm_rank(node_comm, node_rank, ierr)
    call mpi_comm_size(node_comm, node_size, ierr)
    call mpi_comm_free(node_comm, ierr)
    num_devices = acc_get_num_devices(acc_device_nvidia)
    if (num_devices < 1) return
    call acc_set_device_num(node_rank * num_devices / max(node_size, num_devices), &
                            acc_device_nvidia)
    if (nproc == 0) write(*,'(a,i0,a,i0)') 'umwm: env: processes per node: ', node_size, &
      ', GPUs per node: ', min(node_size, num_devices)
  end subroutine select_device
#endif


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
