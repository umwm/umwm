module umwm_mpi
  ! Subroutines to exchange data between processors in parallel mode.

#ifdef MPI

  use mpi
  use umwm_config, only: config_type
  use umwm_grid, only: grid_type
  use umwm_module, only: e, ierr, mpisize, nproc, om, pm
  use umwm_spectrum, only: spectrum_type

  implicit none

  integer :: status(MPI_STATUS_SIZE)

  interface gather_array
    module procedure :: gather_array_1d
    module procedure :: gather_array_2d
    module procedure :: gather_array_3d
  end interface gather_array

contains

  subroutine exchange_halo(config, spectrum, grid)
    ! Exchange halo points between processes.
    ! This version valid only for 1-cell halo width.

    type(config_type), intent(in) :: config
    type(spectrum_type), intent(in) :: spectrum
    type(grid_type), intent(in) :: grid
    integer :: sendcount, recvcount
    integer :: sendtag, recvtag
    integer :: src, dest

    if (nproc < mpisize - 1) then ! communicate with process above

      sendcount = spectrum % num_frequencies * spectrum % num_directions &
                * (grid % iend - grid % iistart_all(nproc + 1) + 1)
      dest = nproc + 1
      sendtag = nproc
      recvcount = spectrum % num_frequencies * spectrum % num_directions &
                * (grid % iiend - grid % iend)
      src = nproc + 1
      recvtag = src

      call mpi_sendrecv(e(:,:,grid % iistart_all(nproc+1):grid % iend), sendcount,&
                        MPI_REAL, dest, sendtag,                 &
                        e(:,:,grid % iend+1:grid % iiend), recvcount,          &
                        MPI_REAL, src, recvtag,                  &
                        MPI_COMM_WORLD, status, ierr)

    end if

    if (nproc > 0) then ! communicate with process below

      sendcount = spectrum % num_frequencies * spectrum % num_directions &
                * (grid % iiend_all(nproc - 1) - grid % istart + 1)
      dest = nproc-1
      sendtag = nproc
      recvcount = spectrum % num_frequencies * spectrum % num_directions &
                * (grid % istart - grid % iistart)
      src  = nproc - 1
      recvtag = src

      call mpi_sendrecv(e(:,:,grid % istart:grid % iiend_all(nproc-1)), sendcount,&
                        MPI_REAL, dest, sendtag,                 &
                        e(:,:,grid % iistart:grid % istart-1), recvcount,      &
                        MPI_REAL, src, recvtag,                  &
                        MPI_COMM_WORLD, status, ierr)
    end if

    if (config % isglobal) then ! if periodic domain

      if (nproc == 0) then

        sendcount = spectrum % num_frequencies * spectrum % num_directions &
                  * grid % first_col_len
        dest = mpisize - 1
        sendtag = nproc
        recvcount = spectrum % num_frequencies * spectrum % num_directions &
                  * grid % last_col_len
        src  = mpisize - 1
        recvtag = src

        call mpi_sendrecv(e(:,:,grid % i_exchange_indices), sendcount,&
                          MPI_REAL, dest, sendtag,             &
                          e(:,:,grid % iistart:grid % istart-1), recvcount,  &
                          MPI_REAL, src, recvtag,              &
                          MPI_COMM_WORLD, status, ierr)

      end if

      if (nproc == mpisize - 1) then

        sendcount = spectrum % num_frequencies * spectrum % num_directions &
                  * grid % last_col_len
        dest = 0
        sendtag = nproc
        recvcount = spectrum % num_frequencies * spectrum % num_directions &
                  * grid % first_col_len
        src = 0
        recvtag = src

        call mpi_sendrecv(e(:,:,grid % i_exchange_indices), sendcount,&
                          MPI_REAL, dest, sendtag,             &
                          e(:,:,grid % iend+1:grid % iiend), recvcount,      &
                          MPI_REAL, src, recvtag,              &
                          MPI_COMM_WORLD, status, ierr)

      end if

    end if

  end subroutine exchange_halo


  subroutine gather_array_1d(srcarray, tgtarray, grid)
    ! Gathers an array of shape (im) to root processor.
    ! Non-blocking implementation.

    type(grid_type), intent(in) :: grid
    real, intent(in) :: srcarray(grid % istart:grid % iend)
    real, intent(out) :: tgtarray(grid % im)

    integer :: nn, requests(mpisize - 1), statuses(MPI_STATUS_SIZE, mpisize - 1)

    if (nproc == 0) then
      tgtarray(grid % istart:grid % iend) = srcarray(grid % istart:grid % iend)
      do nn = 1, mpisize - 1
        call mpi_irecv(tgtarray(grid % istart_all(nn):grid % iend_all(nn)), grid % ilen_all(nn),&
                       MPI_REAL, nn, nn, MPI_COMM_WORLD, requests(nn), ierr)
      end do
    end if

    call mpi_barrier(MPI_COMM_WORLD,ierr)

    if (nproc == 0) then
      call mpi_waitall(mpisize - 1, requests, statuses, ierr)
    else
      call mpi_isend(srcarray(grid % istart:grid % iend), grid % ilen,&
                     MPI_REAL, 0, nproc, MPI_COMM_WORLD, requests(nproc), ierr)
    end if

  end subroutine gather_array_1d


  subroutine gather_array_2d(srcarray, tgtarray, grid)
    ! Gathers an array of shape (om,im) to root processor.
    ! Non-blocking implementation.

    type(grid_type), intent(in) :: grid
    real, intent(in)  :: srcarray(om,grid % istart:grid % iend)
    real, intent(out) :: tgtarray(om,grid % im)

    integer :: nn, requests(mpisize - 1), statuses(MPI_STATUS_SIZE, mpisize - 1)

    if (nproc == 0) then
      tgtarray(:,grid % istart:grid % iend) = srcarray(:,grid % istart:grid % iend)
      do nn = 1, mpisize - 1
        call mpi_irecv(tgtarray(:,grid % istart_all(nn):grid % iend_all(nn)), om * grid % ilen_all(nn),&
                       MPI_REAL, nn, nn, MPI_COMM_WORLD, requests(nn), ierr)
      end do
    end if

    call mpi_barrier(MPI_COMM_WORLD, ierr)

    if (nproc == 0) then
      call mpi_waitall(mpisize - 1, requests, statuses, ierr)
    else
      call mpi_isend(srcarray(:,grid % istart:grid % iend), om * grid % ilen,&
                     MPI_REAL, 0, nproc, MPI_COMM_WORLD, requests(nproc), ierr)
    end if

  end subroutine gather_array_2d


  subroutine gather_array_3d(srcarray, tgtarray, grid)
    ! Gathers an array of shape (om,pm,im) to root processor.
    ! Non-blocking implementation.

    type(grid_type), intent(in) :: grid
    real, intent(in)  :: srcarray(om,pm,grid % istart:grid % iend)
    real, intent(out) :: tgtarray(om,pm,grid % im)

    integer :: nn, requests(mpisize - 1), statuses(MPI_STATUS_SIZE, mpisize - 1)

    if (nproc == 0) then
      tgtarray(:,:,grid % istart:grid % iend) = srcarray(:,:,grid % istart:grid % iend)
      do nn = 1, mpisize - 1
        call mpi_irecv(tgtarray(:,:,grid % istart_all(nn):grid % iend_all(nn)), &
                       om * pm * grid % ilen_all(nn),&
                       MPI_REAL, nn, nn, MPI_COMM_WORLD, requests(nn), ierr)
      end do
    end if

    call mpi_barrier(MPI_COMM_WORLD, ierr)

    if (nproc == 0) then
      call mpi_waitall(mpisize - 1, requests, statuses, ierr)
    else
      call mpi_isend(srcarray(:,:,grid % istart:grid % iend), om * pm * grid % ilen,&
        MPI_REAL, 0, nproc, MPI_COMM_WORLD, requests(nproc), ierr)
    end if

  end subroutine gather_array_3d

#endif

end module umwm_mpi
