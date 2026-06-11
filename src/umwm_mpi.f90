module umwm_mpi
  ! Subroutines to exchange data between processors in parallel mode.

#ifdef MPI

  use mpi
  use umwm_config, only: config_type
  use umwm_module, only: e, first_col_len, i_exchange_indices, iend, &
                         ierr, iiend, iistart, im, istart, &
                         last_col_len, mpisize, nproc, om, pm
  use umwm_spectrum, only: spectrum_type

  implicit none

  integer :: ilen, im_mod, status(MPI_STATUS_SIZE)

  integer, allocatable :: istart_(:), iend_(:), ilen_(:)
  integer, allocatable :: iistart_(:), iiend_(:)

  interface gather_array
    module procedure :: gather_array_1d
    module procedure :: gather_array_2d
    module procedure :: gather_array_3d
  end interface gather_array

contains

  subroutine exchange_halo(config, spectrum)
    ! Exchange halo points between processes.
    ! This version valid only for 1-cell halo width.

    type(config_type), intent(in) :: config
    type(spectrum_type), intent(in) :: spectrum
    integer :: sendcount, recvcount
    integer :: sendtag, recvtag
    integer :: src, dest

    if (nproc < mpisize - 1) then ! communicate with process above

      sendcount = spectrum % num_frequencies * spectrum % num_directions &
                * (iend - iistart_(nproc + 1) + 1)
      dest = nproc + 1
      sendtag = nproc
      recvcount = spectrum % num_frequencies * spectrum % num_directions &
                * (iiend - iend)
      src = nproc + 1
      recvtag = src

      call mpi_sendrecv(e(:,:,iistart_(nproc+1):iend), sendcount,&
                        MPI_REAL, dest, sendtag,                 &
                        e(:,:,iend+1:iiend), recvcount,          &
                        MPI_REAL, src, recvtag,                  &
                        MPI_COMM_WORLD, status, ierr)

    end if

    if (nproc > 0) then ! communicate with process below

      sendcount = spectrum % num_frequencies * spectrum % num_directions &
                * (iiend_(nproc - 1) - istart + 1)
      dest = nproc-1
      sendtag = nproc
      recvcount = spectrum % num_frequencies * spectrum % num_directions &
                * (istart - iistart)
      src  = nproc - 1
      recvtag = src

      call mpi_sendrecv(e(:,:,istart:iiend_(nproc-1)), sendcount,&
                        MPI_REAL, dest, sendtag,                 &
                        e(:,:,iistart:istart-1), recvcount,      &
                        MPI_REAL, src, recvtag,                  &
                        MPI_COMM_WORLD, status, ierr)
    end if

    if (config % isglobal) then ! if periodic domain

      if (nproc == 0) then

        sendcount = spectrum % num_frequencies * spectrum % num_directions &
                  * first_col_len
        dest = mpisize - 1
        sendtag = nproc
        recvcount = spectrum % num_frequencies * spectrum % num_directions &
                  * last_col_len
        src  = mpisize - 1
        recvtag = src

        call mpi_sendrecv(e(:,:,i_exchange_indices), sendcount,&
                          MPI_REAL, dest, sendtag,             &
                          e(:,:,iistart:istart-1), recvcount,  &
                          MPI_REAL, src, recvtag,              &
                          MPI_COMM_WORLD, status, ierr)

      end if

      if (nproc == mpisize - 1) then

        sendcount = spectrum % num_frequencies * spectrum % num_directions &
                  * last_col_len
        dest = 0
        sendtag = nproc
        recvcount = spectrum % num_frequencies * spectrum % num_directions &
                  * first_col_len
        src = 0
        recvtag = src

        call mpi_sendrecv(e(:,:,i_exchange_indices), sendcount,&
                          MPI_REAL, dest, sendtag,             &
                          e(:,:,iend+1:iiend), recvcount,      &
                          MPI_REAL, src, recvtag,              &
                          MPI_COMM_WORLD, status, ierr)

      end if

    end if

  end subroutine exchange_halo


  subroutine gather_array_1d(srcarray, tgtarray)
    ! Gathers an array of shape (im) to root processor.
    ! Non-blocking implementation.

    real, intent(in) :: srcarray(istart:iend)
    real, intent(out) :: tgtarray(im)

    integer :: nn, requests(mpisize - 1), statuses(MPI_STATUS_SIZE, mpisize - 1)

    if (nproc == 0) then
      tgtarray(istart:iend) = srcarray(istart:iend)
      do nn = 1, mpisize - 1
        call mpi_irecv(tgtarray(istart_(nn):iend_(nn)), ilen_(nn),&
                       MPI_REAL, nn, nn, MPI_COMM_WORLD, requests(nn), ierr)
      end do
    end if

    call mpi_barrier(MPI_COMM_WORLD,ierr)

    if (nproc == 0) then
      call mpi_waitall(mpisize - 1, requests, statuses, ierr)
    else
      call mpi_isend(srcarray(istart:iend), ilen,&
                     MPI_REAL, 0, nproc, MPI_COMM_WORLD, requests(nproc), ierr)
    end if

  end subroutine gather_array_1d


  subroutine gather_array_2d(srcarray, tgtarray)
    ! Gathers an array of shape (om,im) to root processor.
    ! Non-blocking implementation.

    real, intent(in)  :: srcarray(om,istart:iend)
    real, intent(out) :: tgtarray(om,im)

    integer :: nn, requests(mpisize - 1), statuses(MPI_STATUS_SIZE, mpisize - 1)

    if (nproc == 0) then
      tgtarray(:,istart:iend) = srcarray(:,istart:iend)
      do nn = 1, mpisize - 1
        call mpi_irecv(tgtarray(:,istart_(nn):iend_(nn)), om * ilen_(nn),&
                       MPI_REAL, nn, nn, MPI_COMM_WORLD, requests(nn), ierr)
      end do
    end if

    call mpi_barrier(MPI_COMM_WORLD, ierr)

    if (nproc == 0) then
      call mpi_waitall(mpisize - 1, requests, statuses, ierr)
    else
      call mpi_isend(srcarray(:,istart:iend), om * ilen,&
                     MPI_REAL, 0, nproc, MPI_COMM_WORLD, requests(nproc), ierr)
    end if

  end subroutine gather_array_2d


  subroutine gather_array_3d(srcarray, tgtarray)
    ! Gathers an array of shape (om,pm,im) to root processor.
    ! Non-blocking implementation.

    real, intent(in)  :: srcarray(om,pm,istart:iend)
    real, intent(out) :: tgtarray(om,pm,im)

    integer :: nn, requests(mpisize - 1), statuses(MPI_STATUS_SIZE, mpisize - 1)

    if (nproc == 0) then
      tgtarray(:,:,istart:iend) = srcarray(:,:,istart:iend)
      do nn = 1, mpisize - 1
        call mpi_irecv(tgtarray(:,:,istart_(nn):iend_(nn)), om * pm * ilen_(nn),&
                       MPI_REAL, nn, nn, MPI_COMM_WORLD, requests(nn), ierr)
      end do
    end if

    call mpi_barrier(MPI_COMM_WORLD, ierr)

    if (nproc == 0) then
      call mpi_waitall(mpisize - 1, requests, statuses, ierr)
    else
      call mpi_isend(srcarray(:,:,istart:iend), om * pm * ilen,&
        MPI_REAL, 0, nproc, MPI_COMM_WORLD, requests(nproc), ierr)
    end if

  end subroutine gather_array_3d

#endif

end module umwm_mpi
