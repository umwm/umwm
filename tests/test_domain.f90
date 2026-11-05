program test_domain
#ifdef MPI
  use mpi
#endif
  use netcdf
  use tuff, only: test, test_result
  use umwm_constants, only: stderr
  use umwm_domain, only: domain_type
  use umwm_env, only: env_init, env_stop, nproc

  implicit none

  type(test_result) :: suite

  call env_init()
  suite = test('test_domain', [ &
    test(rankine_regression), &
    test(waves_initial_state), &
    test(waves_open_boundaries) &
  ])
  call env_stop()
  if (.not. suite % ok) error stop 1

contains

  function rankine_regression() result(res)
    ! Integrates the idealized hurricane case for one hour with gridded and
    ! spectrum output and Stokes drift, and compares the domain mean and
    ! maximum of output fields to reference values. The uniform current over
    ! the flat bottom neither refracts the waves nor shifts their frequency,
    ! so the values are those of UMWM before the domain refactor without the
    ! refraction by its spurious current gradients against the land at the
    ! edges of the domain.
    type(test_result) :: res
    type(domain_type) :: domain
    logical :: ok

    call domain % initialize('namelists/forcing_rankine.nml', output_dir='output')
    call domain % run()

    ! The drag coefficient is ill-conditioned in the calm vortex center and
    ! varies by ~2e-4 between optimized and unoptimized builds.
    associate(n => domain % grid % im)
      ok = matches(domain % ht, n, 3.811970711, 6.314935684, 1e-4) &
           .and. matches(domain % mwp, n, 5.527236938, 6.882996559, 1e-4) &
           .and. matches(domain % mss, n, 5.120557919e-02, 6.678768247e-02, 1e-4) &
           .and. matches(domain % cd, n, 5.015279166e-03, 4.853604138e-01, 1e-3) &
           .and. matches(domain % ustar, n, 9.779370427e-01, 1.610205531, 1e-4) &
           .and. matches(domain % dwp, n, 7.718040943, 10.06773663, 1e-4) &
           .and. matches(domain % taux_form, n, -2.223200165e-02, 2.234684467, 1e-4) &
           .and. matches(domain % taux_ocntop, n, -2.439257503e-02, 2.165755749, 1e-4) &
           .and. matches(domain % epsx_ocn, n, -3.125989065e-02, 4.089962482, 1e-4) &
           .and. matches(domain % epsy_atm, n, 1.487540901e-01, 1.229317474e+01, 1e-4) &
           .and. matches(domain % stokes % us(:,1), n, -5.282592494e-03, 3.113176227e-01, 1e-4) &
           .and. matches(domain % stokes % vs(:,1), n, -5.729752127e-03, 3.947873712e-01, 1e-4)
    end associate

    call domain % finalize()

    res = test('rankine_regression', ok)

  end function rankine_regression


  function waves_initial_state() result(res)
    ! A global domain starts from the spectrum of two steady, uniform wave
    ! partitions, whose significant wave heights add up exactly, and has no
    ! open boundaries.
    type(test_result) :: res
    type(domain_type) :: domain
    logical :: ok

    call write_partitions('output/waves_global/umwmin.nc', 21, 11, &
                          [1., 2.], [1., 2.], [6., 12.], [0.5, -2.])
    call domain % initialize('namelists/waves_global.nml', input_dir='output/waves_global', &
                             output_dir='output/waves_global')
    call domain % diagnostics()
    ok = domain % num_boundary_all == 0 &
         .and. all(abs(domain % ht - sqrt(5.)) <= 1e-4 * sqrt(5.))
    call domain % finalize()

    res = test('waves_initial_state', ok)
  end function waves_initial_state


  function waves_open_boundaries() result(res)
    ! A regional domain starts calm, with 2 m of eastward swell at the open
    ! boundary cells of its western and eastern edges, rising to 4 m an
    ! hour later, and calm at its northern and southern edges. The boundary
    ! cells take the spectrum of the input files, interpolated in time, and
    ! the swell propagates into the domain from its western edge, and not
    ! from its eastern one, from which it travels out of the domain.
    type(test_result) :: res
    type(domain_type) :: domain
    logical :: ok

    call write_partitions('output/waves_regional/umwmin_2026-01-01_00:00:00.nc', 21, 11, &
                          [2.], [0.], [12.], [0.])
    call write_partitions('output/waves_regional/umwmin_2026-01-01_01:00:00.nc', 21, 11, &
                          [4.], [0.], [12.], [0.])
    call domain % initialize('namelists/waves_regional.nml', input_dir='output/waves_regional', &
                             output_dir='output/waves_regional')
    call domain % diagnostics()
    ok = domain % num_boundary_all == 2 * 19 + 2 * 9 - 4 &
         .and. all(abs(height(domain, boundary(domain, .true.)) - 2) <= 1e-4 * 2) &
         .and. all(height(domain, boundary(domain, .false.)) < 1e-3) &
         .and. all(height(domain, interior(domain)) < 1e-3)
    call domain % run()
    ok = ok .and. all(abs(height(domain, boundary(domain, .true.)) - 4) <= 1e-4 * 4) &
         .and. all(height(domain, boundary(domain, .false.)) < 1e-3)
    ! 20 km downwave of the western edge, and 20 km upwave of the eastern one
    ok = ok .and. above(domain, 4, 6, 2.) .and. .not. above(domain, 18, 6, 0.01)
    call domain % finalize()

    res = test('waves_open_boundaries', ok)
  end function waves_open_boundaries


  function boundary(domain, with_swell) result(cells)
    ! The open boundary cells of this process next to the western and
    ! eastern edges, which have swell, or the others, which are calm.
    type(domain_type), intent(in) :: domain
    logical, intent(in) :: with_swell
    integer, allocatable :: cells(:)
    logical :: east_west(domain % num_boundary)
    integer :: j
    east_west = [(any(domain % grid % mi(domain % boundary(j)) == [2, domain % grid % mm - 1]), &
                  j = 1, domain % num_boundary)]
    cells = pack(domain % boundary, east_west .eqv. with_swell)
  end function boundary


  function height(domain, cells) result(values)
    ! The significant wave height at the cells of this process.
    type(domain_type), intent(in) :: domain
    integer, intent(in) :: cells(:)
    real :: values(size(cells))
    values = domain % ht(cells)
  end function height


  function interior(domain) result(cells)
    ! The sea cells of this process that are not open boundary cells.
    type(domain_type), intent(in) :: domain
    integer, allocatable :: cells(:)
    integer :: i
    cells = [(i, i = domain % grid % istart, domain % grid % iend)]
    cells = pack(cells, [(.not. any(domain % boundary == i), i = domain % grid % istart, &
                                                               domain % grid % iend)])
  end function interior


  logical function above(domain, m, n, threshold)
    ! Whether the significant wave height at cell (m,n) is above threshold,
    ! on the process that has the cell.
    type(domain_type), intent(in) :: domain
    integer, intent(in) :: m, n
    real, intent(in) :: threshold
    integer :: i, flag
#ifdef MPI
    integer :: buffer, ierr
#endif
    i = domain % grid % ii(m,n)
    flag = 0
    if (i >= domain % grid % istart .and. i <= domain % grid % iend) then
      flag = merge(1, -1, domain % ht(i) > threshold)
      if (domain % ht(i) <= threshold .eqv. threshold > 0.5) &
        write(stderr, '(a,2(i0,1x),a,f8.4)') 'test_domain: height at ', m, n, ':', domain % ht(i)
    end if
#ifdef MPI
    buffer = flag
    call mpi_allreduce(buffer, flag, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    above = flag > 0
  end function above


  subroutine write_partitions(path, mm, nm, hs_boundary, hs_inside, tm, dir)
    ! Writes an input file of wave partitions on an mm by nm grid, with the
    ! significant wave heights hs_boundary on its western and eastern edges
    ! and the cells next to them, and hs_inside elsewhere; the root process
    ! writes it.
    character(*), intent(in) :: path
    integer, intent(in) :: mm, nm
    real, intent(in) :: hs_boundary(:), hs_inside(:), tm(:), dir(:)
    real :: values(mm,nm,1)
    integer :: ncid, dims(3), j
    character(1) :: p
#ifdef MPI
    integer :: ierr
#endif

    if (nproc == 0) then
      call execute_command_line('mkdir -p ' // path(1:index(path, '/', back=.true.)))
      call nc_check(nf90_create(path, nf90_clobber, ncid), path)
      call nc_check(nf90_def_dim(ncid, 'x', mm, dims(1)), path)
      call nc_check(nf90_def_dim(ncid, 'y', nm, dims(2)), path)
      call nc_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dims(3)), path)
      call nc_check(nf90_enddef(ncid), path)
      do j = 1, size(tm)
        write(p, '(i1)') j - 1
        values = hs_inside(j)
        values([1, 2, mm-1, mm],:,1) = hs_boundary(j)
        call put_field(ncid, dims, 'phs' // p, values, path)
        values = tm(j)
        call put_field(ncid, dims, 'ptm' // p, values, path)
        values = dir(j)
        call put_field(ncid, dims, 'pdir' // p, values, path)
        values = 0.3
        call put_field(ncid, dims, 'pspr' // p, values, path)
      end do
      call nc_check(nf90_close(ncid), path)
    end if
#ifdef MPI
    call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

  end subroutine write_partitions


  subroutine put_field(ncid, dims, name, values, path)
    ! Defines the variable name of dimensions dims in the open file ncid and
    ! writes its values.
    integer, intent(in) :: ncid, dims(3)
    character(*), intent(in) :: name, path
    real, intent(in) :: values(:,:,:)
    integer :: varid
    call nc_check(nf90_redef(ncid), path)
    call nc_check(nf90_def_var(ncid, name, nf90_float, dims, varid), path)
    call nc_check(nf90_enddef(ncid), path)
    call nc_check(nf90_put_var(ncid, varid, values), path)
  end subroutine put_field


  subroutine nc_check(stat, path)
    integer, intent(in) :: stat
    character(*), intent(in) :: path
    if (stat /= nf90_noerr) then
      write(stderr, '(a)') 'test_domain: ' // path // ': ' // trim(nf90_strerror(stat))
      error stop 1
    end if
  end subroutine nc_check


  logical function matches(field, n, mean, maximum, tolerance)
    ! Whether the mean and maximum of field over all n cells of the domain
    ! match within the tolerance relative to the larger of their magnitudes,
    ! because the mean of a field of either sign may nearly cancel.
    real, intent(in) :: field(:)
    integer, intent(in) :: n
    real, intent(in) :: mean, maximum, tolerance
    real :: total, largest, scale
#ifdef MPI
    real :: buffer
    integer :: ierr
#endif
    total = sum(field)
    largest = maxval(field)
#ifdef MPI
    buffer = total
    call mpi_allreduce(buffer, total, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
    buffer = largest
    call mpi_allreduce(buffer, largest, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
    scale = max(abs(mean), abs(maximum))
    matches = abs(total / n - mean) <= tolerance * scale &
              .and. abs(largest - maximum) <= tolerance * scale
    if (.not. matches) write(stderr, '(a,4(es15.7))') 'test_domain: mean, max, expected: ', &
      total / n, largest, mean, maximum
  end function matches

end program test_domain
