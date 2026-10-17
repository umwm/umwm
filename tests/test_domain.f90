program test_domain
#ifdef MPI
  use mpi
#endif
  use tuff, only: test, test_result
  use umwm_constants, only: stderr
  use umwm_domain, only: domain_type
  use umwm_env, only: env_init, env_stop

  implicit none

  type(test_result) :: suite

  call env_init()
  suite = test('test_domain', [test(rankine_regression)])
  call env_stop()
  if (.not. suite % ok) error stop 1

contains

  function rankine_regression() result(res)
    ! Integrates the idealized hurricane case for one hour with gridded and
    ! spectrum output and Stokes drift, and compares the domain mean and
    ! maximum of output fields to those from UMWM before the domain refactor.
    type(test_result) :: res
    type(domain_type) :: domain
    logical :: ok

    call domain % initialize('namelists/forcing_rankine.nml')
    call domain % run()

    ! The drag coefficient is ill-conditioned in the calm vortex center and
    ! varies by ~2e-4 between optimized and unoptimized builds.
    associate(n => domain % grid % im)
      ok = matches(domain % ht, n, 3.812207232, 6.315393925, 1e-4) &
           .and. matches(domain % mwp, n, 5.526470696, 6.882972240, 1e-4) &
           .and. matches(domain % mss, n, 5.120974732e-02, 6.678850949e-02, 1e-4) &
           .and. matches(domain % cd, n, 5.015606267e-03, 4.854240119e-01, 1e-3) &
           .and. matches(domain % ustar, n, 9.779266109e-01, 1.610196948, 1e-4) &
           .and. matches(domain % dwp, n, 7.726010769, 10.06773663, 1e-4) &
           .and. matches(domain % taux_form, n, -2.219776792e-02, 2.233852386, 1e-4) &
           .and. matches(domain % taux_ocntop, n, -2.461001024e-02, 2.166157484, 1e-4) &
           .and. matches(domain % epsx_ocn, n, -3.314319944e-02, 4.064507008, 1e-4) &
           .and. matches(domain % epsy_atm, n, 1.481968960e-01, 1.229311371e+01, 1e-4) &
           .and. matches(domain % stokes % us(:,1), n, -5.597803407e-03, 3.102434874e-01, 1e-4) &
           .and. matches(domain % stokes % vs(:,1), n, -5.990615025e-03, 3.947814107e-01, 1e-4)
    end associate

    call domain % finalize()

    res = test('rankine_regression', ok)

  end function rankine_regression


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
