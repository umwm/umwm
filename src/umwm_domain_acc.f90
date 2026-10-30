submodule (umwm_domain) umwm_domain_acc

  ! The source time step, the diagnostics, and the ocean stress on GPUs, with
  ! OpenACC. Instead of the pipeline over cells of step, a sequence of
  ! kernels each evaluates all cells, with a thread block (gang) for the
  ! spectral column of a cell and its threads (vector) over the spectral
  ! bins. The spectral bins are evaluated by the same elemental functions as
  ! in the column routines of step, and the kernels take the fields as
  ! explicit-shape arrays, which is faster than the components of the domain
  ! in device code. In the tile, i1:i2 are the sea cells and h1:h2 the cells
  ! with the halo and the land cell h1.

  use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc, c_ptr, c_size_t
  use umwm_advection, only: current_rotation_rate, depth_rotation_rate, frequency_outflow, &
                            frequency_rate, propagation_flux, refraction_outflow, upstream_flux
  use umwm_physics, only: equilibrium_bin, integrate_bin, source_rate
  use umwm_source_functions, only: power_bin, power_threshold, s_ice_bin, &
                                   sds_coefficient_bin, sin_d12_bin, snl_d12_bin, snl_gain

  implicit none

  ! largest number of directions, and of frequencies, for the arrays of a
  ! thread in the refraction and frequency shift in more than one substep
  integer, parameter :: max_bins = 128

  ! threads in a thread block of a spectral column
  integer, parameter :: vl = 64

  ! queue of the source integration, which runs during the halo exchange
  integer, parameter :: integrate_queue = 1

  ! The values that kernels return to the CPU, such as the results of
  ! reductions, in GPU memory, from which the CPU copies them
  ! (copy_results). Kernels write host variables in place, and unified
  ! memory shares pages with the small allocations of the CPU; their pages
  ! would migrate between the CPU and the GPU, each time at the cost of a
  ! page fault of 0.2-0.3 ms.
  integer, parameter :: num_results = 6
  real, pointer, contiguous :: results(:) => null()

#ifdef _OPENACC
  ! the routines of the OpenACC C interface, whose device pointers are C
  ! pointers
  interface
    type(c_ptr) function device_malloc(bytes) bind(c, name='acc_malloc')
      import :: c_ptr, c_size_t
      integer(c_size_t), value :: bytes
    end function device_malloc
    subroutine memcpy_from_device(dest, src, bytes) bind(c, name='acc_memcpy_from_device')
      import :: c_ptr, c_size_t
      type(c_ptr), value :: dest, src
      integer(c_size_t), value :: bytes
    end subroutine memcpy_from_device
  end interface
#endif

contains

  module subroutine step_acc(self)
    ! Advances the wave spectrum by one source time step on the GPU, as step.
    class(domain_type), intent(inout) :: self
    real, allocatable :: e_old(:,:,:)
    real :: dts_min(1), current(1)
    logical :: store, currents
#ifdef MPI
    real :: dts_local
    integer :: ierr
#endif

    associate(config => self % config, grid => self % grid, forcing => self % forcing, &
              om => self % spectrum % num_frequencies, pm => self % spectrum % num_directions, &
              i1 => self % grid % istart, i2 => self % grid % iend, &
              h1 => self % grid % iistart - 1, h2 => self % grid % iiend, imm => self % grid % imm)

    if (self % first .and. (om > max_bins .or. pm > max_bins)) &
      error stop 'umwm: step_acc: the GPU build supports at most 128 frequencies and directions'
    call allocate_results()

    call sources_acc(om, pm, i1, i2, h1, h2, imm, self % e, self % f, self % cth, self % sth, &
                     self % cth2, self % sth2, forcing % wspd, forcing % wdir, forcing % uc, &
                     forcing % vc, forcing % rhorat, forcing % fice, self % ustar, &
                     self % logl2overz, self % cp0, self % cg0, self % fkovg, self % k, &
                     self % kdk, self % k3dk, self % k4, self % cothkd, config % snl_fac, &
                     self % snl_bins, self % snl_frac, self % sbf, self % sdv, config % g, &
                     config % fcut_fac, config % fprog, &
                     config % sin_fac, self % fieldscale1, self % fieldscale2, self % dth, &
                     config % mss_fac, self % twopisds_fac, config % sds_power, &
                     config % sdt_fac, config % fice_lth, config % explim, self % oc, &
                     self % shelt, self % mss_moments, self % sds, self % sdt, self % sice, &
                     self % physics_time_step, results)
    call copy_results(1, dts_min)

    ! dynamic time step and in-step model time:
    self % dts = min(dts_min(1), self % dtamin, config % dtg - self % sumt)
#ifdef MPI
    dts_local = self % dts
    call mpi_allreduce(dts_local, self % dts, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif

    ! if first time step, set time step to zero and
    ! integrate only the diagnostic part
    if (.not. config % restart .and. self % first) self % dts = 0

    self % sumt = self % sumt + self % dts

    ! Sin, Snl, and the diagnostic range of Sds are kept for output only, and
    ! the current velocity gradients are updated for varying currents or
    ! first step.
    store = self % first .or. self % sumt >= config % dtg
    if (config % currents .or. self % first) call self % update_current_gradients()

    ! check if currents are non-zero:
    if (.not. config % isglobal) then
      call any_current(grid % iistart, grid % iiend, imm, forcing % uc, forcing % vc, results)
    else
      call any_current(1, imm, imm, forcing % uc, forcing % vc, results)
    end if
    call copy_results(1, current)
    currents = current(1) > 0

#ifdef MPI
    ! Integrate the cells in the halo of other tiles first, and exchange the
    ! halo while the other cells are integrated.
    if (mpisize > 1) then
      call integrate(.true.)
      !$acc wait(integrate_queue)
    end if
#endif
    call integrate(.false.)
#ifdef MPI
    call exchange_halo(self % e, om * pm, grid, config % isglobal)
#endif
    !$acc wait(integrate_queue)

    call advect_acc(om, pm, i1, i2, h1, h2, imm, currents, self % dts, self % dth, &
                    real(self % spectrum % dlnf), self % e, self % ef, self % oc, &
                    forcing % fice, config % fice_uth, forcing % uc, forcing % vc, grid % ie, &
                    grid % iw, grid % in, grid % is, grid % iie, grid % iiw, grid % iin, &
                    grid % iis, grid % dxn, grid % dxs, grid % dye, grid % dyw, grid % oneovar, &
                    grid % cth_curv, grid % sth_curv, self % cg0, self % kovcg, self % invcp0, &
                    self % cgdlncp_dx, self % cgdlncp_dy, self % du_dx, self % du_dy, &
                    self % dv_dx, self % dv_dy)

    call stress_acc(om, pm, i1, i2, h1, h2, imm, self % ef, self % cth, self % sth, &
                    forcing % wspd, forcing % wdir, forcing % uc, forcing % vc, forcing % rhoa, &
                    forcing % rhow, forcing % rhorat, forcing % fice, self % shelt, self % oc, &
                    self % logl2overz, self % cp0, self % fkovg, self % invcp0, self % k, &
                    self % kdk, self % stokes % us(:,1), self % stokes % vs(:,1), &
                    config % sin_fac, self % fieldscale1, self % fieldscale2, self % dthg, &
                    config % z, config % nu_air, config % kappa, self % taux_form, &
                    self % tauy_form, self % tailatmx, self % tailatmy, self % taux_diag, &
                    self % tauy_diag, self % taux_skin, self % tauy_skin, self % taux, &
                    self % tauy, self % cd, self % ustar)

    end associate

    ! the forward spectrum becomes the current one
    call move_alloc(self % e, e_old)
    call move_alloc(self % ef, self % e)
    call move_alloc(e_old, self % ef)

  contains

    subroutine integrate(halo)
      ! Integrates the source functions at the cells in the halo of other
      ! tiles if halo is .true., and at the other cells otherwise.
      logical, intent(in) :: halo
      associate(config => self % config, forcing => self % forcing, &
                om => self % spectrum % num_frequencies, pm => self % spectrum % num_directions)
      call integrate_acc(om, pm, self % grid % istart, self % grid % iend, &
                         self % grid % iistart - 1, self % grid % iiend, self % grid % imm, &
                         halo, self % halo_source, store, self % dts, self % e, self % ef, &
                         self % f, self % cth, self % sth, self % cth2, self % sth2, &
                         forcing % wspd, forcing % wdir, forcing % uc, forcing % vc, &
                         forcing % rhorat, forcing % fice, self % ustar, self % shelt, self % oc, &
                         self % logl2overz, self % cp0, self % fkovg, self % k4, self % cothkd, &
                         self % sbf, self % sdt, self % sdv, self % sice, self % oneoverk4, &
                         self % kdk, self % mss_moments, config % sin_fac, &
                         self % fieldscale1, self % fieldscale2, config % mss_fac, &
                         self % twopisds_fac, config % sds_power, self % inv_sds_power, &
                         config % snl_fac, self % snl_bins, self % snl_frac, self % sds, &
                         self % ssin, self % snl)
      end associate
    end subroutine integrate

  end subroutine step_acc


  subroutine allocate_results()
    ! Allocates results, in the GPU build on the GPU only, unless allocated.
    if (associated(results)) return
#ifdef _OPENACC
    call c_f_pointer(device_malloc(int(num_results, c_size_t) * storage_size(1.) / 8), results, &
                     [num_results])
#else
    allocate(results(num_results))
#endif
  end subroutine allocate_results


  subroutine copy_results(n, res)
    ! Copies the first n values that kernels returned in results to res.
    integer, intent(in) :: n
    real, intent(out), target :: res(n)
#ifdef _OPENACC
    call memcpy_from_device(c_loc(res), c_loc(results), int(n, c_size_t) * storage_size(res) / 8)
#else
    res = results(:n)
#endif
  end subroutine copy_results


  subroutine any_current(n1, n2, imm, uc, vc, res)
    ! Whether the current (uc, vc) is non-zero at any of the cells n1:n2,
    ! returned in res(1) as 1, and as 0 otherwise.
    integer, intent(in) :: n1, n2, imm
    real, intent(in) :: uc(imm), vc(imm)
    real, intent(inout) :: res(1)
    integer :: i
    !$acc serial deviceptr(res)
    res(1) = 0
    !$acc end serial
    !$acc parallel loop reduction(max:res(1)) deviceptr(res)
    do i = n1, n2
      res(1) = max(res(1), merge(1., 0., uc(i) /= 0 .or. vc(i) /= 0))
    end do
  end subroutine any_current


  subroutine sources_acc(om, pm, i1, i2, h1, h2, imm, e, f, cth, sth, cth2, sth2, wspd, &
                         wdir, uc, vc, rhorat, fice, ustar, logl2overz, cp0, cg0, fkovg, k, kdk, &
                         k3dk, k4, cothkd, snl_fac, snl_bins, snl_frac, sbf, sdv, g, fcut_fac, &
                         fprog, sin_fac, fieldscale1, fieldscale2, dth, mss_fac, twopisds_fac, &
                         sds_power, sdt_fac, fice_lth, explim, oc, shelt, moments, sds, sdt, &
                         sice, physics_time_step, res)
    ! Evaluates the source functions, the column-wise sink terms, and the
    ! largest stable source time step at each cell, as sources, and returns
    ! the smallest of these time steps in res(1).
    integer, intent(in) :: om, pm, i1, i2, h1, h2, imm, snl_bins
    real, intent(in) :: e(om,pm,h1:h2), f(om), cth(pm), sth(pm), cth2(pm), sth2(pm)
    real, intent(in), dimension(imm) :: wspd, wdir, uc, vc, rhorat, fice
    real, intent(in) :: ustar(i1:i2)
    real, intent(in), dimension(om,i1:i2) :: logl2overz, fkovg, k, kdk, k3dk, k4, cothkd, &
                                             sbf, sdv
    real, intent(in), dimension(om,h1:h2) :: cp0, cg0
    real, intent(in) :: snl_fac, snl_frac, g, fcut_fac, fprog, sin_fac, fieldscale1, &
                        fieldscale2, dth, mss_fac, twopisds_fac, sds_power, sdt_fac, fice_lth, &
                        explim
    integer, intent(inout) :: oc(i1:i2)
    real, intent(inout) :: shelt(i1:i2), moments(om,3,i1:i2), sds(om,pm,i1:i2), &
                           sdt(om,i1:i2), sice(om,i1:i2), physics_time_step(i1:i2), res(1)

    real :: wspd_i, ustar_i, shelt_i, cos_wdir, sin_wdir, scale, fac, xmin, fcutoff
    real :: d0, dc, ds, partial, running, swh, rate
    integer :: i, o, p, j, no

    fac = max(mss_fac, 0.)
    xmin = power_threshold(sds_power)

    !$acc parallel loop gang vector_length(vl) &
    !$acc   private(wspd_i, ustar_i, shelt_i, cos_wdir, sin_wdir, scale, fcutoff, no, swh, rate)
    do i = i1, i2

      ! protection against low wind speed values
      wspd_i = max(wspd(i), min_wind_speed)
      ustar_i = ustar(i)
      cos_wdir = cos(wdir(i))
      sin_wdir = sin(wdir(i))

      ! cut-off frequency (by default 4*pierson-moskowitz peak frequency) and its
      ! bin, below the bins that the downshifting into it reads
      fcutoff = min(fcut_fac * g / wspd_i, fprog)
      no = oc(i)
      !$acc loop seq
      do o = om - max(2, snl_bins + 1), 1, -1
        no = o
        if (fcutoff > f(o)) exit
      end do
      oc(i) = no

      shelt_i = sheltering_coare35(wspd_i)
      shelt(i) = shelt_i

      ! direction sums of the spectrum, then their sums over the longer
      ! waves: the moments of the mean-squared slope (see mss_moments)
      !$acc loop vector private(d0, dc, ds)
      do o = 1, om
        d0 = 0
        dc = 0
        ds = 0
        !$acc loop seq
        do p = 1, pm
          d0 = d0 + e(o,p,i)
          dc = dc + e(o,p,i) * cth2(p)
          ds = ds + e(o,p,i) * sth2(p)
        end do
        moments(o,1,i) = 0.5 * dth * d0 * k3dk(o,i)
        moments(o,2,i) = 0.5 * dth * dc * k3dk(o,i)
        moments(o,3,i) = 0.5 * dth * ds * k3dk(o,i)
      end do
      !$acc loop vector private(running, partial)
      do j = 1, 3
        running = 0
        !$acc loop seq
        do o = 1, om
          partial = moments(o,j,i)
          moments(o,j,i) = running
          running = running + partial
        end do
      end do

      ! breaking dissipation for the range used by the integration, with the
      ! diagnostic bins that the downshifting reads
      !$acc loop vector collapse(2)
      do p = 1, pm
        do o = 1, no + snl_bins + 1
          sds(o,p,i) = sds_coefficient_bin(f(o), cth2(p), sth2(p), fac, twopisds_fac, &
                                           moments(o,1,i), moments(o,2,i), moments(o,3,i)) &
                     * power_bin(e(o,p,i) * k4(o,i), sds_power, xmin)
        end do
      end do

      ! dissipation due to turbulence
      !$acc loop vector
      do o = 1, om
        sdt(o,i) = sdt_fac * sqrt(rhorat(i)) * ustar_i * k(o,i)
      end do

      ! attenuation by sea ice (see s_ice)
      if (fice(i) > fice_lth) then
        swh = 0
        !$acc loop vector collapse(2) reduction(+:swh)
        do p = 1, pm
          do o = 1, om
            swh = swh + e(o,p,i) * kdk(o,i)
          end do
        end do
        swh = 4 * sqrt(swh * dth)
        !$acc loop vector
        do o = 1, om
          sice(o,i) = s_ice_bin(cg0(o,i), swh)
        end do
      else
        !$acc loop vector
        do o = 1, om
          sice(o,i) = 0
        end do
      end if

      ! largest stable time step, from the largest growth or decay rate
      ! in the prognostic range (see max_source_rate)
      scale = (1 - fice(i)) * twopi * rhorat(i)
      rate = 0
      !$acc loop vector collapse(2) reduction(max:rate)
      do p = 1, pm
        do o = 1, no
          rate = max(rate, abs(source_rate( &
            sin_d12_bin(wspd_i + 2.5 * ustar_i * logl2overz(o,i), &
                        cos_wdir * cth(p) + sin_wdir * sth(p), cp0(o,i), uc(i), vc(i), &
                        cth(p), sth(p), shelt_i, sin_fac, fieldscale1, fieldscale2, scale, &
                        fkovg(o,i), -huge(0.)), &
            sds(o,p,i) * cothkd(o,i) * (1 - snl_gain(kdk(o,i), kdk(o+snl_bins,i), &
                                                     kdk(o+snl_bins+1,i), snl_fac, snl_frac)), &
            sbf(o,i), sdt(o,i), sdv(o,i), sice(o,i))))
        end do
      end do
      physics_time_step(i) = explim / rate

    end do

    ! the smallest of the time steps, by a kernel of its own: a reduction
    ! into res(1) in the kernel over the cells takes more shared memory than
    ! a thread block has
    !$acc serial deviceptr(res)
    res(1) = huge(res)
    !$acc end serial
    !$acc parallel loop reduction(min:res(1)) deviceptr(res)
    do i = i1, i2
      res(1) = min(res(1), physics_time_step(i))
    end do

  end subroutine sources_acc


  subroutine integrate_acc(om, pm, i1, i2, h1, h2, imm, halo, halo_source, store, dts, e, ef, &
                           f, cth, sth, cth2, sth2, wspd, wdir, uc, vc, rhorat, fice, ustar, &
                           shelt, oc, logl2overz, cp0, fkovg, k4, cothkd, sbf, sdt, sdv, sice, &
                           oneoverk4, kdk, moments, sin_fac, fieldscale1, fieldscale2, mss_fac, &
                           twopisds_fac, sds_power, inv_sds_power, snl_fac, snl_bins, snl_frac, &
                           sds, ssin, snl)
    ! Integrates the source functions over dts, as integrate_column, at the
    ! cells in the halo of other tiles if halo is .true., and at the other
    ! cells otherwise: ef is the updated spectrum, and e becomes the mean of
    ! the old and updated spectrum. If store is .true., Sin, Sds and Snl over
    ! the whole spectrum are kept for output. The kernel runs asynchronously
    ! on integrate_queue.
    integer, intent(in) :: om, pm, i1, i2, h1, h2, imm, snl_bins
    logical, intent(in) :: halo, halo_source(i1:i2), store
    real, intent(in) :: dts
    real, intent(inout) :: e(om,pm,h1:h2), ef(om,pm,h1:h2)
    real, intent(in) :: f(om), cth(pm), sth(pm), cth2(pm), sth2(pm)
    real, intent(in), dimension(imm) :: wspd, wdir, uc, vc, rhorat, fice
    real, intent(in) :: ustar(i1:i2), shelt(i1:i2)
    integer, intent(in) :: oc(i1:i2)
    real, intent(in), dimension(om,i1:i2) :: logl2overz, fkovg, k4, cothkd, sbf, sdt, sdv, sice, &
                                             oneoverk4, kdk
    real, intent(in) :: cp0(om,h1:h2), moments(om,3,i1:i2)
    real, intent(in) :: sin_fac, fieldscale1, fieldscale2, mss_fac, twopisds_fac, sds_power, &
                        inv_sds_power, snl_fac, snl_frac
    real, intent(inout) :: sds(om,pm,i1:i2), ssin(om,pm,i1:i2), snl(om,pm,i1:i2)

    real :: wspd_i, ustar_i, shelt_i, cos_wdir, sin_wdir, scale, fac, xmin_sds, xmin_eq
    real :: ssin_bin, snl_bin, coef
    integer :: i, o, p, no, m

    fac = max(mss_fac, 0.)
    xmin_sds = power_threshold(sds_power)
    xmin_eq = power_threshold(inv_sds_power)
    m = snl_bins

    !$acc parallel loop gang vector_length(vl) async(integrate_queue) &
    !$acc   private(wspd_i, ustar_i, shelt_i, cos_wdir, sin_wdir, scale, no)
    do i = i1, i2
      if (halo_source(i) .eqv. halo) then

        no = oc(i)
        wspd_i = max(wspd(i), min_wind_speed)
        ustar_i = ustar(i)
        shelt_i = shelt(i)
        cos_wdir = cos(wdir(i))
        sin_wdir = sin(wdir(i))
        scale = (1 - fice(i)) * twopi * rhorat(i)

        ! the updated spectrum, from the old one, which the
        ! downshifting reads at the higher frequencies
        !$acc loop vector collapse(2) private(ssin_bin, snl_bin, coef)
        do p = 1, pm
          do o = 1, om
            ssin_bin = sin_d12_bin(wspd_i + 2.5 * ustar_i * logl2overz(o,i), &
                                   cos_wdir * cth(p) + sin_wdir * sth(p), cp0(o,i), uc(i), &
                                   vc(i), cth(p), sth(p), shelt_i, sin_fac, fieldscale1, &
                                   fieldscale2, scale, fkovg(o,i), merge(0., -huge(0.), o > no))
            if (o <= no) then
              ! prognostic range: exponential update; the downshifting moves
              ! no energy below the lowest bin (see snl_d12)
              snl_bin = snl_d12_bin(snl_fac * (1 - snl_frac) * kdk(o+m,i) / kdk(o,i), &
                                    snl_fac * snl_frac * kdk(o+m+1,i) / kdk(o,i), &
                                    snl_fac * (merge(1 - snl_frac, 0., o > m) &
                                               + merge(snl_frac, 0., o > m + 1)), &
                                    sds(o,p,i), e(o,p,i), sds(o+m,p,i), e(o+m,p,i), &
                                    sds(o+m+1,p,i), e(o+m+1,p,i))
              ef(o,p,i) = integrate_bin(e(o,p,i), dts, &
                                        source_rate(ssin_bin, sds(o,p,i) * cothkd(o,i), &
                                                    sbf(o,i), sdt(o,i), sdv(o,i), sice(o,i)), &
                                        snl_bin)
            else
              ! diagnostic range: equilibrium of wind input and dissipation
              snl_bin = 0
              coef = sds_coefficient_bin(f(o), cth2(p), sth2(p), fac, twopisds_fac, &
                                         moments(o,1,i), moments(o,2,i), moments(o,3,i))
              ef(o,p,i) = oneoverk4(o,i) &
                * power_bin(equilibrium_bin(ssin_bin, sdt(o,i), sdv(o,i), sice(o,i), coef, &
                                            cothkd(o,i)), inv_sds_power, xmin_eq)
              if (store .and. o > no + m + 1) &
                sds(o,p,i) = coef * power_bin(e(o,p,i) * k4(o,i), sds_power, xmin_sds)
            end if
            if (store) then
              ssin(o,p,i) = ssin_bin
              snl(o,p,i) = snl_bin
            end if
          end do
        end do

        !$acc loop vector collapse(2)
        do p = 1, pm
          do o = 1, om
            e(o,p,i) = 0.5 * (e(o,p,i) + ef(o,p,i))
          end do
        end do

      end if
    end do

  end subroutine integrate_acc


  subroutine advect_acc(om, pm, i1, i2, h1, h2, imm, currents, dts, dth, dlnf, e, ef, oc, &
                        fice, fice_uth, uc, vc, east, west, north, south, iie, iiw, iin, iis, &
                        dxn, dxs, dye, dyw, oneovar, cth, sth, cg0, kovcg, invcp0, cgdlncp_dx, &
                        cgdlncp_dy, du_dx, du_dy, dv_dx, dv_dy)
    ! Propagates and refracts the spectrum e at all cells, and with currents
    ! shifts it in frequency, conserving the wave action, and integrates
    ! into ef, as advect_column. The neighbors east, west, north, and south
    ! are indices in the tile and halo, and iie, iiw, iin, iis true indices.
    integer, intent(in) :: om, pm, i1, i2, h1, h2, imm
    logical, intent(in) :: currents
    real, intent(in) :: dts, dth, dlnf
    real, intent(in) :: e(om,pm,h1:h2)
    real, intent(inout) :: ef(om,pm,h1:h2)
    integer, intent(in) :: oc(i1:i2)
    real, intent(in) :: fice(imm), fice_uth, uc(imm), vc(imm)
    integer, intent(in), dimension(imm) :: east, west, north, south, iie, iiw, iin, iis
    real, intent(in), dimension(i1:i2) :: dxn, dxs, dye, dyw, oneovar, du_dx, du_dy, dv_dx, dv_dy
    real, intent(in) :: cth(pm,i1:i2), sth(pm,i1:i2) ! ray directions adjusted for grid curvature
    real, intent(in), dimension(om,h1:h2) :: cg0, kovcg
    real, intent(in), dimension(om,i1:i2) :: invcp0, cgdlncp_dx, cgdlncp_dy

    real :: en(max_bins)              ! spectrum of a frequency over direction, in refraction substeps
    real :: a(0:max_bins+1), r(0:max_bins+1) ! action and its rate over frequency, in frequency substeps
    real :: dt_over_area, dt_over_dth, dt_over_dlnf, outflow_r, outflow_f
    real :: uc_e, uc_w, vc_n, vc_s, feup, fedn, fwup, fwdn, fnup, fndn, fsup, fsdn, f_out
    real :: rate_minus, rate, rate_plus, a_minus, a_plus, inv_kovcg, re, rw, rn, rs
    real :: rot_minus, rot, rot_plus
    real :: ax_out, ax_in, e_x, ay_out, ay_in, e_y, efop, prev, cur, next, first
    integer :: i, o, p, n, no, nt, ne, nw, nn, ns, pl, pr, nsub_r, nsub_f
    logical :: shift

    !$acc parallel loop gang vector_length(vl) &
    !$acc   private(no, nt, ne, nw, nn, ns, dt_over_area, dt_over_dth, dt_over_dlnf, outflow_r, &
    !$acc           outflow_f, nsub_r, nsub_f, shift, uc_e, uc_w, vc_n, vc_s, feup, fedn, fwup, &
    !$acc           fwdn, fnup, fndn, fsup, fsdn, f_out)
    do i = i1, i2

      no = oc(i)

      if (fice(i) > fice_uth) then

        ! no propagation through ice
        !$acc loop vector collapse(2)
        do p = 1, pm
          do o = 1, no
            ef(o,p,i) = 0
          end do
        end do

      else

        ne = east(i)
        nw = west(i)
        nn = north(i)
        ns = south(i)
        dt_over_area = 0.25 * dts * oneovar(i)

        ! advective action fluxes by the currents through the cell edges
        uc_e = uc(i) + uc(iie(i))
        uc_w = uc(i) + uc(iiw(i))
        vc_n = vc(i) + vc(iin(i))
        vc_s = vc(i) + vc(iis(i))
        feup = (uc_e + abs(uc_e)) * dye(i)
        fedn = (uc_e - abs(uc_e)) * dye(i)
        fwup = (uc_w + abs(uc_w)) * dyw(i)
        fwdn = (uc_w - abs(uc_w)) * dyw(i)
        fnup = (vc_n + abs(vc_n)) * dxn(i)
        fndn = (vc_n - abs(vc_n)) * dxn(i)
        fsup = (vc_s + abs(vc_s)) * dxs(i)
        fsdn = (vc_s - abs(vc_s)) * dxs(i)
        f_out = feup - fwdn + fnup - fsdn

        ! substeps of the refraction and of the frequency shift
        outflow_r = 0
        !$acc loop vector reduction(max:outflow_r)
        do o = 1, no
          outflow_r = max(outflow_r, refraction_outflow(cgdlncp_dx(o,i), cgdlncp_dy(o,i), &
                                                        du_dx(i), du_dy(i), dv_dx(i), &
                                                        dv_dy(i), dth))
        end do
        nsub_r = max(1, ceiling(dts * outflow_r / dth))
        dt_over_dth = dts / (nsub_r * dth)

        nt = min(no + 1, om) ! the highest bin with a rate, the first of the tail if any
        outflow_f = 0
        if (currents) then
          !$acc loop vector reduction(max:outflow_f)
          do o = 1, nt
            outflow_f = max(outflow_f, frequency_outflow(invcp0(o,i), uc(i), vc(i), &
                                                         cgdlncp_dx(o,i), cgdlncp_dy(o,i), &
                                                         cg0(o,i), du_dx(i), du_dy(i), &
                                                         dv_dx(i), dv_dy(i)))
          end do
        end if
        shift = outflow_f > 0
        nsub_f = 1
        if (shift) nsub_f = max(1, ceiling(dts * outflow_f / dlnf))
        dt_over_dlnf = dts / (nsub_f * dlnf)

        ! propagation (see propagation), and the refraction and frequency
        ! shift in a single substep, for each spectral bin
        !$acc loop vector collapse(2) &
        !$acc   private(inv_kovcg, re, rw, rn, rs, ax_out, ax_in, e_x, ay_out, ay_in, e_y, efop, &
        !$acc           pl, pr, rot_minus, rot, rot_plus, rate_minus, rate, rate_plus, a_minus, &
        !$acc           a_plus)
        do p = 1, pm
          do o = 1, no

            ! inflows scaled by the k/cg of the neighbor over that of the cell
            inv_kovcg = 1 / kovcg(o,i)
            re = kovcg(o,ne) * inv_kovcg
            rw = kovcg(o,nw) * inv_kovcg
            rn = kovcg(o,nn) * inv_kovcg
            rs = kovcg(o,ns) * inv_kovcg

            ! out through the downstream edges in x and y, in through the upstream ones
            if (cth(p,i) >= 0) then
              ax_out = (cg0(o,i) + cg0(o,ne)) * dye(i)
              ax_in = (cg0(o,i) + cg0(o,nw)) * dyw(i) * rw
              e_x = e(o,p,nw)
            else
              ax_out = (cg0(o,i) + cg0(o,nw)) * dyw(i)
              ax_in = (cg0(o,i) + cg0(o,ne)) * dye(i) * re
              e_x = e(o,p,ne)
            end if
            if (sth(p,i) >= 0) then
              ay_out = (cg0(o,i) + cg0(o,nn)) * dxn(i)
              ay_in = (cg0(o,i) + cg0(o,ns)) * dxs(i) * rs
              e_y = e(o,p,ns)
            else
              ay_out = (cg0(o,i) + cg0(o,ns)) * dxs(i)
              ay_in = (cg0(o,i) + cg0(o,nn)) * dxn(i) * rn
              e_y = e(o,p,nn)
            end if
            efop = ef(o,p,i) - dt_over_area * propagation_flux(2 * abs(cth(p,i)), ax_out, ax_in, &
                                                               e_x, 2 * abs(sth(p,i)), ay_out, &
                                                               ay_in, e_y, e(o,p,i))
            if (currents) &
              efop = efop - dt_over_area * (f_out * e(o,p,i) + fedn * re * e(o,p,ne) &
                                            - fwup * rw * e(o,p,nw) + fndn * rn * e(o,p,nn) &
                                            - fsup * rs * e(o,p,ns))

            ! refraction in one substep (see refraction)
            if (nsub_r == 1) then
              pl = p + 1 ! anti-clockwise neighbor
              if (pl > pm) pl = 1
              pr = p - 1 ! clockwise neighbor
              if (pr < 1) pr = pm
              rot_minus = depth_rotation_rate(cgdlncp_dx(o,i), cgdlncp_dy(o,i), cth(pr,i), &
                                              sth(pr,i)) &
                + current_rotation_rate(du_dx(i), du_dy(i), dv_dx(i), dv_dy(i), cth(pr,i), sth(pr,i))
              rot = depth_rotation_rate(cgdlncp_dx(o,i), cgdlncp_dy(o,i), cth(p,i), sth(p,i)) &
                + current_rotation_rate(du_dx(i), du_dy(i), dv_dx(i), dv_dy(i), cth(p,i), sth(p,i))
              rot_plus = depth_rotation_rate(cgdlncp_dx(o,i), cgdlncp_dy(o,i), cth(pl,i), &
                                             sth(pl,i)) &
                + current_rotation_rate(du_dx(i), du_dy(i), dv_dx(i), dv_dy(i), cth(pl,i), sth(pl,i))
              efop = efop - dt_over_dth * upstream_flux(rot_minus, rot, rot_plus, &
                                                        e(o,pr,i), e(o,p,i), e(o,pl,i))
            end if

            ! frequency shift in one substep (see frequency_shift): past the
            ! ends, the rate is that of the end bins, and there is no action
            ! below the lowest bin, nor above the highest without a tail
            if (shift .and. nsub_f == 1) then
              rate = frequency_rate(invcp0(o,i), uc(i), vc(i), cgdlncp_dx(o,i), &
                                    cgdlncp_dy(o,i), cg0(o,i), du_dx(i), du_dy(i), dv_dx(i), &
                                    dv_dy(i), cth(p,i), sth(p,i))
              rate_minus = rate
              a_minus = 0
              if (o > 1) then
                rate_minus = frequency_rate(invcp0(o-1,i), uc(i), vc(i), cgdlncp_dx(o-1,i), &
                                            cgdlncp_dy(o-1,i), cg0(o-1,i), du_dx(i), du_dy(i), &
                                            dv_dx(i), dv_dy(i), cth(p,i), sth(p,i))
                a_minus = kovcg(o-1,i) * e(o-1,p,i)
              end if
              rate_plus = rate
              a_plus = 0
              if (o < nt) then
                rate_plus = frequency_rate(invcp0(o+1,i), uc(i), vc(i), cgdlncp_dx(o+1,i), &
                                           cgdlncp_dy(o+1,i), cg0(o+1,i), du_dx(i), du_dy(i), &
                                           dv_dx(i), dv_dy(i), cth(p,i), sth(p,i))
                a_plus = kovcg(o+1,i) * e(o+1,p,i)
              end if
              efop = efop - dt_over_dlnf * inv_kovcg &
                * upstream_flux(rate_minus, rate, rate_plus, a_minus, kovcg(o,i) * e(o,p,i), &
                                a_plus)
            end if

            ef(o,p,i) = efop

          end do
        end do

        ! refraction in more than one substep, by the threads over frequency,
        ! updating the spectrum over direction in place
        if (nsub_r > 1) then
          !$acc loop vector private(en, prev, cur, next, first, pl, pr, rot_minus, rot, rot_plus)
          do o = 1, no
            !$acc loop seq
            do p = 1, pm
              en(p) = e(o,p,i)
            end do
            !$acc loop seq
            do n = 1, nsub_r
              first = en(1)
              prev = en(pm)
              !$acc loop seq
              do p = 1, pm
                pl = p + 1
                if (pl > pm) pl = 1
                pr = p - 1
                if (pr < 1) pr = pm
                rot_minus = depth_rotation_rate(cgdlncp_dx(o,i), cgdlncp_dy(o,i), cth(pr,i), &
                                                sth(pr,i)) &
                  + current_rotation_rate(du_dx(i), du_dy(i), dv_dx(i), dv_dy(i), cth(pr,i), &
                                          sth(pr,i))
                rot = depth_rotation_rate(cgdlncp_dx(o,i), cgdlncp_dy(o,i), cth(p,i), sth(p,i)) &
                  + current_rotation_rate(du_dx(i), du_dy(i), dv_dx(i), dv_dy(i), cth(p,i), sth(p,i))
                rot_plus = depth_rotation_rate(cgdlncp_dx(o,i), cgdlncp_dy(o,i), cth(pl,i), &
                                               sth(pl,i)) &
                  + current_rotation_rate(du_dx(i), du_dy(i), dv_dx(i), dv_dy(i), cth(pl,i), &
                                          sth(pl,i))
                cur = en(p)
                next = first
                if (p < pm) next = en(p+1)
                en(p) = cur - dt_over_dth * upstream_flux(rot_minus, rot, rot_plus, prev, cur, next)
                prev = cur
              end do
            end do
            !$acc loop seq
            do p = 1, pm
              ef(o,p,i) = ef(o,p,i) + en(p) - e(o,p,i)
            end do
          end do
        end if

        ! frequency shift in more than one substep, by the threads over
        ! direction, updating the action over frequency in place
        if (shift .and. nsub_f > 1) then
          !$acc loop vector private(a, r, prev, cur)
          do p = 1, pm
            !$acc loop seq
            do o = 1, nt
              r(o) = frequency_rate(invcp0(o,i), uc(i), vc(i), cgdlncp_dx(o,i), &
                                    cgdlncp_dy(o,i), cg0(o,i), du_dx(i), du_dy(i), dv_dx(i), &
                                    dv_dy(i), cth(p,i), sth(p,i))
              a(o) = kovcg(o,i) * e(o,p,i)
            end do
            r(0) = r(1)
            a(0) = 0
            if (nt == no) then
              r(no+1) = r(no)
              a(no+1) = 0
            end if
            !$acc loop seq
            do n = 1, nsub_f
              prev = a(0)
              !$acc loop seq
              do o = 1, no
                cur = a(o)
                a(o) = cur - dt_over_dlnf * upstream_flux(r(o-1), r(o), r(o+1), prev, cur, a(o+1))
                prev = cur
              end do
            end do
            !$acc loop seq
            do o = 1, no
              ef(o,p,i) = ef(o,p,i) + a(o) / kovcg(o,i) - e(o,p,i)
            end do
          end do
        end if

      end if
    end do

  end subroutine advect_acc


  subroutine stress_acc(om, pm, i1, i2, h1, h2, imm, ef, cth, sth, wspd, wdir, uc, vc, rhoa, &
                        rhow, rhorat, fice, shelt, oc, logl2overz, cp0, fkovg, invcp0, k, kdk, &
                        us, vs, sin_fac, fieldscale1, fieldscale2, dthg, z, nu_air, kappa, &
                        taux_form, tauy_form, tailatmx, tailatmy, taux_diag, tauy_diag, &
                        taux_skin, tauy_skin, taux, tauy, cd, ustar)
    ! Evaluates form and skin drag, drag coefficient, and friction velocity
    ! at all cells from the updated spectrum ef, as wind_stress, with the
    ! surface Stokes drift (us, vs).
    integer, intent(in) :: om, pm, i1, i2, h1, h2, imm
    real, intent(in) :: ef(om,pm,h1:h2), cth(pm), sth(pm)
    real, intent(in), dimension(imm) :: wspd, wdir, uc, vc, rhoa, rhow, rhorat, fice
    real, intent(in) :: shelt(i1:i2)
    integer, intent(in) :: oc(i1:i2)
    real, intent(in), dimension(om,i1:i2) :: logl2overz, fkovg, invcp0, k, kdk
    real, intent(in) :: cp0(om,h1:h2), us(i1:i2), vs(i1:i2)
    real, intent(in) :: sin_fac, fieldscale1, fieldscale2, dthg, z, nu_air, kappa
    real, intent(inout), dimension(i1:i2) :: taux_form, tauy_form, tailatmx, tailatmy, &
                                             taux_diag, tauy_diag, taux_skin, tauy_skin, &
                                             taux, tauy, cd, ustar

    real :: wspd_i, ustar_i, shelt_i, cos_wdir, sin_wdir, scale, flux
    real :: sx, sy, sx_diag, sy_diag, sx_tail, sy_tail
    real :: tail, rhow_dthg, wspdrel, wdirrel, cd_form, cd_skin
    integer :: i, o, p, no

    ! form drag, by a thread block for each cell
    !$acc parallel loop gang vector_length(vl) &
    !$acc   private(no, wspd_i, ustar_i, shelt_i, cos_wdir, sin_wdir, scale, sx, sy, sx_diag, &
    !$acc           sy_diag, sx_tail, sy_tail, tail, rhow_dthg)
    do i = i1, i2

      no = oc(i)
      wspd_i = max(wspd(i), min_wind_speed)
      ustar_i = ustar(i)
      shelt_i = shelt(i)
      cos_wdir = cos(wdir(i))
      sin_wdir = sin(wdir(i))
      scale = (1 - fice(i)) * twopi * rhorat(i)

      ! form drag from atmosphere (positive into waves): the wind input
      ! over direction projected on x and y, over the phase speed, over
      ! all frequencies, over the diagnostic range, and at the highest
      sx = 0
      sy = 0
      sx_diag = 0
      sy_diag = 0
      sx_tail = 0
      sy_tail = 0
      !$acc loop vector collapse(2) private(flux) &
      !$acc   reduction(+:sx, sy, sx_diag, sy_diag, sx_tail, sy_tail)
      do p = 1, pm
        do o = 1, om
          flux = ef(o,p,i) * invcp0(o,i) &
            * sin_d12_bin(wspd_i + 2.5 * ustar_i * logl2overz(o,i), &
                          cos_wdir * cth(p) + sin_wdir * sth(p), cp0(o,i), uc(i), vc(i), &
                          cth(p), sth(p), shelt_i, sin_fac, fieldscale1, fieldscale2, scale, &
                          fkovg(o,i), merge(0., -huge(0.), o > no))
          sx = sx + flux * cth(p) * kdk(o,i)
          sy = sy + flux * sth(p) * kdk(o,i)
          if (o >= no) then
            sx_diag = sx_diag + flux * cth(p) * kdk(o,i)
            sy_diag = sy_diag + flux * sth(p) * kdk(o,i)
          end if
          if (o == om) then
            sx_tail = sx_tail + flux * cth(p)
            sy_tail = sy_tail + flux * sth(p)
          end if
        end do
      end do

      ! wind speed dependent tail
      tail = stress_tail(wspd(i), k(om,i))
      rhow_dthg = rhow(i) * dthg
      tailatmx(i) = sx_tail * k(om,i) * tail * rhow_dthg
      tailatmy(i) = sy_tail * k(om,i) * tail * rhow_dthg
      taux_form(i) = sx * rhow_dthg + tailatmx(i)
      tauy_form(i) = sy * rhow_dthg + tailatmy(i)
      taux_diag(i) = sx_diag * rhow_dthg
      tauy_diag(i) = sy_diag * rhow_dthg

    end do

    ! skin drag, drag coefficient, and friction velocity, by a thread for
    ! each cell, whose variables, unlike those of a thread block, the
    ! subroutine wind_relative can return into
    !$acc parallel loop private(wspdrel, wdirrel, cd_form, cd_skin)
    do i = i1, i2

      ! wind speed and direction relative to surface velocity
      call wind_relative(wspd(i), wdir(i), uc(i) + us(i), vc(i) + vs(i), wspdrel, wdirrel)

      ! form- and skin-induced drag coefficients
      cd_form = drag_coefficient(taux_form(i), tauy_form(i), rhoa(i), wspd(i))
      cd_skin = drag_coefficient_skin(cd_form, wspd(i), wspdrel, z, nu_air, kappa)

      taux_skin(i) = rhoa(i) * cd_skin * wspdrel**2 * cos(wdirrel)
      tauy_skin(i) = rhoa(i) * cd_skin * wspdrel**2 * sin(wdirrel)

      taux(i) = taux_form(i) + taux_skin(i)
      tauy(i) = tauy_form(i) + tauy_skin(i)

      ! total (form + skin) drag coefficient
      cd(i) = drag_coefficient(taux(i), tauy(i), rhoa(i), wspd(i))

      ! friction velocity based on total (form + skin) stress
      ustar(i) = sqrt(sqrt(taux(i)**2 + tauy(i)**2) / rhoa(i))

    end do

  end subroutine stress_acc


  module subroutine diagnostics_acc(self)
    ! Integrated quantities of the wave spectrum for output, as diagnostics.
    class(domain_type), intent(inout) :: self
    associate(om => self % spectrum % num_frequencies, pm => self % spectrum % num_directions, &
              forcing => self % forcing)
    call diagnostics_kernel(om, pm, self % grid % istart, self % grid % iend, &
                            self % grid % iistart - 1, self % grid % iiend, self % grid % imm, &
                            self % e, self % f, self % th, self % cth, self % sth, self % oc, &
                            self % k, self % kdk, self % k3dk, self % invcp0, self % cp0, &
                            self % cg0, forcing % rhow, forcing % uc, forcing % vc, self % dth, &
                            self % dthg, self % mwp, self % momx, self % momy, self % cgmxx, &
                            self % cgmxy, self % cgmyy, self % ht, self % mwd, self % mss, &
                            self % mwl, self % dwd, self % dwl, self % dwp, self % dcp0, &
                            self % dcg0, self % dcp, self % dcg)
    end associate
  end subroutine diagnostics_acc


  subroutine diagnostics_kernel(om, pm, i1, i2, h1, h2, imm, e, f, th, cth, sth, oc, k, kdk, &
                                k3dk, invcp0, cp0, cg0, rhow, uc, vc, dth, dthg, mwp, momx, &
                                momy, cgmxx, cgmxy, cgmyy, ht, mwd, mss, mwl, dwd, dwl, dwp, &
                                dcp0, dcg0, dcp, dcg)
    ! The diagnostics at all cells, by a thread block for each cell.
    integer, intent(in) :: om, pm, i1, i2, h1, h2, imm
    real, intent(in) :: e(om,pm,h1:h2), f(om), th(pm), cth(pm), sth(pm)
    integer, intent(in) :: oc(i1:i2)
    real, intent(in), dimension(om,i1:i2) :: k, kdk, k3dk, invcp0
    real, intent(in), dimension(om,h1:h2) :: cp0, cg0
    real, intent(in), dimension(imm) :: rhow, uc, vc
    real, intent(in) :: dth, dthg
    real, intent(inout), dimension(i1:i2) :: mwp, momx, momy, cgmxx, cgmxy, cgmyy, ht, mwd, &
                                             mss, mwl, dwd, dwl, dwp, dcp0, dcg0, dcp, dcg

    real :: bin, ekdkovcp, m0, m2, mx, my, cxx, cxy, cyy, xcomp, ycomp, mk, peak, rhow_dthg
    integer :: i, o, p, no, first, opeak, ppeak

    !$acc parallel loop gang vector_length(vl) &
    !$acc   private(no, m0, m2, mx, my, cxx, cxy, cyy, xcomp, ycomp, mk, peak, first, opeak, &
    !$acc           ppeak, rhow_dthg)
    do i = i1, i2

      no = oc(i)
      m0 = 0
      m2 = 0
      mx = 0
      my = 0
      cxx = 0
      cxy = 0
      cyy = 0
      xcomp = 0
      ycomp = 0
      mk = 0
      peak = -huge(peak)
      !$acc loop vector collapse(2) private(bin, ekdkovcp) &
      !$acc   reduction(+:m0, m2, mx, my, cxx, cxy, cyy, xcomp, ycomp, mk) reduction(max:peak)
      do p = 1, pm
        do o = 1, om
          bin = e(o,p,i) * kdk(o,i)
          m0 = m0 + bin
          m2 = m2 + f(o)**2 * bin
          xcomp = xcomp + bin * cth(p)
          ycomp = ycomp + bin * sth(p)
          mk = mk + e(o,p,i) * k3dk(o,i)
          peak = max(peak, bin)
          ! total wave momentum and its flux
          if (o <= no) then
            ekdkovcp = bin * invcp0(o,i)
            mx = mx + ekdkovcp * cth(p)
            my = my + ekdkovcp * sth(p)
            cxx = cxx + cg0(o,i) * ekdkovcp * cth(p)**2
            cxy = cxy + cg0(o,i) * ekdkovcp * cth(p) * sth(p)
            cyy = cyy + cg0(o,i) * ekdkovcp * sth(p)**2
          end if
        end do
      end do

      ! the first bin of the peak in array element order, as maxloc
      first = om * pm
      !$acc loop vector collapse(2) reduction(min:first)
      do p = 1, pm
        do o = 1, om
          if (e(o,p,i) * kdk(o,i) == peak) first = min(first, o + om * (p - 1))
        end do
      end do
      opeak = mod(first - 1, om) + 1
      ppeak = (first - 1) / om + 1

      rhow_dthg = rhow(i) * dthg
      mwp(i) = sqrt(m0 / (m2 + tiny(m2)))
      momx(i) = mx * rhow_dthg
      momy(i) = my * rhow_dthg
      cgmxx(i) = cxx * rhow_dthg
      cgmxy(i) = cxy * rhow_dthg
      cgmyy(i) = cyy * rhow_dthg
      ht(i) = 4 * sqrt(m0 * dth)
      mwd(i) = atan2(ycomp, xcomp)
      mss(i) = mk * dth
      mwl(i) = twopi * sqrt(m0 / (mk + tiny(mk)))

      ! dominant (peak) wave quantities
      dwd(i) = th(ppeak)
      dwl(i) = twopi / k(opeak,i)
      dwp(i) = 1. / f(opeak)
      dcp0(i) = cp0(opeak,i)
      dcg0(i) = cg0(opeak,i)
      dcp(i) = cp0(opeak,i) + uc(i) * cth(ppeak) + vc(i) * sth(ppeak)
      dcg(i) = cg0(opeak,i) + uc(i) * cth(ppeak) + vc(i) * sth(ppeak)

    end do

  end subroutine diagnostics_kernel


  module subroutine ocean_stress_acc(self)
    ! Momentum fluxes into the ocean top and bottom, energy fluxes from air
    ! and into the ocean, and the components of form drag, as ocean_stress.
    class(domain_type), intent(inout) :: self
    associate(om => self % spectrum % num_frequencies, pm => self % spectrum % num_directions, &
              forcing => self % forcing)
    call ocean_stress_kernel(om, pm, self % grid % istart, self % grid % iend, &
                             self % grid % iistart - 1, self % grid % iiend, self % grid % imm, &
                             self % e, self % ssin, self % sds, self % snl, self % th, &
                             self % cth, self % sth, self % oc, self % k, self % kdk, &
                             self % invcp0, self % cothkd, self % sdt, self % sdv, self % sbf, &
                             self % cp0, forcing % wspd, forcing % wdir, forcing % rhow, &
                             self % snl_bins, self % snl_frac, self % dthg, self % taux_skin, &
                             self % tauy_skin, self % tailocnx, self % tailocny, &
                             self % taux_ocntop, self % tauy_ocntop, self % taux_ocnbot, &
                             self % tauy_ocnbot, self % taux_snl, self % tauy_snl, &
                             self % epsx_ocn, self % epsy_ocn, self % epsx_atm, self % epsy_atm, &
                             self % taux1, self % tauy1, self % taux2, self % tauy2, &
                             self % taux3, self % tauy3)
    end associate
  end subroutine ocean_stress_acc


  subroutine ocean_stress_kernel(om, pm, i1, i2, h1, h2, imm, e, ssin, sds, snl, th, cth, sth, &
                                 oc, k, kdk, invcp0, cothkd, sdt, sdv, sbf, cp0, wspd, wdir, &
                                 rhow, m, frac, dthg, taux_skin, tauy_skin, tailocnx, tailocny, &
                                 taux_ocntop, tauy_ocntop, taux_ocnbot, tauy_ocnbot, taux_snl, &
                                 tauy_snl, epsx_ocn, epsy_ocn, epsx_atm, epsy_atm, taux1, tauy1, &
                                 taux2, tauy2, taux3, tauy3)
    ! The ocean stress at all cells, by a thread block for each cell. The
    ! sums over frequency and direction of the quantities that ocean_stress
    ! sums over direction and then frequency.
    integer, intent(in) :: om, pm, i1, i2, h1, h2, imm, m
    real, intent(in) :: e(om,pm,h1:h2)
    real, intent(in), dimension(om,pm,i1:i2) :: ssin, sds, snl
    real, intent(in) :: th(pm), cth(pm), sth(pm)
    integer, intent(in) :: oc(i1:i2)
    real, intent(in), dimension(om,i1:i2) :: k, kdk, invcp0, cothkd, sdt, sdv, sbf
    real, intent(in) :: cp0(om,h1:h2)
    real, intent(in), dimension(imm) :: wspd, wdir, rhow
    real, intent(in) :: frac, dthg
    real, intent(in), dimension(i1:i2) :: taux_skin, tauy_skin
    real, intent(inout), dimension(i1:i2) :: tailocnx, tailocny, taux_ocntop, tauy_ocntop, &
                                             taux_ocnbot, tauy_ocnbot, taux_snl, tauy_snl, &
                                             epsx_ocn, epsy_ocn, epsx_atm, epsy_atm, &
                                             taux1, tauy1, taux2, tauy2, taux3, tauy3

    real :: snl_cp, top, form, tail, rhow_dthg
    real :: topx, topy, tailx, taily, botx, boty, snlx, snly, ocnx, ocny, atmx, atmy
    real :: t1x, t1y, t2x, t2y, t3x, t3y
    integer :: i, o, p, no

    !$acc parallel loop gang vector_length(vl) &
    !$acc   private(no, topx, topy, tailx, taily, botx, boty, snlx, snly, ocnx, ocny, atmx, &
    !$acc           atmy, t1x, t1y, t2x, t2y, t3x, t3y, tail, rhow_dthg)
    do i = i1, i2

      no = oc(i)
      topx = 0
      topy = 0
      tailx = 0
      taily = 0
      botx = 0
      boty = 0
      snlx = 0
      snly = 0
      ocnx = 0
      ocny = 0
      atmx = 0
      atmy = 0
      t1x = 0
      t1y = 0
      t2x = 0
      t2y = 0
      t3x = 0
      t3y = 0
      !$acc loop vector collapse(2) private(snl_cp, top, form) &
      !$acc   reduction(+:topx, topy, tailx, taily, botx, boty, snlx, snly, ocnx, ocny, atmx, &
      !$acc               atmy, t1x, t1y, t2x, t2y, t3x, t3y)
      do p = 1, pm
        do o = 1, om

          ! Snl in the diagnostic part, at the phase speed of the bins it moves to
          snl_cp = 0
          if (o >= m + 2) snl_cp = snl(o,p,i) * ((1 - frac) * cp0(o,i) * invcp0(o-m,i) &
                                                + frac * cp0(o,i) * invcp0(o-m-1,i))

          ! stress into ocean top (positive into ocean): dissipation into
          ! currents, corrected for the Snl in the diagnostic part
          top = (e(o,p,i) * (sds(o,p,i) * cothkd(o,i) + sdt(o,i) + sdv(o,i)) - snl_cp) &
              * invcp0(o,i)
          topx = topx + top * cth(p) * kdk(o,i)
          topy = topy + top * sth(p) * kdk(o,i)
          if (o == om) then
            tailx = tailx + top * cth(p)
            taily = taily + top * sth(p)
          end if

          ! stress into ocean bottom (positive into ocean)
          botx = botx + e(o,p,i) * sbf(o,i) * cth(p) * invcp0(o,i) * kdk(o,i)
          boty = boty + e(o,p,i) * sbf(o,i) * sth(p) * invcp0(o,i) * kdk(o,i)

          ! momentum loss due to the non-linear downshifting of energy
          if (o >= m + 2 .and. o <= no) then
            form = snl(o,p,i) * ((1 - frac) * (1 - cp0(o,i) * invcp0(o-m,i)) &
                               + frac * (1 - cp0(o,i) * invcp0(o-m-1,i))) * invcp0(o,i) * kdk(o,i)
            snlx = snlx + form * cth(p)
            snly = snly + form * sth(p)
          end if

          ! energy flux into ocean and from air
          ocnx = ocnx + (e(o,p,i) * (sds(o,p,i) * cothkd(o,i)) - snl_cp) * cth(p) * kdk(o,i)
          ocny = ocny + (e(o,p,i) * (sds(o,p,i) * cothkd(o,i)) - snl_cp) * sth(p) * kdk(o,i)
          atmx = atmx + e(o,p,i) * ssin(o,p,i) * cth(p) * kdk(o,i)
          atmy = atmy + e(o,p,i) * ssin(o,p,i) * sth(p) * kdk(o,i)

          ! components of form drag: wind pushing waves, waves against
          ! wind, and waves overrunning wind
          form = e(o,p,i) * ssin(o,p,i) * invcp0(o,i) * kdk(o,i)
          if (ssin(o,p,i) > 0) then
            t1x = t1x + form * cth(p)
            t1y = t1y + form * sth(p)
          else if (cos(wdir(i) - th(p)) < 0) then
            t2x = t2x + form * cth(p)
            t2y = t2y + form * sth(p)
          else
            t3x = t3x + form * cth(p)
            t3y = t3y + form * sth(p)
          end if

        end do
      end do

      tail = stress_tail(wspd(i), k(om,i))
      rhow_dthg = rhow(i) * dthg
      tailocnx(i) = tailx * k(om,i) * tail * rhow_dthg
      tailocny(i) = taily * k(om,i) * tail * rhow_dthg
      taux_ocntop(i) = topx * rhow_dthg + tailocnx(i) + taux_skin(i)
      tauy_ocntop(i) = topy * rhow_dthg + tailocny(i) + tauy_skin(i)
      taux_ocnbot(i) = botx * rhow_dthg
      tauy_ocnbot(i) = boty * rhow_dthg
      taux_snl(i) = snlx * rhow_dthg
      tauy_snl(i) = snly * rhow_dthg
      epsx_ocn(i) = ocnx * rhow_dthg
      epsy_ocn(i) = ocny * rhow_dthg
      epsx_atm(i) = atmx * rhow_dthg
      epsy_atm(i) = atmy * rhow_dthg
      taux1(i) = t1x * rhow_dthg
      tauy1(i) = t1y * rhow_dthg
      taux2(i) = t2x * rhow_dthg
      tauy2(i) = t2y * rhow_dthg
      taux3(i) = t3x * rhow_dthg
      tauy3(i) = t3y * rhow_dthg

    end do

  end subroutine ocean_stress_kernel


  module function screen_values_acc(self, i) result(res)
    ! Wind speed and direction, significant wave height, mean wave period,
    ! drag coefficient times 10^3, and the cut-off frequency at cell i, for
    ! the screen output of each time step, evaluated on the GPU: the pages
    ! of the fields that the CPU reads migrate to it and back.
    class(domain_type), intent(in) :: self
    integer, intent(in) :: i
    real :: res(6)
    call allocate_results()
    call screen_values_kernel(self % spectrum % num_frequencies, &
                              self % spectrum % num_directions, i, self % grid % istart, &
                              self % grid % iend, self % grid % iistart - 1, self % grid % iiend, &
                              self % grid % imm, self % e, self % f, self % kdk, &
                              self % forcing % wspd, self % forcing % wdir, self % cd, self % oc, &
                              self % dth, results)
    call copy_results(size(res), res)
  end function screen_values_acc


  subroutine screen_values_kernel(om, pm, i, i1, i2, h1, h2, imm, e, f, kdk, wspd, wdir, cd, oc, &
                                  dth, res)
    ! The values of screen_values_acc, returned in res, by a thread block
    ! with threads over direction, summing over frequency first as
    ! sigwaveheight and meanwaveperiod.
    integer, intent(in) :: om, pm, i, i1, i2, h1, h2, imm
    real, intent(in) :: e(om,pm,h1:h2), f(om), kdk(om,i1:i2), wspd(imm), wdir(imm), cd(i1:i2), dth
    integer, intent(in) :: oc(i1:i2)
    real, intent(inout) :: res(6)
    real :: m0, m2, s0, s2
    integer :: o, p
    !$acc parallel num_gangs(1) vector_length(vl) deviceptr(res)
    m0 = 0
    m2 = 0
    !$acc loop vector reduction(+:m0, m2) private(s0, s2)
    do p = 1, pm
      s0 = 0
      s2 = 0
      !$acc loop seq
      do o = 1, om
        s0 = s0 + e(o,p,i) * kdk(o,i)
        s2 = s2 + f(o)**2 * e(o,p,i) * kdk(o,i)
      end do
      m0 = m0 + s0
      m2 = m2 + s2
    end do
    res(1) = wspd(i)
    res(2) = wdir(i)
    res(3) = 4 * sqrt(m0 * dth)
    res(4) = sqrt(m0 / m2)
    res(5) = cd(i) * 1e3
    res(6) = f(oc(i))
    !$acc end parallel
  end subroutine screen_values_kernel

end submodule umwm_domain_acc
