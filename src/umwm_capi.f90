module umwm_capi

  ! C-interoperable interface to the wave model for foreign-language callers,
  ! such as the Python package. A domain is referred to by an opaque handle,
  ! and its fields are exposed as pointers to the memory owned by the domain,
  ! with their type, shape, and lower bounds, so that callers can view them
  ! in place without copying.
  !
  ! The lifecycle of a domain is as in the main program: umwm_env_init,
  ! umwm_domain_new, umwm_domain_initialize, umwm_domain_run (or repeated
  ! umwm_domain_advance), umwm_domain_finalize, umwm_domain_delete, and
  ! umwm_env_stop. Under MPI, initialization, time stepping, and gathering
  ! are collective and must be called by all processes.

  use ieee_arithmetic, only: ieee_get_underflow_mode, ieee_set_underflow_mode, &
                             ieee_support_underflow_control
  use iso_c_binding, only: c_char, c_double, c_f_pointer, c_float, c_int, c_loc, &
                           c_null_char, c_null_ptr, c_ptr
  use iso_fortran_env, only: error_unit, output_unit
#ifdef MPI
  use mpi
#endif
  use umwm_domain, only: domain_type
  use umwm_env, only: env_init, env_stop, mpisize, nproc
  use umwm_io, only: gather_field
  use umwm_restart, only: restart_write

  implicit none

  private

  ! element type codes of umwm_domain_field
  integer(c_int), parameter :: type_real = 1, type_integer = 2, type_logical = 3

  ! status codes
  integer(c_int), parameter :: ok = 0, unknown_name = 1, not_allocated = 2

contains

  subroutine umwm_env_init() bind(c, name='umwm_env_init')
    ! Initializes the parallel environment (MPI if enabled), and flushes
    ! denormal results to zero when the library was compiled with fast math
    ! (or by an Intel compiler, where -ftz is on by default): the compiler
    ! enables that at the startup of an executable but not of a shared
    ! library, and the wave spectrum underflows in much of its range, where
    ! denormal arithmetic is very slow.
    call env_init()
#if defined(__FAST_MATH__) || defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
    call umwm_env_flush_denormals(1_c_int)
#endif
  end subroutine umwm_env_init


  subroutine umwm_env_flush_denormals(enable) bind(c, name='umwm_env_flush_denormals')
    ! Enables (nonzero) or disables (zero) the flushing of denormal
    ! floating-point results to zero in this process.
    integer(c_int), value :: enable
    ! the argument is positional because nvfortran rejects the keyword gradual=
    if (ieee_support_underflow_control(1.0)) &
      call ieee_set_underflow_mode(enable == 0)
  end subroutine umwm_env_flush_denormals


  function umwm_env_denormals_flushed() bind(c, name='umwm_env_denormals_flushed') result(flushed)
    ! 1 if denormal results are flushed to zero in this process, 0 otherwise.
    integer(c_int) :: flushed
    logical :: gradual
    gradual = .true.
    if (ieee_support_underflow_control(1.0)) call ieee_get_underflow_mode(gradual)
    flushed = merge(0, 1, gradual)
  end function umwm_env_denormals_flushed


  subroutine umwm_env_stop() bind(c, name='umwm_env_stop')
    ! Finalizes the parallel environment.
    call env_stop()
  end subroutine umwm_env_stop


  subroutine umwm_env_flush() bind(c, name='umwm_env_flush')
    ! Flushes the standard output and error units of the Fortran runtime,
    ! for callers that redirect them.
    flush(output_unit)
    flush(error_unit)
  end subroutine umwm_env_flush


  subroutine umwm_env_barrier() bind(c, name='umwm_env_barrier')
    ! Waits for all processes; does nothing in serial mode.
#ifdef MPI
    integer :: ierr
    call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
  end subroutine umwm_env_barrier


  function umwm_env_rank() bind(c, name='umwm_env_rank') result(rank)
    ! Rank of this process (0 in serial mode).
    integer(c_int) :: rank
    rank = nproc
  end function umwm_env_rank


  function umwm_env_size() bind(c, name='umwm_env_size') result(size)
    ! Number of processes (1 in serial mode).
    integer(c_int) :: size
    size = mpisize
  end function umwm_env_size


  function umwm_env_mpi() bind(c, name='umwm_env_mpi') result(mpi)
    ! 1 if the library was built with MPI, 0 otherwise.
    integer(c_int) :: mpi
#ifdef MPI
    mpi = 1
#else
    mpi = 0
#endif
  end function umwm_env_mpi


  function umwm_domain_new() bind(c, name='umwm_domain_new') result(handle)
    ! Creates an uninitialized domain and returns its handle.
    type(c_ptr) :: handle
    type(domain_type), pointer :: domain
    allocate(domain)
    handle = c_loc(domain)
  end function umwm_domain_new


  function umwm_domain_initialize(handle, path) bind(c, name='umwm_domain_initialize') result(stat)
    ! Initializes the domain from the namelist file at path (null-terminated),
    ! relative to the current working directory, like the main program.
    ! Returns 0 on success and 1 if the configuration is invalid.
    type(c_ptr), value :: handle
    character(kind=c_char), intent(in) :: path(*)
    integer(c_int) :: stat
    type(domain_type), pointer :: domain
    integer :: istat
    call c_f_pointer(handle, domain)
    call domain % initialize(fortran_string(path), istat)
    stat = istat
  end function umwm_domain_initialize


  subroutine umwm_domain_advance(handle) bind(c, name='umwm_domain_advance')
    ! Integrates the domain over one global time step.
    type(c_ptr), value :: handle
    type(domain_type), pointer :: domain
    call c_f_pointer(handle, domain)
    call domain % advance()
  end subroutine umwm_domain_advance


  subroutine umwm_domain_run(handle) bind(c, name='umwm_domain_run')
    ! Integrates the domain from the current time to the stop time.
    type(c_ptr), value :: handle
    type(domain_type), pointer :: domain
    call c_f_pointer(handle, domain)
    call domain % run()
  end subroutine umwm_domain_run


  function umwm_domain_finished(handle) bind(c, name='umwm_domain_finished') result(finished)
    ! 1 if the domain has reached the stop time, 0 otherwise.
    type(c_ptr), value :: handle
    integer(c_int) :: finished
    type(domain_type), pointer :: domain
    call c_f_pointer(handle, domain)
    finished = merge(1, 0, domain % finished())
  end function umwm_domain_finished


  subroutine umwm_domain_time(handle, buffer) bind(c, name='umwm_domain_time')
    ! Writes the current model time as 'YYYY-MM-DD HH:MM:SS' followed by a
    ! null character into buffer, which must hold at least 20 characters.
    type(c_ptr), value :: handle
    character(kind=c_char), intent(inout) :: buffer(*)
    type(domain_type), pointer :: domain
    character(19) :: timestr
    integer :: i
    call c_f_pointer(handle, domain)
    timestr = domain % current_time % strftime('%Y-%m-%d %H:%M:%S')
    do i = 1, 19
      buffer(i) = timestr(i:i)
    end do
    buffer(20) = c_null_char
  end subroutine umwm_domain_time


  subroutine umwm_domain_diagnostics(handle) bind(c, name='umwm_domain_diagnostics')
    ! Updates the integrated wave diagnostics from the current spectrum.
    type(c_ptr), value :: handle
    type(domain_type), pointer :: domain
    call c_f_pointer(handle, domain)
    call domain % diagnostics()
  end subroutine umwm_domain_diagnostics


  subroutine umwm_domain_ocean_stress(handle) bind(c, name='umwm_domain_ocean_stress')
    ! Updates the ocean-side momentum and energy fluxes and the form drag
    ! components, which the model otherwise evaluates at gridded output times.
    type(c_ptr), value :: handle
    type(domain_type), pointer :: domain
    call c_f_pointer(handle, domain)
    call domain % ocean_stress()
  end subroutine umwm_domain_ocean_stress


  subroutine umwm_domain_write_output(handle) bind(c, name='umwm_domain_write_output')
    ! Writes the gridded output at the current model time, after updating
    ! the ocean stress as the model does before gridded output. Collective.
    type(c_ptr), value :: handle
    type(domain_type), pointer :: domain
    call c_f_pointer(handle, domain)
    call domain % ocean_stress()
    call domain % output_grid_nc(time_string(domain))
  end subroutine umwm_domain_write_output


  subroutine umwm_domain_write_spectrum(handle) bind(c, name='umwm_domain_write_spectrum')
    ! Writes the spectrum output at the configured points at the current
    ! model time, as a record appended to the spectrum output files.
    type(c_ptr), value :: handle
    type(domain_type), pointer :: domain
    call c_f_pointer(handle, domain)
    call domain % output_spectrum_nc(time_string(domain))
  end subroutine umwm_domain_write_spectrum


  subroutine umwm_domain_write_restart(handle) bind(c, name='umwm_domain_write_restart')
    ! Writes a restart file at the current model time. Collective.
    type(c_ptr), value :: handle
    type(domain_type), pointer :: domain
    call c_f_pointer(handle, domain)
    associate(grid => domain % grid)
    call restart_write(time_string(domain), domain % config % reftimestr, domain % spectrum, grid, &
                       domain % e(:,:,grid % istart:grid % iend), domain % k, domain % ustar)
    end associate
  end subroutine umwm_domain_write_restart


  subroutine umwm_domain_finalize(handle) bind(c, name='umwm_domain_finalize')
    ! Releases the memory held by the domain; the handle remains valid.
    type(c_ptr), value :: handle
    type(domain_type), pointer :: domain
    call c_f_pointer(handle, domain)
    call domain % finalize()
  end subroutine umwm_domain_finalize


  subroutine umwm_domain_delete(handle) bind(c, name='umwm_domain_delete')
    ! Destroys the domain; the handle must not be used afterwards.
    type(c_ptr), value :: handle
    type(domain_type), pointer :: domain
    call c_f_pointer(handle, domain)
    deallocate(domain)
  end subroutine umwm_domain_delete


  function umwm_domain_int(handle, name, value) bind(c, name='umwm_domain_int') result(stat)
    ! Integer attributes of the domain: grid and spectrum sizes and the
    ! cell index range of this process. Returns 0 on success and 1 for an
    ! unknown name.
    type(c_ptr), value :: handle
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int), intent(out) :: value
    integer(c_int) :: stat
    type(domain_type), pointer :: domain

    call c_f_pointer(handle, domain)
    stat = ok
    value = 0

    select case (fortran_string(name))
    case ('mm'); value = domain % grid % mm
    case ('nm'); value = domain % grid % nm
    case ('im'); value = domain % grid % im
    case ('imm'); value = domain % grid % imm
    case ('istart'); value = domain % grid % istart
    case ('iend'); value = domain % grid % iend
    case ('iistart'); value = domain % grid % iistart
    case ('iiend'); value = domain % grid % iiend
    case ('om'); value = domain % spectrum % num_frequencies
    case ('pm'); value = domain % spectrum % num_directions
    case ('lm'); value = domain % stokes % lm
    case ('lag'); value = domain % lag
    case ('spectrum_records'); value = domain % spectrum_records
    case default; stat = unknown_name
    end select

  end function umwm_domain_int


  function umwm_domain_real(handle, name, value) bind(c, name='umwm_domain_real') result(stat)
    ! Real attributes of the domain: time steps and spectral spacings.
    ! Returns 0 on success and 1 for an unknown name.
    type(c_ptr), value :: handle
    character(kind=c_char), intent(in) :: name(*)
    real(c_double), intent(out) :: value
    integer(c_int) :: stat
    type(domain_type), pointer :: domain

    call c_f_pointer(handle, domain)
    stat = ok
    value = 0

    select case (fortran_string(name))
    case ('dts'); value = domain % dts
    case ('sumt'); value = domain % sumt
    case ('dtamin'); value = domain % dtamin
    case ('dth'); value = domain % dth
    case ('dlnf'); value = domain % spectrum % dlnf
    case default; stat = unknown_name
    end select

  end function umwm_domain_real


  function umwm_domain_field(handle, name, ptr, dtype, itemsize, rank, shp, lbnd) &
    bind(c, name='umwm_domain_field') result(stat)
    ! Looks up a field of the domain by name and returns the address of its
    ! first element, its element type (1 real, 2 integer, 3 logical) and
    ! size in bytes, its rank, and its shape and lower bounds in Fortran
    ! (column-major) order. Names are those of the gridded and spectrum
    ! output files where they exist, with the Fortran component names as
    ! aliases. Returns 0 on success, 1 for an unknown name, and 2 if the
    ! field is not allocated.
    type(c_ptr), value :: handle
    character(kind=c_char), intent(in) :: name(*)
    type(c_ptr), intent(out) :: ptr
    integer(c_int), intent(out) :: dtype, itemsize, rank, shp(3), lbnd(3)
    integer(c_int) :: stat
    type(domain_type), pointer :: domain

    call c_f_pointer(handle, domain)
    ptr = c_null_ptr
    dtype = 0
    itemsize = 0
    rank = 0
    shp = 0
    lbnd = 0
    stat = ok

    select case (fortran_string(name))

    ! spectral coordinates (o) and (p)
    case ('frequency', 'f'); call real1(domain % f)
    case ('direction', 'theta', 'th'); call real1(domain % th)
    case ('cth'); call real1(domain % cth)
    case ('sth'); call real1(domain % sth)

    ! wave kinematics and spectral weights (o,i)
    case ('wavenumber', 'k'); call real2(domain % k)
    case ('cp0'); call real2(domain % cp0)
    case ('cg0'); call real2(domain % cg0)
    case ('kovcg'); call real2(domain % kovcg)
    case ('kdk'); call real2(domain % kdk)
    case ('k3dk'); call real2(domain % k3dk)
    case ('cothkd'); call real2(domain % cothkd)
    case ('logl2overz'); call real2(domain % logl2overz)
    case ('cgdlncp_dx'); call real2(domain % cgdlncp_dx)
    case ('cgdlncp_dy'); call real2(domain % cgdlncp_dy)
    case ('cgm'); call real2(domain % cgm)
    case ('usl'); call real2(domain % usl)
    case ('vsl'); call real2(domain % vsl)
    case ('mod_moments'); call real3(domain % mod_moments)

    ! wave spectrum and source functions
    case ('F', 'e'); call real3(domain % e)
    case ('ef'); call real3(domain % ef)
    case ('Sin', 'ssin'); call real3(domain % ssin)
    case ('Sds', 'sds'); call real3(domain % sds)
    case ('Snl', 'snl'); call real3(domain % snl)
    case ('Sdt', 'sdt'); call real2(domain % sdt)
    case ('Sdv', 'sdv'); call real2(domain % sdv)
    case ('Sbf', 'sbf'); call real2(domain % sbf)
    case ('Sice', 'sice'); call real2(domain % sice)
    case ('mss_moments'); call real3(domain % mss_moments)
    case ('oc'); call int1(domain % oc)
    case ('du_dx'); call real1(domain % du_dx)
    case ('du_dy'); call real1(domain % du_dy)
    case ('dv_dx'); call real1(domain % dv_dx)
    case ('dv_dy'); call real1(domain % dv_dy)
    case ('shelt'); call real1(domain % shelt)
    case ('physics_time_step'); call real1(domain % physics_time_step)
    case ('psim'); call real1(domain % psim)

    ! stress, drag coefficient and friction velocity
    case ('cd'); call real1(domain % cd)
    case ('ust', 'ustar'); call real1(domain % ustar)
    case ('taux'); call real1(domain % taux)
    case ('tauy'); call real1(domain % tauy)
    case ('taux_form'); call real1(domain % taux_form)
    case ('tauy_form'); call real1(domain % tauy_form)
    case ('taux_form_1', 'taux1'); call real1(domain % taux1)
    case ('tauy_form_1', 'tauy1'); call real1(domain % tauy1)
    case ('taux_form_2', 'taux2'); call real1(domain % taux2)
    case ('tauy_form_2', 'tauy2'); call real1(domain % tauy2)
    case ('taux_form_3', 'taux3'); call real1(domain % taux3)
    case ('tauy_form_3', 'tauy3'); call real1(domain % tauy3)
    case ('taux_skin'); call real1(domain % taux_skin)
    case ('tauy_skin'); call real1(domain % tauy_skin)
    case ('taux_diag'); call real1(domain % taux_diag)
    case ('tauy_diag'); call real1(domain % tauy_diag)
    case ('taux_ocn', 'taux_ocntop'); call real1(domain % taux_ocntop)
    case ('tauy_ocn', 'tauy_ocntop'); call real1(domain % tauy_ocntop)
    case ('taux_bot', 'taux_ocnbot'); call real1(domain % taux_ocnbot)
    case ('tauy_bot', 'tauy_ocnbot'); call real1(domain % tauy_ocnbot)
    case ('taux_snl'); call real1(domain % taux_snl)
    case ('tauy_snl'); call real1(domain % tauy_snl)
    case ('tailatmx'); call real1(domain % tailatmx)
    case ('tailatmy'); call real1(domain % tailatmy)
    case ('tailocnx'); call real1(domain % tailocnx)
    case ('tailocny'); call real1(domain % tailocny)
    case ('epsx_atm'); call real1(domain % epsx_atm)
    case ('epsy_atm'); call real1(domain % epsy_atm)
    case ('epsx_ocn'); call real1(domain % epsx_ocn)
    case ('epsy_ocn'); call real1(domain % epsy_ocn)

    ! integrated diagnostics
    case ('swh', 'ht'); call real1(domain % ht)
    case ('mss'); call real1(domain % mss)
    case ('mwd'); call real1(domain % mwd)
    case ('mwl'); call real1(domain % mwl)
    case ('mwp'); call real1(domain % mwp)
    case ('dwd'); call real1(domain % dwd)
    case ('dwl'); call real1(domain % dwl)
    case ('dwp'); call real1(domain % dwp)
    case ('dcp0'); call real1(domain % dcp0)
    case ('dcp'); call real1(domain % dcp)
    case ('dcg0'); call real1(domain % dcg0)
    case ('dcg'); call real1(domain % dcg)
    case ('momx'); call real1(domain % momx)
    case ('momy'); call real1(domain % momy)
    case ('cgmxx'); call real1(domain % cgmxx)
    case ('cgmxy'); call real1(domain % cgmxy)
    case ('cgmyy'); call real1(domain % cgmyy)
    case ('mod_eps'); call real1(domain % mod_eps)
    case ('mod_us'); call real1(domain % mod_us)
    case ('mod_vs'); call real1(domain % mod_vs)

    ! forcing on the cells (i) and on the native grid (m,n)
    case ('wspd'); call real1(domain % forcing % wspd)
    case ('wdir'); call real1(domain % forcing % wdir)
    case ('uc'); call real1(domain % forcing % uc)
    case ('vc'); call real1(domain % forcing % vc)
    case ('rhoa'); call real1(domain % forcing % rhoa)
    case ('rhow'); call real1(domain % forcing % rhow)
    case ('rhorat'); call real1(domain % forcing % rhorat)
    case ('fice'); call real1(domain % forcing % fice)
    case ('uw'); call real2(domain % forcing % uw)
    case ('vw'); call real2(domain % forcing % vw)
    case ('uwb'); call real2(domain % forcing % uwb)
    case ('vwb'); call real2(domain % forcing % vwb)
    case ('uwf'); call real2(domain % forcing % uwf)
    case ('vwf'); call real2(domain % forcing % vwf)
    case ('wspd_2d'); call real2(domain % forcing % wspd_2d)
    case ('wdir_2d'); call real2(domain % forcing % wdir_2d)
    case ('uc_2d'); call real2(domain % forcing % uc_2d)
    case ('vc_2d'); call real2(domain % forcing % vc_2d)
    case ('ucb'); call real2(domain % forcing % ucb)
    case ('vcb'); call real2(domain % forcing % vcb)
    case ('ucf'); call real2(domain % forcing % ucf)
    case ('vcf'); call real2(domain % forcing % vcf)
    case ('fice_2d'); call real2(domain % forcing % fice_2d)
    case ('ficeb'); call real2(domain % forcing % ficeb)
    case ('ficef'); call real2(domain % forcing % ficef)
    case ('rhoa_2d'); call real2(domain % forcing % rhoa_2d)
    case ('rhow_2d'); call real2(domain % forcing % rhow_2d)

    ! grid geometry (m,n), index maps, and cell geometry (i)
    case ('lon'); call real2(domain % grid % lon)
    case ('lat'); call real2(domain % grid % lat)
    case ('x'); call real2(domain % grid % x)
    case ('y'); call real2(domain % grid % y)
    case ('dlon'); call real2(domain % grid % dlon)
    case ('dlat'); call real2(domain % grid % dlat)
    case ('depth', 'd_2d'); call real2(domain % grid % d_2d)
    case ('dx', 'dx_2d'); call real2(domain % grid % dx_2d)
    case ('dy', 'dy_2d'); call real2(domain % grid % dy_2d)
    case ('area', 'ar_2d'); call real2(domain % grid % ar_2d)
    case ('curvature', 'curv'); call real2(domain % grid % curv)
    case ('seamask', 'mask'); call int2(domain % grid % mask)
    case ('nproc', 'nproc_out'); call int2(domain % grid % nproc_out)
    case ('ii'); call int2(domain % grid % ii)
    case ('mi'); call int1(domain % grid % mi)
    case ('ni'); call int1(domain % grid % ni)
    case ('iw'); call int1(domain % grid % iw)
    case ('ie'); call int1(domain % grid % ie)
    case ('is'); call int1(domain % grid % is)
    case ('in'); call int1(domain % grid % in)
    case ('d'); call real1(domain % grid % d)
    case ('dx_cell'); call real1(domain % grid % dx)
    case ('dy_cell'); call real1(domain % grid % dy)
    case ('ar'); call real1(domain % grid % ar)

    ! Stokes drift
    case ('z', 'stokes_depth'); call real1(domain % stokes % depth)
    case ('u_stokes', 'us'); call real2(domain % stokes % us)
    case ('v_stokes', 'vs'); call real2(domain % stokes % vs)
    case ('d_stokes', 'ds'); call real1(domain % stokes % ds)

    ! spectrum output points
    case ('spectrum_m'); call int1(domain % spectrum_m)
    case ('spectrum_n'); call int1(domain % spectrum_n)
    case ('spectrum_i'); call int1(domain % spectrum_i)

    case default
      stat = unknown_name

    end select

  contains

    subroutine real1(a)
      real, allocatable, intent(in), target :: a(:)
      if (.not. allocated(a)) then
        stat = not_allocated
        return
      end if
      if (size(a) > 0) ptr = c_loc(a)
      dtype = type_real
      itemsize = storage_size(a) / 8
      rank = 1
      shp(1) = size(a, 1)
      lbnd(1) = lbound(a, 1)
    end subroutine real1

    subroutine real2(a)
      real, allocatable, intent(in), target :: a(:,:)
      if (.not. allocated(a)) then
        stat = not_allocated
        return
      end if
      if (size(a) > 0) ptr = c_loc(a)
      dtype = type_real
      itemsize = storage_size(a) / 8
      rank = 2
      shp(1:2) = shape(a)
      lbnd(1:2) = lbound(a)
    end subroutine real2

    subroutine real3(a)
      real, allocatable, intent(in), target :: a(:,:,:)
      if (.not. allocated(a)) then
        stat = not_allocated
        return
      end if
      if (size(a) > 0) ptr = c_loc(a)
      dtype = type_real
      itemsize = storage_size(a) / 8
      rank = 3
      shp(1:3) = shape(a)
      lbnd(1:3) = lbound(a)
    end subroutine real3

    subroutine int1(a)
      integer, allocatable, intent(in), target :: a(:)
      if (.not. allocated(a)) then
        stat = not_allocated
        return
      end if
      if (size(a) > 0) ptr = c_loc(a)
      dtype = type_integer
      itemsize = storage_size(a) / 8
      rank = 1
      shp(1) = size(a, 1)
      lbnd(1) = lbound(a, 1)
    end subroutine int1

    subroutine int2(a)
      integer, allocatable, intent(in), target :: a(:,:)
      if (.not. allocated(a)) then
        stat = not_allocated
        return
      end if
      if (size(a) > 0) ptr = c_loc(a)
      dtype = type_integer
      itemsize = storage_size(a) / 8
      rank = 2
      shp(1:2) = shape(a)
      lbnd(1:2) = lbound(a)
    end subroutine int2

  end function umwm_domain_field


  function umwm_domain_gather(handle, name, level, values) &
    bind(c, name='umwm_domain_gather') result(stat)
    ! Gathers a cell field (one value per sea cell, over all processes)
    ! on the root process and remaps it to the native grid, into values of
    ! size mm*nm in column-major order, with land cells set to NaN. For the
    ! Stokes drift components, level selects the depth level (1-based).
    ! Only the root process receives the field. Collective under MPI.
    ! Returns 0 on success and 1 for an unknown name.
    type(c_ptr), value :: handle
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int), value :: level
    real(c_float), intent(inout) :: values(*)
    integer(c_int) :: stat
    type(domain_type), pointer :: domain
    real, pointer :: field(:)
    real, allocatable :: field_mn(:,:)
    logical :: all_cells

    call c_f_pointer(handle, domain)
    stat = ok
    all_cells = .false.

    select case (fortran_string(name))
    case ('shelt'); field => domain % shelt
    case ('physics_time_step'); field => domain % physics_time_step
    case ('du_dx'); field => domain % du_dx
    case ('du_dy'); field => domain % du_dy
    case ('dv_dx'); field => domain % dv_dx
    case ('dv_dy'); field => domain % dv_dy
    case ('cd'); field => domain % cd
    case ('ust', 'ustar'); field => domain % ustar
    case ('taux'); field => domain % taux
    case ('tauy'); field => domain % tauy
    case ('taux_form'); field => domain % taux_form
    case ('tauy_form'); field => domain % tauy_form
    case ('taux_form_1', 'taux1'); field => domain % taux1
    case ('tauy_form_1', 'tauy1'); field => domain % tauy1
    case ('taux_form_2', 'taux2'); field => domain % taux2
    case ('tauy_form_2', 'tauy2'); field => domain % tauy2
    case ('taux_form_3', 'taux3'); field => domain % taux3
    case ('tauy_form_3', 'tauy3'); field => domain % tauy3
    case ('taux_skin'); field => domain % taux_skin
    case ('tauy_skin'); field => domain % tauy_skin
    case ('taux_diag'); field => domain % taux_diag
    case ('tauy_diag'); field => domain % tauy_diag
    case ('taux_ocn', 'taux_ocntop'); field => domain % taux_ocntop
    case ('tauy_ocn', 'tauy_ocntop'); field => domain % tauy_ocntop
    case ('taux_bot', 'taux_ocnbot'); field => domain % taux_ocnbot
    case ('tauy_bot', 'tauy_ocnbot'); field => domain % tauy_ocnbot
    case ('taux_snl'); field => domain % taux_snl
    case ('tauy_snl'); field => domain % tauy_snl
    case ('tailatmx'); field => domain % tailatmx
    case ('tailatmy'); field => domain % tailatmy
    case ('tailocnx'); field => domain % tailocnx
    case ('tailocny'); field => domain % tailocny
    case ('epsx_atm'); field => domain % epsx_atm
    case ('epsy_atm'); field => domain % epsy_atm
    case ('epsx_ocn'); field => domain % epsx_ocn
    case ('epsy_ocn'); field => domain % epsy_ocn
    case ('swh', 'ht'); field => domain % ht
    case ('mss'); field => domain % mss
    case ('mwd'); field => domain % mwd
    case ('mwl'); field => domain % mwl
    case ('mwp'); field => domain % mwp
    case ('dwd'); field => domain % dwd
    case ('dwl'); field => domain % dwl
    case ('dwp'); field => domain % dwp
    case ('dcp0'); field => domain % dcp0
    case ('dcp'); field => domain % dcp
    case ('dcg0'); field => domain % dcg0
    case ('dcg'); field => domain % dcg
    case ('momx'); field => domain % momx
    case ('momy'); field => domain % momy
    case ('cgmxx'); field => domain % cgmxx
    case ('cgmxy'); field => domain % cgmxy
    case ('cgmyy'); field => domain % cgmyy
    case ('mod_eps'); field => domain % mod_eps
    case ('mod_us'); field => domain % mod_us
    case ('mod_vs'); field => domain % mod_vs
    case ('wspd'); field => domain % forcing % wspd; all_cells = .true.
    case ('wdir'); field => domain % forcing % wdir; all_cells = .true.
    case ('uc'); field => domain % forcing % uc; all_cells = .true.
    case ('vc'); field => domain % forcing % vc; all_cells = .true.
    case ('rhoa'); field => domain % forcing % rhoa; all_cells = .true.
    case ('rhow'); field => domain % forcing % rhow; all_cells = .true.
    case ('rhorat'); field => domain % forcing % rhorat; all_cells = .true.
    case ('fice'); field => domain % forcing % fice; all_cells = .true.
    case ('psim'); field => domain % psim; all_cells = .true.
    case ('d'); field => domain % grid % d; all_cells = .true.
    case ('u_stokes', 'us'); field => domain % stokes % us(:,level)
    case ('v_stokes', 'vs'); field => domain % stokes % vs(:,level)
    case ('d_stokes', 'ds'); field => domain % stokes % ds
    case default
      stat = unknown_name
      return
    end select

    associate(grid => domain % grid)
    allocate(field_mn(grid % mm, grid % nm))
    if (all_cells) then
      call gather_field(field(grid % istart:grid % iend), field_mn, grid)
    else
      call gather_field(field, field_mn, grid)
    end if
    if (nproc == 0) values(1:grid % mm * grid % nm) = reshape(field_mn, [grid % mm * grid % nm])
    end associate

  end function umwm_domain_gather


  function time_string(domain) result(timestr)
    ! The current model time in the form of the output file names.
    type(domain_type), intent(in) :: domain
    character(19) :: timestr
    timestr = domain % current_time % strftime('%Y-%m-%d_%H:%M:%S')
  end function time_string


  function fortran_string(cstring) result(string)
    ! Converts a null-terminated C string to a Fortran string.
    character(kind=c_char), intent(in) :: cstring(*)
    character(:), allocatable :: string
    integer :: i, n

    n = 0
    do while (cstring(n + 1) /= c_null_char)
      n = n + 1
    end do

    allocate(character(n) :: string)
    do i = 1, n
      string(i:i) = cstring(i)
    end do

  end function fortran_string

end module umwm_capi
