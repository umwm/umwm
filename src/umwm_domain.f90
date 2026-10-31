module umwm_domain

  ! The wave model domain: configuration, spectral and horizontal grids,
  ! forcing, and the wave spectrum together with the fields derived from it.

#ifdef MPI
  use mpi
#endif
  use datetime_module, only: datetime, strptime, timedelta
  use umwm_advection, only: frequency_shift, propagation, refraction
  use umwm_config, only: config_type
  use umwm_constants, only: pi, rk, twopi
  use umwm_dispersion, only: group_speed, wavenumber
  use umwm_env, only: mpisize, nproc
  use umwm_forcing, only: forcing_type, min_wind_speed
  use umwm_grid, only: grid_type
  use umwm_io, only: output_grid
#ifdef MPI
  use umwm_mpi, only: exchange_halo
#endif
  use umwm_physics, only: integrate_sources, max_source_rate
  use umwm_restart, only: restart_read, restart_write
  use umwm_sheltering, only: sheltering_coare35
  use umwm_source_functions, only: mss_moments, s_ice, sds_coefficient, sds_d12, sin_d12, snl_d12, &
                                   snl_shift
  use umwm_spectrum, only: spectrum_type
  use umwm_stokes, only: stokes_type
  use umwm_stress, only: directional_integral, drag_coefficient, &
                         drag_coefficient_skin, stress_tail, wind_relative
  use umwm_util, only: raiseexception

  implicit none

  private
  public :: domain_type

  type :: domain_type

    type(config_type) :: config
    type(spectrum_type) :: spectrum
    type(grid_type) :: grid
    type(forcing_type) :: forcing
    type(stokes_type) :: stokes

    ! model time
    type(datetime) :: current_time, stop_time

    ! time stepping
    logical :: first = .true. ! .true. during the first time step
    real :: sumt = 0          ! time integrated in the current global time step [s]
    real :: dts = 0           ! source (and advection) time step [s]
    real :: dtamin = 0        ! advection stability limit [s]

    ! grid cell and process for the screen output
    integer :: iip = 0, nproc_plot = 0

    ! derived constants
    real :: dth, dthg
    integer :: snl_bins = 0 ! the shift of the downshifting, snl_shift, in whole bins
    real :: snl_frac = 0    ! and the fraction of a bin beyond them
    real :: twopisds_fac, fieldscale1, fieldscale2, inv_sds_power

    ! spectral coordinates (o) and (p)
    real, allocatable :: f(:), th(:), cth(:), sth(:)
    real, allocatable :: cth2(:), sth2(:) ! cos(2 th), sin(2 th)

    ! wave kinematics and derived factors (o,i)
    real, allocatable :: k(:,:), cp0(:,:), cg0(:,:) ! cp0 and cg0 include halo
    real, allocatable :: kovcg(:,:) ! k/cg, action per unit ln(f) over e; includes halo
    real, allocatable :: invcp0(:,:), kdk(:,:), k3dk(:,:), k4(:,:), oneoverk4(:,:)
    real, allocatable :: fkovg(:,:), cothkd(:,:), logl2overz(:,:), psiml2(:,:)
    real, allocatable :: sbf(:,:), sdv(:,:)
    real, allocatable :: cgdlncp_dx(:,:), cgdlncp_dy(:,:) ! (cg/cp) grad(cp), for refraction

    ! wave spectrum (o,p,i); ef is the spectrum forward in time
    real, allocatable :: e(:,:,:), ef(:,:,:)

    ! source functions (o,p,i); sds excludes the plunging breaker factor cothkd.
    ! ssin and snl are kept only for output; the time step evaluates them per cell.
    real, allocatable :: ssin(:,:,:), sds(:,:,:), snl(:,:,:)
    real, allocatable :: sdt(:,:), sice(:,:)     ! (o,i)
    real, allocatable :: mss_moments(:,:,:)      ! (o,3,i) mean-squared slope moments
    integer, allocatable :: oc(:)            ! highest prognostic frequency bin
    real, allocatable :: shelt(:), physics_time_step(:), psim(:)

    ! current velocity gradient, for refraction
    real, allocatable :: du_dx(:), du_dy(:), dv_dx(:), dv_dy(:)

    ! pipelining of integration and advection over cells
    integer :: lag = 0                    ! largest distance to a neighbor in the tile
    logical, allocatable :: halo_source(:) ! cells in the halo of other tiles

    ! stress, drag coefficient and friction velocity
    real, allocatable :: cd(:), ustar(:), taux(:), tauy(:)
    real, allocatable :: taux_form(:), tauy_form(:), taux_skin(:), tauy_skin(:)
    real, allocatable :: taux_diag(:), tauy_diag(:), tailatmx(:), tailatmy(:)
    real, allocatable :: taux_ocntop(:), tauy_ocntop(:), taux_ocnbot(:), tauy_ocnbot(:)
    real, allocatable :: taux_snl(:), tauy_snl(:), tailocnx(:), tailocny(:)
    real, allocatable :: taux1(:), tauy1(:), taux2(:), tauy2(:), taux3(:), tauy3(:)
    real, allocatable :: epsx_atm(:), epsy_atm(:), epsx_ocn(:), epsy_ocn(:)

    ! integrated diagnostics
    real, allocatable :: ht(:), mss(:), mwd(:), mwl(:), mwp(:)
    real, allocatable :: dwd(:), dwl(:), dwp(:), dcp0(:), dcp(:), dcg0(:), dcg(:)
    real, allocatable :: momx(:), momy(:), cgmxx(:), cgmxy(:), cgmyy(:)

    ! spectrum output points
    integer :: spectrum_records = 0
    character(19) :: spectrum_timestr = ''
    integer, allocatable :: spectrum_m(:), spectrum_n(:), spectrum_i(:)
    character(40), allocatable :: spectrum_id(:)

  contains

    procedure :: initialize
    procedure :: run
    procedure :: advance
    procedure :: finished
    procedure :: diagnostics
    procedure :: ocean_stress
    procedure :: finalize
    procedure, private :: allocate_fields
    procedure, private :: initialize_physics
    procedure, private :: initialize_pipeline
    procedure, private :: step
    procedure, private :: step_acc
    procedure, private :: sources
    procedure, private :: wind_input
    procedure, private :: integrate_column
    procedure, private :: advect_column
    procedure, private :: update_current_gradients
    procedure, private :: wind_stress
    procedure :: output_grid_nc
    procedure :: output_spectrum_nc

  end type domain_type

  interface

    module subroutine diagnostics(self)
      class(domain_type), intent(inout) :: self
    end subroutine diagnostics

    module subroutine ocean_stress(self)
      class(domain_type), intent(inout) :: self
    end subroutine ocean_stress

    module subroutine output_grid_nc(self, timestr)
      class(domain_type), intent(in) :: self
      character(19), intent(in) :: timestr
    end subroutine output_grid_nc

    module subroutine output_spectrum_nc(self, timestr)
      class(domain_type), intent(inout) :: self
      character(19), intent(in) :: timestr
    end subroutine output_spectrum_nc

    module subroutine step_acc(self)
      class(domain_type), intent(inout) :: self
    end subroutine step_acc

    module subroutine diagnostics_acc(self)
      class(domain_type), intent(inout) :: self
    end subroutine diagnostics_acc

    module subroutine ocean_stress_acc(self)
      class(domain_type), intent(inout) :: self
    end subroutine ocean_stress_acc

    module function screen_values_acc(self, i) result(res)
      class(domain_type), intent(in) :: self
      integer, intent(in) :: i
      real :: res(6)
    end function screen_values_acc

  end interface

contains

  subroutine initialize(self, path, stat, input_dir)
    ! Reads the configuration from the namelist file at path, and
    ! constructs the grids, forcing, and the initial wave state.
    ! If stat is present, it is set to 0 on success and to 1 if the
    ! configuration is invalid, in which case the domain is left
    ! uninitialized; otherwise an invalid configuration stops the program.
    ! The input files are read from input_dir, input/ by default.
    class(domain_type), intent(inout) :: self
    character(*), intent(in) :: path
    integer, intent(out), optional :: stat
    character(*), intent(in), optional :: input_dir

    self % config = config_type(path)
    if (present(input_dir)) self % config % input_dir = input_dir
    if (.not. self % config % validate(stop_on_error=.not. present(stat), rank=nproc)) then
      if (present(stat)) stat = 1
      return
    end if
    if (present(stat)) stat = 0

    associate(config => self % config, grid => self % grid)

    self % current_time = strptime(config % starttimestr, '%Y-%m-%d %H:%M:%S')
    self % stop_time = strptime(config % stoptimestr, '%Y-%m-%d %H:%M:%S')

    self % spectrum = spectrum_type(config % om, config % pm, &
      real(config % fmin, rk), real(config % fmax, rk))

    call grid % initialize(config)
    call grid % print_diagnostics()
    call grid % initialize_direction_projection(self % spectrum)
    call self % forcing % initialize(config, grid)

    call self % allocate_fields()
    call output_grid(grid)
    call self % forcing % load(config, config % starttimestr, grid)
    call self % initialize_physics()
    call self % initialize_pipeline()

    end associate

  end subroutine initialize


  subroutine allocate_fields(self)
    class(domain_type), intent(inout) :: self

    associate(om => self % spectrum % num_frequencies, &
              pm => self % spectrum % num_directions, &
              is => self % grid % istart, ie => self % grid % iend, &
              iis => self % grid % iistart - 1, iie => self % grid % iiend)

    allocate(self % f(om), self % th(pm), self % cth(pm), self % sth(pm), &
             self % cth2(pm), self % sth2(pm))

    allocate(self % k(om,is:ie), self % invcp0(om,is:ie), self % kdk(om,is:ie), &
             self % k3dk(om,is:ie), self % k4(om,is:ie), self % oneoverk4(om,is:ie), &
             self % fkovg(om,is:ie), self % cothkd(om,is:ie), self % logl2overz(om,is:ie), &
             self % sbf(om,is:ie), self % sdv(om,is:ie), &
             self % cgdlncp_dx(om,is:ie), self % cgdlncp_dy(om,is:ie), &
             self % psiml2(om,is:ie), self % sdt(om,is:ie), self % sice(om,is:ie), &
             self % mss_moments(om,3,is:ie), source=0.)
    allocate(self % cp0(om,iis:iie), self % cg0(om,iis:iie), self % kovcg(om,iis:iie), &
             source=tiny(1.))

    ! the spectral fields are filled by a loop that runs on the GPU in the GPU
    ! build, which places their memory there
    allocate(self % e(om,pm,iis:iie), self % ef(om,pm,iis:iie))
    allocate(self % ssin(om,pm,is:ie), self % sds(om,pm,is:ie), self % snl(om,pm,is:ie))
    call fill(size(self % e), self % e, tiny(1.))
    call fill(size(self % ef), self % ef, tiny(1.))
    call fill(size(self % ssin), self % ssin, 0.)
    call fill(size(self % sds), self % sds, 0.)
    call fill(size(self % snl), self % snl, 0.)

    allocate(self % oc(is:ie), source=0)
    allocate(self % psim(self % grid % imm), source=0.)

    allocate(self % shelt(is:ie), self % physics_time_step(is:ie), &
             self % du_dx(is:ie), self % du_dy(is:ie), self % dv_dx(is:ie), self % dv_dy(is:ie), &
             self % cd(is:ie), self % ustar(is:ie), self % taux(is:ie), self % tauy(is:ie), &
             self % taux_form(is:ie), self % tauy_form(is:ie), &
             self % taux_skin(is:ie), self % tauy_skin(is:ie), &
             self % taux_diag(is:ie), self % tauy_diag(is:ie), &
             self % tailatmx(is:ie), self % tailatmy(is:ie), &
             self % taux_ocntop(is:ie), self % tauy_ocntop(is:ie), &
             self % taux_ocnbot(is:ie), self % tauy_ocnbot(is:ie), &
             self % taux_snl(is:ie), self % tauy_snl(is:ie), &
             self % tailocnx(is:ie), self % tailocny(is:ie), &
             self % taux1(is:ie), self % tauy1(is:ie), self % taux2(is:ie), &
             self % tauy2(is:ie), self % taux3(is:ie), self % tauy3(is:ie), &
             self % epsx_atm(is:ie), self % epsy_atm(is:ie), &
             self % epsx_ocn(is:ie), self % epsy_ocn(is:ie), &
             self % ht(is:ie), self % mss(is:ie), self % mwd(is:ie), &
             self % mwl(is:ie), self % mwp(is:ie), self % dwd(is:ie), &
             self % dwl(is:ie), self % dwp(is:ie), self % dcp0(is:ie), &
             self % dcp(is:ie), self % dcg0(is:ie), self % dcg(is:ie), &
             self % momx(is:ie), self % momy(is:ie), self % cgmxx(is:ie), &
             self % cgmxy(is:ie), self % cgmyy(is:ie), source=0.)

    end associate

  end subroutine allocate_fields


  subroutine fill(n, field, value)
    ! Sets all n elements of field to value.
    integer, intent(in) :: n
    real, intent(out) :: field(n)
    real, intent(in) :: value
    integer :: i
    !$acc parallel loop
    do i = 1, n
      field(i) = value
    end do
  end subroutine fill


  subroutine initialize_physics(self)
    ! Initializes spectral coordinates, wave kinematics and derived factors,
    ! time step limits, drag, friction velocity, and the Stokes drift.
    class(domain_type), intent(inout) :: self

    real, allocatable :: kd(:,:), l2(:,:), dwn(:,:)
    real :: dlnf, cfllim, mindelx
    integer :: i, o
#ifdef MPI
    integer :: n, ierr
#endif

    associate(config => self % config, grid => self % grid, forcing => self % forcing, &
              om => self % spectrum % num_frequencies, &
              pm => self % spectrum % num_directions, &
              is => self % grid % istart, ie => self % grid % iend, &
              f => self % f, k => self % k, cp0 => self % cp0, cg0 => self % cg0)

    dlnf = real(self % spectrum % dlnf, kind(dlnf))
    self % dth = real(self % spectrum % dth, kind(self % dth))
    f = real(self % spectrum % frequency, kind(f))
    self % th = real(self % spectrum % direction, kind(self % th))
    self % cth = cos(self % th)
    self % sth = sin(self % th)
    self % cth2 = cos(2 * self % th)
    self % sth2 = sin(2 * self % th)

    self % dthg = self % dth * config % g
    self % twopisds_fac = twopi * config % sds_fac
    self % fieldscale1 = config % sin_diss1 / config % sin_fac
    self % fieldscale2 = config % sin_diss2 / config % sin_diss1
    self % inv_sds_power = 1. / config % sds_power

    ! this limits the Courant number to its theoretical value slightly
    ! larger than 1/sqrt(2), depending on the number of directional bins;
    ! this also limits the number of directions to be divisible by 8,
    ! not 4! (more isotropic in Cartesian projection)
    if (mod(pm, 8) == 0) then
      cfllim = cos(0.25 * pi - 0.5 * self % dth)
    else
      cfllim = 1. / sqrt(2.)
    end if

    ! the downshifting shift in bins, whole and fraction:
    self % snl_bins = int(snl_shift / dlnf)
    self % snl_frac = snl_shift / dlnf - self % snl_bins

    ! wave numbers, phase speeds, and group velocities:
    if (nproc == 0) write(*,'(a)') 'umwm: dispersion: solving for dispersion relationship;'

    !$acc parallel loop collapse(2)
    do i = is, ie
      do o = 1, om
        k(o,i) = real(wavenumber(real(f(o), rk), real(grid % d(i), rk), real(config % rhow0, rk), &
                                 real(config % g, rk), real(config % sfct, rk)), kind(k))
      end do
    end do

    if (nproc == 0) write(*,'(a)') 'umwm: dispersion: dispersion relationship done;'

    allocate(kd(om,is:ie), l2(om,is:ie), dwn(om,is:ie))

    ! limit kd to avoid floating overflow in transcendental functions:
    !$acc parallel loop collapse(2)
    do i = is, ie
      do o = 1, om
        kd(o,i) = min(k(o,i) * grid % d(i), 20.)
      end do
    end do

    !$acc parallel loop collapse(2)
    do i = is, ie
      do o = 1, om
        cp0(o,i) = twopi * f(o) / k(o,i)
        cg0(o,i) = real(group_speed(real(k(o,i), rk), real(grid % d(i), rk), &
                                    real(config % rhow0, rk), real(config % g, rk), &
                                    real(config % sfct, rk)), kind(cg0))
      end do
    end do

    ! frequently used factors:
    !$acc parallel loop collapse(2)
    do i = is, ie
      do o = 1, om
        dwn(o,i)       = twopi * dlnf * f(o) / abs(cg0(o,i))           ! dk
        l2(o,i)        = 0.5 * abs(cp0(o,i)) / f(o)                    ! lambda/2
        self % k4(o,i) = k(o,i)**4.                                    ! k^4
        self % oneoverk4(o,i) = 1. / self % k4(o,i)                    ! k^-4
        self % kdk(o,i)  = k(o,i) * dwn(o,i)                           ! k*dk
        self % k3dk(o,i) = k(o,i)**3. * dwn(o,i)                       ! k*k*k*dk
        self % fkovg(o,i) = f(o) * k(o,i) / config % g                 ! f*k/g
        self % cothkd(o,i) = cosh(0.2 * kd(o,i)) / sinh(0.2 * kd(o,i)) ! coth(0.2*kd)
        self % invcp0(o,i) = 1. / cp0(o,i)                             ! 1/cp
        self % sbf(o,i) = config % sbf_fac * k(o,i) / (sinh(2. * kd(o,i))) &      ! bottom friction
                        + config % sbp_fac * k(o,i) / (cosh(kd(o,i)) * cosh(kd(o,i))) ! bottom percolation
        self % sdv(o,i) = 4. * config % nu_water * k(o,i)**2.          ! viscosity
      end do
    end do

    ! half-wavelength over z, limited to 10 m for l/2 > 10 m:
    self % logl2overz = log(l2 / config % z)
    where (l2 > 20.) self % logl2overz = log(20. / config % z)

#ifdef MPI
    call exchange_halo(cp0, om, grid, config % isglobal)
    call exchange_halo(cg0, om, grid, config % isglobal)
#endif

    ! land points for cp and cg (needed for advection/refraction), the
    ! smallest over all sea points:
    do o = 1, om
      cp0(o,grid % iistart-1) = minval(cp0(o,is:ie))
      cg0(o,grid % iistart-1) = minval(cg0(o,is:ie))
    end do
#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE, cp0(:,grid % iistart-1), om, MPI_REAL, MPI_MIN, &
                       MPI_COMM_WORLD, ierr)
    call mpi_allreduce(MPI_IN_PLACE, cg0(:,grid % iistart-1), om, MPI_REAL, MPI_MIN, &
                       MPI_COMM_WORLD, ierr)
#endif

    ! k/cg, the Jacobian from the variance spectrum to the action per unit
    ! ln(f), on the sea cells, the halo, and land:
    !$acc parallel loop collapse(2)
    do i = grid % iistart-1, grid % iiend
      do o = 1, om
        self % kovcg(o,i) = twopi * f(o) / (cp0(o,i) * cg0(o,i))
      end do
    end do

    ! phase speed gradients at fixed frequency, times cg/cp, for the
    ! bottom-induced refraction, over the sea neighbors of each cell:
    do i = is, ie
      self % cgdlncp_dx(:,i) = cg0(:,i) * self % invcp0(:,i) * grid % oneovdx(i) &
        * sea_difference(cp0(:,i), cp0(:,grid % ie(i)), cp0(:,grid % iw(i)), &
                         grid % iie(i) <= grid % im, grid % iiw(i) <= grid % im)
      self % cgdlncp_dy(:,i) = cg0(:,i) * self % invcp0(:,i) * grid % oneovdy(i) &
        * sea_difference(cp0(:,i), cp0(:,grid % in(i)), cp0(:,grid % is(i)), &
                         grid % iin(i) <= grid % im, grid % iis(i) <= grid % im)
    end do

    mindelx = min(minval(grid % dx_2d, grid % mask == 1), minval(grid % dy_2d, grid % mask == 1))
    self % dtamin = 0.98 * cfllim * mindelx / maxval(cg0(:,is:ie))

    ! if sea ice from file, update the fice field
    if (config % seaice) forcing % fice = grid % remap_mn2i(forcing % ficef)

    ! if forcing from file, update the wspd field for ustar first guess
    if (config % winds) forcing % wspd = grid % remap_mn2i(sqrt(forcing % uwf**2 + forcing % vwf**2))
    call forcing % apply_wind_speed_floor()

    ! initialize drag coefficient (Large and Pond, 1981) and friction velocity:
    self % cd = 1.2e-3
    where (forcing % wspd(is:ie) > 11) self % cd = (0.49 + 0.065 * forcing % wspd(is:ie)) * 1e-3
    self % ustar = sqrt(self % cd) * forcing % wspd(is:ie)

    ! figure out which process will print to screen:
    self % iip = grid % ii(config % xpl, config % ypl)
#ifdef MPI
    if (grid % mask(config % xpl, config % ypl) == 0) then
      if (nproc == 0) then
        write(0,*) config % xpl, config % ypl, grid % mask(config % xpl, config % ypl)
        call raiseexception('warning', 'init', &
                            'land-point chosen for stdout, may go out of bounds')
        call raiseexception('warning', 'init', &
                            'check xpl and ypl in the output namelist in ' // trim(config % path))
        stop
      end if
    end if

    if (nproc == 0) then
      do n = 0, mpisize - 1
        if (self % iip >= grid % istart_all(n) .and. self % iip <= grid % iend_all(n)) &
          self % nproc_plot = n
      end do
    end if

    call mpi_bcast(self % nproc_plot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(self % iip, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#else
    write(*,fmt=102)
    write(*,*) 'initialization summary:'
    write(*,fmt=103)
    write(*,*) 'bin  f[hz]    t[s]    min(k)   max(k)   min(c)   max(c)   min(cg)   max(cg)'
    write(*,fmt=103)
    do o = 1, om
      write(*,fmt=101) o, f(o), 1. / f(o), minval(k(o,:), dim=1), maxval(k(o,:), dim=1), &
                       minval(cp0(o,1:grid % im), dim=1), maxval(cp0(o,1:grid % im), dim=1), &
                       minval(cg0(o,1:grid % im), dim=1), maxval(cg0(o,1:grid % im), dim=1)
    end do
    write(*,fmt=102)
101 format(1x,i2,8(2x,f7.4))
102 format('!',77('='),'!')
103 format('!',77('-'),'!')
#endif

    call self % stokes % initialize(config % stokes_depths, f, k, dwn, grid % d(is:ie), &
                                    self % dth, is, ie)
    if (nproc == 0) write(*, '(a)') 'umwm: stokes_drift: initialized'
    call self % stokes % compute(self % e(:,:,is:ie), self % cth, self % sth)

    end associate

  end subroutine initialize_physics


  subroutine initialize_pipeline(self)
    ! Finds the cells sent to the halo of other tiles, and the largest
    ! distance to a neighbor within the tile, for pipelining the integration
    ! and advection over cells (see step).
    class(domain_type), intent(inout) :: self
    integer :: i

    associate(grid => self % grid, is => self % grid % istart, ie => self % grid % iend)

    allocate(self % halo_source(is:ie), source=.false.)
#ifdef MPI
    if (nproc > 0) self % halo_source(is:min(grid % iiend_all(nproc-1), ie)) = .true.
    if (nproc < mpisize - 1) self % halo_source(max(grid % iistart_all(nproc+1), is):ie) = .true.
    if (self % config % isglobal .and. mpisize > 1 .and. (nproc == 0 .or. nproc == mpisize - 1)) &
      self % halo_source(grid % i_exchange_indices) = .true.
#endif

    self % lag = 0
    do i = is, ie
      associate(neighbors => [grid % ie(i), grid % iw(i), grid % in(i), grid % is(i)])
        self % lag = max(self % lag, maxval(neighbors - i, mask=neighbors >= is .and. neighbors <= ie))
      end associate
    end do

    end associate

  end subroutine initialize_pipeline


  subroutine run(self)
    ! Integrates the model from the current time to the stop time.
    class(domain_type), intent(inout) :: self

    do while (.not. self % finished())
      call self % advance()
    end do

  end subroutine run


  pure logical function finished(self)
    ! Whether the model has reached the stop time.
    class(domain_type), intent(in) :: self
    finished = .not. (self % current_time < self % stop_time)
  end function finished


  subroutine advance(self)
    ! Integrates the model over one global time step: updates the forcing,
    ! integrates the source time steps, and writes gridded, spectrum, and
    ! restart output at their output times.
    class(domain_type), intent(inout) :: self

    character(19) :: timestr
    logical :: fullhour
#ifdef MPI
    integer :: ierr
#endif

    associate(config => self % config, grid => self % grid, forcing => self % forcing, &
              current_time => self % current_time, &
              is => self % grid % istart, ie => self % grid % iend)

    timestr = trim(current_time % strftime('%Y-%m-%d_%H:%M:%S'))

    ! read wave spectrum field from a restart file if necessary:
    if (self % first .and. config % restart) &
      call restart_read(config % starttimestr, self % spectrum, grid, &
                        self % e(:,:,is:ie), self % ustar)

    ! report current and next checkpoint time:
    if (nproc == 0) write(*,fmt='(a)') 'umwm: solver: current time is:     ' // timestr

    current_time = current_time + timedelta(seconds=nint(config % dtg))
    timestr = trim(current_time % strftime('%Y-%m-%d_%H:%M:%S'))

    if (nproc == 0) write(*,fmt='(a)') 'umwm: solver: integrating to time: ' // timestr

    ! read atmosphere and ocean input fields from file.
    ! in case of esmf coupling, fields are assumed to be updated
    ! externally, and this call is not used.
#ifndef ESMF
    call forcing % update(config, timestr, grid)
#endif

    ! inner time loop: global time step
    self % sumt = 0
    do while (self % sumt < config % dtg)

#ifndef ESMF
      call forcing % interpolate(config, grid, self % sumt) ! interpolate force fields in time
#else
      call forcing % apply_wind_speed_floor()
#endif

#if defined(_OPENACC) && !defined(ESMF)
      call self % step_acc()
#else
      call self % step()
#endif

      if (self % first) then

        ! diagnostic calculations before output
        call self % diagnostics()

        if (config % outgrid > 0) call self % output_grid_nc(config % starttimestr)
        if (config % outspec > 0) call self % output_spectrum_nc(config % starttimestr)

#ifdef MPI
        call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

        if (nproc == 0) then
          ! diagnostic output header on screen
          write(*,fmt='(a)') 'dtg frac   tstep [s]  wspd [m/s] wdir [rad]   ' &
                             // 'swh [m]    mwp [s]    cd*10^3   fprog [hz]'
        end if

        self % first = .false.

      end if

      ! diagnostic output on screen
      if (nproc == self % nproc_plot) then
        associate(i => self % iip)
#ifdef _OPENACC
          write(*,fmt=100) self % sumt / config % dtg, self % dts, screen_values_acc(self, i)
#else
          write(*,fmt=100) self % sumt / config % dtg, self % dts, forcing % wspd(i), &
                           forcing % wdir(i), sigwaveheight(self, i), &
                           meanwaveperiod(self, i), self % cd(i) * 1e3, self % f(self % oc(i))
#endif
        end associate
      end if

    end do ! end while(sumt<dtg) loop

    if (config % stokes) call self % stokes % compute(self % e(:,:,is:ie), self % cth, self % sth)

    call self % diagnostics() ! model diagnostics for output

    fullhour = current_time % getminute() == 0 .and. current_time % getsecond() == 0

    ! gridded output
    if (config % outgrid > 0) then
      if (mod(current_time % gethour(), config % outgrid) == 0 .and. fullhour) then
#ifndef ESMF
        call self % ocean_stress()
#endif
        call self % output_grid_nc(timestr)
      end if
    else if (config % outgrid == -1) then
#ifndef ESMF
      call self % ocean_stress()
#endif
      call self % output_grid_nc(timestr)
    end if

    ! spectrum output
    if (config % outspec > 0) then
      if (mod(current_time % gethour(), config % outspec) == 0 .and. fullhour) &
        call self % output_spectrum_nc(timestr)
    else if (config % outspec == -1) then
      call self % output_spectrum_nc(timestr)
    end if

    ! restart output
    if (config % outrst > 0) then
      if (mod(current_time % gethour(), config % outrst) == 0 .and. fullhour) &
        call restart_write(timestr, config % reftimestr, self % spectrum, grid, &
                           self % e(:,:,is:ie), self % k, self % ustar)
    else if (config % outrst == -1) then
      call restart_write(timestr, config % reftimestr, self % spectrum, grid, &
                         self % e(:,:,is:ie), self % k, self % ustar)
    end if

    end associate

100 format(2x, f4.2, 2x, f8.3, 6(1x, f9.6))

  end subroutine advance


  subroutine finalize(self)
    ! Releases all memory held by the domain; intent(out) deallocates
    ! all allocatable components and restores default values.
    class(domain_type), intent(out) :: self
  end subroutine finalize


  subroutine step(self)
    ! Advances the wave spectrum by one source time step: evaluates and
    ! integrates the source functions, then propagates and refracts the
    ! spectrum, and evaluates the wind stress.
    class(domain_type), intent(inout) :: self
    real, allocatable :: e_old(:,:,:)
    logical :: store, currents, update_gradients
    integer :: i
#ifdef MPI
    real :: dts_local
    integer :: ierr
#endif

    associate(config => self % config, grid => self % grid, forcing => self % forcing)

    call self % sources()

    ! dynamic time step and in-step model time:
    self % dts = min(minval(self % physics_time_step), self % dtamin, config % dtg - self % sumt)
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
    ! first step; both always in coupled mode.
#ifdef ESMF
    store = .true.
    update_gradients = .true.
#else
    store = self % first .or. self % sumt >= config % dtg
    update_gradients = config % currents .or. self % first
#endif
    if (update_gradients) call self % update_current_gradients()

    ! check if currents are non-zero:
    if (.not. config % isglobal) then
      currents = any(forcing % uc(grid % iistart:grid % iiend) /= 0) &
            .or. any(forcing % vc(grid % iistart:grid % iiend) /= 0)
    else
      currents = any(forcing % uc /= 0) .or. any(forcing % vc /= 0)
    end if

    ! Integrate the cells in the halo of other tiles first, and exchange the halo.
    do i = grid % istart, grid % iend
      if (self % halo_source(i)) call self % integrate_column(i, store)
    end do
#ifdef MPI
    call exchange_halo(self % e, size(self % e, 1) * size(self % e, 2), grid, config % isglobal)
#endif

    ! Integrate and advect the remaining cells, with advection lagging behind
    ! so that the spectra of all neighbors are integrated and still in cache.
    do i = grid % istart, grid % iend + self % lag
      if (i <= grid % iend) then
        if (.not. self % halo_source(i)) call self % integrate_column(i, store)
      end if
      if (i - self % lag >= grid % istart) call self % advect_column(i - self % lag, currents)
    end do

    end associate

    ! the forward spectrum becomes the current one
    call move_alloc(self % e, e_old)
    call move_alloc(self % ef, self % e)
    call move_alloc(e_old, self % ef)

#ifdef ESMF
    call self % ocean_stress() ! stress into ocean top and bottom
#endif

  end subroutine step


  subroutine sources(self)
    ! Evaluates the source functions, the column-wise sink terms, and the
    ! largest stable source time step at each cell.
    class(domain_type), intent(inout) :: self
    real :: ssin(self % spectrum % num_frequencies, self % spectrum % num_directions)
    real :: wspd, fcutoff
    integer :: i, o

    associate(config => self % config, forcing => self % forcing, &
              om => self % spectrum % num_frequencies, pm => self % spectrum % num_directions)

    do i = self % grid % istart, self % grid % iend

      ! protection against low wind speed values
      wspd = max(forcing % wspd(i), min_wind_speed)

      ! cut-off frequency (by default 4*pierson-moskowitz peak frequency) and its
      ! bin, below the bins that the downshifting into it reads
      fcutoff = min(config % fcut_fac * config % g / wspd, config % fprog)
      do o = om - max(2, self % snl_bins + 1), 1, -1
        self % oc(i) = o
        if (fcutoff > self % f(o)) exit
      end do

      self % shelt(i) = sheltering_coare35(wspd)
      call self % wind_input(i, self % oc(i), ssin)

      ! breaking dissipation for the range used by the integration, with the
      ! diagnostic bins that the downshifting reads; the rest of the
      ! diagnostic range is evaluated for output only
      call mss_moments(om, pm, self % e(:,:,i), self % k3dk(:,i), self % cth2, self % sth2, &
                       self % dth, self % mss_moments(:,:,i))
      call sds_d12(om, pm, 1, self % oc(i) + self % snl_bins + 1, self % e(:,:,i), &
                   self % k4(:,i), self % f, self % cth2, self % sth2, config % mss_fac, &
                   self % twopisds_fac, self % mss_moments(:,:,i), config % sds_power, &
                   self % sds(:,:,i))

      ! dissipation due to turbulence
      self % sdt(:,i) = config % sdt_fac * sqrt(forcing % rhorat(i)) * self % ustar(i) * self % k(:,i)

      call s_ice(om, pm, self % e(:,:,i), self % kdk(:,i), self % cg0(:,i), self % dth, &
                 forcing % fice(i), config % fice_lth, self % sice(:,i))

      self % physics_time_step(i) = config % explim &
        / max_source_rate(om, pm, self % oc(i), ssin, self % sds(:,:,i), &
                          self % cothkd(:,i), self % kdk(:,i), config % snl_fac, &
                          self % snl_bins, self % snl_frac, self % sbf(:,i), &
                          self % sdt(:,i), self % sdv(:,i), self % sice(:,i))

    end do

    end associate

  end subroutine sources


  subroutine wind_input(self, i, no, ssin)
    ! Evaluates the wind input at cell i for frequencies 1:no.
    class(domain_type), intent(in) :: self
    integer, intent(in) :: i, no
    real, intent(inout) :: ssin(self % spectrum % num_frequencies, self % spectrum % num_directions)
    real :: wspd, wind(self % spectrum % num_frequencies), cosd(self % spectrum % num_directions)

    associate(forcing => self % forcing)

    ! wind speed at height of half wavelength, protected against low wind speed
    wspd = max(forcing % wspd(i), min_wind_speed)
#ifdef ESMF
    wind = wspd + 2.5 * self % ustar(i) * (self % logl2overz(:,i) + self % psim(i) - self % psiml2(:,i))
#else
    wind = wspd + 2.5 * self % ustar(i) * self % logl2overz(:,i)
#endif
    cosd = cos(forcing % wdir(i) - self % th)

    call sin_d12(self % spectrum % num_frequencies, self % spectrum % num_directions, &
                 no, self % oc(i), wind, cosd, self % cp0(:,i), forcing % uc(i), forcing % vc(i), &
                 self % cth, self % sth, self % shelt(i), self % config % sin_fac, &
                 self % fieldscale1, self % fieldscale2, &
                 (1 - forcing % fice(i)) * twopi * forcing % rhorat(i), self % fkovg(:,i), ssin)

    end associate

  end subroutine wind_input


  subroutine integrate_column(self, i, store)
    ! Integrates the source functions at cell i over the time step dts.
    ! If store is .true., Sin, Sds and Snl over the whole spectrum are kept
    ! for output.
    class(domain_type), intent(inout) :: self
    integer, intent(in) :: i
    logical, intent(in) :: store
    real, dimension(self % spectrum % num_frequencies, self % spectrum % num_directions) :: &
      ssin, snl, sds_coef

    associate(config => self % config, &
              om => self % spectrum % num_frequencies, pm => self % spectrum % num_directions)

    call self % wind_input(i, om, ssin)
    if (store) &
      call sds_d12(om, pm, self % oc(i) + self % snl_bins + 2, om, self % e(:,:,i), &
                   self % k4(:,i), self % f, self % cth2, self % sth2, config % mss_fac, &
                   self % twopisds_fac, self % mss_moments(:,:,i), config % sds_power, &
                   self % sds(:,:,i))
    call snl_d12(om, pm, self % oc(i), self % e(:,:,i), self % kdk(:,i), config % snl_fac, &
                 self % snl_bins, self % snl_frac, self % sds(:,:,i), snl)
    ! dissipation coefficient for the diagnostic range
    call sds_coefficient(om, pm, self % oc(i) + 1, om, self % f, self % cth2, self % sth2, &
                         config % mss_fac, self % twopisds_fac, self % mss_moments(:,:,i), sds_coef)
    call integrate_sources(om, pm, self % oc(i), self % dts, self % e(:,:,i), ssin, &
                           self % sds(:,:,i), snl, sds_coef, self % sbf(:,i), self % sdt(:,i), &
                           self % sdv(:,i), self % sice(:,i), self % cothkd(:,i), &
                           self % oneoverk4(:,i), self % inv_sds_power, self % ef(:,:,i))
    if (store) then
      self % ssin(:,:,i) = ssin
      self % snl(:,:,i) = snl
    end if

    end associate

  end subroutine integrate_column


  subroutine advect_column(self, i, currents)
    ! Propagates and refracts the spectrum at cell i, and with currents
    ! shifts it in frequency, conserving the wave action, and integrates
    ! into ef; evaluates the wind stress from the updated spectrum. Advection
    ! by ocean currents is included if currents is .true.
    class(domain_type), intent(inout) :: self
    integer, intent(in) :: i
    logical, intent(in) :: currents
#ifdef ESMF
    real :: e_propagated(self % spectrum % num_frequencies, self % spectrum % num_directions)
#endif

    associate(config => self % config, grid => self % grid, &
              om => self % spectrum % num_frequencies, pm => self % spectrum % num_directions, &
              e => self % e, ef => self % ef, cg0 => self % cg0, kovcg => self % kovcg, &
              uc => self % forcing % uc, vc => self % forcing % vc, &
              ie => self % grid % ie(i), iw => self % grid % iw(i), &
              in => self % grid % in(i), is => self % grid % is(i))

    if (self % forcing % fice(i) > config % fice_uth) then
      ef(1:self % oc(i),:,i) = 0 ! no propagation through ice
    else

      call propagation(om, pm, self % oc(i), e(:,:,i), e(:,:,ie), e(:,:,iw), e(:,:,in), &
                       e(:,:,is), cg0(:,i), cg0(:,ie), cg0(:,iw), cg0(:,in), cg0(:,is), &
                       kovcg(:,i), kovcg(:,ie), kovcg(:,iw), kovcg(:,in), kovcg(:,is), &
                       grid % cth_curv(:,i), grid % sth_curv(:,i), grid % dxn(i), grid % dxs(i), &
                       grid % dye(i), grid % dyw(i), currents, &
                       uc(i) + uc(grid % iie(i)), uc(i) + uc(grid % iiw(i)), &
                       vc(i) + vc(grid % iin(i)), vc(i) + vc(grid % iis(i)), &
                       0.25 * self % dts * grid % oneovar(i), ef(:,:,i))

#ifdef ESMF
      ! in coupled mode, refraction and the frequency shift act on the
      ! propagated spectrum
      e_propagated = ef(:,:,i)
      call refraction(om, pm, self % oc(i), e_propagated, &
                      self % cgdlncp_dx(:,i), self % cgdlncp_dy(:,i), &
                      self % du_dx(i), self % du_dy(i), self % dv_dx(i), self % dv_dy(i), &
                      grid % cth_curv(:,i), grid % sth_curv(:,i), self % dts, self % dth, &
                      ef(:,:,i))
      if (currents) &
        call frequency_shift(om, pm, self % oc(i), e_propagated, kovcg(:,i), &
                             self % invcp0(:,i), cg0(:,i), &
                             self % cgdlncp_dx(:,i), self % cgdlncp_dy(:,i), uc(i), vc(i), &
                             self % du_dx(i), self % du_dy(i), self % dv_dx(i), self % dv_dy(i), &
                             grid % cth_curv(:,i), grid % sth_curv(:,i), self % dts, &
                             real(self % spectrum % dlnf), ef(:,:,i))
#else
      call refraction(om, pm, self % oc(i), e(:,:,i), &
                      self % cgdlncp_dx(:,i), self % cgdlncp_dy(:,i), &
                      self % du_dx(i), self % du_dy(i), self % dv_dx(i), self % dv_dy(i), &
                      grid % cth_curv(:,i), grid % sth_curv(:,i), self % dts, self % dth, &
                      ef(:,:,i))
      if (currents) &
        call frequency_shift(om, pm, self % oc(i), e(:,:,i), kovcg(:,i), &
                             self % invcp0(:,i), cg0(:,i), &
                             self % cgdlncp_dx(:,i), self % cgdlncp_dy(:,i), uc(i), vc(i), &
                             self % du_dx(i), self % du_dy(i), self % dv_dx(i), self % dv_dy(i), &
                             grid % cth_curv(:,i), grid % sth_curv(:,i), self % dts, &
                             real(self % spectrum % dlnf), ef(:,:,i))
#endif

    end if

    end associate

    call self % wind_stress(i)

  end subroutine advect_column


  subroutine update_current_gradients(self)
    ! Updates the current velocity gradient for the current-induced
    ! refraction and frequency shift, over the sea neighbors of each cell.
    class(domain_type), intent(inout) :: self
    associate(grid => self % grid)
    call current_gradients(grid % istart, grid % iend, grid % im, grid % imm, &
                           self % forcing % uc, self % forcing % vc, grid % iie, grid % iiw, &
                           grid % iin, grid % iis, grid % oneovdx, grid % oneovdy, &
                           self % du_dx, self % du_dy, self % dv_dx, self % dv_dy)
    end associate
  end subroutine update_current_gradients


  subroutine current_gradients(i1, i2, im, imm, uc, vc, iie, iiw, iin, iis, oneovdx, oneovdy, &
                               du_dx, du_dy, dv_dx, dv_dy)
    ! Gradient of the current velocity (uc, vc) at the cells i1:i2, over
    ! their sea neighbors, with true indices iie, iiw, iin, and iis; on the
    ! GPU in the GPU build.
    integer, intent(in) :: i1, i2, im, imm
    real, intent(in) :: uc(imm), vc(imm)
    integer, intent(in) :: iie(imm), iiw(imm), iin(imm), iis(imm)
    real, intent(in) :: oneovdx(i1:i2), oneovdy(i1:i2)
    real, intent(out) :: du_dx(i1:i2), du_dy(i1:i2), dv_dx(i1:i2), dv_dy(i1:i2)
    logical :: sea_e, sea_w, sea_n, sea_s
    integer :: i

    !$acc parallel loop private(sea_e, sea_w, sea_n, sea_s)
    do i = i1, i2
      sea_e = iie(i) <= im
      sea_w = iiw(i) <= im
      sea_n = iin(i) <= im
      sea_s = iis(i) <= im
      du_dx(i) = sea_difference(uc(i), uc(iie(i)), uc(iiw(i)), sea_e, sea_w) * oneovdx(i)
      du_dy(i) = sea_difference(uc(i), uc(iin(i)), uc(iis(i)), sea_n, sea_s) * oneovdy(i)
      dv_dx(i) = sea_difference(vc(i), vc(iie(i)), vc(iiw(i)), sea_e, sea_w) * oneovdx(i)
      dv_dy(i) = sea_difference(vc(i), vc(iin(i)), vc(iis(i)), sea_n, sea_s) * oneovdy(i)
    end do

  end subroutine current_gradients


  pure elemental real function sea_difference(f, f_plus, f_minus, sea_plus, sea_minus) result(res)
    ! Difference per cell width of a field across a cell, from its value f
    ! at the cell and f_plus and f_minus at the neighbors in the positive
    ! and negative directions: centered between two sea neighbors, one-sided
    ! next to land, and zero between two land neighbors. Land values, such
    ! as a current set to zero or the phase speed of land cells, do not
    ! enter.
    !$acc routine seq
    real, intent(in) :: f, f_plus, f_minus
    logical, intent(in) :: sea_plus, sea_minus
    if (sea_plus .and. sea_minus) then
      res = 0.5 * (f_plus - f_minus)
    else if (sea_plus) then
      res = f_plus - f
    else if (sea_minus) then
      res = f - f_minus
    else
      res = 0
    end if
  end function sea_difference


  subroutine wind_stress(self, i)
    ! Evaluates form and skin drag, drag coefficient, and friction velocity
    ! at cell i from the updated spectrum ef.
    class(domain_type), intent(inout) :: self
    integer, intent(in) :: i
    real :: ssin(self % spectrum % num_frequencies, self % spectrum % num_directions)
    real, dimension(self % spectrum % num_frequencies) :: taux_util, tauy_util
    real :: tail, rhow_dthg, wspdrel, wdirrel, cd_form, cd_skin

    associate(config => self % config, forcing => self % forcing, &
              om => self % spectrum % num_frequencies, oc => self % oc(i))

    ! form drag from atmosphere (positive into waves)
    call self % wind_input(i, om, ssin)
    call directional_integral(om, self % spectrum % num_directions, self % ef(:,:,i), &
                              ssin, self % cth, self % sth, taux_util, tauy_util)
    taux_util = taux_util * self % invcp0(:,i)
    tauy_util = tauy_util * self % invcp0(:,i)

    ! wind speed dependent tail
    tail = stress_tail(forcing % wspd(i), self % k(om,i))
    rhow_dthg = forcing % rhow(i) * self % dthg

    self % tailatmx(i) = taux_util(om) * self % k(om,i) * tail * rhow_dthg
    self % tailatmy(i) = tauy_util(om) * self % k(om,i) * tail * rhow_dthg

    self % taux_form(i) = sum(taux_util * self % kdk(:,i)) * rhow_dthg + self % tailatmx(i)
    self % tauy_form(i) = sum(tauy_util * self % kdk(:,i)) * rhow_dthg + self % tailatmy(i)

    self % taux_diag(i) = sum(taux_util(oc:om) * self % kdk(oc:om,i)) * rhow_dthg
    self % tauy_diag(i) = sum(tauy_util(oc:om) * self % kdk(oc:om,i)) * rhow_dthg

    ! wind speed and direction relative to surface velocity
    call wind_relative(forcing % wspd(i), forcing % wdir(i), &
                       forcing % uc(i) + self % stokes % us(i,1), &
                       forcing % vc(i) + self % stokes % vs(i,1), wspdrel, wdirrel)

    ! form- and skin-induced drag coefficients
    cd_form = drag_coefficient(self % taux_form(i), self % tauy_form(i), &
                               forcing % rhoa(i), forcing % wspd(i))
    cd_skin = drag_coefficient_skin(cd_form, forcing % wspd(i), wspdrel, &
                                    config % z, config % nu_air, config % kappa)

    self % taux_skin(i) = forcing % rhoa(i) * cd_skin * wspdrel**2 * cos(wdirrel)
    self % tauy_skin(i) = forcing % rhoa(i) * cd_skin * wspdrel**2 * sin(wdirrel)

    self % taux(i) = self % taux_form(i) + self % taux_skin(i)
    self % tauy(i) = self % tauy_form(i) + self % tauy_skin(i)

    ! total (form + skin) drag coefficient
    self % cd(i) = drag_coefficient(self % taux(i), self % tauy(i), forcing % rhoa(i), forcing % wspd(i))

    ! friction velocity based on total (form + skin) stress
    self % ustar(i) = sqrt(sqrt(self % taux(i)**2 + self % tauy(i)**2) / forcing % rhoa(i))

    end associate

  end subroutine wind_stress


  pure real function sigwaveheight(self, i) result(swh)
    ! Significant wave height at cell i.
    class(domain_type), intent(in) :: self
    integer, intent(in) :: i
    integer :: p
    swh = 0
    do p = 1, self % spectrum % num_directions
      swh = swh + sum(self % e(:,p,i) * self % kdk(:,i))
    end do
    swh = 4. * sqrt(swh * self % dth)
  end function sigwaveheight


  pure real function meanwaveperiod(self, i) result(mwp)
    ! Mean wave period at cell i.
    class(domain_type), intent(in) :: self
    integer, intent(in) :: i
    real :: m0, m2
    integer :: p
    m0 = 0
    m2 = 0
    do p = 1, self % spectrum % num_directions
      m0 = m0 + sum(self % e(:,p,i) * self % kdk(:,i))
      m2 = m2 + sum(self % f**2 * self % e(:,p,i) * self % kdk(:,i))
    end do
    mwp = sqrt(m0 / m2)
  end function meanwaveperiod

end module umwm_domain
