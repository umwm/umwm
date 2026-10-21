module umwm_domain

  ! The wave model domain: configuration, spectral and horizontal grids,
  ! forcing, and the wave spectrum together with the fields derived from it.

#ifdef MPI
  use mpi
#endif
  use datetime_module, only: datetime, strptime, timedelta
  use umwm_advection, only: max_rotation_rate, propagation, refraction
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
  use umwm_source_functions, only: mss_moments, s_ice, sds_coefficient, sds_d12, sin_d12, snl_d12
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
    real :: dtr = 0           ! refraction time step [s]
    real :: dtr_max = 0       ! refraction stability limit [s]
    real :: dtamin = 0        ! advection stability limit [s]

    ! grid cell and process for the screen output
    integer :: iip = 0, nproc_plot = 0

    ! derived constants
    real :: dth, dthg, oneovdth, bf1, bf2
    real :: twopisds_fac, fieldscale1, fieldscale2, inv_sds_power

    ! spectral coordinates (o) and (p)
    real, allocatable :: f(:), th(:), cth(:), sth(:)
    real, allocatable :: cth2(:), sth2(:) ! cos(2 th), sin(2 th)

    ! wave kinematics and derived factors (o,i)
    real, allocatable :: k(:,:), cp0(:,:), cg0(:,:) ! cp0 and cg0 include halo
    real, allocatable :: invcp0(:,:), kdk(:,:), k3dk(:,:), k4(:,:), oneoverk4(:,:)
    real, allocatable :: fkovg(:,:), cothkd(:,:), logl2overz(:,:), psiml2(:,:)
    real, allocatable :: sbf(:,:), sdv(:,:), snl_arg(:,:), bf1_renorm(:,:), bf2_renorm(:,:)
    real, allocatable :: dcp_dx(:,:), dcp_dy(:,:) ! halved phase speed gradients

    ! wave spectrum (o,p,i); ef is the spectrum forward in time
    real, allocatable :: e(:,:,:), ef(:,:,:)

    ! source functions (o,p,i); sds excludes the plunging breaker factor cothkd.
    ! ssin and snl are kept only for output; the time step evaluates them per cell.
    real, allocatable :: ssin(:,:,:), sds(:,:,:), snl(:,:,:)
    real, allocatable :: sdt(:,:), sice(:,:)     ! (o,i)
    real, allocatable :: mss_moments(:,:,:)      ! (o,3,i) mean-squared slope moments
    integer, allocatable :: oc(:)            ! highest prognostic frequency bin
    real, allocatable :: shelt(:), physics_time_step(:), psim(:)

    ! refraction
    integer, allocatable :: oc_rot(:) ! oc at the last rotation rate update
    real, allocatable :: curl(:)      ! halved current vorticity

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
    procedure, private :: sources
    procedure, private :: wind_input
    procedure, private :: integrate_column
    procedure, private :: advect_column
    procedure, private :: update_rotation
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
             self % sbf(om,is:ie), self % sdv(om,is:ie), self % snl_arg(om,is:ie), &
             self % bf1_renorm(om,is:ie), self % bf2_renorm(om,is:ie), &
             self % dcp_dx(om,is:ie), self % dcp_dy(om,is:ie), &
             self % psiml2(om,is:ie), self % sdt(om,is:ie), self % sice(om,is:ie), &
             self % mss_moments(om,3,is:ie), source=0.)
    allocate(self % cp0(om,iis:iie), self % cg0(om,iis:iie), source=tiny(1.))

    allocate(self % e(om,pm,iis:iie), self % ef(om,pm,iis:iie), source=tiny(1.))
    allocate(self % ssin(om,pm,is:ie), self % sds(om,pm,is:ie), self % snl(om,pm,is:ie), &
             source=0.)

    allocate(self % oc(is:ie), self % oc_rot(is:ie), source=0)
    allocate(self % psim(self % grid % imm), source=0.)

    allocate(self % shelt(is:ie), self % physics_time_step(is:ie), self % curl(is:ie), &
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


  subroutine initialize_physics(self)
    ! Initializes spectral coordinates, wave kinematics and derived factors,
    ! time step limits, drag, friction velocity, and the Stokes drift.
    class(domain_type), intent(inout) :: self

    real, allocatable :: kd(:,:), l2(:,:), dwn(:,:)
    real :: dlnf, bf1, bf2, cfllim, mindelx
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
    self % oneovdth = 1. / self % dth
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

    ! diffusion values in 2 frequencies:
    bf1 = exp(-16 * dlnf * dlnf)
    bf2 = exp(-64 * dlnf * dlnf)
    self % bf1 = bf1 / (bf1 + bf2)
    self % bf2 = bf2 / (bf1 + bf2)

    ! wave numbers, phase speeds, and group velocities:
    if (nproc == 0) write(*,'(a)') 'umwm: dispersion: solving for dispersion relationship;'

    do concurrent (o = 1:om, i = is:ie)
      k(o,i) = real(wavenumber(real(f(o), rk), real(grid % d(i), rk), real(config % rhow0, rk), &
                               real(config % g, rk), real(config % sfct, rk)), kind(k))
    end do

    if (nproc == 0) write(*,'(a)') 'umwm: dispersion: dispersion relationship done;'

    allocate(kd(om,is:ie), l2(om,is:ie), dwn(om,is:ie))

    ! limit kd to avoid floating overflow in transcendental functions:
    do concurrent (o = 1:om, i = is:ie)
      kd(o,i) = min(k(o,i) * grid % d(i), 20.)
    end do

    do concurrent (o = 1:om, i = is:ie)
      cp0(o,i) = twopi * f(o) / k(o,i)
      cg0(o,i) = real(group_speed(real(k(o,i), rk), real(grid % d(i), rk), &
                                  real(config % rhow0, rk), real(config % g, rk), &
                                  real(config % sfct, rk)), kind(cg0))
    end do

    ! frequently used factors:
    do concurrent (o = 1:om, i = is:ie)
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

    ! renormalization factors for snl:
    do concurrent (o = 1:om - 2, i = is:ie)
      self % bf1_renorm(o,i) = config % snl_fac * self % bf1 * self % kdk(o+1,i) / self % kdk(o,i)
      self % bf2_renorm(o,i) = config % snl_fac * self % bf2 * self % kdk(o+2,i) / self % kdk(o,i)
      self % snl_arg(o,i) = 1. - (self % bf1_renorm(o,i) + self % bf2_renorm(o,i))
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

    ! phase speed gradients for the bottom-induced refraction:
    do i = is, ie
      self % dcp_dx(:,i) = 0.5 * (cp0(:,grid % ie(i)) - cp0(:,grid % iw(i))) * grid % oneovdx(i)
      self % dcp_dy(:,i) = 0.5 * (cp0(:,grid % in(i)) - cp0(:,grid % is(i))) * grid % oneovdy(i)
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
                            'check xpl and ypl in the output namelist in namelists/main.nml')
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

      call self % step()

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
          write(*,fmt=100) self % sumt / config % dtg, self % dts, forcing % wspd(i), &
                           forcing % wdir(i), sigwaveheight(self, i), &
                           meanwaveperiod(self, i), self % cd(i) * 1e3, self % f(self % oc(i))
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
    logical :: store, currents, update_rotation
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
    ! the rotation rate of rays is updated for varying currents or first step;
    ! both always in coupled mode.
#ifdef ESMF
    store = .true.
    update_rotation = .true.
#else
    store = self % first .or. self % sumt >= config % dtg
    update_rotation = config % currents .or. self % first
#endif
    if (update_rotation) call self % update_rotation()
    self % dtr = min(self % dtr_max, self % dts)

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

      ! cut-off frequency (by default 4*pierson-moskowitz peak frequency) and its bin
      fcutoff = min(config % fcut_fac * config % g / wspd, config % fprog)
      do o = om - 2, 1, -1
        self % oc(i) = o
        if (fcutoff > self % f(o)) exit
      end do

      self % shelt(i) = sheltering_coare35(wspd)
      call self % wind_input(i, self % oc(i), ssin)

      ! breaking dissipation for the range used by the integration;
      ! the rest of the diagnostic range is evaluated for output only
      call mss_moments(om, pm, self % e(:,:,i), self % k3dk(:,i), self % cth2, self % sth2, &
                       self % dth, self % mss_moments(:,:,i))
      call sds_d12(om, pm, 1, min(self % oc(i) + 2, om), self % e(:,:,i), self % k4(:,i), &
                   self % f, self % cth2, self % sth2, config % mss_fac, self % twopisds_fac, &
                   self % mss_moments(:,:,i), config % sds_power, self % sds(:,:,i))

      ! dissipation due to turbulence
      self % sdt(:,i) = config % sdt_fac * sqrt(forcing % rhorat(i)) * self % ustar(i) * self % k(:,i)

      call s_ice(om, pm, self % e(:,:,i), self % kdk(:,i), self % cg0(:,i), self % dth, &
                 forcing % fice(i), config % fice_lth, self % sice(:,i))

      self % physics_time_step(i) = config % explim &
        / max_source_rate(om, pm, self % oc(i), ssin, self % sds(:,:,i), &
                          self % cothkd(:,i), self % snl_arg(:,i), self % sbf(:,i), &
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
      call sds_d12(om, pm, self % oc(i) + 3, om, self % e(:,:,i), self % k4(:,i), &
                   self % f, self % cth2, self % sth2, config % mss_fac, self % twopisds_fac, &
                   self % mss_moments(:,:,i), config % sds_power, self % sds(:,:,i))
    call snl_d12(om, pm, self % oc(i), self % e(:,:,i), self % bf1_renorm(:,i), &
                 self % bf2_renorm(:,i), config % snl_fac, self % sds(:,:,i), snl)
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
    ! Propagates and refracts the spectrum at cell i and integrates into ef,
    ! and evaluates the wind stress from the updated spectrum. Advection by
    ! ocean currents is included if currents is .true.
    class(domain_type), intent(inout) :: self
    integer, intent(in) :: i
    logical, intent(in) :: currents
#ifdef ESMF
    real :: e_propagated(self % spectrum % num_frequencies, self % spectrum % num_directions)
#endif

    associate(config => self % config, grid => self % grid, &
              om => self % spectrum % num_frequencies, pm => self % spectrum % num_directions, &
              e => self % e, ef => self % ef, cg0 => self % cg0, uc => self % forcing % uc, &
              vc => self % forcing % vc, ie => self % grid % ie(i), iw => self % grid % iw(i), &
              in => self % grid % in(i), is => self % grid % is(i))

    if (self % forcing % fice(i) > config % fice_uth) then
      ef(1:self % oc(i),:,i) = 0 ! no propagation through ice
    else

      call propagation(om, pm, self % oc(i), e(:,:,i), e(:,:,ie), e(:,:,iw), e(:,:,in), &
                       e(:,:,is), cg0(:,i), cg0(:,ie), cg0(:,iw), cg0(:,in), cg0(:,is), &
                       grid % cth_curv(:,i), grid % sth_curv(:,i), grid % dxn(i), grid % dxs(i), &
                       grid % dye(i), grid % dyw(i), currents, &
                       uc(i) + uc(grid % iie(i)), uc(i) + uc(grid % iiw(i)), &
                       vc(i) + vc(grid % iin(i)), vc(i) + vc(grid % iis(i)), &
                       0.25 * self % dts * grid % oneovar(i), ef(:,:,i))

#ifdef ESMF
      ! in coupled mode, refraction acts on the propagated spectrum
      e_propagated = ef(:,:,i)
      call refraction(om, pm, min(self % oc(i), self % oc_rot(i)), e_propagated, &
                      self % dcp_dx(:,i), self % dcp_dy(:,i), self % curl(i), &
                      self % cth, self % sth, self % dtr, self % oneovdth, ef(:,:,i))
#else
      call refraction(om, pm, min(self % oc(i), self % oc_rot(i)), e(:,:,i), &
                      self % dcp_dx(:,i), self % dcp_dy(:,i), self % curl(i), &
                      self % cth, self % sth, self % dtr, self % oneovdth, ef(:,:,i))
#endif

    end if

    end associate

    call self % wind_stress(i)

  end subroutine advect_column


  subroutine update_rotation(self)
    ! Updates the current-induced ray rotation and the refraction time step limit.
    class(domain_type), intent(inout) :: self
    real :: rmax
    integer :: i
#ifdef MPI
    real :: rmax_local
    integer :: ierr
#endif

    associate(grid => self % grid, uc => self % forcing % uc, vc => self % forcing % vc)

    rmax = 0
    do i = grid % istart, grid % iend
      self % curl(i) = 0.5 * ((vc(grid % iie(i)) - vc(grid % iiw(i))) * grid % oneovdx(i) &
                            - (uc(grid % iin(i)) - uc(grid % iis(i))) * grid % oneovdy(i))
      self % oc_rot(i) = self % oc(i)
      rmax = max(rmax, max_rotation_rate(self % spectrum % num_frequencies, &
                                         self % spectrum % num_directions, self % oc(i), &
                                         self % dcp_dx(:,i), self % dcp_dy(:,i), self % curl(i), &
                                         self % cth, self % sth))
    end do

#ifdef MPI
    rmax_local = rmax
    call mpi_allreduce(rmax_local, rmax, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif

    if (rmax <= tiny(rmax)) then
      self % dtr_max = self % dts
    else
      self % dtr_max = self % dth / rmax
    end if

    end associate

  end subroutine update_rotation


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
