module umwm_config

  use umwm_constants, only: rk

  implicit none

  private
  public :: config_type

  integer, parameter :: max_stokes_depths = 100
  integer, parameter :: allowedoutputtimes(10) = [-1, 0, 1, 2, 3, 4, 6, 8, 12, 24]

  type :: config_type
    logical :: read_ok = .false.
    character(len=256) :: path = ''
    character(len=256) :: read_error = ''

    logical :: isglobal = .false.
    integer :: mm = 0
    integer :: nm = 0
    integer :: om = 0
    integer :: pm = 0
    real :: fmin = 0.
    real :: fmax = 0.
    real :: fprog = 0.
    character(len=19) :: starttimestr = ''
    character(len=19) :: stoptimestr = ''
    character(len=19) :: reftimestr = '1970-01-01 00:00:00'
    real :: dtg = 0.
    logical :: restart = .false.

    real :: g = 0.
    real :: nu_air = 0.
    real :: nu_water = 0.
    real :: sfct = 0.
    real :: kappa = 0.
    real :: z = 0.
    real :: gustiness = 0.
    real :: dmin = 0.
    real :: explim = 0.
    real :: sin_fac = 0.
    real :: sin_diss1 = 0.
    real :: sin_diss2 = 0.
    real :: sds_fac = 0.
    real :: sds_power = 0.
    real :: mss_fac = 0.
    real :: snl_fac = 0.
    real :: sdt_fac = 0.
    real :: sbf_fac = 0.
    real :: sbp_fac = 0.

    logical :: gridfromfile = .false.
    real :: delx = 0.
    real :: dely = 0.
    logical :: topofromfile = .false.
    real :: dpt = 0.
    logical :: fillestuaries = .false.
    logical :: filllakes = .false.

    logical :: winds = .false.
    logical :: currents = .false.
    logical :: air_density = .false.
    logical :: water_density = .false.
    logical :: seaice = .false.

    real :: wspd0 = 0.
    real :: wdir0 = 0.
    real :: uc0 = 0.
    real :: vc0 = 0.
    real :: rhoa0 = 0.
    real :: rhow0 = 0.
    real :: fice0 = 0.
    real :: fice_lth = 0.
    real :: fice_uth = 0.

    integer :: outgrid = 0
    integer :: outspec = 0
    integer :: outrst = 0
    integer :: xpl = 0
    integer :: ypl = 0
    logical :: stokes = .false.

    real, allocatable :: stokes_depths(:)

  contains
    procedure :: validate
  end type config_type

  interface config_type
    module procedure :: config_type_cons
  end interface config_type

contains

  function config_type_cons(path) result(res)
    character(len=*), intent(in) :: path
    type(config_type) :: res

    integer :: unit, stat
    real :: depths(max_stokes_depths)
    character(len=256) :: iomsg

    res%path = path
    depths = -1.

    open(newunit=unit, file=path, status='old', form='formatted', &
         access='sequential', action='read', iostat=stat, iomsg=iomsg)
    if (stat /= 0) then
      res%read_error = 'could not open ' // trim(path) // ': ' // trim(iomsg)
      return
    end if

    call read_domain(unit, res%isglobal, res%mm, res%nm, res%om, res%pm, &
      res%fmin, res%fmax, res%fprog, res%starttimestr, res%stoptimestr, &
      res%reftimestr, res%dtg, res%restart, stat, iomsg)
    if (stat == 0) call read_physics(unit, res%g, res%nu_air, res%nu_water, &
      res%sfct, res%kappa, res%z, res%gustiness, res%dmin, res%explim, &
      res%sin_fac, res%sin_diss1, res%sin_diss2, res%sds_fac, &
      res%sds_power, res%mss_fac, res%snl_fac, res%sdt_fac, res%sbf_fac, &
      res%sbp_fac, stat, iomsg)
    if (stat == 0) call read_grid(unit, res%gridfromfile, res%delx, res%dely, &
      res%topofromfile, res%dpt, res%fillestuaries, res%filllakes, stat, iomsg)
    if (stat == 0) call read_forcing(unit, res%winds, res%currents, &
      res%air_density, res%water_density, res%seaice, stat, iomsg)
    if (stat == 0) call read_forcing_constant(unit, res%wspd0, res%wdir0, &
      res%uc0, res%vc0, res%rhoa0, res%rhow0, res%fice0, res%fice_lth, &
      res%fice_uth, stat, iomsg)
    if (stat == 0) call read_output(unit, res%outgrid, res%outspec, res%outrst, &
      res%xpl, res%ypl, res%stokes, stat, iomsg)
    close(unit)

    if (stat /= 0) then
      res%read_error = 'could not read all namelists from ' // trim(path) // &
                       ': ' // trim(iomsg)
      return
    end if

    call read_stokes_depths(path, depths, stat, iomsg)
    if (stat /= 0) then
      res%read_error = 'could not read &STOKES from ' // trim(path) // &
                       ': ' // trim(iomsg)
      return
    end if

    call set_stokes_depths(res, depths)
    res%read_ok = .true.

  end function config_type_cons


  subroutine read_domain(unit, isglobal, mm, nm, om, pm, fmin, fmax, fprog, &
                         starttimestr, stoptimestr, reftimestr, dtg, restart, &
                         stat, iomsg)
    integer, intent(in) :: unit
    logical, intent(inout) :: isglobal, restart
    integer, intent(inout) :: mm, nm, om, pm
    real, intent(inout) :: fmin, fmax, fprog, dtg
    character(len=*), intent(inout) :: starttimestr, stoptimestr, reftimestr
    integer, intent(out) :: stat
    character(len=*), intent(out) :: iomsg

    namelist /domain/ isglobal, mm, nm, om, pm, fmin, fmax, fprog, &
      starttimestr, stoptimestr, reftimestr, dtg, restart

    read(unit, nml=domain, iostat=stat, iomsg=iomsg)
  end subroutine read_domain


  subroutine read_physics(unit, g, nu_air, nu_water, sfct, kappa, z, &
                          gustiness, dmin, explim, sin_fac, sin_diss1, &
                          sin_diss2, sds_fac, sds_power, mss_fac, snl_fac, &
                          sdt_fac, sbf_fac, sbp_fac, stat, iomsg)
    integer, intent(in) :: unit
    real, intent(inout) :: g, nu_air, nu_water, sfct, kappa, z
    real, intent(inout) :: gustiness, dmin, explim, sin_fac, sin_diss1
    real, intent(inout) :: sin_diss2, sds_fac, sds_power, mss_fac, snl_fac
    real, intent(inout) :: sdt_fac, sbf_fac, sbp_fac
    integer, intent(out) :: stat
    character(len=*), intent(out) :: iomsg

    namelist /physics/ g, nu_air, nu_water, sfct, kappa, z, &
      gustiness, dmin, explim, sin_fac, sin_diss1, sin_diss2, &
      sds_fac, sds_power, mss_fac, snl_fac, sdt_fac, sbf_fac, sbp_fac

    read(unit, nml=physics, iostat=stat, iomsg=iomsg)
  end subroutine read_physics


  subroutine read_grid(unit, gridfromfile, delx, dely, topofromfile, dpt, &
                       fillestuaries, filllakes, stat, iomsg)
    integer, intent(in) :: unit
    logical, intent(inout) :: gridfromfile, topofromfile
    logical, intent(inout) :: fillestuaries, filllakes
    real, intent(inout) :: delx, dely, dpt
    integer, intent(out) :: stat
    character(len=*), intent(out) :: iomsg

    namelist /grid/ gridfromfile, delx, dely, topofromfile, dpt, &
      fillestuaries, filllakes

    read(unit, nml=grid, iostat=stat, iomsg=iomsg)
  end subroutine read_grid


  subroutine read_forcing(unit, winds, currents, air_density, water_density, &
                          seaice, stat, iomsg)
    integer, intent(in) :: unit
    logical, intent(inout) :: winds, currents, air_density, water_density, seaice
    integer, intent(out) :: stat
    character(len=*), intent(out) :: iomsg

    namelist /forcing/ winds, currents, air_density, water_density, seaice

    read(unit, nml=forcing, iostat=stat, iomsg=iomsg)
  end subroutine read_forcing


  subroutine read_forcing_constant(unit, wspd0, wdir0, uc0, vc0, rhoa0, rhow0, &
                                   fice0, fice_lth, fice_uth, stat, iomsg)
    integer, intent(in) :: unit
    real, intent(inout) :: wspd0, wdir0, uc0, vc0, rhoa0, rhow0
    real, intent(inout) :: fice0, fice_lth, fice_uth
    integer, intent(out) :: stat
    character(len=*), intent(out) :: iomsg

    namelist /forcing_constant/ wspd0, wdir0, uc0, vc0, rhoa0, rhow0, &
      fice0, fice_lth, fice_uth

    read(unit, nml=forcing_constant, iostat=stat, iomsg=iomsg)
  end subroutine read_forcing_constant


  subroutine read_output(unit, outgrid, outspec, outrst, xpl, ypl, stokes, &
                         stat, iomsg)
    integer, intent(in) :: unit
    integer, intent(inout) :: outgrid, outspec, outrst, xpl, ypl
    logical, intent(inout) :: stokes
    integer, intent(out) :: stat
    character(len=*), intent(out) :: iomsg

    namelist /output/ outgrid, outspec, outrst, xpl, ypl, stokes

    read(unit, nml=output, iostat=stat, iomsg=iomsg)
  end subroutine read_output


  subroutine read_stokes_depths(path, depths, stat, iomsg)
    character(len=*), intent(in) :: path
    real, intent(out) :: depths(:)
    integer, intent(out) :: stat
    character(len=*), intent(out) :: iomsg

    integer :: unit

    namelist /stokes/ depths

    depths = -1.
    open(newunit=unit, file=path, status='old', form='formatted', &
         access='sequential', action='read', iostat=stat, iomsg=iomsg)
    if (stat /= 0) return

    read(unit, nml=stokes, iostat=stat, iomsg=iomsg)
    close(unit)
  end subroutine read_stokes_depths


  subroutine set_stokes_depths(config, depths)
    type(config_type), intent(inout) :: config
    real, intent(in) :: depths(:)
    integer :: n, count_positive

    count_positive = count(depths > 0.)
    if (allocated(config%stokes_depths)) deallocate(config%stokes_depths)
    allocate(config%stokes_depths(count_positive))

    count_positive = 0
    do n = 1, size(depths)
      if (depths(n) > 0.) then
        count_positive = count_positive + 1
        config%stokes_depths(count_positive) = depths(n)
      end if
    end do
  end subroutine set_stokes_depths


  function validate(self, stop_on_error, rank) result(ok)
    class(config_type), intent(inout) :: self
    logical, intent(in), optional :: stop_on_error
    integer, intent(in), optional :: rank
    logical :: ok

    logical :: should_stop
    integer :: rank_value
    integer :: year, month, day, hour, minute, second

    ok = .true.
    should_stop = .true.
    if (present(stop_on_error)) should_stop = stop_on_error
    rank_value = 0
    if (present(rank)) rank_value = rank

    if (.not. self%read_ok) then
      call diagnostic('error', trim(self%read_error), ok, rank_value)
    end if

    if (self%mm < 3 .or. self%nm < 3) &
      call diagnostic('error', &
        'bad value in main.nml: mm and nm must be > 2', ok, rank_value)

    if (self%om < 3) &
      call diagnostic('error', &
        'bad value in main.nml: om must be >= 3', ok, rank_value)

    if (mod(self%pm, 4) /= 0) then
      call diagnostic('error', &
        'bad value in main.nml: pm must be divisible by 4', ok, rank_value)
    elseif (mod(self%pm, 8) /= 0) then
      call diagnostic('warning', &
        'pm should be divisible by 8 for optimal propagation properties', &
        ok, rank_value)
    end if

    if (self%fmin <= 0 .or. self%fmax <= 0 .or. self%fprog <= 0) &
      call diagnostic('error', &
        'bad value in main.nml: fmin, fmax and fprog must be > 0', ok, rank_value)

    if (self%fmin >= self%fmax) &
      call diagnostic('error', &
        'bad value in main.nml: fmin must be < fmax', ok, rank_value)

    if (self%fprog > self%fmax) then
      call diagnostic('warning', &
        'bad value in main.nml: fprog must be <= fmax; using highest allowed value', &
        ok, rank_value)
      self%fprog = self%fmax
    end if

    if (self%dtg <= 0) &
      call diagnostic('error', &
        'bad value in main.nml: dtg must be > 0', ok, rank_value)

    if (.not. valid_timestamp(self%starttimestr, year, month, day, hour, minute, second)) &
      call diagnostic('error', &
        'bad value in main.nml: startTimeStr must be YYYY-MM-DD HH:MM:SS', &
        ok, rank_value)

    if (.not. valid_timestamp(self%stoptimestr, year, month, day, hour, minute, second)) &
      call diagnostic('error', &
        'bad value in main.nml: stopTimeStr must be YYYY-MM-DD HH:MM:SS', &
        ok, rank_value)

    if (.not. valid_timestamp(self%reftimestr, year, month, day, hour, minute, second)) &
      call diagnostic('error', &
        'bad value in main.nml: refTimeStr must be YYYY-MM-DD HH:MM:SS', &
        ok, rank_value)

    if (valid_timestamp(self%starttimestr, year, month, day, hour, minute, second) .and. &
        valid_timestamp(self%stoptimestr, year, month, day, hour, minute, second)) then
      if (self%stoptimestr <= self%starttimestr) &
        call diagnostic('error', &
          'bad value in main.nml: stopTimeStr must be after startTimeStr', &
          ok, rank_value)
    end if

    if (self%z <= 0) &
      call diagnostic('error', &
        'bad value in main.nml: z must be > 0', ok, rank_value)

    if (self%gustiness < 0) then
      call diagnostic('error', &
        'bad value in main.nml: gustiness must be positive', ok, rank_value)
    elseif (self%gustiness > 0.2) then
      call diagnostic('warning', &
        '(bad) value in main.nml: gustiness > 0.2; proceed with caution', &
        ok, rank_value)
    end if

    if (self%dmin <= 0) &
      call diagnostic('error', &
        'bad value in main.nml: dmin must be > 0', ok, rank_value)

    if (.not. self%gridfromfile) then
      if (self%delx <= 0 .or. self%dely <= 0) &
        call diagnostic('error', &
          'bad value in main.nml: delx and dely must be > 0', ok, rank_value)
    end if

    if (.not. self%topofromfile) then
      if (self%dpt <= 0) &
        call diagnostic('error', &
          'bad value in main.nml: dpt must be > 0', ok, rank_value)
    end if

    if (self%rhoa0 <= 0 .or. self%rhow0 <= 0) &
      call diagnostic('error', &
        'bad value in main.nml: rhoa0 and rhow0 must be > 0', ok, rank_value)

    if (self%fice0 < 0 .or. self%fice0 > 1 .or. &
        self%fice_lth < 0 .or. self%fice_lth > 1 .or. &
        self%fice_uth < 0 .or. self%fice_uth > 1) &
      call diagnostic('error', &
        'bad value in main.nml: sea-ice fractions must be in [0, 1]', ok, rank_value)

    if (self%fice_lth > self%fice_uth) &
      call diagnostic('error', &
        'bad value in main.nml: fice_lth must be <= fice_uth', ok, rank_value)

    if (.not. any(self%outgrid == allowedoutputtimes)) &
      call diagnostic('error', &
        'bad value in main.nml: outgrid must be -1, 0, 1, 2, 3, 4, 6, 8, 12 or 24', &
        ok, rank_value)

    if (.not. any(self%outspec == allowedoutputtimes)) &
      call diagnostic('error', &
        'bad value in main.nml: outspec must be -1, 0, 1, 2, 3, 4, 6, 8, 12 or 24', &
        ok, rank_value)

    if (.not. any(self%outrst == allowedoutputtimes)) &
      call diagnostic('error', &
        'bad value in main.nml: outrst must be -1, 0, 1, 2, 3, 4, 6, 8, 12 or 24', &
        ok, rank_value)

    if (self%xpl < 1 .or. self%xpl > self%mm .or. &
        self%ypl < 1 .or. self%ypl > self%nm) &
      call diagnostic('error', &
        'bad value in main.nml: xpl and ypl must be inside the grid', ok, rank_value)

    if (.not. allocated(self%stokes_depths)) then
      call diagnostic('error', &
        'bad value in main.nml: at least one positive Stokes depth is required', &
        ok, rank_value)
    elseif (size(self%stokes_depths) < 1) then
      call diagnostic('error', &
        'bad value in main.nml: at least one positive Stokes depth is required', &
        ok, rank_value)
    end if

    if (.not. ok .and. should_stop) error stop 1

  end function validate


  subroutine diagnostic(level, message, ok, rank)
    character(len=*), intent(in) :: level
    character(len=*), intent(in) :: message
    logical, intent(inout) :: ok
    integer, intent(in) :: rank

    if (level == 'error') ok = .false.
    if (rank == 0) write(*, '(a)') 'umwm: config: ' // trim(level) // ': ' // trim(message)
  end subroutine diagnostic


  logical function valid_timestamp(value, year, month, day, hour, minute, second) result(ok)
    character(len=*), intent(in) :: value
    integer, intent(out) :: year, month, day, hour, minute, second
    integer :: stat

    ok = .false.
    year = 0
    month = 0
    day = 0
    hour = 0
    minute = 0
    second = 0

    if (len_trim(value) /= 19) return
    if (value(5:5) /= '-' .or. value(8:8) /= '-' .or. &
        value(11:11) /= ' ' .or. value(14:14) /= ':' .or. &
        value(17:17) /= ':') return

    read(value(1:4), *, iostat=stat) year
    if (stat /= 0) return
    read(value(6:7), *, iostat=stat) month
    if (stat /= 0) return
    read(value(9:10), *, iostat=stat) day
    if (stat /= 0) return
    read(value(12:13), *, iostat=stat) hour
    if (stat /= 0) return
    read(value(15:16), *, iostat=stat) minute
    if (stat /= 0) return
    read(value(18:19), *, iostat=stat) second
    if (stat /= 0) return

    if (year < 1) return
    if (month < 1 .or. month > 12) return
    if (day < 1 .or. day > days_in_month(year, month)) return
    if (hour < 0 .or. hour > 23) return
    if (minute < 0 .or. minute > 59) return
    if (second < 0 .or. second > 59) return

    ok = .true.
  end function valid_timestamp


  integer function days_in_month(year, month) result(days)
    integer, intent(in) :: year
    integer, intent(in) :: month
    integer, parameter :: days_per_month(12) = &
      [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    days = days_per_month(month)
    if (month == 2 .and. is_leap_year(year)) days = 29
  end function days_in_month


  logical function is_leap_year(year) result(res)
    integer, intent(in) :: year

    res = (mod(year, 4) == 0 .and. mod(year, 100) /= 0) .or. mod(year, 400) == 0
  end function is_leap_year

end module umwm_config
