module umwm_forcing

  use netcdf
  use umwm_config, only: config_type
  use umwm_grid, only: grid_type
  use umwm_module, only: sumt

  implicit none

  private
  public :: forcing_type

  type :: forcing_type
    logical :: readfile = .false.

    ! Native grid snapshots and interpolation scratch fields.
    real, dimension(:,:), allocatable :: gustu, gustv
    real, dimension(:,:), allocatable :: rhoa_2d, rhow_2d
    real, dimension(:,:), allocatable :: wspd_2d, wdir_2d
    real, dimension(:,:), allocatable :: fice_2d, ficeb, ficef
    real, dimension(:,:), allocatable :: uwb, vwb, uw, vw, uwf, vwf
    real, dimension(:,:), allocatable :: ucb, uc_2d, ucf, vcb, vc_2d, vcf

    ! Remapped forcing fields used by the model.
    real, dimension(:), allocatable :: rhoab, rhoa, rhoaf
    real, dimension(:), allocatable :: rhowb, rhow, rhowf
    real, dimension(:), allocatable :: rhorat
    real, dimension(:), allocatable :: uc, vc
    real, dimension(:), allocatable :: wspd, wdir
    real, dimension(:), allocatable :: fice

  contains
    procedure :: initialize => forcing_initialize
    procedure :: finalize => forcing_finalize
    procedure :: load => forcing_load
    procedure :: update => forcing_update
    procedure :: interpolate => forcing_interpolate
  end type forcing_type

contains

  subroutine forcing_initialize(self, config, grid)
    class(forcing_type), intent(inout) :: self
    type(config_type), intent(in) :: config
    type(grid_type), intent(in) :: grid

    real :: uw0, vw0

    call self % finalize()

    self % readfile = any([config % winds, config % currents, config % air_density, &
                           config % water_density, config % seaice])

    allocate(self % gustu(grid % mm, grid % nm), self % gustv(grid % mm, grid % nm))
    allocate(self % rhoa_2d(grid % mm, grid % nm), self % rhow_2d(grid % mm, grid % nm))
    allocate(self % wspd_2d(grid % mm, grid % nm), self % wdir_2d(grid % mm, grid % nm))
    allocate(self % fice_2d(grid % mm, grid % nm), self % ficeb(grid % mm, grid % nm), &
             self % ficef(grid % mm, grid % nm))
    allocate(self % uw(grid % mm, grid % nm), self % uwb(grid % mm, grid % nm), &
             self % uwf(grid % mm, grid % nm))
    allocate(self % vw(grid % mm, grid % nm), self % vwb(grid % mm, grid % nm), &
             self % vwf(grid % mm, grid % nm))
    allocate(self % uc_2d(grid % mm, grid % nm), self % ucb(grid % mm, grid % nm), &
             self % ucf(grid % mm, grid % nm))
    allocate(self % vc_2d(grid % mm, grid % nm), self % vcb(grid % mm, grid % nm), &
             self % vcf(grid % mm, grid % nm))

    allocate(self % wspd(grid % imm), self % wdir(grid % imm))
    allocate(self % fice(grid % imm))
    allocate(self % uc(grid % imm), self % vc(grid % imm))
    allocate(self % rhoa(grid % imm), self % rhoab(grid % imm), self % rhoaf(grid % imm))
    allocate(self % rhow(grid % imm), self % rhowb(grid % imm), self % rhowf(grid % imm))
    allocate(self % rhorat(grid % imm))

    self % gustu = 0.0
    self % gustv = 0.0

    uw0 = config % wspd0 * cos(config % wdir0)
    vw0 = config % wspd0 * sin(config % wdir0)

    self % uw = uw0
    self % uwb = uw0
    self % uwf = uw0
    self % vw = vw0
    self % vwb = vw0
    self % vwf = vw0
    self % wspd_2d = config % wspd0
    self % wdir_2d = config % wdir0
    self % wspd = config % wspd0
    self % wdir = config % wdir0

    self % uc_2d = config % uc0
    self % ucb = config % uc0
    self % ucf = config % uc0
    self % vc_2d = config % vc0
    self % vcb = config % vc0
    self % vcf = config % vc0
    self % uc = config % uc0
    self % vc = config % vc0

    self % fice_2d = config % fice0
    self % ficeb = config % fice0
    self % ficef = config % fice0
    self % fice = config % fice0

    self % rhoa_2d = config % rhoa0
    self % rhoab = config % rhoa0
    self % rhoaf = config % rhoa0
    self % rhoa = config % rhoa0

    self % rhow_2d = config % rhow0
    self % rhowb = config % rhow0
    self % rhowf = config % rhow0
    self % rhow = config % rhow0

    self % rhorat = self % rhoa / self % rhow

  end subroutine forcing_initialize


  subroutine forcing_finalize(self)
    class(forcing_type), intent(inout) :: self

    if (allocated(self % gustu)) deallocate(self % gustu, self % gustv)
    if (allocated(self % rhoa_2d)) deallocate(self % rhoa_2d, self % rhow_2d)
    if (allocated(self % wspd_2d)) deallocate(self % wspd_2d, self % wdir_2d)
    if (allocated(self % fice_2d)) deallocate(self % fice_2d, self % ficeb, self % ficef)
    if (allocated(self % uw)) deallocate(self % uw, self % uwb, self % uwf)
    if (allocated(self % vw)) deallocate(self % vw, self % vwb, self % vwf)
    if (allocated(self % uc_2d)) deallocate(self % uc_2d, self % ucb, self % ucf)
    if (allocated(self % vc_2d)) deallocate(self % vc_2d, self % vcb, self % vcf)

    if (allocated(self % wspd)) deallocate(self % wspd, self % wdir)
    if (allocated(self % fice)) deallocate(self % fice)
    if (allocated(self % uc)) deallocate(self % uc, self % vc)
    if (allocated(self % rhoa)) deallocate(self % rhoa, self % rhoab, self % rhoaf)
    if (allocated(self % rhow)) deallocate(self % rhow, self % rhowb, self % rhowf)
    if (allocated(self % rhorat)) deallocate(self % rhorat)

    self % readfile = .false.

  end subroutine forcing_finalize


  subroutine forcing_update(self, config, timestr, grid)
    ! Advance forcing snapshots and load the next file-backed snapshot.
    class(forcing_type), intent(inout) :: self
    type(config_type), intent(in) :: config
    character(len=19), intent(in) :: timestr
    type(grid_type), intent(in) :: grid

    if (config % winds) then
      self % uwb = self % uwf
      self % vwb = self % vwf
    end if

    if (config % currents) then
      self % ucb = self % ucf
      self % vcb = self % vcf
    end if

    if (config % air_density) self % rhoab = self % rhoaf
    if (config % water_density) self % rhowb = self % rhowf
    if (config % seaice) self % ficeb = self % ficef

    self % readfile = any([config % winds, config % currents, config % air_density, &
                           config % water_density, config % seaice])
    if (self % readfile) call self % load(config, timestr, grid)

  end subroutine forcing_update


  subroutine forcing_load(self, config, timestr, grid)
    ! Load atmospheric and oceanic input fields for one forcing time level.
    class(forcing_type), intent(inout) :: self
    type(config_type), intent(in) :: config
    character(len=19), intent(in) :: timestr
    type(grid_type), intent(in) :: grid

    character(999) :: nc_infile
    character(19) :: readstr
    integer :: ncid, varid

    readstr = timestr
    readstr(11:11) = '_'

    nc_infile = 'input/umwmin_' // readstr // '.nc'
    self % readfile = any([config % winds, config % currents, config % air_density, &
                           config % water_density, config % seaice])

    if (self % readfile) call forcing_nc_check(nf90_open(trim(nc_infile), nf90_nowrite, ncid))

    if (config % winds) then
      call forcing_nc_check(nf90_inq_varid(ncid, 'uw', varid))
      call forcing_nc_check(nf90_get_var(ncid, varid, self % uwf))
      call forcing_nc_check(nf90_inq_varid(ncid, 'vw', varid))
      call forcing_nc_check(nf90_get_var(ncid, varid, self % vwf))
    else
      self % wspd_2d = config % wspd0
      self % wdir_2d = config % wdir0
      self % wspd = config % wspd0
      self % wdir = config % wdir0
    end if

    if (config % currents) then
      call forcing_nc_check(nf90_inq_varid(ncid, 'uc', varid))
      call forcing_nc_check(nf90_get_var(ncid, varid, self % ucf))
      call forcing_nc_check(nf90_inq_varid(ncid, 'vc', varid))
      call forcing_nc_check(nf90_get_var(ncid, varid, self % vcf))
    else
      self % ucf = config % uc0
      self % vcf = config % vc0
      self % uc = config % uc0
      self % vc = config % vc0
    end if

    if (config % seaice) then
      call forcing_nc_check(nf90_inq_varid(ncid, 'fice', varid))
      call forcing_nc_check(nf90_get_var(ncid, varid, self % ficef))
    else
      self % fice_2d = config % fice0
      self % fice = config % fice0
    end if

    where (grid % mask == 0)
      self % ucf = 0
      self % vcf = 0
    end where

    if (config % air_density) then
      call forcing_nc_check(nf90_inq_varid(ncid, 'rhoa', varid))
      call forcing_nc_check(nf90_get_var(ncid, varid, self % rhoa_2d))
    else
      self % rhoa_2d = config % rhoa0
      self % rhoa = config % rhoa0
    end if

    if (config % water_density) then
      call forcing_nc_check(nf90_inq_varid(ncid, 'rhow', varid))
      call forcing_nc_check(nf90_get_var(ncid, varid, self % rhow_2d))
    else
      self % rhow_2d = config % rhow0
      self % rhow = config % rhow0
    end if

    if (self % readfile) call forcing_nc_check(nf90_close(ncid))

    self % rhoaf = grid % remap_mn2i(self % rhoa_2d)
    self % rhowf = grid % remap_mn2i(self % rhow_2d)
    self % rhorat = self % rhoa / self % rhow

  end subroutine forcing_load


  subroutine forcing_interpolate(self, config, grid)
    ! Interpolate atmospheric and oceanic forcing fields in time.
    class(forcing_type), intent(inout) :: self
    type(config_type), intent(in) :: config
    type(grid_type), intent(in) :: grid

    real :: alpha

    alpha = sumt / config % dtg

    if (config % winds) then

      self % uw = self % uwb * (1 - alpha) + self % uwf * alpha
      self % vw = self % vwb * (1 - alpha) + self % vwf * alpha

      if (config % gustiness > 0) then
        call random_number(self % gustu)
        call random_number(self % gustv)

        self % gustu = config % gustiness * (2 * self % gustu - 1)
        self % gustv = config % gustiness * (2 * self % gustv - 1)

        self % uw = self % uw * (1 + self % gustu)
        self % vw = self % vw * (1 + self % gustv)
      end if

      self % wspd_2d = sqrt(self % uw**2 + self % vw**2)
      self % wdir_2d = atan2(self % vw, self % uw)

      self % wspd = grid % remap_mn2i(self % wspd_2d)
      self % wdir = grid % remap_mn2i(self % wdir_2d)

    end if

    if (config % seaice) then
      self % fice_2d = self % ficeb * (1 - alpha) + self % ficef * alpha
      self % fice = grid % remap_mn2i(self % fice_2d)
    end if

    if (config % currents) then
      self % uc_2d = self % ucb * (1 - alpha) + self % ucf * alpha
      self % vc_2d = self % vcb * (1 - alpha) + self % vcf * alpha

      self % uc = grid % remap_mn2i(self % uc_2d)
      self % vc = grid % remap_mn2i(self % vc_2d)
    end if

    if (config % air_density) self % rhoa = self % rhoab * (1 - alpha) + self % rhoaf * alpha
    if (config % water_density) self % rhow = self % rhowb * (1 - alpha) + self % rhowf * alpha

    self % rhorat = self % rhoa / self % rhow

  end subroutine forcing_interpolate


  subroutine forcing_nc_check(stat)
    integer, intent(in) :: stat

    if (stat /= nf90_noerr) then
      write(*,*) 'error in netcdf forcing input'
      write(*,*) trim(nf90_strerror(stat))
      stop
    end if

  end subroutine forcing_nc_check

end module umwm_forcing
