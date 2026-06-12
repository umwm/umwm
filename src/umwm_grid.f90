module umwm_grid
! Horizontal grid geometry, masks, remapping, and partition metadata.

#ifdef MPI
use mpi
#endif

use umwm_config, only: config_type
use umwm_constants, only: dr, pi, r_earth, twopi
use umwm_spectrum, only: spectrum_type

implicit none

private
public :: grid_type

type :: grid_type
  integer :: mm = 0
  integer :: nm = 0
  integer :: im = 0
  integer :: imm = 0

  integer :: istart = 0
  integer :: iend = 0
  integer :: iistart = 0
  integer :: iiend = 0

  integer :: first_col_len = 0
  integer :: last_col_len = 0
  integer :: ilen = 0
  integer :: im_mod = 0
  character :: remap_dir = 'v'
  integer :: filled_estuary_cells = -1
  integer :: lake_fill_events = 0

  integer, allocatable :: iw(:), ie(:), is(:), in(:)
  integer, allocatable :: iiw(:), iie(:), iis(:), iin(:)
  integer, allocatable :: i_exchange_indices(:)
  integer, allocatable :: mi(:), ni(:)
  integer, allocatable :: ii(:,:), mask(:,:), nproc_out(:,:)
  integer, allocatable :: lake_fill_m(:), lake_fill_n(:), lake_fill_count(:)

  integer, allocatable :: istart_all(:), iend_all(:), ilen_all(:)
  integer, allocatable :: iistart_all(:), iiend_all(:)

  real, allocatable :: cth_curv(:,:), sth_curv(:,:)

  real, allocatable :: ar_2d(:,:), curv(:,:)
  real, allocatable :: d_2d(:,:), dlon(:,:), dlat(:,:)
  real, allocatable :: dx_2d(:,:), dy_2d(:,:)
  real, allocatable :: lat(:,:), lon(:,:), x(:,:), y(:,:)

  real, allocatable :: ar(:), d(:), dx(:), dy(:)
  real, allocatable :: dxn(:), dxs(:), dyw(:), dye(:)
  real, allocatable :: oneovar(:), oneovdx(:), oneovdy(:)

contains
  procedure :: initialize
  procedure :: initialize_direction_projection
  procedure :: print_diagnostics
  procedure :: finalize
  procedure :: remap_i2mn
  procedure :: remap_mn2i
  procedure, private :: allocate_native
  procedure, private :: define_grid
  procedure, private :: define_masks
  procedure, private :: fill
  procedure, private :: partition
  procedure, private :: allocate_remapped
  procedure, private :: remap
end type grid_type

contains

subroutine initialize(self, config)
  class(grid_type), intent(inout) :: self
  type(config_type), intent(in) :: config

  call self % allocate_native(config)
  call self % define_grid(config)
  call self % define_masks(config)
  call self % partition(config)
  call self % allocate_remapped()
  call self % remap(config)

end subroutine initialize


subroutine print_diagnostics(self)
  use umwm_module, only: nproc
#ifdef MPI
  use umwm_module, only: mpisize
#endif

  class(grid_type), intent(in) :: self
  integer :: event
#ifdef MPI
  integer :: nn
#endif

  if (nproc /= 0) return

  write(*,fmt=101) 'umwm: grid: dx min/max/mean [m]:     ', &
                  minval(self % dx_2d), maxval(self % dx_2d), &
                  sum(self % dx_2d) / (self % mm * self % nm)
  write(*,fmt=101) 'umwm: grid: dy min/max/mean [m]:     ', &
                  minval(self % dy_2d), maxval(self % dy_2d), &
                  sum(self % dy_2d) / (self % mm * self % nm)
  write(*,fmt=101) 'umwm: grid: area min/max/mean [m^2]: ', &
                  minval(self % ar_2d), maxval(self % ar_2d), &
                  sum(self % ar_2d) / (self % mm * self % nm)

  if (self % filled_estuary_cells >= 0) then
    write(*,fmt=201) 'umwm: masks: filled cells with 3-land neighbours,', &
                    self % filled_estuary_cells, ' cells total.'
  end if

  do event = 1, self % lake_fill_events
    write(*,fmt=202) 'umwm: masks: filled closed sea at i,j:', &
                    self % lake_fill_m(event), self % lake_fill_n(event), &
                    ', ', self % lake_fill_count(event), ' cells total.'
  end do

  write(unit=*,fmt=302) self % imm
  write(unit=*,fmt=303) self % imm - self % im, &
                        float(self % imm - self % im) / float(self % imm) * 100.
  write(unit=*,fmt=304) self % im, float(self % im) / float(self % imm) * 100.
  write(unit=*,fmt=305) minval(self % d_2d, self % mask == 1)
  write(unit=*,fmt=306) maxval(self % d_2d, self % mask == 1)

#ifdef MPI
  if (allocated(self % ilen_all)) then
    write(*,fmt='(a)') 'umwm: partition: tiling summary:'
    write(*,fmt='(a)') '+---------+------------+----------------+----------------+'
    write(*,fmt='(a)') '|  nproc  |    ilen    |     istart     |      iend      |'
    write(*,fmt='(a)') '+---------+------------+----------------+----------------+'
    do nn = 0, mpisize - 1
      write(*,fmt=307) nn, self % ilen_all(nn), self % istart_all(nn), self % iend_all(nn)
    end do

    write(*,fmt='(a)') '+---------+------------+----------------+----------------+'
  end if

  if (allocated(self % iistart_all)) then
    write(*,fmt='(a)') 'umwm: remap: tiling with halo summary:'
    write(*,fmt='(a)') '+---------+------------+----------------+----------------+'
    write(*,fmt='(a)') '|  nproc  |   iilen    |    iistart     |     iiend      |'
    write(*,fmt='(a)') '+---------+------------+----------------+----------------+'
    do nn = 0, mpisize - 1
      write(*,fmt=310) nn, self % iiend_all(nn) - self % iistart_all(nn) + 1, &
                        self % iistart_all(nn), self % iiend_all(nn)
    end do

    write(*,fmt='(a)') '+---------+------------+----------------+----------------+'
  end if
#endif

  101 format(a,3(f15.2,1x))
  201 format(a,i8,a)
  202 format(a,2(i5,1x),a,i8,a)
  302 format('umwm: masks: total number of grid points: ',i9)
  303 format('umwm: masks: number of land points:       ',i9,', ',f5.1,'%')
  304 format('umwm: masks: number of sea points:        ',i9,', ',f5.1,'%')
  305 format('umwm: masks: shallowest point:            ',f9.3,' meters')
  306 format('umwm: masks: deepest point:               ',f9.3,' meters')
#ifdef MPI
  307 format('| ',i7,' | ',i10,' | ',2(i14,' | '))
  310 format('| ',i7,' | ',i10,' | ',2(i14,' | '))
#endif

end subroutine print_diagnostics


subroutine allocate_native(self, config)
  class(grid_type), intent(inout) :: self
  type(config_type), intent(in) :: config

  self % mm = config % mm
  self % nm = config % nm
  self % filled_estuary_cells = -1
  self % lake_fill_events = 0

  allocate(self % ar_2d(self % mm,self % nm))
  allocate(self % curv(self % mm,self % nm))
  allocate(self % d_2d(self % mm,self % nm), self % dx_2d(self % mm,self % nm), &
           self % dy_2d(self % mm,self % nm))
  allocate(self % dlon(self % mm,self % nm), self % dlat(self % mm,self % nm))
  allocate(self % ii(self % mm,self % nm))
  allocate(self % lat(self % mm,self % nm), self % lon(self % mm,self % nm))
  allocate(self % x(self % mm,self % nm), self % y(self % mm,self % nm))
  allocate(self % mask(self % mm,self % nm))
  allocate(self % nproc_out(self % mm,self % nm))
  if (allocated(self % lake_fill_m)) deallocate(self % lake_fill_m)
  if (allocated(self % lake_fill_n)) deallocate(self % lake_fill_n)
  if (allocated(self % lake_fill_count)) deallocate(self % lake_fill_count)
  allocate(self % lake_fill_m(0), self % lake_fill_n(0), self % lake_fill_count(0))
  self % nproc_out = 0

end subroutine allocate_native


subroutine allocate_remapped(self)
  class(grid_type), intent(inout) :: self

  allocate(self % iw(self % imm), self % ie(self % imm), &
           self % is(self % imm), self % in(self % imm))
  allocate(self % iiw(self % imm), self % iie(self % imm), &
           self % iis(self % imm), self % iin(self % imm))
  allocate(self % mi(self % imm), self % ni(self % imm))

  allocate(self % ar(self % istart:self % iend))
  allocate(self % d(self % imm), self % dx(self % imm), self % dy(self % imm))
  allocate(self % dxs(self % istart:self % iend), self % dxn(self % istart:self % iend))
  allocate(self % dyw(self % istart:self % iend), self % dye(self % istart:self % iend))
  allocate(self % oneovar(self % istart:self % iend), &
           self % oneovdx(self % istart:self % iend), &
           self % oneovdy(self % istart:self % iend))

end subroutine allocate_remapped


subroutine define_grid(self, config)
! Defines grid spacing and grid cell areas.
  use netcdf
  use umwm_module, only: nproc
  use umwm_util, only: distance_haversine, raiseexception

  class(grid_type), intent(inout) :: self
  type(config_type), intent(in) :: config
  logical :: loniscontinuous = .true.

  integer :: m, n
  integer :: ncid, varid, stat

  real, allocatable :: abscoslat(:,:), lon_tmp(:,:), rotx(:,:), roty(:,:)
  real, allocatable :: rlon(:,:), rlat(:,:)

  if (config % gridfromfile) then

    stat = nf90_open('input/umwm.gridtopo', nf90_nowrite, ncid)

    if (stat /= 0) then

      if (nproc == 0) then
        call raiseexception('warning', 'grid', &
                            'input/umwm.gridtopo not found; trying input/umwm.grid')
      end if

      stat = nf90_open('input/umwm.grid', nf90_nowrite, ncid)

      if (stat /= 0) then

        if (nproc == 0) then
          call raiseexception('abort', 'grid', &
                              'input/umwm.grid not found. make sure input files are in place')
        end if

        stop

      end if

    end if

    call grid_nc_check(nf90_inq_varid(ncid, 'lon', varid))
    call grid_nc_check(nf90_get_var(ncid, varid, self % lon))
    call grid_nc_check(nf90_inq_varid(ncid, 'lat', varid))
    call grid_nc_check(nf90_get_var(ncid, varid, self % lat))
    call grid_nc_check(nf90_close(ncid))

    allocate(abscoslat(self % mm,self % nm), lon_tmp(self % mm,self % nm), &
             rotx(self % mm,self % nm), roty(self % mm,self % nm))
    allocate(rlon(self % mm,self % nm), rlat(self % mm,self % nm))

    ! figure out if longitude field is continuous:
    if (minval(self % lon) < -175 .and. maxval(self % lon) > 175) loniscontinuous = .false.

    lon_tmp = self % lon

    ! store original lon array before modifying:
    if (.not. loniscontinuous) where(self % lon < 0) self % lon = self % lon + 360

    do n = 1, self % nm
      do m = 2, self % mm - 1
        self % dlon(m,n) = 0.5 * (self % lon(m+1,n) - self % lon(m-1,n))
      end do
    end do

    self % dlon(1,:) = 2 * self % dlon(2,:) - self % dlon(3,:)
    self % dlon(self % mm,:) = 2 * self % dlon(self % mm-1,:) - self % dlon(self % mm-2,:)

    do n = 2, self % nm - 1
      do m = 1, self % mm
        self % dlat(m,n) = 0.5 * (self % lat(m,n+1) - self % lat(m,n-1))
      end do
    end do
    self % dlat(:,1) = 2 * self % dlat(:,2) - self % dlat(:,3)
    self % dlat(:,self % nm) = 2 * self % dlat(:,self % nm-1) - self % dlat(:,self % nm-2)

    abscoslat = abs(cos(dr * self % lat))

    ! revert to original lon array:
    self % lon = lon_tmp

    rlon = self % lon * twopi / 360.
    rlat = self % lat * twopi / 360.

    do n = 1, self % nm
      do m = 2, self % mm - 1
        self % dx_2d(m,n) = r_earth * distance_haversine(0.5 * (rlon(m-1,n) + rlon(m  ,n)), &
                                                         0.5 * (rlon(m  ,n) + rlon(m+1,n)), &
                                                         0.5 * (rlat(m-1,n) + rlat(m  ,n)), &
                                                         0.5 * (rlat(m  ,n) + rlat(m+1,n)))
      end do
    end do

    self % dx_2d(1,:) = 2 * self % dx_2d(2,:) - self % dx_2d(3,:)
    self % dx_2d(self % mm,:) = 2 * self % dx_2d(self % mm-1,:) - self % dx_2d(self % mm-2,:)

    do n = 2, self % nm - 1
      do m = 1, self % mm
        self % dy_2d(m,n) = r_earth * distance_haversine(0.5 * (rlon(m,n-1) + rlon(m,n  )), &
                                                         0.5 * (rlon(m,n  ) + rlon(m,n+1)), &
                                                         0.5 * (rlat(m,n-1) + rlat(m,n  )), &
                                                         0.5 * (rlat(m,n  ) + rlat(m,n+1)))
      end do
    end do

    self % dy_2d(:,1) = 2 * self % dy_2d(:,2) - self % dy_2d(:,3)
    self % dy_2d(:,self % nm) = 2 * self % dy_2d(:,self % nm-1) - self % dy_2d(:,self % nm-2)

    ! compute grid rotation for great circle propagation
    self % curv = 0
    do n = 1, self % nm
      do m = 2, self % mm - 1
        self % curv(m,n) = atan2(sin(rlon(m+1,n) - rlon(m-1,n)) * cos(rlat(m+1,n)), &
                                 cos(rlat(m-1,n)) * sin(rlat(m+1,n)) &
                                -sin(rlat(m-1,n)) * cos(rlat(m+1,n)) &
                                * cos(rlon(m+1,n) - rlon(m-1,n)))
      end do
    end do

    if (config % isglobal) then
      m = 1
      self % curv(m,:) = atan2(sin(rlon(m+1,:) - rlon(self % mm,:)) * cos(rlat(m+1,:)), &
                               cos(rlat(self % mm,:)) * sin(rlat(m+1,:)) &
                              -sin(rlat(self % mm,:)) * cos(rlat(m+1,:)) &
                              * cos(rlon(m+1,:) - rlon(self % mm,:)))
      m = self % mm
      self % curv(m,:) = atan2(sin(rlon(1,:) - rlon(m-1,:)) * cos(rlat(1,:)), &
                               cos(rlat(m-1,:)) * sin(rlat(1,:)) &
                              -sin(rlat(m-1,:)) * cos(rlat(1,:)) &
                              * cos(rlon(1,:) - rlon(m-1,:)))
    else
      self % curv(1,:) = self % curv(2,:)
      self % curv(self % mm,:) = self % curv(self % mm-1,:)
    end if

    self % curv = self % curv - 0.5 * pi

    deallocate(abscoslat, lon_tmp, rotx, roty, rlon, rlat)

  else ! use constant value from namelist

    self % dx_2d = config % delx
    self % dy_2d = config % dely

    self % curv = 0

    self % lon = 0
    self % lat = 0
    self % dlon = 0
    self % dlat = 0

    self % x(1,:) = 0
    do m = 2, self % mm
      self % x(m,:) = self % x(m-1,:) + 0.5 * (self % dx_2d(m-1,:) + self % dx_2d(m,:))
    end do

    self % y(:,1) = 0
    do n = 2, self % nm
      self % y(:,n) = self % y(:,n-1) + 0.5 * (self % dy_2d(:,n-1) + self % dy_2d(:,n))
    end do

  end if

  if (config % topofromfile) then ! read depth field from file

    call grid_nc_check(nf90_open('input/umwm.gridtopo', nf90_nowrite, ncid))
    call grid_nc_check(nf90_inq_varid(ncid, 'z', varid))
    call grid_nc_check(nf90_get_var(ncid, varid, self % d_2d))
    call grid_nc_check(nf90_close(ncid))

  else ! use constant value from namelist

    self % d_2d = config % dpt

  end if

  ! compute cell areas and reciprocals:
  self % ar_2d = self % dx_2d * self % dy_2d

end subroutine define_grid


subroutine define_masks(self, config)
! Defines landmasks and optionally removes one-cell wide estuaries and lakes.

  class(grid_type), intent(inout) :: self
  type(config_type), intent(in) :: config
  logical :: iterate

  integer :: m, n
  integer :: exm, exn
  integer :: cnt, fillcount

  ! set masks:
  if (config % topofromfile) then

    ! set initial seamask everywhere:
    self % mask = 1

    ! set up boundary points:
    self % mask(:,1) = 0
    self % mask(:,self % nm) = 0

    ! close e and w edges if limited area:
    if (.not. config % isglobal) then
      self % mask(1,:) = 0
      self % mask(self % mm,:) = 0
    end if

    ! set land mask where depth is non-negative, and then set depth to dmin:
    where(self % d_2d >= 0)
      self % mask = 0
      self % d_2d = config % dmin
    endwhere

    ! make depths positive and limit to dmin:
    self % d_2d = abs(self % d_2d)
    where(self % d_2d < config % dmin) self % d_2d = config % dmin

    ! fill estuaries and isolated sea points:
    if (config % fillestuaries) then
      fillcount = 0
      iterate = .true.
      do while (iterate)
        iterate = .false.
        do n = 2, self % nm - 1
          do m = 2, self % mm - 1

            if (self % mask(m,n) == 1) then

              cnt = 0
              if (self % mask(m-1,n) == 1) cnt = cnt + 1
              if (self % mask(m+1,n) == 1) cnt = cnt + 1
              if (self % mask(m,n-1) == 1) cnt = cnt + 1
              if (self % mask(m,n+1) == 1) cnt = cnt + 1

              if (cnt <= 1) then
                self % mask(m,n) = 0
                fillcount = fillcount + 1
                iterate = .true.
              end if

            end if

          end do
        end do
      end do
      self % filled_estuary_cells = fillcount
    end if

    ! discard lakes or unwanted closed basin from the domain:
    if (config % filllakes) then
      open(unit=24, file='namelists/exclude.nml')
      do
        read(unit=24, fmt=*, end=107) exm, exn
        fillcount = 0
        call self % fill(exm, exn, fillcount)
        self % lake_fill_events = self % lake_fill_events + 1
        self % lake_fill_m = [self % lake_fill_m, exm]
        self % lake_fill_n = [self % lake_fill_n, exn]
        self % lake_fill_count = [self % lake_fill_count, fillcount]
      end do
    end if

    107 close(unit=24)

  else

    ! set initial mask everywhere:
    self % mask = 1

    ! set up boundary points
    self % mask(:,1) = 0
    self % mask(:,self % nm) = 0

    ! close e and w edges if limited area:
    if (.not. config % isglobal) then
      self % mask(1,:) = 0
      self % mask(self % mm,:) = 0
    end if

    where(self % mask == 0) self % d_2d = config % dmin

  end if

  ! calculate the upper index for 1-d arrays:
  self % im = count(self % mask == 1)
  self % imm = self % mm * self % nm

end subroutine define_masks


recursive subroutine fill(self, m, n, fillcount)
! Mask out enclosed seas chosen by user.
  class(grid_type), intent(inout) :: self
  integer, intent(in) :: m, n
  integer, intent(in out) :: fillcount

  ! return if reached domain edge:
  if (any(m == [1, self % mm]) .or. any(n == [1, self % nm])) return

  ! fill:
  self % mask(m,n) = 0

  fillcount = fillcount + 1

  ! recurse in all directions:
  if (self % mask(m-1,n) == 1) call self % fill(m-1, n, fillcount)
  if (self % mask(m+1,n) == 1) call self % fill(m+1, n, fillcount)
  if (self % mask(m,n-1) == 1) call self % fill(m, n-1, fillcount)
  if (self % mask(m,n+1) == 1) call self % fill(m, n+1, fillcount)

end subroutine fill


subroutine partition(self, config)
! Partitions the domain for parallel computation.
#ifdef MPI
  use umwm_module, only: ierr, mpisize, nproc
#endif

  class(grid_type), intent(inout) :: self
  type(config_type), intent(in) :: config

#ifdef MPI
  integer :: nn
#ifdef ESMF
  integer :: i, m, n
#endif
#endif

  if (config % isglobal) continue

#ifndef MPI
  self % istart = 1
  self % iend = self % im
  self % iistart = self % istart
  self % iiend = self % iend
  self % ilen = self % iend - self % istart + 1
#else

  allocate(self % istart_all(0:mpisize-1), self % iend_all(0:mpisize-1), &
           self % ilen_all(0:mpisize-1))

  ! find out what is the length of my part:
  self % im_mod = mod(self % im, mpisize)
  if (self % im_mod == 0) then
    self % ilen = self % im / mpisize
    self % istart = nproc * self % ilen + 1
    self % iend = self % istart + self % ilen - 1
  else
    self % ilen = floor(float(self % im) / float(mpisize))
    self % istart = nproc * self % ilen + 1
    self % iend = self % istart + self % ilen - 1
    do nn = 1, self % im_mod
      if (nproc == nn - 1) self % ilen = self % ilen + 1
    end do
  end if

  ! adjust start/end boundaries:
  if (nproc < self % im_mod) then
    if (nproc == 0) then
      self % iend = self % iend + 1
    else
      self % istart = self % istart + nproc
      self % iend = self % iend + nproc + 1
    end if
  else
    self % istart = self % istart + self % im_mod
    self % iend = self % istart + self % ilen - 1
  end if

  ! Which direction for remapping? (matters only in parallel or global mode)
  if (self % mm >= self % nm) then
    self % remap_dir = 'v'
  else
    self % remap_dir = 'h'
  end if

  if (config % isglobal) self % remap_dir = 'v'

  ! The code below adjusts the start and end indices of each tile
  ! because currently ESMF DEBlockList accepts only regular rectangular
  ! domains.

#ifdef ESMF
  if (self % remap_dir == 'h') then ! row-major remapping

    ! adjust ends first:
    i = 0
    outer1: do n = 1, self % nm
      inner1: do m = 1, self % mm
        if (self % mask(m,n) == 1) i = i + 1
        if (i == self % iend .and. m /= self % mm) then
          self % iend = self % iend + count(self % mask(m+1:self % mm,n) == 1)
          exit outer1
        end if
      end do inner1
    end do outer1

    ! now adjust beginnings (all but proc 0):
    if (nproc /= 0) then
      i = 0
      outer2: do n = 1, self % nm
        inner2: do m = 1, self % mm
          if (self % mask(m,n) == 1) i = i + 1
          if (i == self % istart .and. m /= 1 .and. count(self % mask(1:m-1,n) == 1) > 0) then
            self % istart = self % istart + count(self % mask(m:self % mm,n) == 1)
            exit outer2
          end if
        end do inner2
      end do outer2
    end if

  else if (self % remap_dir == 'v') then ! column-major remapping

    ! adjust ends first:
    i = 0
    outer3: do m = 1, self % mm
      inner3: do n = 1, self % nm
        if (self % mask(m,n) == 1) i = i + 1
        if (i == self % iend .and. n /= self % nm) then
          self % iend = self % iend + count(self % mask(m,n+1:self % nm) == 1)
          exit outer3
        end if
      end do inner3
    end do outer3

    ! now adjust beginnings (all but proc 0):
    if (nproc /= 0) then
      i = 0
      outer4: do m = 1, self % mm
        inner4: do n = 1, self % nm
          if (self % mask(m,n) == 1) i = i + 1
          if (i == self % istart .and. n /= 1 .and. count(self % mask(m,1:n-1) == 1) > 0) then
            self % istart = self % istart + count(self % mask(m,n:self % nm) == 1)
            exit outer4
          end if
        end do inner4
      end do outer4
    end if

  end if

  ! adjust tile length:
  self % ilen = self % iend - self % istart + 1

#endif

  ! gather tile mpisize information to root process:
  call mpi_gather(self % istart, 1, MPI_INTEGER, self % istart_all, 1, MPI_INTEGER, &
                  0, MPI_COMM_WORLD, ierr)
  call mpi_gather(self % iend, 1, MPI_INTEGER, self % iend_all, 1, MPI_INTEGER, &
                  0, MPI_COMM_WORLD, ierr)
  call mpi_gather(self % ilen, 1, MPI_INTEGER, self % ilen_all, 1, MPI_INTEGER, &
                  0, MPI_COMM_WORLD, ierr)

#endif

end subroutine partition


subroutine remap(self, config)
! Remaps two-dimensional arrays into one-dimensional arrays while leaving
! out land points, assigns neighboring points, and builds halo metadata.
#ifdef MPI
  use umwm_module, only: ierr, mpisize, nproc
#endif

  class(grid_type), intent(inout) :: self
  type(config_type), intent(in) :: config
  integer :: i, m, n

#ifdef MPI
  integer :: nn
  integer :: counter, itemp
  integer, dimension(:), allocatable :: n_exchange_indices
#endif

  ! in serial mode this does not matter, but we must pick one:
#ifndef MPI
  self % remap_dir = 'v'
#endif

  ! first construct (m,n)->(i) transformation:
  self % ii = 0
  i = 0
  if (self % remap_dir == 'h') then ! column-major

    ! sea points:
    do n = 1, self % nm
      do m = 1, self % mm
        if (self % mask(m,n) == 1) then
          i = i + 1
          self % ii(m,n) = i
          self % dx(i) = self % dx_2d(m,n)
          self % dy(i) = self % dy_2d(m,n)
          self % d(i) = self % d_2d(m,n)
        end if
      end do
    end do

    ! land points:
    do n = 1, self % nm
      do m = 1, self % mm
        if (self % mask(m,n) == 0) then
          i = i + 1
          self % ii(m,n) = i
          self % dx(i) = self % dx_2d(m,n)
          self % dy(i) = self % dy_2d(m,n)
          self % d(i) = self % d_2d(m,n)
        end if
      end do
    end do

  else if (self % remap_dir == 'v') then ! row-major

    ! sea points:
    do m = 1, self % mm
      do n = 1, self % nm
        if (self % mask(m,n) == 1) then
          i = i + 1
          self % ii(m,n) = i
          self % dx(i) = self % dx_2d(m,n)
          self % dy(i) = self % dy_2d(m,n)
          self % d(i) = self % d_2d(m,n)
        end if
      end do
    end do

    ! land points:
    do m = 1, self % mm
      do n = 1, self % nm
        if (self % mask(m,n) == 0) then
          i = i + 1
          self % ii(m,n) = i
          self % dx(i) = self % dx_2d(m,n)
          self % dy(i) = self % dy_2d(m,n)
          self % d(i) = self % d_2d(m,n)
        end if
      end do
    end do

  else

    stop 'umwm: remap: error - remap_dir must be ''h'' or ''v'''

  end if

  ! now construct (i)->(m,n) transformation:
  self % mi = 0
  self % ni = 0
  do n = 1, self % nm
    do m = 1, self % mm
      i = self % ii(m,n)
      self % mi(i) = m
      self % ni(i) = n
    end do
  end do

  ! neighboring point indices for advection and refraction:
  self % is = 0
  self % in = 0
  self % iw = 0
  self % ie = 0
  do n = 2, self % nm - 1
    do m = 2, self % mm - 1
      i = self % ii(m,n)
      self % is(i) = self % ii(m,n-1)
      self % in(i) = self % ii(m,n+1)
      self % iw(i) = self % ii(m-1,n)
      self % ie(i) = self % ii(m+1,n)
    end do
  end do

  ! adjust periodic boundary:
  if (config % isglobal) then
    do n = 2, self % nm - 1

      self % is(self % ii(1,n)) = self % ii(1,n-1)
      self % in(self % ii(1,n)) = self % ii(1,n+1)
      self % iw(self % ii(1,n)) = self % ii(self % mm,n)
      self % ie(self % ii(1,n)) = self % ii(2,n)

      self % is(self % ii(self % mm,n)) = self % ii(self % mm,n-1)
      self % in(self % ii(self % mm,n)) = self % ii(self % mm,n+1)
      self % iw(self % ii(self % mm,n)) = self % ii(self % mm-1,n)
      self % ie(self % ii(self % mm,n)) = self % ii(1,n)

    end do
  end if

  ! true indices:
  ! (is,in,iw,ie will be aliased for land points)
  self % iis = self % is
  self % iin = self % in
  self % iiw = self % iw
  self % iie = self % ie

  ! cell edges in x and y:
  do i = self % istart, self % iend
    self % dxn(i) = 0.5 * (self % dx(i) + self % dx(self % iin(i)))
    self % dxs(i) = 0.5 * (self % dx(i) + self % dx(self % iis(i)))
    self % dye(i) = 0.5 * (self % dy(i) + self % dy(self % iie(i)))
    self % dyw(i) = 0.5 * (self % dy(i) + self % dy(self % iiw(i)))
  end do

#ifdef MPI
  if (nproc == 0) then ! root process

    self % iistart = self % istart
    itemp = self % iend

    do
      if (self % remap_dir == 'h') self % iiend = self % in(itemp)
      if (self % remap_dir == 'v') self % iiend = self % ie(itemp)
      if (self % iiend > self % im) then
        itemp = itemp - 1
        cycle
      else
        exit
      end if
    end do

  else if (nproc == mpisize - 1) then

    itemp = self % istart
    self % iiend = self % iend

    do
      if (self % remap_dir == 'h') self % iistart = self % is(itemp)
      if (self % remap_dir == 'v') self % iistart = self % iw(itemp)
      if (self % iistart > self % im) then
        itemp = itemp + 1
        cycle
      else
        exit
      end if
    end do

  else

    itemp = self % istart

    do
      if (self % remap_dir == 'h') self % iistart = self % is(itemp)
      if (self % remap_dir == 'v') self % iistart = self % iw(itemp)
      if (self % iistart > self % im) then
        itemp = itemp + 1
        cycle
      else
        exit
      end if
    end do

    itemp = self % iend

    do
      if (self % remap_dir == 'h') self % iiend = self % in(itemp)
      if (self % remap_dir == 'v') self % iiend = self % ie(itemp)
      if (self % iiend > self % im) then
        itemp = itemp - 1
        cycle
      else
        exit
      end if
    end do

  end if

  if (config % isglobal) then

    allocate(n_exchange_indices(0))

    self % first_col_len = 0
    do n = 2, self % nm - 1
      if (self % mask(1,n) == 1 .and. self % mask(self % mm,n) == 1) then
        n_exchange_indices = [n_exchange_indices, n]
        self % first_col_len = self % first_col_len + 1
      end if
    end do
    self % last_col_len = self % first_col_len

    allocate(self % i_exchange_indices(self % first_col_len))

    ! find west neighbor indices on first processor
    if (nproc == 0) then

      do n = 1, self % first_col_len
        self % i_exchange_indices(n) = self % ii(1,n_exchange_indices(n))
      end do

      self % iistart = self % iistart - self % last_col_len

      counter = 0
      m = 1
      do n = 2, self % nm - 1
        if (self % mask(m,n) == 1) then
          i = self % ii(m,n)
          if (self % mask(self % mi(self % iw(i)), self % ni(self % iw(i))) == 0) then
            self % iw(i) = self % iistart - 1
          else
            self % iw(i) = self % iistart + counter
            counter = counter + 1
          end if
        end if
      end do

    end if

    ! find east neighbor indices on last processor
    if (nproc == mpisize - 1) then

      do n = 1, self % first_col_len
        self % i_exchange_indices(n) = self % ii(self % mm,n_exchange_indices(n))
      end do

      self % iiend = self % iiend + self % first_col_len

      counter = 1
      m = self % mm
      do n = 2, self % nm - 1
        if (self % mask(m,n) == 1) then
          i = self % ii(m,n)
          if (self % mask(self % mi(self % ie(i)), self % ni(self % ie(i))) == 0) then
            self % ie(i) = self % iistart - 1
          else
            self % ie(i) = self % iend + counter
            counter = counter + 1
          end if
        end if
      end do

    end if

    deallocate(n_exchange_indices)

  end if
#endif

#ifndef MPI
  self % iistart = self % istart
  self % iiend = self % iend
#endif

#ifdef MPI
  ! distribute iistart and iiend to everyone:
  allocate(self % iistart_all(0:mpisize-1), self % iiend_all(0:mpisize-1))
  call mpi_allgather(self % iistart, 1, MPI_INTEGER, self % iistart_all, 1, &
                     MPI_INTEGER, MPI_COMM_WORLD, ierr)
  call mpi_allgather(self % iiend, 1, MPI_INTEGER, self % iiend_all, 1, &
                     MPI_INTEGER, MPI_COMM_WORLD, ierr)
#endif

  ! treat land points for halo-adjacent arrays:
  do i = self % istart, self % iend
    if (self % is(i) < self % iistart .or. self % is(i) > self % iiend) self % is(i) = self % iistart - 1
    if (self % in(i) < self % iistart .or. self % in(i) > self % iiend) self % in(i) = self % iistart - 1
    if (self % ie(i) < self % iistart .or. self % ie(i) > self % iiend) self % ie(i) = self % iistart - 1
    if (self % iw(i) < self % iistart .or. self % iw(i) > self % iiend) self % iw(i) = self % iistart - 1
  end do

#ifdef MPI
  ! create the domain partitioning field for output:
  if (nproc == 0) then
    self % nproc_out = -1
    do n = 1, self % nm
      do m = 1, self % mm
        do nn = 0, mpisize - 1
          if (self % ii(m,n) >= self % istart_all(nn) .and. self % ii(m,n) <= self % iend_all(nn)) then
            self % nproc_out(m,n) = nn
            exit
          end if
        end do
      end do
    end do
  end if
#endif

  ! compute grid cell areas and reciprocals:
  self % ar = self % dx(self % istart:self % iend) * self % dy(self % istart:self % iend)
  self % oneovdx = 1. / self % dx(self % istart:self % iend)
  self % oneovdy = 1. / self % dy(self % istart:self % iend)
  self % oneovar = 1. / self % ar

end subroutine remap


subroutine initialize_direction_projection(self, spectrum)
  class(grid_type), intent(inout) :: self
  type(spectrum_type), intent(in) :: spectrum

  integer :: i, p

  if (allocated(self % cth_curv)) deallocate(self % cth_curv)
  if (allocated(self % sth_curv)) deallocate(self % sth_curv)

  allocate(self % cth_curv(spectrum % num_directions,self % istart:self % iend))
  allocate(self % sth_curv(spectrum % num_directions,self % istart:self % iend))

  ! calculate wave ray directions adjusted for grid curvature:
  do i = self % istart, self % iend
    do p = 1, spectrum % num_directions
      self % cth_curv(p,i) = cos(real(spectrum % direction(p), kind(self % cth_curv)) &
                            + self % curv(self % mi(i), self % ni(i)))
      self % sth_curv(p,i) = sin(real(spectrum % direction(p), kind(self % sth_curv)) &
                            + self % curv(self % mi(i), self % ni(i)))
    end do
  end do

end subroutine initialize_direction_projection


pure function remap_i2mn(self, field_i) result(field_mn)
! Remaps an (i) indexed array to (m,n).
  class(grid_type), intent(in) :: self
  real, dimension(self % imm), intent(in) :: field_i
  real, dimension(self % mm,self % nm) :: field_mn

  integer :: m, n

  do n = 1, self % nm
    do m = 1, self % mm
      field_mn(m,n) = field_i(self % ii(m,n))
    end do
  end do

end function remap_i2mn


pure function remap_mn2i(self, field_mn) result(field_i)
! Remaps an (m,n) indexed array to (i).
  class(grid_type), intent(in) :: self
  real, dimension(self % mm,self % nm), intent(in) :: field_mn
  real, dimension(self % imm) :: field_i

  integer :: m, n

  do n = 1, self % nm
    do m = 1, self % mm
      field_i(self % ii(m,n)) = field_mn(m,n)
    end do
  end do

end function remap_mn2i


subroutine finalize(self)
  class(grid_type), intent(inout) :: self

  if (allocated(self % iw)) deallocate(self % iw, self % ie, self % is, self % in)
  if (allocated(self % iiw)) deallocate(self % iiw, self % iie, self % iis, self % iin)
  if (allocated(self % i_exchange_indices)) deallocate(self % i_exchange_indices)
  if (allocated(self % mi)) deallocate(self % mi, self % ni)
  if (allocated(self % ii)) deallocate(self % ii)
  if (allocated(self % mask)) deallocate(self % mask)
  if (allocated(self % nproc_out)) deallocate(self % nproc_out)
  if (allocated(self % lake_fill_m)) deallocate(self % lake_fill_m)
  if (allocated(self % lake_fill_n)) deallocate(self % lake_fill_n)
  if (allocated(self % lake_fill_count)) deallocate(self % lake_fill_count)
  if (allocated(self % istart_all)) deallocate(self % istart_all, self % iend_all, self % ilen_all)
  if (allocated(self % iistart_all)) deallocate(self % iistart_all, self % iiend_all)
  if (allocated(self % cth_curv)) deallocate(self % cth_curv)
  if (allocated(self % sth_curv)) deallocate(self % sth_curv)
  if (allocated(self % ar_2d)) deallocate(self % ar_2d)
  if (allocated(self % curv)) deallocate(self % curv)
  if (allocated(self % d_2d)) deallocate(self % d_2d, self % dlon, self % dlat, self % dx_2d, self % dy_2d)
  if (allocated(self % lat)) deallocate(self % lat, self % lon)
  if (allocated(self % x)) deallocate(self % x, self % y)
  if (allocated(self % ar)) deallocate(self % ar)
  if (allocated(self % d)) deallocate(self % d, self % dx, self % dy)
  if (allocated(self % dxn)) deallocate(self % dxn, self % dxs, self % dyw, self % dye)
  if (allocated(self % oneovar)) deallocate(self % oneovar, self % oneovdx, self % oneovdy)

end subroutine finalize


subroutine grid_nc_check(stat)
  use netcdf

  integer, intent(in) :: stat

  if (stat /= nf90_noerr) then
    write(*,*) 'error in netcdf i/o'
    write(*,*) trim(nf90_strerror(stat))
    stop
  end if

end subroutine grid_nc_check

end module umwm_grid
