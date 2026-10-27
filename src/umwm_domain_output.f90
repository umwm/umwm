submodule (umwm_domain) umwm_domain_output

  ! Integrated diagnostics, ocean stress, and gridded and spectrum output.

  use, intrinsic :: iso_fortran_env, only: real64
  use netcdf
  use umwm_io, only: gather_field, put_description, put_time_metadata, &
                     seconds_since_reference

  implicit none

contains

  module subroutine diagnostics(self)
    ! Integrated quantities of the wave spectrum for output.
    class(domain_type), intent(inout) :: self
    real :: spectrumbin(self % spectrum % num_frequencies, self % spectrum % num_directions)
    real :: m0, m2, xcomp, ycomp, mag, ekdkovcp, rhow_dthg
    integer :: i, o, p, opeak, ppeak, peak(2)

#ifdef _OPENACC
    call diagnostics_acc(self)
    return
#endif

    associate(e => self % e, kdk => self % kdk, f => self % f, cth => self % cth, &
              sth => self % sth, forcing => self % forcing, &
              om => self % spectrum % num_frequencies, pm => self % spectrum % num_directions)

    do i = self % grid % istart, self % grid % iend

      do p = 1, pm
        spectrumbin(:,p) = e(:,p,i) * kdk(:,i)
      end do

      ! mean wave period
      m0 = 0
      m2 = 0
      do p = 1, pm
        do o = 1, om
          m0 = m0 + spectrumbin(o,p)
          m2 = m2 + f(o)**2 * spectrumbin(o,p)
        end do
      end do
      self % mwp(i) = sqrt(m0 / (m2 + tiny(m2)))

      ! total wave momentum and its flux
      self % momx(i) = 0
      self % momy(i) = 0
      self % cgmxx(i) = 0
      self % cgmxy(i) = 0
      self % cgmyy(i) = 0
      do p = 1, pm
        do o = 1, self % oc(i)
          ekdkovcp = spectrumbin(o,p) * self % invcp0(o,i)
          self % momx(i) = self % momx(i) + ekdkovcp * cth(p)
          self % momy(i) = self % momy(i) + ekdkovcp * sth(p)
          self % cgmxx(i) = self % cgmxx(i) + self % cg0(o,i) * ekdkovcp * cth(p)**2
          self % cgmxy(i) = self % cgmxy(i) + self % cg0(o,i) * ekdkovcp * cth(p) * sth(p)
          self % cgmyy(i) = self % cgmyy(i) + self % cg0(o,i) * ekdkovcp * sth(p)**2
        end do
      end do
      rhow_dthg = forcing % rhow(i) * self % dthg
      self % momx(i) = self % momx(i) * rhow_dthg
      self % momy(i) = self % momy(i) * rhow_dthg
      self % cgmxx(i) = self % cgmxx(i) * rhow_dthg
      self % cgmxy(i) = self % cgmxy(i) * rhow_dthg
      self % cgmyy(i) = self % cgmyy(i) * rhow_dthg

      ! significant wave height
      self % ht(i) = 4 * sqrt(sum(spectrumbin) * self % dth)

      ! mean wave direction
      xcomp = 0
      ycomp = 0
      do p = 1, pm
        mag = sum(spectrumbin(:,p))
        xcomp = xcomp + mag * cth(p)
        ycomp = ycomp + mag * sth(p)
      end do
      self % mwd(i) = atan2(ycomp, xcomp)

      ! mean-squared slope and mean wavelength from wavenumber spectrum moments
      m2 = 0
      do p = 1, pm
        m2 = m2 + sum(e(:,p,i) * self % k3dk(:,i))
      end do
      self % mss(i) = m2 * self % dth
      self % mwl(i) = twopi * sqrt(m0 / (m2 + tiny(m2)))

      ! dominant (peak) wave quantities
      peak = maxloc(spectrumbin)
      opeak = peak(1)
      ppeak = peak(2)

      self % dwd(i) = self % th(ppeak)
      self % dwl(i) = twopi / self % k(opeak,i)
      self % dwp(i) = 1. / f(opeak)
      self % dcp0(i) = self % cp0(opeak,i)
      self % dcg0(i) = self % cg0(opeak,i)
      self % dcp(i) = self % dcp0(i) + forcing % uc(i) * cth(ppeak) + forcing % vc(i) * sth(ppeak)
      self % dcg(i) = self % dcg0(i) + forcing % uc(i) * cth(ppeak) + forcing % vc(i) * sth(ppeak)

    end do

    end associate

  end subroutine diagnostics


  module subroutine ocean_stress(self)
    ! Momentum fluxes into the ocean top and bottom, energy fluxes from
    ! air and into the ocean, and the components of form drag.
    class(domain_type), intent(inout) :: self
    real, dimension(self % spectrum % num_frequencies) :: taux_util, tauy_util, snl_cp
    real :: tail, rhow_dthg, form
    integer :: i, o, p

#ifdef _OPENACC
    call ocean_stress_acc(self)
    return
#endif

    associate(e => self % e, ssin => self % ssin, sds => self % sds, snl => self % snl, &
              cothkd => self % cothkd, cp0 => self % cp0, invcp0 => self % invcp0, kdk => self % kdk, &
              cth => self % cth, sth => self % sth, forcing => self % forcing, &
              om => self % spectrum % num_frequencies, pm => self % spectrum % num_directions, &
              bf1 => self % bf1, bf2 => self % bf2)

    do i = self % grid % istart, self % grid % iend

      tail = stress_tail(forcing % wspd(i), self % k(om,i))
      rhow_dthg = forcing % rhow(i) * self % dthg

      ! stress into ocean top (positive into ocean): dissipation into currents,
      ! corrected for the Snl in the diagnostic part
      taux_util = 0
      tauy_util = 0
      do p = 1, pm
        do o = 1, om
          taux_util(o) = taux_util(o) + e(o,p,i) &
                       * (sds(o,p,i) * cothkd(o,i) + self % sdt(o,i) + self % sdv(o,i)) &
                       * cth(p) * invcp0(o,i)
          tauy_util(o) = tauy_util(o) + e(o,p,i) &
                       * (sds(o,p,i) * cothkd(o,i) + self % sdt(o,i) + self % sdv(o,i)) &
                       * sth(p) * invcp0(o,i)
        end do
        do o = 3, om
          snl_cp(o) = snl(o,p,i) * (bf1 * cp0(o,i) * invcp0(o-1,i) + bf2 * cp0(o,i) * invcp0(o-2,i))
          taux_util(o) = taux_util(o) - snl_cp(o) * cth(p) * invcp0(o,i)
          tauy_util(o) = tauy_util(o) - snl_cp(o) * sth(p) * invcp0(o,i)
        end do
      end do

      self % tailocnx(i) = taux_util(om) * self % k(om,i) * tail * rhow_dthg
      self % tailocny(i) = tauy_util(om) * self % k(om,i) * tail * rhow_dthg
      self % taux_ocntop(i) = sum(taux_util * kdk(:,i)) * rhow_dthg + self % tailocnx(i) + self % taux_skin(i)
      self % tauy_ocntop(i) = sum(tauy_util * kdk(:,i)) * rhow_dthg + self % tailocny(i) + self % tauy_skin(i)

      ! stress into ocean bottom (positive into ocean)
      taux_util = 0
      tauy_util = 0
      do p = 1, pm
        taux_util = taux_util + e(:,p,i) * self % sbf(:,i) * cth(p) * invcp0(:,i)
        tauy_util = tauy_util + e(:,p,i) * self % sbf(:,i) * sth(p) * invcp0(:,i)
      end do
      self % taux_ocnbot(i) = sum(taux_util * kdk(:,i)) * rhow_dthg
      self % tauy_ocnbot(i) = sum(tauy_util * kdk(:,i)) * rhow_dthg

      ! Snl conserves energy, but not momentum. the following term is
      ! the momentum loss due to non-linear downshifting of energy.
      taux_util = 0
      tauy_util = 0
      do p = 1, pm
        do o = 3, self % oc(i)
          snl_cp(o) = snl(o,p,i) * (bf1 * (1 - cp0(o,i) * invcp0(o-1,i)) &
                                  + bf2 * (1 - cp0(o,i) * invcp0(o-2,i)))
          taux_util(o) = taux_util(o) + snl_cp(o) * cth(p) * invcp0(o,i)
          tauy_util(o) = tauy_util(o) + snl_cp(o) * sth(p) * invcp0(o,i)
        end do
      end do
      self % taux_snl(i) = sum(taux_util * kdk(:,i)) * rhow_dthg
      self % tauy_snl(i) = sum(tauy_util * kdk(:,i)) * rhow_dthg

      ! energy flux into ocean
      taux_util = 0
      tauy_util = 0
      do p = 1, pm
        taux_util = taux_util + e(:,p,i) * (sds(:,p,i) * cothkd(:,i)) * cth(p)
        tauy_util = tauy_util + e(:,p,i) * (sds(:,p,i) * cothkd(:,i)) * sth(p)
        do o = 3, om
          snl_cp(o) = snl(o,p,i) * (bf1 * cp0(o,i) * invcp0(o-1,i) + bf2 * cp0(o,i) * invcp0(o-2,i))
          taux_util(o) = taux_util(o) - snl_cp(o) * cth(p)
          tauy_util(o) = tauy_util(o) - snl_cp(o) * sth(p)
        end do
      end do
      self % epsx_ocn(i) = sum(taux_util * kdk(:,i)) * rhow_dthg
      self % epsy_ocn(i) = sum(tauy_util * kdk(:,i)) * rhow_dthg

      ! energy flux from air
      call directional_integral(om, pm, e(:,:,i), ssin(:,:,i), cth, sth, taux_util, tauy_util)
      self % epsx_atm(i) = sum(taux_util * kdk(:,i)) * rhow_dthg
      self % epsy_atm(i) = sum(tauy_util * kdk(:,i)) * rhow_dthg

      ! components of form drag; not momentum fluxes into the ocean,
      ! but computed here because they are needed only for output.
      self % taux1(i) = 0
      self % tauy1(i) = 0
      self % taux2(i) = 0
      self % tauy2(i) = 0
      self % taux3(i) = 0
      self % tauy3(i) = 0
      do p = 1, pm
        do o = 1, om
          form = e(o,p,i) * ssin(o,p,i) * invcp0(o,i) * kdk(o,i)
          if (ssin(o,p,i) > 0) then ! wind pushing waves
            self % taux1(i) = self % taux1(i) + form * cth(p)
            self % tauy1(i) = self % tauy1(i) + form * sth(p)
          else if (cos(forcing % wdir(i) - self % th(p)) < 0) then ! waves against wind
            self % taux2(i) = self % taux2(i) + form * cth(p)
            self % tauy2(i) = self % tauy2(i) + form * sth(p)
          else ! waves overrunning wind
            self % taux3(i) = self % taux3(i) + form * cth(p)
            self % tauy3(i) = self % tauy3(i) + form * sth(p)
          end if
        end do
      end do
      self % taux1(i) = self % taux1(i) * rhow_dthg
      self % tauy1(i) = self % tauy1(i) * rhow_dthg
      self % taux2(i) = self % taux2(i) * rhow_dthg
      self % tauy2(i) = self % tauy2(i) * rhow_dthg
      self % taux3(i) = self % taux3(i) * rhow_dthg
      self % tauy3(i) = self % tauy3(i) * rhow_dthg

    end do

    end associate

  end subroutine ocean_stress


  module subroutine output_grid_nc(self, timestr)
    ! Writes out model gridded output in a netcdf format.
    class(domain_type), intent(in) :: self
    character(19), intent(in) :: timestr

    character(19) :: timestrnew
    integer :: stat, ncid, timeid, l
    integer :: xdimid, ydimid, zdimid, fdimid, thdimid, tdimid
    logical :: define
    real :: output_field(self % grid % mm, self % grid % nm)
    real(real64) :: time_value(1)

    timestrnew = timestr
    timestrnew(11:11) = '_'

    associate(grid => self % grid, is => self % grid % istart, ie => self % grid % iend, &
              om => self % spectrum % num_frequencies, pm => self % spectrum % num_directions, &
              lm => self % stokes % lm, stokes => self % config % stokes)

    if (nproc == 0) then

      stat = nf90_create('output/umwmout_' // timestrnew // '.nc', nf90_clobber, ncid)

      stat = nf90_def_dim(ncid, 'x', grid % mm, xdimid)
      stat = nf90_def_dim(ncid, 'y', grid % nm, ydimid)
      stat = nf90_def_dim(ncid, 'f', om, fdimid)
      stat = nf90_def_dim(ncid, 'th', pm, thdimid)
      stat = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, tdimid)

      stat = nf90_def_var(ncid, 'time', NF90_DOUBLE, [tdimid], timeid)
      call put_time_metadata(ncid, timeid, self % config % reftimestr)

      if (stokes) then
        stat = nf90_def_dim(ncid, 'z', lm, zdimid)
        call define_var('z', 'depth', 'm', [zdimid])
        call define_var('u_stokes', 'stokes drift x-component', 'm/s', [xdimid, ydimid, zdimid, tdimid])
        call define_var('v_stokes', 'stokes drift y-component', 'm/s', [xdimid, ydimid, zdimid, tdimid])
        call define_var('d_stokes', 'stokes drift e-folding depth', 'm', [xdimid, ydimid, tdimid])
      end if

      call define_var('frequency', 'frequency', 'hz', [fdimid])
      call define_var('theta', 'directions', 'rad', [thdimid])
      call define_var('lon', 'longitude', 'degrees east', [xdimid, ydimid, tdimid])
      call define_var('lat', 'latitude', 'degrees north', [xdimid, ydimid, tdimid])
      call define_var('seamask', 'seamask', 'non-dimensional', [xdimid, ydimid, tdimid], nf90_int)
      call define_var('depth', 'ocean depth', 'm', [xdimid, ydimid, tdimid])

    end if

    define = .true.
    call fields()

    if (nproc == 0) then
      stat = nf90_enddef(ncid)
      time_value(1) = seconds_since_reference(timestr, self % config % reftimestr)
      stat = nf90_put_var(ncid, varid('time'), time_value, start=[1], count=[1])
      stat = nf90_put_var(ncid, varid('frequency'), self % spectrum % frequency, start=[1], count=[om])
      stat = nf90_put_var(ncid, varid('theta'), self % spectrum % direction, start=[1], count=[pm])
      stat = nf90_put_var(ncid, varid('lon'), grid % lon, start=[1,1,1], count=[grid % mm, grid % nm, 1])
      stat = nf90_put_var(ncid, varid('lat'), grid % lat, start=[1,1,1], count=[grid % mm, grid % nm, 1])
      stat = nf90_put_var(ncid, varid('seamask'), grid % mask, start=[1,1,1], count=[grid % mm, grid % nm, 1])
      stat = nf90_put_var(ncid, varid('depth'), grid % d_2d, start=[1,1,1], count=[grid % mm, grid % nm, 1])
    end if

    define = .false.
    call fields()

    if (stokes) then
      if (nproc == 0) stat = nf90_put_var(ncid, varid('z'), self % stokes % depth, start=[1], count=[lm])
      do l = 1, lm
        call write_field('u_stokes', self % stokes % us(is:ie,l), [1,1,l,1], [grid % mm, grid % nm, 1, 1])
        call write_field('v_stokes', self % stokes % vs(is:ie,l), [1,1,l,1], [grid % mm, grid % nm, 1, 1])
      end do
      call write_field('d_stokes', self % stokes % ds(is:ie), [1,1,1], [grid % mm, grid % nm, 1])
    end if

    if (nproc == 0) then
      stat = nf90_close(ncid)
      write(unit=*,fmt='(a)') 'umwm: output_nc: output written to output/umwmout_' // timestrnew // '.nc'
    end if

    end associate

  contains

    subroutine fields()
      ! Visits all gridded output fields, in the order of definition.
      associate(is => self % grid % istart, ie => self % grid % iend, forcing => self % forcing)
      call field('wspd', 'wind speed', 'm/s', forcing % wspd(is:ie))
      call field('wdir', 'wind direction', 'rad', forcing % wdir(is:ie))
      call field('uc', 'ocean current, x-component', 'm/s', forcing % uc(is:ie))
      call field('vc', 'ocean current, y-component', 'm/s', forcing % vc(is:ie))
      call field('rhoa', 'air density', 'kg/m^3', forcing % rhoa(is:ie))
      call field('rhow', 'water density', 'kg/m^3', forcing % rhow(is:ie))
      call field('fice', 'seaice fraction', 'non-dimensional', forcing % fice(is:ie))
      call field('psim', 'universal stability function for momentum', 'non-dimensional', &
                 self % psim(is:ie))
      call field('momx', 'momentum, x-component', 'kgm/s', self % momx)
      call field('momy', 'momentum, y-component', 'kgm/s', self % momy)
      call field('cgmxx', 'cg*momentum, xx-component', 'kgm^2/s^2', self % cgmxx)
      call field('cgmxy', 'cg*momentum, xy-component', 'kgm^2/s^2', self % cgmxy)
      call field('cgmyy', 'cg*momentum, yy-component', 'kgm^2/s^2', self % cgmyy)
      call field('shelt', 'sheltering coefficient', 'non-dimensional', self % shelt)
      call field('epsx_atm', 'wave energy growth flux, x-component', 'kg/s^3', self % epsx_atm)
      call field('epsy_atm', 'wave energy growth flux, y-component', 'kg/s^3', self % epsy_atm)
      call field('epsx_ocn', 'wave energy dissipation flux, x-component', 'kg/s^3', self % epsx_ocn)
      call field('epsy_ocn', 'wave energy dissipation flux, y-component', 'kg/s^3', self % epsy_ocn)
      call field('taux_form', 'form drag, x-component', 'n/m^2', self % taux_form)
      call field('tauy_form', 'form drag, y-component', 'n/m^2', self % tauy_form)
      call field('taux_form_1', 'form drag, part 1, x-component', 'n/m^2', self % taux1)
      call field('tauy_form_1', 'form drag, part 1, y-component', 'n/m^2', self % tauy1)
      call field('taux_form_2', 'form drag, part 2, x-component', 'n/m^2', self % taux2)
      call field('tauy_form_2', 'form drag, part 2, y-component', 'n/m^2', self % tauy2)
      call field('taux_form_3', 'form drag, part 3, x-component', 'n/m^2', self % taux3)
      call field('tauy_form_3', 'form drag, part 3, y-component', 'n/m^2', self % tauy3)
      call field('taux_skin', 'skin drag, x-component', 'n/m^2', self % taux_skin)
      call field('tauy_skin', 'skin drag, y-component', 'n/m^2', self % tauy_skin)
      call field('taux_diag', 'diagnostic form drag, x-component', 'n/m^2', self % taux_diag)
      call field('tauy_diag', 'diagnostic form drag, y-component', 'n/m^2', self % tauy_diag)
      call field('taux_ocn', 'momentum flux from breaking waves to ocean top, x-component', &
                 'n/m^2', self % taux_ocntop)
      call field('tauy_ocn', 'momentum flux from breaking waves to ocean top, y-component', &
                 'n/m^2', self % tauy_ocntop)
      call field('taux_bot', 'momentum flux from waves to ocean bottom, x-component', &
                 'n/m^2', self % taux_ocnbot)
      call field('tauy_bot', 'momentum flux from waves to ocean bottom, y-component', &
                 'n/m^2', self % tauy_ocnbot)
      call field('taux_snl', 'momentum flux due to snl, x-component', 'n/m^2', self % taux_snl)
      call field('tauy_snl', 'momentum flux due to snl, y-component', 'n/m^2', self % tauy_snl)
      call field('tailatmx', 'atmosphere tail stress part, x-component', 'n/m^2', self % tailatmx)
      call field('tailatmy', 'atmosphere tail stress part, y-component', 'n/m^2', self % tailatmy)
      call field('tailocnx', 'ocean tail stress part, x-component', 'n/m^2', self % tailocnx)
      call field('tailocny', 'ocean tail stress part, y-component', 'n/m^2', self % tailocny)
      call field('cd', 'drag coefficient of air', 'non-dimensional', self % cd)
      call field('ust', 'friction velocity of air', 'm/s', self % ustar)
      call field('swh', 'significant wave height', 'm', self % ht)
      call field('mss', 'mean-squared slope', 'non-dimensional', self % mss)
      call field('mwp', 'mean wave period', 's', self % mwp)
      call field('mwl', 'mean wavelength', 'm', self % mwl)
      call field('mwd', 'mean wave direction', 'rad', self % mwd)
      call field('dwp', 'dominant wave period', 's', self % dwp)
      call field('dwl', 'dominant wavelength', 'm', self % dwl)
      call field('dwd', 'dominant wave direction', 'rad', self % dwd)
      call field('dcp0', 'dominant phase speed, intrinsic', 'm/s', self % dcp0)
      call field('dcg0', 'dominant group speed, intrinsic', 'm/s', self % dcg0)
      call field('dcp', 'dominant phase speed', 'm/s', self % dcp)
      call field('dcg', 'dominant group speed', 'm/s', self % dcg)
      call field('physics_time_step', 'Physics time step', 's', self % physics_time_step)
      end associate
    end subroutine fields

    subroutine field(name, description, units, values)
      ! Defines or writes a gridded output field, depending on define.
      character(*), intent(in) :: name, description, units
      real, intent(in) :: values(:)
      if (define) then
        if (nproc == 0) call define_var(name, description, units, [xdimid, ydimid, tdimid])
      else
        call write_field(name, values, [1,1,1], [self % grid % mm, self % grid % nm, 1])
      end if
    end subroutine field

    subroutine define_var(name, description, units, dimids, xtype)
      character(*), intent(in) :: name, description, units
      integer, intent(in) :: dimids(:)
      integer, intent(in), optional :: xtype
      integer :: id
      if (present(xtype)) then
        stat = nf90_def_var(ncid, name, xtype, dimids, id)
      else
        stat = nf90_def_var(ncid, name, NF90_FLOAT, dimids, id)
      end if
      call put_description(ncid, id, description)
      stat = nf90_put_att(ncid, id, name='units', values=units)
    end subroutine define_var

    subroutine write_field(name, values, start, count)
      ! Gathers a tile field on the root process and writes it.
      character(*), intent(in) :: name
      real, intent(in) :: values(:)
      integer, intent(in) :: start(:), count(:)
      call gather_field(values, output_field, self % grid)
      if (nproc == 0) stat = nf90_put_var(ncid, varid(name), output_field, start=start, count=count)
    end subroutine write_field

    integer function varid(name)
      character(*), intent(in) :: name
      integer :: status
      status = nf90_inq_varid(ncid, name, varid)
    end function varid

  end subroutine output_grid_nc


  module subroutine output_spectrum_nc(self, timestr)
    ! Writes out model spectrum output in a netcdf format.
    class(domain_type), intent(inout) :: self
    character(19), intent(in) :: timestr

    character(9999) :: filename
    integer :: stat, ncid, nn, record, i
    integer :: fdimid, thdimid, tdimid, scalarid
    integer :: lon_scalarid, lat_scalarid, timeid, wspdid, wdirid
    integer :: specid, sinid, sdsid, snlid, freqid, thetaid, wlid, sdtid, sdvid, sbfid
    real(real64) :: time_value(1)

    if (.not. allocated(self % spectrum_i)) call read_points()

    record = self % spectrum_records + 1

    associate(om => self % spectrum % num_frequencies, pm => self % spectrum % num_directions)

    do nn = 1, size(self % spectrum_i)

      i = self % spectrum_i(nn)
      if (i < self % grid % istart .or. i > self % grid % iend) cycle

      filename = 'output/umwmspc_' // trim(self % spectrum_id(nn)) // '_' &
                 // self % spectrum_timestr // '.nc'

      if (record == 1) then

        ! create output file:
        stat = nf90_create(trim(filename), nf90_clobber, ncid)

        ! define dimensions:
        stat = nf90_def_dim(ncid, 'scalar', 1, scalarid)
        stat = nf90_def_dim(ncid, 'frequency', om, fdimid)
        stat = nf90_def_dim(ncid, 'direction', pm, thdimid)
        stat = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, tdimid)

        ! define variables:
        stat = nf90_def_var(ncid, 'time', NF90_DOUBLE, [tdimid], timeid)
        stat = nf90_def_var(ncid, 'frequency', NF90_FLOAT, [fdimid], freqid)
        stat = nf90_def_var(ncid, 'wavenumber', NF90_FLOAT, [fdimid], wlid)
        stat = nf90_def_var(ncid, 'direction', NF90_FLOAT, [thdimid], thetaid)
        stat = nf90_def_var(ncid, 'longitude', NF90_FLOAT, [scalarid], lon_scalarid)
        stat = nf90_def_var(ncid, 'latitude', NF90_FLOAT, [scalarid], lat_scalarid)
        stat = nf90_def_var(ncid, 'wspd', NF90_FLOAT, [scalarid,tdimid], wspdid)
        stat = nf90_def_var(ncid, 'wdir', NF90_FLOAT, [scalarid,tdimid], wdirid)
        stat = nf90_def_var(ncid, 'F', NF90_FLOAT, [fdimid,thdimid,tdimid], specid)
        stat = nf90_def_var(ncid, 'Sin', NF90_FLOAT, [fdimid,thdimid,tdimid], sinid)
        stat = nf90_def_var(ncid, 'Sds', NF90_FLOAT, [fdimid,thdimid,tdimid], sdsid)
        stat = nf90_def_var(ncid, 'Sdt', NF90_FLOAT, [fdimid,tdimid], sdtid)
        stat = nf90_def_var(ncid, 'Sdv', NF90_FLOAT, [fdimid], sdvid)
        stat = nf90_def_var(ncid, 'Sbf', NF90_FLOAT, [fdimid], sbfid)
        stat = nf90_def_var(ncid, 'Snl', NF90_FLOAT, [fdimid,thdimid,tdimid], snlid)

        call put_time_metadata(ncid, timeid, self % config % reftimestr)
        call put_description(ncid, freqid, 'Frequency')
        call put_description(ncid, wlid, 'Wavenumber')
        call put_description(ncid, thetaid, 'Direction')
        call put_description(ncid, lon_scalarid, 'Longitude')
        call put_description(ncid, lat_scalarid, 'Latitude')
        call put_description(ncid, wspdid, 'Wind speed')
        call put_description(ncid, wdirid, 'Wind direction')
        call put_description(ncid, specid, 'Surface elevation variance spectrum')
        call put_description(ncid, sinid, 'Wind input/export rate')
        call put_description(ncid, sdsid, 'Wave breaking decay rate')
        call put_description(ncid, sdtid, 'Turbulent dissipation rate')
        call put_description(ncid, sdvid, 'Viscous dissipation rate')
        call put_description(ncid, sbfid, 'Combined bottom friction and bottom percolation rate')
        call put_description(ncid, snlid, 'Nonlinear wave-wave interaction absolute source/sink tendency')

        ! end of definition mode:
        stat = nf90_enddef(ncid)

        ! fill in static fields:
        stat = nf90_put_var(ncid, freqid, self % spectrum % frequency, start=[1], count=[om])
        stat = nf90_put_var(ncid, wlid, self % k(:,i), start=[1], count=[om])
        stat = nf90_put_var(ncid, thetaid, self % spectrum % direction, start=[1], count=[pm])
        stat = nf90_put_var(ncid, lon_scalarid, &
                            self % grid % lon(self % spectrum_m(nn), self % spectrum_n(nn)))
        stat = nf90_put_var(ncid, lat_scalarid, &
                            self % grid % lat(self % spectrum_m(nn), self % spectrum_n(nn)))
        stat = nf90_put_var(ncid, sdvid, self % sdv(:,i), start=[1], count=[om])
        stat = nf90_put_var(ncid, sbfid, self % sbf(:,i), start=[1], count=[om])

      else

        ! open file for writing:
        stat = nf90_open(trim(filename), nf90_write, ncid)

        stat = nf90_inq_varid(ncid, 'F', specid)
        stat = nf90_inq_varid(ncid, 'Sin', sinid)
        stat = nf90_inq_varid(ncid, 'Sds', sdsid)
        stat = nf90_inq_varid(ncid, 'Sdt', sdtid)
        stat = nf90_inq_varid(ncid, 'Snl', snlid)
        stat = nf90_inq_varid(ncid, 'time', timeid)
        stat = nf90_inq_varid(ncid, 'wspd', wspdid)
        stat = nf90_inq_varid(ncid, 'wdir', wdirid)

      end if

      ! fill in variables:
      time_value(1) = seconds_since_reference(timestr, self % config % reftimestr)
      stat = nf90_put_var(ncid, timeid, time_value, start=[record], count=[1])
      stat = nf90_put_var(ncid, specid, self % e(:,:,i), start=[1,1,record], count=[om,pm,1])
      stat = nf90_put_var(ncid, sinid, self % ssin(:,:,i), start=[1,1,record], count=[om,pm,1])
      stat = nf90_put_var(ncid, sdsid, self % sds(:,:,i) * spread(self % cothkd(:,i), 2, pm), &
                          start=[1,1,record], count=[om,pm,1])
      stat = nf90_put_var(ncid, sdtid, self % sdt(:,i), start=[1,record], count=[om,1])
      stat = nf90_put_var(ncid, snlid, self % snl(:,:,i), start=[1,1,record], count=[om,pm,1])
      stat = nf90_put_var(ncid, wspdid, self % forcing % wspd(i), start=[1,record])
      stat = nf90_put_var(ncid, wdirid, self % forcing % wdir(i), start=[1,record])

      ! close file:
      stat = nf90_close(ncid)

      write(*,'(a)') 'umwm: output_spectrum_nc: spectrum written to ' // trim(filename)

    end do

    end associate

    self % spectrum_records = record

  contains

    subroutine read_points()
      ! Reads the spectrum output points from the spectrum.nml next to the
      ! main namelist (namelists/spectrum.nml by default).
      integer, parameter :: max_points = 999
      integer :: m(max_points), n(max_points), npts, unit, xy(2)
      character(40) :: id(max_points)
      character(2) :: coord_input
      character(:), allocatable :: path
      real :: lon, lat

      self % spectrum_timestr = timestr
      self % spectrum_timestr(11:11) = '_'
      npts = 0

      ! open spectrum points list file:
      path = self % config % sibling_path('spectrum.nml')
      open(newunit=unit, file=path, status='old', &
           form='formatted', access='sequential', iostat=stat)

      if (stat == 0) then

        ! are points given in lat/lon or grid indices?
        read(unit, fmt='(a)') coord_input

        ! loop over points in the list:
        do

          if (any(coord_input == ['xy', 'XY'])) then

            read(unit, *, iostat=stat) m(npts+1), n(npts+1), id(npts+1)
            if (is_iostat_end(stat)) exit
            if (stat /= 0) stop 'umwm: output_nc: error: cannot read a point in ' // path

            if (m(npts+1) < 2 .or. m(npts+1) > self % config % mm - 1 .or. &
                n(npts+1) < 2 .or. n(npts+1) > self % config % nm - 1) &
              stop 'umwm: output_nc: error: a requested point in ' // path &
                // ' is out of bounds'

          else if (any(coord_input == ['ll', 'LL'])) then

            read(unit, *, iostat=stat) lon, lat, id(npts+1)
            if (is_iostat_end(stat)) exit
            if (stat /= 0) stop 'umwm: output_nc: error: cannot read a point in ' // path
            xy = minloc((lon - self % grid % lon)**2 + (lat - self % grid % lat)**2)
            m(npts+1) = xy(1)
            n(npts+1) = xy(2)

          else

            stop 'umwm: output_nc: error: first line in ' // path // ' must ' &
              // 'contain "xy" (or "XY") or "ll" (or "LL").'

          end if

          npts = npts + 1

        end do

        close(unit)

      end if

      self % spectrum_m = m(1:npts)
      self % spectrum_n = n(1:npts)
      self % spectrum_id = id(1:npts)
      allocate(self % spectrum_i(npts))
      do nn = 1, npts
        self % spectrum_i(nn) = self % grid % ii(m(nn), n(nn))
      end do

    end subroutine read_points

  end subroutine output_spectrum_nc

end submodule umwm_domain_output
