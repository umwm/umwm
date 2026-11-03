module umwm_partitions

  ! Wave partitions: the integrated parameters of the wind sea and swell
  ! systems of a wave field, read from the input files, and the parametric
  ! spectra that rebuild the wave spectrum from them, for the initial state
  ! and the open boundaries of the model.
  !
  ! Each partition p (0 to max_partitions - 1; by convention 0 is the wind
  ! sea) of an input file is a set of variables on the grid named after it:
  !
  !   phs<p>   significant wave height [m]; negative where there are no wave
  !            data, such as on land, and 0 where the partition is absent
  !   ptm<p>   mean period Tm-1,0 = m-1/m0 [s], or instead
  !   ptp<p>   peak period [s]
  !   pdir<p>  mean direction of propagation, counterclockwise from east (x),
  !            in the direction convention of the model [rad]
  !   pspr<p>  directional spread sqrt(2 (1 - M1)) [rad] (optional)
  !   pgam<p>  peak enhancement factor of the JONSWAP spectrum (optional)
  !
  ! The spectrum of a partition is a JONSWAP spectrum in frequency, whose
  ! peak frequency is found so that the mean period of the discrete
  ! spectrum on the frequencies of the model is the given one, times a
  ! cos^2s((theta - mean direction) / 2) distribution in direction, whose
  ! s = 2 / spread**2 - 1 has the given spread (Kuik et al., 1988). Its
  ! variance is distributed over the bins in the weights k dk dtheta of the
  ! model, so that the significant wave height of the discrete spectrum is
  ! the given one exactly, whatever the depth and the spectral resolution.

  use netcdf

  implicit none

  private
  public :: partitions_type, max_partitions, partition_spectrum, &
            period_mean, period_peak, default_spread, default_gamma

  integer, parameter :: max_partitions = 10
  integer, parameter :: period_mean = 1 ! mean period Tm-1,0
  integer, parameter :: period_peak = 2 ! peak period
  real, parameter :: default_spread = 0.5 ! [rad], s = 7
  real, parameter :: default_gamma = 3.3
  real, parameter :: min_spread = 0.01 ! [rad], s = 20000

  type :: partitions_type
    ! The partitions of one input file, on the grid (m,n).
    integer :: num = 0                        ! number of partitions
    integer :: id(max_partitions) = -1        ! their p in the variable names
    integer :: period_kind(max_partitions) = 0
    real, allocatable :: hs(:,:,:), period(:,:,:), dir(:,:,:) ! (m,n,partition)
    real, allocatable :: spread(:,:,:), gamma(:,:,:)
  contains
    procedure :: read => partitions_read
    procedure :: valid => partitions_valid
    procedure :: spectrum => partitions_spectrum
    procedure :: finalize => partitions_finalize
  end type partitions_type

contains

  subroutine partitions_read(self, ncid, path, mm, nm)
    ! Reads the partitions in the open input file ncid, of path, on the grid
    ! of mm by nm cells. num is 0 if the file has none.
    class(partitions_type), intent(inout) :: self
    integer, intent(in) :: ncid
    character(*), intent(in) :: path
    integer, intent(in) :: mm, nm
    integer :: p, j, varid

    call self % finalize()

    do p = 0, max_partitions - 1
      if (nf90_inq_varid(ncid, 'phs' // digit(p), varid) == nf90_noerr) then
        self % num = self % num + 1
        self % id(self % num) = p
      end if
    end do
    if (self % num == 0) return

    allocate(self % hs(mm,nm,self % num), self % period(mm,nm,self % num), &
             self % dir(mm,nm,self % num), self % spread(mm,nm,self % num), &
             self % gamma(mm,nm,self % num))

    do j = 1, self % num
      p = self % id(j)
      call read_field('phs' // digit(p), self % hs(:,:,j))
      if (has('ptm' // digit(p))) then
        self % period_kind(j) = period_mean
        call read_field('ptm' // digit(p), self % period(:,:,j))
      else if (has('ptp' // digit(p))) then
        self % period_kind(j) = period_peak
        call read_field('ptp' // digit(p), self % period(:,:,j))
      else
        call fail('no period ptm' // digit(p) // ' or ptp' // digit(p) // ' for phs' // digit(p))
      end if
      if (.not. has('pdir' // digit(p))) call fail('no direction pdir' // digit(p) // &
                                                   ' for phs' // digit(p))
      call read_field('pdir' // digit(p), self % dir(:,:,j))
      if (has('pspr' // digit(p))) then
        call read_field('pspr' // digit(p), self % spread(:,:,j))
      else
        self % spread(:,:,j) = default_spread
      end if
      if (has('pgam' // digit(p))) then
        call read_field('pgam' // digit(p), self % gamma(:,:,j))
      else
        self % gamma(:,:,j) = default_gamma
      end if
    end do

  contains

    logical function has(name)
      character(*), intent(in) :: name
      integer :: id
      has = nf90_inq_varid(ncid, name, id) == nf90_noerr
    end function has

    subroutine read_field(name, field)
      character(*), intent(in) :: name
      real, intent(out) :: field(mm,nm)
      integer :: id, stat
      stat = nf90_inq_varid(ncid, name, id)
      if (stat == nf90_noerr) stat = nf90_get_var(ncid, id, field, start=[1, 1, 1], &
                                                  count=[mm, nm, 1])
      if (stat /= nf90_noerr) call fail('cannot read ' // name // ': ' // trim(nf90_strerror(stat)))
    end subroutine read_field

    subroutine fail(message)
      character(*), intent(in) :: message
      write(*, '(a)') 'umwm: forcing: error: ' // trim(path) // ': ' // message
      stop 1
    end subroutine fail

  end subroutine partitions_read


  pure character(1) function digit(p)
    integer, intent(in) :: p
    digit = achar(iachar('0') + p)
  end function digit


  pure logical function partitions_valid(self, m, n) result(valid)
    ! Whether the file has wave data at cell (m,n): a partition there with
    ! a significant wave height that is not negative.
    class(partitions_type), intent(in) :: self
    integer, intent(in) :: m, n
    valid = .false.
    if (self % num > 0) valid = any(self % hs(m,n,:) >= 0)
  end function partitions_valid


  pure subroutine partitions_spectrum(self, m, n, om, pm, f, th, dth, kdk, e)
    ! The spectrum e of the partitions at cell (m,n), the sum of their
    ! parametric spectra, in the bins of frequencies f, directions th,
    ! directional spacing dth, and weights kdk = k dk of the cell; tiny where
    ! it has no partitions.
    class(partitions_type), intent(in) :: self
    integer, intent(in) :: m, n, om, pm
    real, intent(in) :: f(om), th(pm), dth, kdk(om)
    real, intent(out) :: e(om,pm)
    integer :: j

    e = 0
    do j = 1, self % num
      call partition_spectrum(om, pm, f, th, dth, kdk, self % hs(m,n,j), self % period(m,n,j), &
                              self % period_kind(j), self % dir(m,n,j), self % spread(m,n,j), &
                              self % gamma(m,n,j), e)
    end do
    e = max(e, tiny(e))

  end subroutine partitions_spectrum


  pure subroutine partition_spectrum(om, pm, f, th, dth, kdk, hs, period, period_kind, dir, &
                                     spread, gamma, e)
    ! Adds to e the parametric spectrum of one partition of significant wave
    ! height hs, mean (period_kind = period_mean, Tm-1,0) or peak
    ! (period_peak) period, mean direction dir, directional spread, and
    ! JONSWAP peak enhancement gamma, in the bins of logarithmically spaced
    ! frequencies f, directions th of spacing dth, and weights kdk = k dk.
    ! A partition without energy or period adds nothing; a spread or gamma
    ! that is not positive is replaced by its default.
    integer, intent(in) :: om, pm
    real, intent(in) :: f(om), th(pm), dth, kdk(om)
    real, intent(in) :: hs, period, dir, spread, gamma
    integer, intent(in) :: period_kind
    real, intent(inout) :: e(om,pm)
    real :: shape(om), weight(om), spreading(pm), fp, gam, s, variance, total, c_max
    real :: c(pm)
    integer :: p

    if (.not. (hs > 0 .and. period > 0)) return

    gam = default_gamma
    if (gamma > 0) gam = max(gamma, 1.)
    fp = peak_frequency(om, f, period, period_kind, gam)
    call jonswap(om, f, fp, gam, shape)

    ! fraction of the variance in each frequency bin, of width f dlnf
    weight = f * shape
    total = sum(weight)
    if (.not. total > 0) return
    weight = weight / total

    ! cos^2s((theta - dir) / 2), relative to its largest value on the bins
    ! so that it does not underflow there
    s = default_spread
    if (spread > 0) s = max(spread, min_spread)
    s = max(2 / s**2 - 1, 0.)
    c = abs(cos(0.5 * (th - dir)))
    c_max = maxval(c)
    do p = 1, pm
      if (s == 0) then
        spreading(p) = 1
      else if (c(p) > 0) then
        spreading(p) = exp(2 * s * log(c(p) / c_max))
      else
        spreading(p) = 0
      end if
    end do
    spreading = spreading / sum(spreading)

    variance = hs**2 / 16
    do p = 1, pm
      e(:,p) = e(:,p) + variance * weight * spreading(p) / (kdk * dth)
    end do

  end subroutine partition_spectrum


  pure real function peak_frequency(om, f, period, period_kind, gamma) result(fp)
    ! The peak frequency of the JONSWAP spectrum of peak enhancement gamma
    ! with the given peak period, or with the given mean period Tm-1,0 on
    ! the discrete frequencies f; within [f(1)/2, f(om)] if the mean period
    ! cannot be reached there.
    integer, intent(in) :: om, period_kind
    real, intent(in) :: f(om), period, gamma
    real :: shape(om), ratio, m0, fp_new
    integer :: iteration

    if (period_kind == period_peak) then
      fp = 1 / period
      return
    end if

    ! Tm-1,0 is near 0.86/fp for PM and 0.90/fp for JONSWAP; the mean period
    ! of the discrete spectrum scales nearly as 1/fp, which the iteration
    ! corrects for.
    fp = min(max(0.9 / period, 0.5 * f(1)), f(om))
    do iteration = 1, 50
      call jonswap(om, f, fp, gamma, shape)
      ! with bins of width f dlnf, m-1 ~ sum(shape) and m0 ~ sum(f * shape)
      m0 = sum(f * shape)
      if (.not. m0 > 0) exit
      ratio = sum(shape) / m0 / period
      fp_new = min(max(fp * ratio, 0.5 * f(1)), f(om))
      if (abs(ratio - 1) < 1e-5 .or. fp_new == fp) exit
      fp = fp_new
    end do

  end function peak_frequency


  pure subroutine jonswap(om, f, fp, gamma, shape)
    ! The JONSWAP spectrum of peak frequency fp and peak enhancement gamma,
    ! not normalized, at the frequencies f; the Pierson-Moskowitz spectrum
    ! for gamma = 1.
    integer, intent(in) :: om
    real, intent(in) :: f(om), fp, gamma
    real, intent(out) :: shape(om)
    real :: sigma
    integer :: o

    do o = 1, om
      if (f(o) <= fp) then
        sigma = 0.07
      else
        sigma = 0.09
      end if
      shape(o) = exp(-1.25 * (fp / f(o))**4 - 5 * log(f(o) / fp)) &
               * gamma**exp(-0.5 * ((f(o) - fp) / (sigma * fp))**2)
    end do

  end subroutine jonswap


  subroutine partitions_finalize(self)
    class(partitions_type), intent(inout) :: self
    if (allocated(self % hs)) deallocate(self % hs, self % period, self % dir, &
                                         self % spread, self % gamma)
    self % num = 0
    self % id = -1
    self % period_kind = 0
  end subroutine partitions_finalize

end module umwm_partitions
