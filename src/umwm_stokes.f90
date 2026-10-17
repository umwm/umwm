module umwm_stokes
  ! Wave-induced Stokes drift.

  use umwm_constants, only: eulerinv, twopi

  implicit none

  private
  public :: stokes_type

  type :: stokes_type
    integer :: lm = 0                     ! number of depth levels
    real, allocatable :: depth(:)         ! depth levels (negative down) [m]
    real, allocatable :: weight(:,:,:)    ! (om,lm,i) spectral weights of the drift
    real, allocatable :: us(:,:), vs(:,:) ! (i,lm) drift velocity components [m/s]
    real, allocatable :: ds(:)            ! e-folding depth [m]
  contains
    procedure :: initialize
    procedure :: compute
  end type stokes_type

contains

  subroutine initialize(self, depths, f, k, dwn, d, dth, istart, iend)
    ! Precomputes the depth-dependent spectral weights of the Stokes drift.
    class(stokes_type), intent(inout) :: self
    real, intent(in) :: depths(:)                ! positive depth levels [m]
    real, intent(in) :: f(:)                     ! frequency [Hz]
    integer, intent(in) :: istart, iend
    real, intent(in) :: k(:,istart:)             ! wavenumber [rad/m]
    real, intent(in) :: dwn(:,istart:)           ! wavenumber bin width [rad/m]
    real, intent(in) :: d(istart:)               ! water depth [m]
    real, intent(in) :: dth                      ! direction bin width [rad]
    integer :: i, l, o
    real :: arg, kd

    self % lm = size(depths)
    self % depth = - depths
    allocate(self % weight(size(f), self % lm, istart:iend))
    allocate(self % us(istart:iend, self % lm), self % vs(istart:iend, self % lm))
    allocate(self % ds(istart:iend))

    do i = istart, iend
      do l = 1, self % lm
        do o = 1, size(f)
          arg = 2 * k(o,i) * (self % depth(l) + d(i))
          kd = k(o,i) * d(i)
          if (abs(arg) > 50 .or. kd > 50) then
            ! hyperbolic trig. functions would overflow;
            ! use deep water approximation instead
            self % weight(o,l,i) = twopi * f(o) * 2 * k(o,i)**2 &
                                 * exp(2 * k(o,i) * self % depth(l)) * dwn(o,i) * dth
          else
            ! first order approximation for arbitrary depth
            self % weight(o,l,i) = twopi * f(o) * k(o,i)**2 &
                                 * cosh(2 * k(o,i) * (self % depth(l) + d(i))) &
                                 / sinh(kd)**2 * dwn(o,i) * dth
          end if
        end do
        if (abs(self % depth(l)) > d(i)) self % weight(:,l,i) = 0
      end do
    end do

  end subroutine initialize


  subroutine compute(self, e, cth, sth)
    ! Computes Stokes drift velocities and e-folding depth from the spectrum.
    class(stokes_type), intent(inout) :: self
    real, intent(in) :: e(:,:,lbound(self % us, 1):) ! (om,pm,i)
    real, intent(in) :: cth(:), sth(:)
    real :: ex(size(e, 1)), ey(size(e, 1)), usmag(self % lm), ust_efolding
    integer :: i, l, p

    do i = lbound(self % us, 1), ubound(self % us, 1)

      ! direction-integrated spectrum, projected on x and y
      ex = 0
      ey = 0
      do p = 1, size(e, 2)
        ex = ex + e(:,p,i) * cth(p)
        ey = ey + e(:,p,i) * sth(p)
      end do

      do l = 1, self % lm
        self % us(i,l) = sum(self % weight(:,l,i) * ex)
        self % vs(i,l) = sum(self % weight(:,l,i) * ey)
      end do
      usmag = sqrt(self % us(i,:)**2 + self % vs(i,:)**2)

      ! Stokes e-folding depth
      self % ds(i) = 0
      if (usmag(1) == 0) cycle
      ust_efolding = usmag(1) * eulerinv
      do l = 2, self % lm
        if (usmag(l) < ust_efolding) then
          self % ds(i) = - (abs(usmag(l-1) - ust_efolding) * self % depth(l)  &
                          + abs(usmag(l) - ust_efolding) * self % depth(l-1)) &
                          / (usmag(l-1) - usmag(l))
          exit
        end if
      end do

    end do

  end subroutine compute

end module umwm_stokes
