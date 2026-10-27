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

    !$acc parallel loop collapse(3) private(arg, kd)
    do i = istart, iend
      do l = 1, self % lm
        do o = 1, size(f)
          arg = 2 * k(o,i) * (self % depth(l) + d(i))
          kd = k(o,i) * d(i)
          if (abs(self % depth(l)) > d(i)) then
            ! below the bottom
            self % weight(o,l,i) = 0
          else if (abs(arg) > 50 .or. kd > 50) then
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

#ifdef _OPENACC
    call compute_acc(size(e, 1), size(e, 2), self % lm, lbound(self % us, 1), &
                     ubound(self % us, 1), e, cth, sth, self % weight, self % depth, &
                     self % us, self % vs, self % ds)
    return
#endif

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


  subroutine compute_acc(om, pm, lm, i1, i2, e, cth, sth, weight, depth, us, vs, ds)
    ! Computes the Stokes drift as compute, on the GPU, by a thread block for
    ! each cell, with its threads over the depth levels.
    integer, intent(in) :: om, pm, lm, i1, i2
    real, intent(in) :: e(om,pm,i1:i2), cth(pm), sth(pm), weight(om,lm,i1:i2), depth(lm)
    real, intent(inout) :: us(i1:i2,lm), vs(i1:i2,lm), ds(i1:i2)
    real :: ex, ey, sx, sy, usmag, usmag_prev, ust_efolding
    integer :: i, l, o, p

    !$acc parallel loop gang vector_length(64) private(usmag, usmag_prev, ust_efolding)
    do i = i1, i2

      ! the drift at each level, from the direction-integrated spectrum
      ! projected on x and y
      !$acc loop vector private(sx, sy, ex, ey)
      do l = 1, lm
        sx = 0
        sy = 0
        !$acc loop seq
        do o = 1, om
          ex = 0
          ey = 0
          !$acc loop seq
          do p = 1, pm
            ex = ex + e(o,p,i) * cth(p)
            ey = ey + e(o,p,i) * sth(p)
          end do
          sx = sx + weight(o,l,i) * ex
          sy = sy + weight(o,l,i) * ey
        end do
        us(i,l) = sx
        vs(i,l) = sy
      end do

      ! Stokes e-folding depth
      ds(i) = 0
      usmag_prev = sqrt(us(i,1)**2 + vs(i,1)**2)
      if (usmag_prev /= 0) then
        ust_efolding = usmag_prev * eulerinv
        !$acc loop seq
        do l = 2, lm
          usmag = sqrt(us(i,l)**2 + vs(i,l)**2)
          if (usmag < ust_efolding) then
            ds(i) = - (abs(usmag_prev - ust_efolding) * depth(l) &
                     + abs(usmag - ust_efolding) * depth(l-1)) / (usmag_prev - usmag)
            exit
          end if
          usmag_prev = usmag
        end do
      end if

    end do

  end subroutine compute_acc

end module umwm_stokes
