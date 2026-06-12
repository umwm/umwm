module umwm_stokes

  implicit none

  integer :: lm
  real, allocatable :: depth(:), ds(:)
  real, allocatable :: us(:,:), vs(:,:), usmag(:,:)
  real, allocatable :: util(:,:,:), arg(:,:,:)
contains

  subroutine stokes_drift(spectrum, config, grid, option)
    ! Computes wave-induced Stokes drift

    use umwm_constants, only: eulerinv
    use umwm_config, only: config_type
    use umwm_grid, only: grid_type
    use umwm_module, only: twopi, e, f, k, dwn, dth, cth, sth,nproc
    use umwm_spectrum, only: spectrum_type

    type(spectrum_type), intent(in) :: spectrum
    type(config_type), intent(in) :: config
    type(grid_type), intent(in) :: grid
    character(4), intent(in), optional :: option
    integer :: i, l, o, p
    real :: ust_efolding
    real, allocatable :: kd(:,:)

    if (present(option)) then
      if (option == 'init') then

        ! get size of depth array
        lm = size(config % stokes_depths)
        allocate(depth(lm))
        depth = - config % stokes_depths

        ! allocate stokes velocities and utility array
        allocate(us(grid % istart:grid % iend,lm))
        allocate(vs(grid % istart:grid % iend,lm))
        allocate(usmag(grid % istart:grid % iend,lm))
        allocate(ds(grid % istart:grid % iend))
        allocate(util(spectrum % num_frequencies,grid % istart:grid % iend,lm))
        allocate(arg(spectrum % num_frequencies,grid % istart:grid % iend,lm))
        allocate(kd(spectrum % num_frequencies,grid % istart:grid % iend))

        do concurrent (o=1:spectrum % num_frequencies, i=grid % istart:grid % iend)
          kd(o,i) = k(o,i) * grid % d(i)
        end do

        ! compute exponent
        do l = 1, lm
          do i = grid % istart, grid % iend
            do o = 1, spectrum % num_frequencies

              arg(o,i,l) = 2 * k(o,i) * (depth(l) + grid % d(i))

              if(abs(arg(o,i,l)) > 50 .or. kd(o,i) > 50) then
                ! hyperbolic trig. functions would overflow;
                ! use deep water approximation instead
                util(o,i,l) = twopi * f(o) * 2 * k(o,i)**2 &
                            * exp(2 * k(o,i) * depth(l)) * dwn(o,i) * dth
              else
                ! first order approximation for arbitrary depth
                util(o,i,l) = twopi * f(o) * k(o,i)**2 &
                            * cosh(2 * k(o,i) * (depth(l) + grid % d(i)))&
                            / sinh(kd(o,i))**2 * dwn(o,i) * dth
              end if

            end do
          
            if (abs(depth(l)) > grid % d(i)) util(:,i,l) = 0

          end do
        end do

        deallocate(arg, kd)

        if (nproc == 0) write(*, fmt=101) 'umwm: stokes_drift: initialized'

      end if ! option == 'init'
    end if ! present(option)
    
    us = 0
    vs = 0
    ds = 0

    ! stokes velocities
    do l = 1, lm
      do i = grid % istart, grid % iend
        do p = 1, spectrum % num_directions
          do o = 1, spectrum % num_frequencies
            us(i,l) = us(i,l) + util(o,i,l) * e(o,p,i) * cth(p)
            vs(i,l) = vs(i,l) + util(o,i,l) * e(o,p,i) * sth(p)
          end do
        end do
      end do
    end do

    ! Stokes drift magnitude
    usmag = sqrt(us**2 + vs**2)

    ! Stokes e-folding depth
    do i = grid % istart, grid % iend

      if (usmag(i,1) == 0) then
        ds(i) = 0
        cycle
      end if

      ust_efolding = usmag(i,1) * eulerinv

      depth_loop: do l = 2, lm
        if(usmag(i,l) < ust_efolding)then
          ds(i) = (abs(usmag(i,l-1) - ust_efolding) * depth(l)  &
                 + abs(usmag(i,l) - ust_efolding) * depth(l-1)) &
                 / (usmag(i,l-1) - usmag(i,l))

          exit depth_loop
        end if
      end do depth_loop

    end do

    ds = -ds

    101 format(a)

  end subroutine stokes_drift

end module umwm_stokes
