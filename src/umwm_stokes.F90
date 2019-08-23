module umwm_stokes

  implicit none

  integer :: lm
  real, allocatable :: depth(:), ds(:)
  real, allocatable :: us(:,:), vs(:,:), usmag(:,:)
  real, allocatable :: util(:,:,:), arg(:,:,:)
  real :: depths(100) = -1.

contains

  subroutine stokes_drift(option)
    ! Computes wave-induced Stokes drift

    use umwm_constants, only: eulerinv
    use umwm_module, only: twopi, d, e, f, k, dwn, dth, istart, iend, om, pm,&
                           cth, sth,nproc

    character(len=4), intent(in), optional :: option
    integer :: i, l, o, p
    integer :: coords(2)
    real :: ust_efolding
    real, allocatable :: kd(:,:)

    namelist /stokes/ depths

    if (present(option)) then
      if (option == 'init') then

        ! read depth levels from namelist
        open(unit=21, file='namelists/main.nml', status='old',&
             form='formatted', access='sequential', action='read')
        read(unit=21, nml=stokes)
        close(unit=21)

        ! get mpisize of depth array
        lm = count(depths >= 0)
        allocate(depth(lm))
        depth = - depths(1:lm)

        ! allocate stokes velocities and utility array
        allocate(us(istart:iend,lm))
        allocate(vs(istart:iend,lm))
        allocate(usmag(istart:iend,lm))
        allocate(ds(istart:iend))
        allocate(util(om,istart:iend,lm))
        allocate(arg(om,istart:iend,lm))
        allocate(kd(om,istart:iend))

        do concurrent (o=1:om, i=istart:iend) 
          kd(o,i) = k(o,i) * d(i)
        end do

        ! compute exponent
        do l = 1,lm
          do i = istart,iend
            do o = 1,om

              arg(o,i,l) = 2 * k(o,i) * (depth(l) + d(i))

              if(abs(arg(o,i,l)) > 50 .or. kd(o,i) > 50) then
                ! hyperbolic trig. functions would overflow;
                ! use deep water approximation instead
                util(o,i,l) = twopi * f(o) * 2 * k(o,i)**2 &
                            * exp(2 * k(o,i) * depth(l)) * dwn(o,i) * dth
              else
                ! first order approximation for arbitrary depth
                util(o,i,l) = twopi * f(o) * k(o,i)**2 &
                            * cosh(2 * k(o,i) * (depth(l) + d(i)))&
                            / sinh(kd(o,i))**2 * dwn(o,i) * dth
              end if

            end do
          
            if (abs(depth(l)) > d(i)) util(:,i,l) = 0

          end do
        end do

        deallocate(arg,kd)

        if (nproc == 0) write(*, fmt=101) 'umwm: stokes_drift: initialized'

      end if ! option == 'init'
    end if ! present(option)
    
    us = 0
    vs = 0
    ds = 0

    ! stokes velocities
    do l = 1, lm
      do i = istart, iend
        do p = 1, pm
          do o = 1, om
            us(i,l) = us(i,l) + util(o,i,l) * e(o,p,i) * cth(p)
            vs(i,l) = vs(i,l) + util(o,i,l) * e(o,p,i) * sth(p)
          end do
        end do
      end do
    end do

    ! Stokes drift magnitude
    usmag = sqrt(us**2 + vs**2)

    ! Stokes e-folding depth
    do i = istart, iend

      if (usmag(i,1) == 0) then
        ds(i) = 0
        continue
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
