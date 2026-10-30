module umwm_physics
  ! Integration of the source functions for one spectral column. The
  ! elemental functions evaluate one spectral bin, for the column routines
  ! and for the GPU kernels.

  use umwm_source_functions, only: power

  implicit none

  private
  public :: equilibrium_bin, integrate_bin, integrate_sources, max_source_rate, source_rate

contains

  pure real function max_source_rate(nf, nd, oc, ssin, sds, cothkd, snl_arg, sbf, sdt, sdv, sice) &
    result(res)
    ! Largest magnitude of the exponential growth/decay rate in the
    ! prognostic range, which limits the source integration time step.
    ! The breaking dissipation sds excludes the plunging breaker factor cothkd.
    integer, intent(in) :: nf, nd, oc
    real, intent(in) :: ssin(nf,nd), sds(nf,nd)
    real, intent(in) :: cothkd(nf), snl_arg(nf), sbf(nf), sdt(nf), sdv(nf), sice(nf)
    integer :: o, p

    res = 0
    do p = 1, nd
      do o = 1, oc
        res = max(res, abs(source_rate(ssin(o,p), sds(o,p) * cothkd(o) * snl_arg(o), &
                                       sbf(o), sdt(o), sdv(o), sice(o))))
      end do
    end do

  end function max_source_rate


  pure elemental real function source_rate(ssin, sds, sbf, sdt, sdv, sice) result(res)
    ! Exponential growth (positive) or decay rate of one spectral bin, from
    ! the wind input ssin, and the breaking dissipation sds, bottom friction
    ! and percolation sbf, turbulent and viscous dissipation sdt and sdv, and
    ! the sea ice attenuation sice (not positive).
    !$acc routine seq
    real, intent(in) :: ssin, sds, sbf, sdt, sdv, sice
    res = ssin - sds - sbf - sdt - sdv + sice
  end function source_rate


  pure elemental real function integrate_bin(e, dts, rate, snl) result(res)
    ! The spectrum e of one prognostic bin integrated over dts, growing or
    ! decaying exponentially at rate, with the downshifting snl.
    !$acc routine seq
    real, intent(in) :: e, dts, rate, snl
    res = e * exp(dts * rate) + dts * snl
  end function integrate_bin


  pure elemental real function equilibrium_bin(ssin, sdt, sdv, sice, sds_coef, cothkd) &
    result(res)
    ! The saturation spectrum e k^4 of one diagnostic bin in equilibrium of
    ! wind input and dissipation, raised to the power sds_power.
    !$acc routine seq
    real, intent(in) :: ssin, sdt, sdv, sice, sds_coef, cothkd
    res = max(ssin - sdt - sdv + sice, 0.) / (sds_coef * cothkd)
  end function equilibrium_bin


  pure subroutine integrate_sources(nf, nd, oc, dts, e, ssin, sds, snl, sds_coef, &
                                    sbf, sdt, sdv, sice, cothkd, oneoverk4, inv_sds_power, ef)
    ! Integrates the source functions over dts with an exponential update in
    ! the prognostic range (o <= oc), and sets the diagnostic range to the
    ! equilibrium of wind input and dissipation. The breaking dissipation sds
    ! and its coefficient sds_coef exclude the plunging breaker factor cothkd.
    ! On return, ef is the updated spectrum and e is the mean of the old and
    ! updated spectrum.
    integer, intent(in) :: nf, nd, oc
    real, intent(in) :: dts
    real, intent(inout) :: e(nf,nd)
    real, intent(in) :: ssin(nf,nd), sds(nf,nd), snl(nf,nd), sds_coef(nf,nd)
    real, intent(in) :: sbf(nf), sdt(nf), sdv(nf), sice(nf), cothkd(nf), oneoverk4(nf)
    real, intent(in) :: inv_sds_power
    real, intent(out) :: ef(nf,nd)
    real :: equilibrium(oc+1:nf,nd) ! diagnostic range in equilibrium of input and dissipation
    real :: eforward
    integer :: o, p

    do p = 1, nd
      do o = 1, oc
        eforward = integrate_bin(e(o,p), dts, &
                                 source_rate(ssin(o,p), sds(o,p) * cothkd(o), sbf(o), sdt(o), &
                                             sdv(o), sice(o)), snl(o,p))
        ef(o,p) = eforward
        e(o,p) = 0.5 * (e(o,p) + eforward)
      end do
      do o = oc + 1, nf
        equilibrium(o,p) = equilibrium_bin(ssin(o,p), sdt(o), sdv(o), sice(o), sds_coef(o,p), &
                                           cothkd(o))
      end do
    end do

    call power(size(equilibrium), equilibrium, inv_sds_power)

    do p = 1, nd
      do o = oc + 1, nf
        eforward = oneoverk4(o) * equilibrium(o,p)
        ef(o,p) = eforward
        e(o,p) = 0.5 * (e(o,p) + eforward)
      end do
    end do

  end subroutine integrate_sources

end module umwm_physics
