module umwm_advection

  ! Wave energy advection in geographical (x, y) and directional (theta)
  ! space for one spectral column. 1st order upstream finite differences.

  implicit none

  private
  public :: max_rotation_rate, propagation, refraction

contains

  pure subroutine propagation(nf, nd, oc, e, e_e, e_w, e_n, e_s, cg, cg_e, cg_w, cg_n, cg_s, &
                              cth, sth, dxn, dxs, dye, dyw, currents, uc_e, uc_w, vc_n, vc_s, &
                              dt_over_area, ef)
    ! Advects the spectrum e of a cell with neighbors to the east, west,
    ! north, and south (e_e, e_w, e_n, e_s) by the group velocity and,
    ! if currents is .true., by the ocean currents; integrates into ef.
    integer, intent(in) :: nf, nd, oc
    real, intent(in) :: e(nf,nd), e_e(nf,nd), e_w(nf,nd), e_n(nf,nd), e_s(nf,nd)
    real, intent(in) :: cg(nf), cg_e(nf), cg_w(nf), cg_n(nf), cg_s(nf) ! group speed
    real, intent(in) :: cth(nd), sth(nd) ! ray directions adjusted for grid curvature
    real, intent(in) :: dxn, dxs, dye, dyw ! cell edge lengths
    logical, intent(in) :: currents
    real, intent(in) :: uc_e, uc_w, vc_n, vc_s ! double current speed at cell edges
    real, intent(in) :: dt_over_area ! 0.25 * dt / cell area
    real, intent(inout) :: ef(nf,nd)

    integer :: o, p
    real :: ae(nf), aw(nf), an(nf), as(nf)
    real :: cxp, cxm, cyp, cym, flux
    real :: feup, fedn, fwup, fwdn, fnup, fndn, fsup, fsdn

    ! The double group velocity at the east cell edge is (cg + cg_e) * cth,
    ! and so on. Group speed is positive, so the upstream direction is given
    ! by the signs of cth and sth alone.
    ae = (cg + cg_e) * dye
    aw = (cg + cg_w) * dyw
    an = (cg + cg_n) * dxn
    as = (cg + cg_s) * dxs

    feup = (uc_e + abs(uc_e)) * dye
    fedn = (uc_e - abs(uc_e)) * dye
    fwup = (uc_w + abs(uc_w)) * dyw
    fwdn = (uc_w - abs(uc_w)) * dyw
    fnup = (vc_n + abs(vc_n)) * dxn
    fndn = (vc_n - abs(vc_n)) * dxn
    fsup = (vc_s + abs(vc_s)) * dxs
    fsdn = (vc_s - abs(vc_s)) * dxs

    do p = 1, nd

      cxp = cth(p) + abs(cth(p))
      cxm = cth(p) - abs(cth(p))
      cyp = sth(p) + abs(sth(p))
      cym = sth(p) - abs(sth(p))

      do o = 1, oc

        ! advective energy flux in and out of all 4 directions
        flux = ae(o) * (cxp * e(o,p) + cxm * e_e(o,p)) & ! east, out and in
             - aw(o) * (cxp * e_w(o,p) + cxm * e(o,p)) & ! west, in and out
             + an(o) * (cyp * e(o,p) + cym * e_n(o,p)) & ! north, out and in
             - as(o) * (cyp * e_s(o,p) + cym * e(o,p))   ! south, in and out

        if (currents) &
          flux = flux + (feup * e(o,p)   + fedn * e_e(o,p) &
                      -  fwup * e_w(o,p) - fwdn * e(o,p)   &
                      +  fnup * e(o,p)   + fndn * e_n(o,p) &
                      -  fsup * e_s(o,p) - fsdn * e(o,p))

        ef(o,p) = ef(o,p) - dt_over_area * flux

      end do
    end do

  end subroutine propagation


  pure real function rotation_rate(dcp_dx, dcp_dy, curl, cth, sth) result(res)
    ! Rate of rotation of wave rays [rad/s] due to bottom- and
    ! current-induced refraction.
    real, intent(in) :: dcp_dx, dcp_dy ! phase speed gradient, halved
    real, intent(in) :: curl ! current vorticity, halved
    real, intent(in) :: cth, sth
    res = dcp_dx * sth - dcp_dy * cth + curl
  end function rotation_rate


  pure real function max_rotation_rate(nf, nd, no, dcp_dx, dcp_dy, curl, cth, sth) result(res)
    ! Largest anti-clockwise ray rotation rate over frequencies 1:no, or 0.
    integer, intent(in) :: nf, nd, no
    real, intent(in) :: dcp_dx(nf), dcp_dy(nf), curl
    real, intent(in) :: cth(nd), sth(nd)
    integer :: o, p
    res = 0
    do p = 1, nd
      do o = 1, no
        res = max(res, rotation_rate(dcp_dx(o), dcp_dy(o), curl, cth(p), sth(p)))
      end do
    end do
  end function max_rotation_rate


  pure subroutine refraction(nf, nd, no, e, dcp_dx, dcp_dy, curl, cth, sth, dtr, oneovdth, ef)
    ! Advects the spectrum e of a cell in directional space for frequencies
    ! 1:no and integrates the tendency into ef.
    integer, intent(in) :: nf, nd, no
    real, intent(in) :: e(nf,nd)
    real, intent(in) :: dcp_dx(nf), dcp_dy(nf), curl
    real, intent(in) :: cth(nd), sth(nd)
    real, intent(in) :: dtr, oneovdth
    real, intent(inout) :: ef(nf,nd)

    real :: rot(no,0:nd+1) ! ray rotation rate, periodic in direction
    real :: rotl, rotr, flux
    integer :: o, p, pl, pr

    do p = 1, nd
      do o = 1, no
        rot(o,p) = rotation_rate(dcp_dx(o), dcp_dy(o), curl, cth(p), sth(p))
      end do
    end do
    rot(:,0) = rot(:,nd)
    rot(:,nd+1) = rot(:,1)

    do p = 1, nd
      pl = modulo(p, nd) + 1     ! anti-clockwise neighbor
      pr = modulo(p - 2, nd) + 1 ! clockwise neighbor
      do o = 1, no
        ! rotation at cell edges
        rotl = 0.5 * (rot(o,p) + rot(o,p+1))
        rotr = 0.5 * (rot(o,p) + rot(o,p-1))
        flux = 0.5 * ((rotl + abs(rotl)) * e(o,p)  &
                    + (rotl - abs(rotl)) * e(o,pl) &
                    - (rotr + abs(rotr)) * e(o,pr) &
                    - (rotr - abs(rotr)) * e(o,p)) &
             * oneovdth
        ef(o,p) = ef(o,p) - dtr * flux
      end do
    end do

  end subroutine refraction

end module umwm_advection
