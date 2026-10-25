module umwm_advection

  ! Wave energy advection in geographical (x, y) and directional (theta)
  ! space for one spectral column. 1st order upstream finite differences.

  implicit none

  private
  public :: current_rotation_rate, depth_rotation_rate, propagation, refraction

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


  pure real function depth_rotation_rate(cgdlncp_dx, cgdlncp_dy, cth, sth) result(res)
    ! Rate of anti-clockwise rotation [rad/s] of the wavenumber direction th
    ! due to depth-induced refraction, -(cg/cp) dcp/dn, where n is the
    ! coordinate normal to th, to its left, and cp the phase speed at fixed
    ! frequency: Snell's law, dth/ds = -(1/cp) dcp/dn, along the ray s that
    ! the energy travels at the group speed cg.
    real, intent(in) :: cgdlncp_dx, cgdlncp_dy ! (cg/cp) grad(cp)
    real, intent(in) :: cth, sth ! cos(th), sin(th)
    res = cgdlncp_dx * sth - cgdlncp_dy * cth
  end function depth_rotation_rate


  pure real function current_rotation_rate(du_dx, du_dy, dv_dx, dv_dy, cth, sth) result(res)
    ! Rate of anti-clockwise rotation [rad/s] of the wavenumber direction th
    ! due to current-induced refraction, -(k/|k|).dU/dn, from dk/dt = -grad(k.U),
    ! where n is the coordinate normal to th, to its left. It is half the
    ! vorticity, (dv/dx - du/dy)/2, plus a part from the strain of the
    ! current that varies with direction.
    real, intent(in) :: du_dx, du_dy, dv_dx, dv_dy ! current velocity gradient
    real, intent(in) :: cth, sth ! cos(th), sin(th)
    res = sth**2 * dv_dx - cth**2 * du_dy + sth * cth * (du_dx - dv_dy)
  end function current_rotation_rate


  pure subroutine refraction(nf, nd, no, e, cgdlncp_dx, cgdlncp_dy, du_dx, du_dy, dv_dx, dv_dy, &
                             cth, sth, dt, dth, ef)
    ! Advects the spectrum e of a cell in directional space for frequencies
    ! 1:no over the time step dt and integrates the change into ef. The
    ! time step is split into as many equal substeps as the upstream scheme
    ! needs in this cell to stay stable and positive.
    integer, intent(in) :: nf, nd, no
    real, intent(in) :: e(nf,nd)
    real, intent(in) :: cgdlncp_dx(nf), cgdlncp_dy(nf) ! (cg/cp) grad(cp)
    real, intent(in) :: du_dx, du_dy, dv_dx, dv_dy ! current velocity gradient
    real, intent(in) :: cth(nd), sth(nd) ! ray directions in the grid frame
    real, intent(in) :: dt, dth
    real, intent(inout) :: ef(nf,nd)

    real :: rot(no,0:nd+1) ! ray rotation rate, periodic in direction
    real :: en(no,nd), en_next(no,nd) ! spectrum at the start and end of a substep
    real :: rot_current, rot_depth, rot_vorticity, rot_strain, outflow
    integer :: o, p, n, nsub

    do p = 1, nd
      rot_current = current_rotation_rate(du_dx, du_dy, dv_dx, dv_dy, cth(p), sth(p))
      do o = 1, no
        rot(o,p) = depth_rotation_rate(cgdlncp_dx(o), cgdlncp_dy(o), cth(p), sth(p)) + rot_current
      end do
    end do
    rot(:,0) = rot(:,nd)
    rot(:,nd+1) = rot(:,1)

    ! The upstream scheme is stable and positive if the fraction of a bin
    ! that leaves it through both faces in a substep does not exceed 1. Over
    ! direction, the rotation rate is a sinusoid in th from the depth, of
    ! amplitude rot_depth, plus half the vorticity, plus a sinusoid in 2 th
    ! from the strain of the current, of amplitude rot_strain. Per unit time,
    ! a bin loses at most the largest rotation rate through one face, and,
    ! where the rotation diverges, at most dth times the largest derivative
    ! of the rotation rate over th through both faces.
    rot_vorticity = 0.5 * abs(dv_dx - du_dy)
    rot_strain = 0.5 * sqrt((du_dx - dv_dy)**2 + (dv_dx + du_dy)**2)
    outflow = 0
    do o = 1, no
      rot_depth = sqrt(cgdlncp_dx(o)**2 + cgdlncp_dy(o)**2)
      outflow = max(outflow, rot_depth + rot_vorticity + rot_strain, &
                    dth * (rot_depth + 2 * rot_strain))
    end do
    nsub = max(1, ceiling(dt * outflow / dth))

    if (nsub == 1) then
      call upstream(nf, nd, no, rot, e, dt / dth, ef)
    else
      en = e(1:no,:)
      do n = 1, nsub
        en_next = en
        call upstream(no, nd, no, rot, en, dt / (nsub * dth), en_next)
        en = en_next
      end do
      ef(1:no,:) = ef(1:no,:) + en - e(1:no,:)
    end if

  end subroutine refraction


  pure subroutine upstream(n1, nd, no, rot, e, dt_over_dth, ef)
    ! Integrates the upstream advection of e in directional space for
    ! frequencies 1:no over dt_over_dth into ef, given the rotation rate rot.
    integer, intent(in) :: n1, nd, no
    real, intent(in) :: rot(no,0:nd+1), e(n1,nd), dt_over_dth
    real, intent(inout) :: ef(n1,nd)
    real :: rotl, rotr
    integer :: o, p, pl, pr
    do p = 1, nd
      pl = modulo(p, nd) + 1     ! anti-clockwise neighbor
      pr = modulo(p - 2, nd) + 1 ! clockwise neighbor
      do o = 1, no
        ! rotation at cell edges
        rotl = 0.5 * (rot(o,p) + rot(o,p+1))
        rotr = 0.5 * (rot(o,p) + rot(o,p-1))
        ef(o,p) = ef(o,p) - dt_over_dth * (max(rotl, 0.) * e(o,p) + min(rotl, 0.) * e(o,pl) &
                                         - max(rotr, 0.) * e(o,pr) - min(rotr, 0.) * e(o,p))
      end do
    end do
  end subroutine upstream

end module umwm_advection
