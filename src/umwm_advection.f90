module umwm_advection

  ! Wave action advection in geographical (x, y), directional (theta), and
  ! frequency space for one spectral column. 1st order upstream finite
  ! differences.

  implicit none

  private
  public :: current_rotation_rate, depth_rotation_rate, frequency_shift, propagation, &
            refraction

contains

  pure subroutine propagation(nf, nd, oc, e, e_e, e_w, e_n, e_s, cg, cg_e, cg_w, cg_n, cg_s, &
                              kovcg, kovcg_e, kovcg_w, kovcg_n, kovcg_s, &
                              cth, sth, dxn, dxs, dye, dyw, currents, uc_e, uc_w, vc_n, vc_s, &
                              dt_over_area, ef)
    ! Advects the wave action of a cell with neighbors to the east, west,
    ! north, and south (e_e, e_w, e_n, e_s) by the group velocity and,
    ! if currents is .true., by the ocean currents, and integrates the change
    ! of the variance spectrum e into ef. The action per unit ln(frequency)
    ! is e k/cg (kovcg), so the action flowing in from a neighbor is its
    ! spectrum times the ratio of its k/cg to that of the cell.
    integer, intent(in) :: nf, nd, oc
    real, intent(in) :: e(nf,nd), e_e(nf,nd), e_w(nf,nd), e_n(nf,nd), e_s(nf,nd)
    real, intent(in) :: cg(nf), cg_e(nf), cg_w(nf), cg_n(nf), cg_s(nf) ! group speed
    real, intent(in) :: kovcg(nf), kovcg_e(nf), kovcg_w(nf), kovcg_n(nf), kovcg_s(nf) ! k/cg
    real, intent(in) :: cth(nd), sth(nd) ! ray directions adjusted for grid curvature
    real, intent(in) :: dxn, dxs, dye, dyw ! cell edge lengths
    logical, intent(in) :: currents
    real, intent(in) :: uc_e, uc_w, vc_n, vc_s ! double current speed at cell edges
    real, intent(in) :: dt_over_area ! 0.25 * dt / cell area
    real, intent(inout) :: ef(nf,nd)

    integer :: o, p
    real :: ae(nf), aw(nf), an(nf), as(nf) ! out of the cell through each edge
    real :: ae_in(nf), aw_in(nf), an_in(nf), as_in(nf) ! into the cell through each edge
    real :: re(nf), rw(nf), rn(nf), rs(nf) ! k/cg of the neighbors over that of the cell
    real :: fe_in(nf), fw_in(nf), fn_in(nf), fs_in(nf) ! into the cell by the currents
    real :: feup, fedn, fwup, fwdn, fnup, fndn, fsup, fsdn, f_out

    ! The double group velocity at the east cell edge is (cg + cg_e) * cth,
    ! and so on.
    ae = (cg + cg_e) * dye
    aw = (cg + cg_w) * dyw
    an = (cg + cg_n) * dxn
    as = (cg + cg_s) * dxs

    ! inflows scaled by the k/cg of the neighbor over that of the cell
    do o = 1, oc
      re(o) = kovcg_e(o) / kovcg(o)
      rw(o) = kovcg_w(o) / kovcg(o)
      rn(o) = kovcg_n(o) / kovcg(o)
      rs(o) = kovcg_s(o) / kovcg(o)
    end do
    ae_in(1:oc) = ae(1:oc) * re(1:oc)
    aw_in(1:oc) = aw(1:oc) * rw(1:oc)
    an_in(1:oc) = an(1:oc) * rn(1:oc)
    as_in(1:oc) = as(1:oc) * rs(1:oc)

    ! Group speed is positive, so in each direction the action flows out
    ! through the downstream edges in x and y, and in through the upstream
    ! ones, given by the signs of cth and sth alone.
    do p = 1, nd
      if (cth(p) >= 0) then
        if (sth(p) >= 0) then
          call propagate_direction(oc, 2 * abs(cth(p)), ae, aw_in, e_w(:,p), &
                                   2 * abs(sth(p)), an, as_in, e_s(:,p), e(:,p), dt_over_area, ef(:,p))
        else
          call propagate_direction(oc, 2 * abs(cth(p)), ae, aw_in, e_w(:,p), &
                                   2 * abs(sth(p)), as, an_in, e_n(:,p), e(:,p), dt_over_area, ef(:,p))
        end if
      else
        if (sth(p) >= 0) then
          call propagate_direction(oc, 2 * abs(cth(p)), aw, ae_in, e_e(:,p), &
                                   2 * abs(sth(p)), an, as_in, e_s(:,p), e(:,p), dt_over_area, ef(:,p))
        else
          call propagate_direction(oc, 2 * abs(cth(p)), aw, ae_in, e_e(:,p), &
                                   2 * abs(sth(p)), as, an_in, e_n(:,p), e(:,p), dt_over_area, ef(:,p))
        end if
      end if
    end do

    if (.not. currents) return

    feup = (uc_e + abs(uc_e)) * dye
    fedn = (uc_e - abs(uc_e)) * dye
    fwup = (uc_w + abs(uc_w)) * dyw
    fwdn = (uc_w - abs(uc_w)) * dyw
    fnup = (vc_n + abs(vc_n)) * dxn
    fndn = (vc_n - abs(vc_n)) * dxn
    fsup = (vc_s + abs(vc_s)) * dxs
    fsdn = (vc_s - abs(vc_s)) * dxs
    f_out = feup - fwdn + fnup - fsdn
    fe_in(1:oc) = fedn * re(1:oc)
    fw_in(1:oc) = fwup * rw(1:oc)
    fn_in(1:oc) = fndn * rn(1:oc)
    fs_in(1:oc) = fsup * rs(1:oc)

    ! advective action flux by the currents in and out of all 4 directions
    do p = 1, nd
      do o = 1, oc
        ef(o,p) = ef(o,p) - dt_over_area * (f_out * e(o,p) + fe_in(o) * e_e(o,p) &
                                          - fw_in(o) * e_w(o,p) + fn_in(o) * e_n(o,p) &
                                          - fs_in(o) * e_s(o,p))
      end do
    end do

  end subroutine propagation


  pure subroutine propagate_direction(oc, cx, ax_out, ax_in, e_x, cy, ay_out, ay_in, e_y, &
                                      e, dt_over_area, ef)
    ! Integrates into ef the action of one direction flowing out of a cell
    ! through its downstream edges in x and y, and in from its upstream
    ! neighbors in x and y, with spectra e_x and e_y; cx and cy are twice
    ! the absolute direction cosines, and the a are the double group speed
    ! times the edge length, for the inflows scaled by the ratio of k/cg.
    integer, intent(in) :: oc
    real, intent(in) :: cx, ax_out(*), ax_in(*), e_x(*), cy, ay_out(*), ay_in(*), e_y(*)
    real, intent(in) :: e(*), dt_over_area
    real, intent(inout) :: ef(*)
    integer :: o
    do o = 1, oc
      ef(o) = ef(o) - dt_over_area * (cx * (ax_out(o) * e(o) - ax_in(o) * e_x(o)) &
                                    + cy * (ay_out(o) * e(o) - ay_in(o) * e_y(o)))
    end do
  end subroutine propagate_direction


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

  pure subroutine frequency_shift(nf, nd, no, e, kovcg, invcp, cg, cgdlncp_dx, cgdlncp_dy, &
                                  uc, vc, du_dx, du_dy, dv_dx, dv_dy, cth, sth, dt, dlnf, ef)
    ! Advects the wave action of a cell in frequency for frequencies 1:no
    ! over the time step dt and integrates the change of the variance
    ! spectrum e into ef. Along the rays, the intrinsic frequency sigma of
    ! the waves changes where the current flows across the depth contours,
    ! and where it varies along the waves:
    !
    !   d(ln sigma)/dt = (1/sigma) (dsigma/dd) U.grad(d) - (cg/cp) (k/|k|).grad(U).(k/|k|)
    !
    ! where (dsigma/dd) grad(d) = k (cg/cp) grad(cp) at fixed sigma. The
    ! action per unit ln(sigma) is e k/cg (kovcg). Action leaves the lowest
    ! bin through its lower face, and flows between the highest bin and the
    ! diagnostic tail above it. The time step is split into as many equal
    ! substeps as the upstream scheme needs in this cell to stay stable and
    ! positive.
    integer, intent(in) :: nf, nd, no
    real, intent(in) :: e(nf,nd)
    real, intent(in) :: kovcg(nf), invcp(nf), cg(nf) ! k/cg, 1/cp, cg
    real, intent(in) :: cgdlncp_dx(nf), cgdlncp_dy(nf) ! (cg/cp) grad(cp)
    real, intent(in) :: uc, vc ! current velocity
    real, intent(in) :: du_dx, du_dy, dv_dx, dv_dy ! current velocity gradient
    real, intent(in) :: cth(nd), sth(nd) ! ray directions in the grid frame
    real, intent(in) :: dt, dlnf
    real, intent(inout) :: ef(nf,nd)

    real :: rate(0:no+1,nd) ! rate of change of ln(sigma), extended past both ends
    real :: a(0:no+1,nd)    ! action per unit ln(sigma), 0 below the lowest bin
    real :: a_next(no,nd)
    real :: dt_over_dlnf_kovcg(no) ! the substep over the bin width, over k/cg
    real :: strain, strain_max, outflow, dt_over_dlnf
    integer :: o, p, n, nsub, nt

    nt = min(no + 1, nf) ! the highest bin with a rate, the first of the tail if any

    ! The upstream scheme is stable and positive if the fraction of a bin
    ! that leaves it through both faces in a substep does not exceed 1. Per
    ! unit time, a bin loses at most the largest rate over frequency and
    ! direction: the strain rate (k/|k|).grad(U).(k/|k|) is half the
    ! divergence plus a sinusoid in 2 th of amplitude strain_max - |div|/2.
    strain_max = 0.5 * abs(du_dx + dv_dy) + 0.5 * sqrt((du_dx - dv_dy)**2 + (du_dy + dv_dx)**2)
    outflow = 0
    do o = 1, nt
      outflow = max(outflow, invcp(o) * (abs(uc * cgdlncp_dx(o) + vc * cgdlncp_dy(o)) &
                                         + cg(o) * strain_max))
    end do
    if (outflow <= 0) return
    nsub = max(1, ceiling(dt * outflow / dlnf))
    dt_over_dlnf = dt / (nsub * dlnf)

    ! Past the ends, the rate is that of the end bins, so that the rate at
    ! the faces is the mean of the rates of the bins on either side; there is
    ! no action below the lowest bin, nor above the highest without a tail.
    a(0,:) = 0
    a(no+1,:) = 0
    do p = 1, nd
      strain = cth(p)**2 * du_dx + sth(p) * cth(p) * (du_dy + dv_dx) + sth(p)**2 * dv_dy
      do o = 1, nt
        rate(o,p) = invcp(o) * (uc * cgdlncp_dx(o) + vc * cgdlncp_dy(o) - cg(o) * strain)
        a(o,p) = kovcg(o) * e(o,p)
      end do
      rate(0,p) = rate(1,p)
      if (nt == no) rate(no+1,p) = rate(no,p)
    end do

    if (nsub == 1) then
      dt_over_dlnf_kovcg = dt_over_dlnf / kovcg(1:no)
      do p = 1, nd
        do o = 1, no
          ef(o,p) = ef(o,p) - dt_over_dlnf_kovcg(o) * face_fluxes(o, p)
        end do
      end do
    else
      do n = 1, nsub
        do p = 1, nd
          do o = 1, no
            a_next(o,p) = a(o,p) - dt_over_dlnf * face_fluxes(o, p)
          end do
        end do
        a(1:no,:) = a_next
      end do
      do p = 1, nd
        do o = 1, no
          ef(o,p) = ef(o,p) + a(o,p) / kovcg(o) - e(o,p)
        end do
      end do
    end if

  contains

    pure real function face_fluxes(o, p) result(res)
      ! Upstream action flux out of bin o through its upper face, less that
      ! in through its lower face.
      integer, intent(in) :: o, p
      real :: rate_up, rate_dn
      rate_up = 0.5 * (rate(o,p) + rate(o+1,p))
      rate_dn = 0.5 * (rate(o-1,p) + rate(o,p))
      res = max(rate_up, 0.) * a(o,p) + min(rate_up, 0.) * a(o+1,p) &
          - max(rate_dn, 0.) * a(o-1,p) - min(rate_dn, 0.) * a(o,p)
    end function face_fluxes

  end subroutine frequency_shift

end module umwm_advection
