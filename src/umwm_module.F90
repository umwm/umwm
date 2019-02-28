module umwm_module
!=======================================================================
!
! description: a global definitions module and main memory pool.
!
!=======================================================================
use datetime_module,only:datetime

implicit none

character(len=5),parameter :: version = '2.0.0'

! use blocking mpi routines?
logical,parameter :: mpiisblocking = .false.

! time objects
type(datetime) :: starttime,stoptime,currenttime

character(len=19) :: starttimestr_nml,stoptimestr_nml

! domain dimensions
integer :: mm  ! domain mpisize in x
integer :: nm  ! domain mpisize in y
integer :: om  ! number of frequency/wavenumber bins
integer :: pm  ! number of direction bins

integer :: im  ! unrolled domain length (sea-points only)
integer :: imm ! unrolled domain length (all, equals mm*nm)

integer :: istart,iend   ! tile exclusive range (no halo)
integer :: iistart,iiend ! tile computational range (exclusive+halo)

integer :: nproc ! process rank
integer :: mpisize  ! mpi pool mpisize
integer :: ierr  ! mpi error return status

! lengths of leftmost and rightmost columns
integer :: first_col_len,last_col_len

! stdout related variables
integer :: iip,nproc_plot,xpl,ypl

! time stepping
integer :: step,timestep

! .true. during first time step
logical :: first,firstdtg

! main control switches
logical :: isglobal,restart

! grid and bathymetry related switches
logical :: gridfromfile,topofromfile,islandsFromFile,filllakes,fillestuaries

! output related switches
integer :: outgrid,outspec,outrst
logical :: stokes

! time steps:
real :: dta,dtr,dtamin

! miscellaneous variables
real :: bf1,bf1a,bf2
real :: cgmax,cfllim
real :: delx,dely
real :: dpt,alphax,alphay,dlnf,dmin,dtg,dts,dth,dthg
real :: explim
real :: fmin,fmax,fprog
real :: fieldscale1,fieldscale2
real :: g,gustiness
real :: inv_sds_power
real :: kappa
real :: log10overz
real :: mindelx,mss_fac
real :: nu_air,nu_water
real :: oneovdth
real :: rhoa0,rhow0
real :: sbf_fac,sbp_fac,sds_fac,sds_power,sdt_fac,sfct,sin_diss1
real :: sin_diss2,sin_fac,snl_fac,sumt
real :: temp0,twopisds_fac,twonu
real :: wspd0,wdir0,uc0,vc0,z
real :: fice0,fice_lth,fice_uth

!=======================================================================

! global constants
real,parameter :: pi     = 3.141592653589793 ! pi
real,parameter :: invpi  = 1/pi              ! 1/pi
real,parameter :: twopi  = 2*pi              ! 2*pi
real,parameter :: twopi2 = twopi**2          ! 4*pi^2
real,parameter :: dr     = pi/180.           ! deg -> rad

integer,dimension(10),parameter :: &
allowedoutputtimes = [-1,0,1,2,3,4,6,8,12,24]

!=======================================================================

! allocatable arrays

! 1-dimensional arrays:

! neighbor grid indices, aliased
integer,dimension(:),allocatable :: iw,ie,is,in

! neighbor grid indices, true
integer,dimension(:),allocatable :: iiw,iie,iis,iin

! exchange indices for periodic bc
integer,dimension(:),allocatable :: i_exchange_indices

! indices m and n as functions of i
integer,dimension(:),allocatable :: mi,ni

! cut-off frequency index (maximum prognostic)
integer,dimension(:),allocatable :: oc

! directional indices, anti-clockwise and clockwise
integer,dimension(:),allocatable :: pl,pr

real,dimension(:),allocatable :: th,cth,sth
real,dimension(:),allocatable :: cth2
real,dimension(:),allocatable :: dom
real,dimension(:),allocatable :: f

! wave ray directions with grid curvature correction:
real,dimension(:,:),allocatable :: cth_curv,sth_curv

integer,dimension(:,:),allocatable :: ii
integer,dimension(:,:),allocatable :: mask
integer,dimension(:,:),allocatable :: nproc_out

! 2-dimensional, unrolled arrays:
real,dimension(:,:),allocatable :: ar_2d
real,dimension(:,:),allocatable :: curv
real,dimension(:,:),allocatable :: d_2d,dlon,dlat,dx_2d,dy_2d,alphax_2d,alphay_2d
real,dimension(:,:),allocatable :: gustu,gustv
real,dimension(:,:),allocatable :: lat,lon
real,dimension(:,:),allocatable :: x,y
real,dimension(:,:),allocatable :: rhoa_2d,rhow_2d
real,dimension(:,:),allocatable :: wspd_2d,wdir_2d
real,dimension(:,:),allocatable :: fice_2d,ficeb,ficef
real,dimension(:,:),allocatable :: uwb,vwb,uw,vw,uwf,vwf
real,dimension(:,:),allocatable :: ucb,uc_2d,ucf,vcb,vc_2d,vcf

real,dimension(:),allocatable :: ar,cd,d,dx,dy,dwd,dwl,dwp,fcutoff,mwf,pwf
real,dimension(:),allocatable :: dxn,dxs,dyw,dye
real,dimension(:),allocatable :: dcp0,dcp,dcg0,dcg
real,dimension(:),allocatable :: ht,mss,mwd,mwl,mwp,shelt
real,dimension(:),allocatable :: oneovar,oneovdx,oneovdy

real,dimension(:),allocatable :: momx,momy ! momentum in x- and y-direction
real,dimension(:),allocatable :: cgmxx,cgmxy,cgmyy ! horizontal momentum fluxes

! air and water density:
real,dimension(:),allocatable :: rhoab,rhoa,rhoaf
real,dimension(:),allocatable :: rhowb,rhow,rhowf
real,dimension(:),allocatable :: rhorat

! stability function
real,dimension(:),allocatable :: psim

! stress (momentum flux) arrays [n/m^2]
real,dimension(:),allocatable :: taux,tauy
real,dimension(:),allocatable :: taux_form,tauy_form
real,dimension(:),allocatable :: taux_skin,tauy_skin
real,dimension(:),allocatable :: taux_diag,tauy_diag
real,dimension(:),allocatable :: taux_ocntop,tauy_ocntop
real,dimension(:),allocatable :: taux_ocnbot,tauy_ocnbot
real,dimension(:),allocatable :: taux_snl,tauy_snl

! wave energy growth flux [kg/s^3]
real,dimension(:),allocatable :: epsx_atm, epsy_atm

! wave energy dissipation flux [kg/s^3]
real,dimension(:),allocatable :: epsx_ocn, epsy_ocn

! form drag components:
real,dimension(:),allocatable :: taux1,tauy1
real,dimension(:),allocatable :: taux2,tauy2
real,dimension(:),allocatable :: taux3,tauy3

! tail stress components:
real,dimension(:),allocatable :: tailatmx,tailatmy
real,dimension(:),allocatable :: tailocnx,tailocny

real,dimension(:),allocatable :: uc,vc,ustar
real,dimension(:),allocatable :: wspd,wdir
real,dimension(:),allocatable :: fice

! snl downshifting weights, used in snl routine:
real,dimension(:,:),allocatable :: bf1_renorm,bf2_renorm

! utility array used for mss in sds routine:
real,dimension(:,:),allocatable :: cth2pp

! group and phase velocities:
real,dimension(:,:),allocatable :: cg0,cp0

real,dimension(:,:),allocatable :: cothkd
real,dimension(:,:),allocatable :: dwn,invcp0
real,dimension(:,:),allocatable :: fkovg
real,dimension(:,:),allocatable :: k,k4,kdk,k3dk
real,dimension(:,:),allocatable :: l2,logl2overz,oneoverk4,psiml2
real,dimension(:,:),allocatable :: sbf,sdv,sdt,snl_arg

real,dimension(:,:,:),allocatable :: dummy,e,ef,rotl,rotr,sds,snl,ssin,sice

!=======================================================================
end module umwm_module
