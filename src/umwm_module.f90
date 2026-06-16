module umwm_module
  ! Module containing global variables and parameters
  use datetime_module, only: datetime
  use umwm_constants

  implicit none

  ! use blocking mpi routines?
  logical, parameter :: mpiisblocking = .false.

  ! time objects
  type(datetime) :: starttime, stoptime, currenttime

  ! spectrum dimensions
  integer :: om  ! number of frequency/wavenumber bins
  integer :: pm  ! number of direction bins

  integer :: nproc ! process rank
  integer :: mpisize  ! mpi pool mpisize
  integer :: ierr  ! mpi error return status

  ! stdout related variables
  integer :: iip, nproc_plot, xpl, ypl

  ! .true. during first time step
  logical :: first, firstdtg

  ! main control switches
  logical :: restart

  ! output related switches
  integer :: outgrid, outspec, outrst
  logical :: stokes

  ! time steps:
  real :: dta, dtr, dtamin

  ! miscellaneous variables
  real :: bf1,bf1a,bf2
  real :: cgmax,cfllim
  real :: delx,dely
  real :: dpt,dlnf,dmin,dtg,dts,dth,dthg
  real :: explim
  real :: fmin,fmax,fprog
  real :: fieldscale1,fieldscale2
  real :: g,gustiness
  real :: inv_sds_power
  real :: kappa
  real :: oneovdth
  real :: sbf_fac,sbp_fac,sds_fac,sds_power,sdt_fac,sfct,sin_diss1
  real :: sin_diss2,sin_fac,snl_fac,sumt
  real :: twopisds_fac
  real :: wspd0,wdir0,uc0,vc0,z
  real :: fice0,fice_lth,fice_uth

  integer, parameter :: allowedoutputtimes(10) = [-1, 0, 1, 2, 3, 4, 6, 8, 12, 24]

  ! 1-D allocatable arrays:

  ! cut-off frequency index (maximum prognostic)
  integer,dimension(:),allocatable :: oc

  ! directional indices, anti-clockwise and clockwise
  integer,dimension(:),allocatable :: pl,pr

  real,dimension(:),allocatable :: th,cth,sth
  real,dimension(:),allocatable :: cth2
  real,dimension(:),allocatable :: dom
  real,dimension(:),allocatable :: f

  real,dimension(:),allocatable :: cd,dwd,dwl,dwp,fcutoff
  real,dimension(:),allocatable :: dcp0,dcp,dcg0,dcg
  real,dimension(:),allocatable :: ht,mss,mwd,mwl,mwp,shelt

  real,dimension(:),allocatable :: momx,momy ! momentum in x- and y-direction
  real,dimension(:),allocatable :: cgmxx,cgmxy,cgmyy ! horizontal momentum fluxes
  real,dimension(:),allocatable :: physics_time_step

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

  real,dimension(:),allocatable :: ustar

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
  real,dimension(:,:),allocatable :: sbf,sdv,sdt,snl_arg,sice

  real,dimension(:,:,:),allocatable :: dummy,e,ef,rotl,rotr,sds,snl,ssin

end module umwm_module
