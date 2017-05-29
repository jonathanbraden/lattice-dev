!
! When everything is done, I only want user adjustable parameters in here !!! (everything else should be hidden)
!

module params

#ifndef L_VALUE
#define L_VALUE 13.*rinit   !lenrat*rinit   !400.0   !10.0
#endif

#ifndef ALPH
#define ALPH  25.   !5.   !400./(128.*0.005)  ! AV   !1280.0
#endif

! These are all stored here to prevent me from needing many copies of each

!#ifndef 
!#define 
!#endif

#ifndef CPOW
#define CPOW CP  ! 6
#endif

use modparam
use mpi

implicit none

include "fftw3.f"

integer, parameter :: dl = kind(1.d0)

! some useful constants
real(dl), parameter :: twopi = 6.2831853071795864769252867665590_dl
real(dl), parameter :: sqrt3 = 1.7320508075688772935274463415059_dl


! solver control parameters
integer, parameter :: cubepow = 10
integer, parameter :: n = 2**cubepow            ! sampled grid size (simulation cube is n^3 pts
integer, parameter :: p = n+2                   ! padded grid size (>n, adjust for cache efficiency)
integer, parameter :: tt = 2**14                ! total number of time steps to take (i.e. runtime)
integer, parameter :: nn = n/2+1                ! Nyquist frequency (calculated, leave it alone)
integer, parameter :: ns = sqrt3*(n/2) + 2      ! highest wavenumber on 3D grid (leave it alone)

integer, parameter :: ntser = 2**6              ! number of time steps to store output for later time series analysis

!real, parameter :: alpha = 80.0                 ! dx/dt (be careful not to violate Courant condition)
!real, parameter :: alpha = 10.
real(dl), parameter :: alpha = ALPH
real(dl), parameter :: len = L_VALUE
real(dl), parameter :: dx = len/n                  ! grid spacing   (physical grid size is n*dx)
real(dl), parameter :: dt = dx/alpha !0.0001_dl !dx/alpha          ! time step size (simulated timespan is tt*dt)
real(dl), parameter :: dk = twopi/(n*dx)            ! frequency domain grid spacing (leave it alone)

integer, parameter, dimension(1:3) :: nside = (/ n, n, n /)
integer, parameter, dimension(1:3) :: nnside = (/ nside(1)/2+1, nside(2)/2+1, nside(3)/2+1 /)
integer, parameter :: ns_2 = 2.**0.5*(nside(1)/2) + 2  ! only useful if x and y directions are the same
integer, parameter :: nvol = nside(1)*nside(2)*nside(3)
real(dl), parameter :: rnvol = dble(nside(1))*dble(nside(2))*dble(nside(3))
real(dl), parameter, dimension(1:3) :: lside = (/ len, len, len /)
real(dl), parameter :: lenvol = lside(1)*lside(2)*lside(3)
real(dl), parameter :: dk2_eff = twopi/(nside(1)*dx)
! output control parameters
!integer, parameter :: nx = n                  ! spatial grid is downsampled to nx^3 pts for output
integer, parameter :: ds = 4			! the down-sample factor
integer, parameter :: nt=2**4  !nt = 2**6                 ! simulation will be logged every nt time steps


!real, parameter :: c3 = 0.0, c2 = 0.0, c1 = 1.0, c0 = -6.0, cc = 1.0
!real, parameter :: c3 = 0.0, c2 = 1.0, c1 = 2.0, c0 = -24.0, cc = 6.0
!real, parameter :: c3 = 1.0, c2 = 0.0, c1 = 8.0, c0 = -56.0, cc = 12.0
real(dl), parameter :: c3 = 1.0_dl, c2 = 3.0_dl, c1 = 14.0_dl, c0 = -128.0_dl, cc = 30.0_dl


! Output flags (set 1)
logical, parameter :: output = .true.           ! set this to false to disable all file output at once
logical, parameter :: oscale = .true.  !.true.          ! scale output variables to counter-act expansion

logical, parameter :: output$fvar = .true.      ! output the field variances


logical, parameter :: output$gnu = .true.       ! output curves in gnuplot format (sinle file)
logical, parameter :: output$vis = .false.       ! output curves in VisIt X-Y format (per frame)

logical, parameter :: write$checkpt = .false.
logical, parameter :: resume$run = .false.



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! fields are referred to by their symbolic aliases
!integer, parameter :: fields = 2                ! total number of scalar fields being evolved (now in model file)
!integer, parameter :: phi = 1, psi = 2          ! symbolic aliases for scalar fields


integer :: outframenum

complex, allocatable, dimension(:,:,:) :: Fk, Fk2


!integer :: mpisize
integer :: mpistatus(MPI_STATUS_SIZE)
integer :: mpierror

!
! Storage and parameters for entropy calculations
!
integer, parameter :: nument = 6
real, dimension(nument) :: entropy

real, dimension(fields) :: avef, avefd, varf, varfd, davef, f2
real, dimension(fields) :: omege_eff


!!!!!!
!
! Stuff to move into FRWmetric class
!
!!!!!!

! Needs to be changed accordingly.

! initial conditions for homogeneous field component
!real, parameter ::  phi0 =  SQRT((m2phi/mpl**2)*(3.0/5.0)*(ll/gg**2))
!real, parameter :: dphi0 =  -SQRT((1./2.0)*phi0**2-(ll/4.0)*(mpl**2/m2phi)*phi0**4+(mpl**4/m2phi**2)*(gg**2/6.0)*phi0**6)
!real, parameter ::    H0 =   -(1.0/SQRT(2.0))*dphi0 
!real, parameter ::   dH0 =   -H0**2
!real, parameter :: ddphi0 =  -3.*H0*dphi0-(phi0-ll*(mpl**2/m2phi)*phi0**3+(gg**2)*(mpl**4/m2phi**2)*phi0**5)
!real, parameter ::  ddH0 = -(1./2.)*dphi0**2-H0**3

!real, parameter :: phi0 =   1.0093430384226378929425913902459
!real, parameter :: dphi0 = -0.7137133070120812430962278466136
!real, parameter :: H0 = 0.5046715192113189464712956951230
!real, parameter :: dH0 = -H0**2
!real, parameter :: ddphi0 = -(3.0*H0*dphi0 + m2phi*phi0)
!real, parameter :: ddH0 = -dphi0*ddphi0

! scale factor and horizon size (sampled on two subsequent time slices)
!real :: LA(2) = (/ 1.0 - H0*dt, 1.0 /)
!real :: LA(2) = (/ 1.0 - H0*dt - dt**2*(dphi0**2+lambda*phi0**4)/12., 1.0 /)
!real :: LA(2) = (/ 1.0 - H0*dt, 1.0 /)
!real :: LH(2) = (/ 1.0, 1.0/)    ! for no expansion
!real :: LH(2) = 1.0/(/ H0 - dH0*dt + ddH0*dt**2/2.0, H0 /)
!real :: LH(2) = (/ H0-dH0*dt+ddH0*dt**2/2., H0 /)




!!!!!!!
!
! Stuff to move into the lattice parameters class
!
!!!!!!

real :: tbegin
real :: tconf   ! add the initial one to output


!!!!!!!!!
!
! Stuff to move into some sort of output class
!
!!!!!!!!!


real, parameter :: numstat = 2*fields+15  ! number of statistics to store in CDF and PSD, this should probably be overestimated if anything since it's cheap storage (problem is that I open numstat PDF files!!!)

real PD(n+1,numstat)



integer, parameter :: numcrossstat = 15
real PSD_cr(ns, numcrossstat)
character(12) :: DVAR_cr(numcrossstat)

! Storage for the various files.  Above they are hardcoded, currently I'm doing this better in a loop
integer, parameter :: numtraj = 16
integer, parameter :: numtstat = 7

integer, parameter :: numblock = cubepow
integer, parameter :: numrgstat = 10

real, allocatable, dimension(:,:,:) :: tmp, tmp2

real :: PSD_EM(ns,10), CDF_EM(n+1,10), PDF_EM(n+1, 10)
character(12), dimension(10), parameter :: DVAR_EM = (/"T00-rho_hom ","T01    ","T02    ","T03    ","T11-P_iso ","T12    ","T13    ","T22-P_iso  ","T23   ","T33-P_iso "/)

! Parameters controlling trajectories
real, dimension(numtstat, numtraj) :: TVAR

! edit plotTrajectories subroutine if you change this!!!
integer, allocatable :: trajloc(:,:)
integer :: ntlocal

real :: rhoave, prsave, rhozero, prszero, rhozero2, prszero2

real :: npart(ns,fields), denspart(fields)
real, dimension(fields) :: effm
real, dimension(ns,fields) :: w2eff
real, allocatable :: fldspec(:,:), flddspec(:,:)

integer, parameter :: numpent = 4
real, dimension(numpent, fields) :: ent

!integer, parameter :: numblock = cubepow

real, dimension(0:numblock, numtraj) :: rgtvar

real, dimension(:,:), allocatable :: rgentropy
real, dimension(0:numblock,0:1) :: rgent2
real, dimension(0:numblock,0:1) :: rgent1

real, dimension(0:numblock,1:n+1,1:2) :: rgpivots

integer :: rgcomm(numblock)  ! Can make this smaller, but why bother
integer :: rgstart
integer :: xstart(numblock), xend(numblock), kend(numblock)  ! labels are x,y coord, data type, rg level

integer :: rgxs(2,3,numblock), rgxe(2,3,numblock)
integer :: ncount(numblock)

! Arrays to store stuff for rgblock sorting
integer :: rgsortcomm(1:numblock)
!integer, parameter :: numrgstat = 10
real, dimension(1:numrgstat, 1:n+1, 0:numblock) :: CDF_RG !, PD_RG
real, dimension(1:n+1, 0:numblock) :: PD_RG !,CDF_RG
!integer :: idrg

real, allocatable :: tmpsort(:), tmpsort2(:)
integer :: npiv(0:numblock)

! Arrays for pdf's of differences and ratios
!real, dimension(1:n+1, 0:numblock-1) :: CDF_DIFF, CDF_RAT
real, dimension(numrgstat, 1:n+1, 0:numblock-1) :: CDF_DIFF, CDF_RAT
real :: PDF_DIFF(n+1,0:numblock-1)
real :: PDF_RAT(n+1,0:numblock-1)
character(12), dimension(1:numrgstat) :: DVAR_RG

real, dimension(0:numblock,2) :: sigrg, meanrg  ! 2nd index is for rho vs. lnrho
integer, parameter :: maxmom = 6
real, dimension(1:numrgstat,1:maxmom,0:numblock,1:2) :: rgmom

real, dimension(0:numblock,3) :: trajrg

integer :: seednum

integer :: slicerank

! Very temporary array to store the various RG levels before I output them
integer, parameter :: winlevs=4
real, dimension(1:winlevs), parameter :: winscl = (/ 0.1*dx, 2.*dx, 5.*dx, 10.*dx/)
real :: rgslice(1:1, 1:winlevs, 1:n, 1:n)  ! first column stores different stats


end module params
