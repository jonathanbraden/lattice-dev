module modparam
  
  implicit none

  integer, parameter :: phi=1
  integer, parameter :: fields=1
  real, parameter :: lenrat = 10


! This is only used to set the scale of fluctuations
  real, parameter :: mpl = 100.

! Potential parameters
! Laura's model
!  real, parameter :: phimin = 10.
!  real, parameter :: m2 =    ! this sets the overall scale of the potential
!  real, parameter :: gam2 = 1.   ! height of bump
!  real, parameter :: beta2 = 1.  ! width of bump

  real, parameter :: del = 1./10.  !1./30.   ! check
  real, parameter :: rhocorr = 0.5*del**2*(1+del/4.)**2 - del**2  ! accurate to leading order only
!  real, parameter :: rhovac = 0.
!  real, parameter :: rhocorr = -del**2 + del**3 + del**4/4.

! Since I'm scaling the field to phi0, I need to add in it's relation to mpl
  real, parameter :: phi0 = 1.0  ! assume it's just the Plank mass

! Parameters for initializing file values
! These are computed from the given potential
! Eventually, I should get these out of a numerical simulation
  real, parameter :: phifalse = -0.94564927392359022 !-1. + del/2.
  real, parameter :: phitrue = 1.0466805318046015 !1. + del/2.

! Giblin parameter
!  real, parameter :: rinit = 3./(32.**0.5*del)

! Kosowsky parameter
!  real, parameter :: rinit = 1./del

! Copeland Parameter
  real, parameter :: rinit = 2.**0.5/del


! Parameters controlling the expansion basically
! For nonhomogeneous fields, discard the phi0 and dphi0
!  real, parameter :: H0 = ((del*2.+ sig)/3.)**0.5
  real, parameter :: H0 = 0.   !0.1/rinit !currently, this is just setting the expansion rate
  real, parameter :: dH0 = 0.0
  real, parameter :: ddH0 = 0.0

! Now add in the initial cond. on the Hubble and scale factor
! a0 is assumed to be 1
!  real, parameter :: aprev = 1.0 - H0*dt
!  real, parameter :: hlprev = 1.0 / (H0-dH0*dt + ddH0*dt**2/2.0)

end module modparam
