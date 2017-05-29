module modparam
  
  implicit none

  integer, parameter :: phi=1
  integer, parameter :: fields=1
  real, parameter :: lenrat = 10


! This is only used to set the scale of fluctuations
  real, parameter :: mpl = 1.0e7/3.0

! Potential parameters
! Laura's model
!  real, parameter :: phimin = 10.
!  real, parameter :: m2 =    ! this sets the overall scale of the potential
!  real, parameter :: gam2 = 1.   ! height of bump
!  real, parameter :: beta2 = 1.  ! width of bump

  real, parameter :: del = 1./30.   ! check
  real, parameter :: rhovac = 0.

! Since I'm scaling the field to phi0, I need to add in it's relation to mpl
  real, parameter :: phi0 = 1.0  ! assume it's just the Plank mass

! Parameters for initializing file values
! These are computed from the given potential
! Eventually, I should get these out of a numerical simulation
  real, parameter :: phifalse = 1. - del
  real, parameter :: phitrue = -1. - del

! Giblin parameter
!  real, parameter :: rinit = 3./(32.**0.5*del)

! Kosowsky parameter
  real, parameter :: rinit = 1/del


! Parameters controlling the expansion basically
! For nonhomogeneous fields, discard the phi0 and dphi0
!  real, parameter :: H0 = ((del*2.+ sig)/3.)**0.5
  real, parameter :: H0 = 1.e-3 !currently, this is just setting the expansion rate
  real, parameter :: dH0 = 0.0
  real, parameter :: ddH0 = 0.0

! Now add in the initial cond. on the Hubble and scale factor
! a0 is assumed to be 1
!  real, parameter :: aprev = 1.0 - H0*dt
!  real, parameter :: hlprev = 1.0 / (H0-dH0*dt + ddH0*dt**2/2.0)

end module modparam
