module modparam
  
  implicit none

  integer, parameter :: phi=1
  integer, parameter :: fields=1

! This is only used to set the scale of fluctuations
  real, parameter :: mpl = 1.0e7/3.0

! Parameters for initializing file values
! These are computed from the given potential
! Eventually, I should get these out of a numerical simulation
  real, parameter :: phifalse = 1.
  real, parameter :: phitrue = -1.

! Parameters controlling the expansion basically
! For nonhomogeneous fields, discard the phi0 and dphi0
!  real, parameter :: H0 = ((del*2.+ sig)/3.)**0.5
  real, parameter :: H0 = 1.e-3 !currently, this is just setting the expansion rate
  real, parameter :: dH0 = 0.0
  real, parameter :: ddH0 = 0.0

end module modparam
