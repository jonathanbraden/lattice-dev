module modparam
  
  implicit none

  integer, parameter :: fields=1
  integer, parameter :: phi=1

! Parameters for potential 
!  real, parameter :: m2phi = 1.0
!  real, parameter :: m2psi = 0.0
  real, parameter :: mpl = 1.0e7/3.0
  real, parameter :: lambda = 1.0

! Initial conditions for homogeneous field
  real, parameter :: phi0 =   2.3393837654714997732962993666073
  real, parameter :: dphi0 = -2.7363582010758065274616992909302
  real, parameter :: H0 = 1.9348974397391251388968698880012
  real, parameter :: dH0 = 0.0
  real, parameter :: ddphi0 = -(2.0*H0*dphi0 + lambda*phi0**3)
  real, parameter :: ddH0 = -dphi0*ddphi0

! These need to be adjusted for every model
  real, parameter, dimension(fields) :: fld0 = (/ phi0 /)
  real, parameter, dimension(fields) :: dfld0 = (/ dphi0 /)
  real, parameter, dimension(fields) :: ddfld0 = (/ ddphi0 /)


end module modparam
