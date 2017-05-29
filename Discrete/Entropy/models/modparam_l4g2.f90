module modparam
  
  implicit none

  integer, parameter :: fields=2
  integer, parameter :: phi=1, psi=2

! Parameters for potential 
!  real, parameter :: m2phi = 1.0
!  real, parameter :: m2psi = 0.0
  real, parameter :: mpl = 1.0e7/3.0
  real, parameter :: lambda = 1.0
  real, parameter :: g2 = 200.  ! this is g^2/lambda

! Initial conditions for homogeneous field
  real, parameter :: phi0 =   2.3393837654714997732962993666073
  real, parameter :: dphi0 = -2.7363582010758065274616992909302
  real, parameter :: H0 = 1.9348974397391251388968698880012
  real, parameter :: dH0 = 0.0
  real, parameter :: ddphi0 = -(2.0*H0*dphi0 + lambda*phi0**3)
  real, parameter :: ddH0 = -dphi0*ddphi0

  real, parameter :: psi0 = 3.9e-7

! These need to be adjusted for every model
  real, parameter, dimension(fields) :: fld0 = (/ phi0, psi0 /)
  real, parameter, dimension(fields) :: dfld0 = (/ dphi0, 0.0/)
  real, parameter, dimension(fields) :: ddfld0 = (/ ddphi0, 0.0/)


end module modparam
