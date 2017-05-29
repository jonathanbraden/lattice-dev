module modparam
  
  implicit none

  integer, parameter :: fields=3
  integer, parameter :: phi=1, psi=2, sigma=3

! Parameters for potential 
  real, parameter :: m2phi = 1.0
  real, parameter :: m2psi = 0.0
  real, parameter :: m2sig = 0.
  real, parameter :: mpl = 2.0e5
  real, parameter :: mpl2 = mpl**2
  real, parameter :: g2 = 100.0**2

! Initial conditions for homogeneous field
  real, parameter :: phi0 =   1.0093430384226378929425913902459
  real, parameter :: dphi0 = -0.7137133070120812430962278466136
  real, parameter :: H0 = 0.5046715192113189464712956951230
  real, parameter :: dH0 = -H0**2
  real, parameter :: ddphi0 = -(3.0*H0*dphi0 + m2phi*phi0)
  real, parameter :: ddH0 = -dphi0*ddphi0

  real, parameter :: sig0 = 1.  ! adjust as necessary

! These need to be adjusted for every model
  real, parameter, dimension(fields) :: fld0 = (/ phi0, 0.0, sig0/)
  real, parameter, dimension(fields) :: dfld0 = (/ dphi0, 0.0, 0./)
  real, parameter, dimension(fields) :: ddfld0 = (/ ddphi0, 0.0, 0./)


end module modparam
