module modparam

#include "macros.h"  
  implicit none

  integer, parameter :: fields=1
  integer, parameter :: phi=1

! Parameters for potential 
  real, parameter :: mpl = 2.0e5
  real, parameter :: mpl2 = mpl**2
  real, parameter :: m2phi = 1.0
  real, parameter :: ll = 2.8125e-6*(mpl2/m2phi)
  real, parameter :: gg2 = 10.*ll**2*(mpl2/m2phi)**2  !10.*lambda**2

! Initial conditions for homogeneous field
  real, parameter :: phi0 = (3.*ll/5.*gg2)**0.5
  real, parameter :: pot0 = -(0.5*phi0**2-0.25*ll*phi0**4 + (gg2/6.)*phi0**6)**0.5
  real, parameter :: dphi0 = 0.
  real, parameter :: H0 = -(1./3.**0.5)*pot0
  real, parameter :: dH0 = -H0**2
  real, parameter :: ddphi0 = -3.0*H0*dphi0 -(phi0-ll*phi0**3+gg2*phi0**5)
  real, parameter :: ddH0 = -(1./2.)*dphi0**2-H0**3

! These need to be adjusted for every model
  real, parameter, dimension(fields) :: fld0 = (/ phi0/)
  real, parameter, dimension(fields) :: dfld0 = (/ 0./)
  real, parameter, dimension(fields) :: ddfld0 = (/ ddphi0/)


end module modparam
