!
!
! Define the model to run in this module
! Need to include:
!   - potential
!   - derivative of potential w.r.t each field
!   - effective masses
!   - summary of model to print in files
!
!


module model

  use modparam
  use params
  use latpar

#include "macros.h"

  implicit none
 
  contains

    subroutine modelsummary(fd)
      integer fd

      write (fd,'(a,4(f0.5,a))') "# V(phi,psi) = ", &
           lambda, "phi^4/4"      

    end subroutine modelsummary

    ! calculates the value of the potential at each point

    function modelv(hr, i, j, k)
      real modelv
      real, dimension(fields,SIRANGE) :: hr; integer i, j, k
      !		 modelv = (hr(phi,i,j,k)**2-(lambdapf*ll/2.0)*(mpl**2/m2phi)*hr(phi,i,j,k)**4+(gpf*gg**2/3.0)*(mpl**4/m2phi**2)*hr(phi,i,j,k)**6) ! This is 2*V. a factor of mpl**2 is ignored

      modelv = 0.5*lambda*hr(phi,i,j,k)**4

    end function modelv

    ! calculates the potential derivatives at each point wrt the fields. modeldv is a 2 by 1 array. First row is (dV/d\phi)/phi and second row is (dV/dpsi)/psi

    function modeldv(hr, i, j, k)
      real, dimension(fields) :: modeldv
      real, dimension(fields,SIRANGE) :: hr; integer i, j, k

      ! basis vectors for vectorizing potential derivative
      real, dimension(fields), parameter :: V1 = 1
      !real, dimension(fields), parameter :: V1 = (/1,0/), V2 = (/0,1/)

      modeldv =  lambda*hr(phi,i,j,k)**3*V1

    end function modeldv

    ! Contribution to the effective mass from the model
    function modelm2eff()
      real, dimension(fields) :: modelm2eff
      !        modelm2eff = 0.
!      modelm2eff = (/ 0.0, g2*phi0**2 /)
!      modelm2eff = (/ m2phi, m2psi + g2*phi0**2/)

      modelm2eff = (/ 3.*lambda*phi0**2 /)

      return
    end function modelm2eff


    ! This is a far more useful way to output this
    !
    ! As input, put in say the mean values of the fields (or their variances
    ! or whatever is required for the particular model.
    !
    ! It might be best to try many different definitions and see how they compare
    !
    function m2eff(fldave, fldvar)
      real, dimension(fields) :: m2eff
      real, dimension(fields) :: fldave, fldvar

      m2eff(:) = 0.

!      m2eff(phi) = m2phi + g2*fldave(psi)**2 + g2*fldvar(psi)**2
!      m2eff(psi) = g2*fldave(phi)**2 + g2*fldvar(phi)**2

      m2eff(phi) = 3.*lambda*(fldave(phi)**2+fldvar(phi)**2)
     
    end function m2eff

!
!   returns the potential for the mean field
!   Fix this later to include contribution from mean chi field
!    
    function potave(avefld)
      real :: potave
      real, dimension(fields) :: avefld

!      potave = 0.5*m2phi*avefld(phi)**2
! This one is probably a better definition
!      potave = 0.5*m2phi*avefld(phi)**2 + 0.5*g2*avefld(phi)**2*avefld(psi)**2

      potave = 0.25*lambda*avefld(phi)**4

    end function potave

    function potavefull(avefld)
      real :: potavefull
      real, dimension(fields) :: avefld

!      potavefull = 0.5*m2phi*avefld(phi)**2 + 0.5*g2*avefld(phi)**2*avefld(psi)**2

      potavefull = 0.25*lambda*avefld(phi)**4 
    end function potavefull

    function fvprime(fval)
      real, dimension(fields,fields) :: fvprime
      real, dimension(fields) :: fval
      
!      fvprime(1,1) = m2phi*fval(1)**2 + g2*fval(1)**2*fval(2)**2
!      fvprime(1,2) = g2*fval(1)**3*fval(2)
!      fvprime(2,1) = m2phi*fval(1)*fval(2) + g2*fval(1)*fval(2)**3
!      fvprime(2,2) = g2*fval(1)**2*fval(2)**2

      fvprime(1,1) = lambda*fval(1)**4
    end function fvprime

  end module model
