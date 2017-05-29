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

      write (fd,'(a,4(f0.5,a))') "# V(phi) = bubble run"

    end subroutine modelsummary

! Subroutines specifying the background
#ifdef FIXEDBG
    real(dl) function get_scale_factor()
      get_scale_factor = 1.
    end function get_scale_factor

    real(dl) function get_hubble()
      get_hubble = 0.
    end function get_hubble
#endif

    ! calculates 2* the value of the potential at each point
    function modelv(hr, i, j, k)
      real modelv
      real, dimension(fields,SIRANGE) :: hr
      integer i, j, k

      real :: tmp

      tmp = hr(phi,i,j,k)

      modelv = 0.25*(tmp**2-1.)**2
    end function modelv

    ! calculates the potential derivatives at each point wrt the fields. modeldv is a 2 by 1 array. First row is (dV/d\phi)/phi and second row is (dV/dpsi)/psi

    function modeldv(hr, i, j, k)
      real, dimension(fields) :: modeldv
      real, dimension(fields,SIRANGE) :: hr; integer i, j, k

      ! basis vectors for vectorizing potential derivative
      !        real, dimension(fields), parameter :: V1 = 1
      real, dimension(fields), parameter :: V1 = (/1./)
      real :: tmp, tmp2
      tmp = hr(phi,i,j,k)

      modeldv = 0.5*(tmp**2-1)*tmp *V1

    end function modeldv

    function m2eff(fldave, fldvar)
      real, dimension(fields) :: m2eff
      real, dimension(fields) :: fldave, fldvar

      m2eff(:) = 0.

    end function m2eff

  end module model
