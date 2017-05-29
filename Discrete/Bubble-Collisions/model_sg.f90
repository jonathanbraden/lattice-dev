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

      write (fd,'(a,4(f0.5,a))') "# V(phi) = sine-gordon"

    end subroutine modelsummary

! Subroutines specifying the background
#ifdef FIXEDBG
    real(dl) function get_hubble()
      get_hubble = H0
    end function get_hubble
#ifdef CONFORMAL
    real(dl) function get_scale_factor()
      get_scale_factor = 1._dl / (1._dl - H0*time)
    end function get_scale_factor

    real(dl) function integral_a2(tstart, deltat)
      real(dl) :: tstart, deltat
      integral_a2 = deltat
    end function integral_a2

    real(dl) function integral_a4(tstart, deltat)
      real(dl) :: tstart,deltat
      integral_a4 = deltat
    end function integral_a4

    real(dl) function integral_aneg2(tstart, deltat)
      real(dl) :: tstart, deltat
      integral_aneg2 = deltat
    end function integral_aneg2
#else
    real(dl) function get_scale_factor()
      get_scale_factor = exp(H0*time)
    end function get_scale_factor

    real(dl) function integral_a1(tstart, deltat)
      real(dl) :: tstart, deltat

      real(dl) :: hub
      hub = get_hubble()

      integral_a1 = exp(hub*tstart)/hub * (exp(hub*deltat)-1._dl)
!      print*,"a1",integral_a1, tstart 
!      integral_a1 = deltat
    end function integral_a1

    real(dl) function integral_a3(tstart,deltat)
      real(dl) :: tstart, deltat

      real(dl) :: hub
      hub = get_hubble()

      integral_a3 = exp(3.*hub*tstart)/(3._dl*hub)*(exp(3._dl*hub*deltat)-1._dl)
!      print*,"a3",integral_a3, tstart
!      integral_a3 = deltat
    end function integral_a3

    real(dl) function integral_aneg3(tstart,deltat)
      real(dl) :: tstart, deltat

      real(dl) :: hub
      hub = get_hubble()
      integral_aneg3 = exp(-3._dl*hub*tstart)/(3._dl*hub)*(1._dl - exp(-3._dl*hub*deltat))
!      print*,integral_aneg3, tstart, deltat
!      integral_aneg3 = dt
    end function integral_aneg3
#endif    
#endif

    ! calculates 2* the value of the potential at each point
    function modelv(hr, i, j, k)
      real modelv
      real, dimension(fields,SIRANGE) :: hr
      integer i, j, k

      real :: tmp

      tmp = hr(phi,i,j,k)

      modelv = 2.*(1-cos(tmp))
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

      modeldv = (sin(tmp)) *V1

    end function modeldv

    function m2eff(fldave, fldvar)
      real, dimension(fields) :: m2eff
      real, dimension(fields) :: fldave, fldvar

      m2eff(:) = 0.

    end function m2eff

  end module model
