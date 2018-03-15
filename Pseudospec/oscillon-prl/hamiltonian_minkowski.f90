module Hamiltonian
#include "macros.h"
  use fftw3
  use params

  implicit none

  real(dl), parameter :: V0 = 1._dl/72._dl  ! Time in units of V''(v), field in units of v
  real(dl), parameter :: mpl = 1.e7/3._dl  ! only need for fluctuations
  integer, parameter :: nfld=1

  integer, parameter :: n_Hamiltonian_terms = 2

! Evolution variables
  real(dl), dimension(nfld, SIRANGE) :: fld
  real(dl), dimension(nfld, IRANGE) :: fldp
  real(dl) :: yscl, ysclp

  type(C_PTR) :: planf, planb
#ifdef THREEDIM
  real(C_DOUBLE), pointer :: laplace(:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)
#endif
#ifdef TWODIM
  real(C_DOUBLE), pointer :: laplace(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:)
#endif
#ifdef ONEDIM
  real(C_DOUBLE), pointer :: laplace(:)
  complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
#endif

#ifdef THREEDIM
  real(dl), parameter :: c3=1.0_dl, c2=3.0_dl, c1=14.0_dl, c0=-128.0_dl, cc=30.0_dl
#endif
#ifdef TWODIM
!  real(dl), parameter :: c2=0._dl, c1=1._dl, c0=-4._dl
  real(dl), parameter :: c2=1._dl/6._dl, c1=2._dl/3._dl, c0=-10._dl/3._dl
#endif

contains

  subroutine symp_o2step(dt, w1, w2)
    real*8 :: dt, w1, w2

    integer :: i

    do i=2,n_Hamiltonian_terms-1
       call Hamiltonian_Split(w1*dt/2._dl, i)
    enddo
    call Hamiltonian_split(w1*dt, n_Hamiltonian_terms)
    do i=n_Hamiltonian_terms-1,2,-1
       call Hamiltonian_Split(w1*dt/2._dl,i)
    enddo
    call Hamiltonian_Split((w1+w2)*dt/2._dl,1)
    return
  end subroutine symp_o2step

  subroutine Hamiltonian_Split(dt, hamind)
    real(dl) :: dt
    integer :: hamind

    select case(hamind)
    case(1)
       call Hamiltonian_fields_kin(dt)
    case(2)
       call Hamiltonian_fields_grad_pot(dt)
    case default
       print*,"Undefined Hamiltonian term"
       stop
    end select
  end subroutine Hamiltonian_Split
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now define the individual Hamiltonian pieces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>@brief Kinetic energy of scalar fields
  subroutine Hamiltonian_fields_kin(dt)
    real(dl) :: dt
    integer :: i,j,k

#ifdef VECTORIZE
    fld(:,IRANGE) = fld(:,IRANGE) + dt*fldp(:,IRANGE)
#endif
#ifdef LOOPEVOLVE
!$OMP PARALLEL DO FIRSTPRIVATE(dt)
    FLOOP
      fld(:,LATIND) = fld(:,LATIND) + dt*fldp(:,LATIND)
    FLOOPEND
!$OMP END PARALLEL DO
#endif
  end subroutine Hamiltonian_fields_kin

  !>@brief Gradient and potential energy of the fields
  subroutine Hamiltonian_fields_grad_pot(dt)
    real(dl) :: dt
    real(dl) :: GE2, PE
    
    integer :: l,i,j,k
    real(dl) :: elap, lap(nfld), m2(nfld)
    real(dl) :: b0,b1,b2,b3

#ifdef SPECTRAL
    do l=1,nfld
       laplace(IRANGE) = fld(l,IRANGE)
#ifdef THREEDIM
       call laplacian_3d(nx, ny, nz, laplace, Fk, dk, planf, planb)
#endif
#ifdef TWODIM
       call laplacian_2d(nx, ny, laplace, Fk, dk, planf, planb)
#endif
#ifdef ONEDIM
       call laplacian_1d(nx, laplace, Fk, dk, planf, planb)
#endif

#ifdef VECTORIZE
       fldp(l,IRANGE) = fldp(l,IRANGE) +  &
            dt * (laplace(IRANGE) - modeldv(fld(1,IRANGE),l)) 
#endif
#ifdef LOOPEVOLVE
!$OMP PARALLEL DO PRIVATE(l)
       FLOOP
         fldp(l,LATIND) = fldp(l,LATIND) + &
              dt * (laplace(LATIND) - modeldv(fld(1,LATIND),l)
       FLOOPEND
!$OMP END PARALLEL DO
#endif
    enddo ! end of loop over fields
#endif

#ifdef DISCRETE
    elap = 1./(dx)**2/cc
! Scale laplacian coefficients to reduce number of multiplies
    b0=c0*elap; b1=c1*elap; b2=c2*elap; b3=c3*elap  ! Don't need this now, since these are actual constants
    call wrap_fields()
!$OMP PARALLEL DO PRIVATE(lap,m2,l) FIRSTPRIVATE(b0,b1,b2,b3,dt)
    FLOOP
      lap(:) = STENCIL(b,LAPLACIAN)
      do l=1,nfld
         m2(l) = modeldv(fld(1,LATIND),l)
      enddo
      fldp(:,LATIND) = fldp(:,LATIND) + dt*(lap(:) - m2(:))
    FLOOPEND
!$OMP END PARALLEL DO
#endif
  end subroutine Hamiltonian_fields_grad_pot

  elemental function potential(f1)
    real(dl) :: potential
    real(dl), intent(in) :: f1

    potential = V0*(1._dl-f1**6)**2
  end function potential

  elemental function modeldv(f1,ind)
    real(dl) :: modeldv
    real(dl), intent(in) :: f1
    integer, intent(in) :: ind

    modeldv = -12._dl*V0*(1._dl-f1**6)*f1**5
  end function modeldv

  function get_scale_factor()
    real(dl) :: get_scale_factor
    get_scale_factor = 1._dl
  end function get_scale_factor

  function get_hubble()
    real(dl) :: get_hubble
    get_hubble = 0._dl
  end function get_hubble

  function kinetic_norm()
    real(dl) :: kinetic_norm
    kinetic_norm = 1._dl
  end function kinetic_norm

! This is generic and can be moved elsewhere
  function grav_energy()
    real(dl) :: grav_energy
    grav_energy = 0._dl
  end function grav_energy

  subroutine calc_metric()
    ysclp = 0._dl
  end subroutine calc_metric
                                                                      
  !>@brief
  !> Wrap the fields to impose periodic boundary conditions
  !>
  !>@todo Modularise this and take the field we want to wrap as an input.  Write for MPI distribution.
    subroutine wrap_fields()
      integer :: i,j,k, l

#ifdef THREEDIM
      do l=1,nfld
         fld(:,0,:,:) = fld(:,nx,:,:); fld(:,nx+1,:,:) = fld(:,1,:,:)
         fld(:,:,0,:) = fld(:,:,ny,:); fld(:,:,ny+1,:) = fld(:,:,1,:)
         fld(:,:,:,0) = fld(:,:,:,nz); fld(:,:,:,nz+1) = fld(:,:,:,1)
      enddo
#endif
#ifdef TWODIM
      do l=1,nfld
         fld(l,0,:) = fld(l,nx,:); fld(l,nx+1,:) = fld(l,1,:)
         fld(l,:,0) = fld(l,:,ny); fld(l,:,ny+1) = fld(l,:,1)
      enddo
#endif
#ifdef ONEDIM
      do l=1,nfld
         fld(l,0) = fld(l,nx); fld(l,nx+1) = fld(l,1)
      enddo
#endif
    end subroutine wrap_fields

end module Hamiltonian
