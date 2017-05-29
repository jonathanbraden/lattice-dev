#define VPRIME_MACRO fld(l,IRANGE)**3

#define VPRIME_MACRO_VEC fld(:,i,j,k)**3

#define POTENTIAL_MACRO_SUM 0.25*fld(1,IRANGE)**4
#define POTENTIAL_MACRO_NOSUM 0.25*fld(1,LATIND)**4

module Hamiltonian
#include "macros.h"
  use mpi
  use fftw3
  use params

  implicit none

  real(dl), parameter :: g2=2.
  real(dl), parameter :: mpl=1.e7/3.
  integer, parameter :: nfld=1

  real(dl), parameter :: phi0 = 2.3393837654714997732962993666073
  real(dl), parameter :: dphi0 = -2.7363582010758065274616992909302
  real(dl), parameter, dimension(nfld) :: fld0 = (/phi0/)
  real(dl), parameter, dimension(nfld) :: dfld0 = (/dphi0/)

  integer, parameter :: n_Hamiltonian_terms = 2

! These are only here for compatibility with the evolve_spectral.f90 file
! In the future I should be able to eliminate it
  real(dl):: yscl, ysclp

  type(C_PTR) :: planf, planb
#ifdef THREEDIM
  real(C_DOUBLE), allocatable, dimension(:,:,:,:) :: fld, fldp
  real(C_DOUBLE), pointer :: laplace(:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)
#endif
#ifdef TWODIM
  real(C_DOUBLE), allocatable, dimension(:,:,:) :: fld, fldp
  real(C_DOUBLE), pointer :: laplace(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:)
#endif
#ifdef ONEDIM
  real(C_DOUBLE), allocatable, dimension(:,:) :: fld, fldp
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

  real(dl), private :: kecur, gecur, pecur, rhocur

contains

  subroutine symp_o2step(dt, w1, w2)
    real(dl) :: dt, w1, w2

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
  
!
! Now define the individual Hamiltonian pieces
!
  subroutine Hamiltonian_fields_kin(dt)
    real(dl) :: dt
    real(dl) :: KE2
    integer :: i,j,k

#ifdef VECTORIZE
    fld(:,IRANGE) = fld(:,IRANGE) + dt*fldp(:,IRANGE)
#endif
#ifdef LOOPEVOLVE
    FLDLOOP
      fld(:,i,j,k) = fld(:,i,j,k) + dt*fldp(:,i,j,k)
    FLDLOOPEND
#endif
  end subroutine Hamiltonian_fields_kin

  subroutine Hamiltonian_fields_grad_pot(dt)
    real(dl) :: dt

    integer :: l,i,j,k
    real(dl) :: elap, lap(nfld)

#ifdef SPECTRAL
!#ifdef VECTORIZE
    do l=1,nfld
       laplace(IRANGE) = fld(l,IRANGE)
#ifdef THREEDIM
!       call laplacian_3d(nx, ny, nz, laplace, Fk, dk, planf, planb)
       call laplacian_spectral(nx,ny,nz,laplace,Fk,dk,planf,planb,k_size,k_start_index)
#endif
#ifdef TWODIM
       call laplacian_2d(nx, ny, laplace, Fk, dk, planf, planb)
#endif
#ifdef ONEDIM
       call laplacian_1d(nx, laplace, Fk, dk, planf, planb)
#endif
       fldp(l,IRANGE) = fldp(l,IRANGE) +   &
            dt * (laplace(IRANGE) -  VPRIME_MACRO ) !modeldv(fld(1,IRANGE),fld(2,IRANGE),l) )
    enddo
!#endif
#endif
! Have to fill this in for performance testing still
#ifdef LOOPEVOLVE
    FLDLOOP
    FLDLOOPEND
#endif

#ifdef DISCRETE
    elap = 1./dx**2/cc
    call wrap_fields()
    FLOOP
      lap(:) = STENCIL(c,LAPLACIAN)
      lap = lap*elap
!      do l=1,nfld
!         m2(l) = modeldv(fld(1,i,j,k),fld(2,i,j,k),l)
!      enddo
      fldp(:,i,j,k) = fldp(:,i,j,k) + dt*(lap - VPRIME_MACRO_VEC)
    FLOOPEND
#endif

  end subroutine Hamiltonian_fields_grad_pot

! For speed, might be better to inline these above
    pure function modeldv(f)
      real(dl), dimension(1:nfld) :: modeldv
      real(dl), intent(in), dimension(1:nfld) :: f

      modeldv(1) = f(1)**3
!      modeldv(1) = f(1)**3 + g2*f(1)*f(2)**2
!      modeldv(2) = g2*f(1)**2*f(2)
    end function modeldv

    pure function potential(f)
      real(dl) :: potential
      real(dl), dimension(1:nfld), intent(in) :: f

      potential = 0.25*f(1)**4
!      potential = 0.25*f(1)**4 + 0.5*g2*f(1)**2*f(2)**2
    end function potential

    function kinetic_energy()
      real(dl) :: kinetic_energy
      kinetic_energy = kecur
    end function kinetic_energy
    function grad_energy()
      real(dl) :: grad_energy
      grad_energy = gecur
    end function grad_energy
    function potential_energy()
      real(dl) :: potential_energy
      potential_energy = pecur
    end function potential_energy
    

    function get_scale_factor()
      real(dl) :: get_scale_factor
      get_scale_factor = 1._dl
    end function get_scale_factor

    function get_hubble()
      real(dl) :: get_hubble

      get_hubble = 0.
    end function get_hubble

    function kinetic_norm()
      real(dl) :: kinetic_norm
      kinetic_norm = 1._dl
    end function kinetic_norm

    function grav_energy()
      real(dl) :: grav_energy
      grav_energy = 0.
    end function grav_energy

    subroutine calc_metric()
      ysclp = 0.
    end subroutine calc_metric

    function scalar_rho()
      real(dl) :: scalar_rho
      real(dl) :: GE, PE, KE, rho
      real(dl) :: elap, lap(1:nfld)
      integer :: i,j,k, l

#ifdef SPECTRAL
      GE = 0._dl
      do l=1,nfld
         laplace(IRANGE) = fld(l,IRANGE)
#ifdef THREEDIM
         call laplacian_spectral(nx,ny,nz,laplace,Fk,dk,planf,planb,k_size,k_start_index)
#endif
#ifdef TWODIM
         call laplacian_2d(nx,ny,laplace,Fk,dk,planf,planb)
#endif
#ifdef ONEDIM
         call laplacian_1d(nx,laplace,Fk,dk,planf,planb)
#endif

#ifdef VECTORIZE
         GE = GE - sum(fld(l,IRANGE)*laplace(IRANGE))
#endif
#ifdef LOOPEVOLVE
!$OMP PARALLEL DO REDUCTION(+:GE)
         FLOOP
           GE = GE - fld(l,LATIND)*laplace(LATIND)
         FLOOPEND
!$OMP END PARALLEL DO
#endif
      enddo
      GE = 0.5_dl*GE / nvol
#endif
#ifdef DISCRETE
      call wrap_fields()
      elap = 0.5 / dx**2 / cc
      GE = 0._dl
!$OMP PARALLEL DO PRIVATE(lap) REDUCTION(+:GE)
      FLOOP
        lap = STENCIL(c,LAPLACIAN)
        GE = GE - sum(fld(:,LATIND)*lap(:))
      FLOOPEND
!$OMP END PARALLEL DO
      GE = elap*GE / nvol
#endif

#ifdef VECTORIZE
      PE = sum(POTENTIAL_MACRO_SUM)
      KE = sum(fldp(:,IRANGE)**2)
#endif
#ifdef LOOPEVOLVE
      KE = 0.
!$OMP PARALLEL DO REDUCTION(+:KE)
      FLOOP
        KE = KE + sum(fldp(:,LATIND)**2)
      FLOOPEND
!$OMP END PARALLEL DO

      PE = 0.
!$OMP PARALLEL DO REDUCTION(+:PE)
      FLOOP
        PE = PE + POTENTIAL_MACRO_NOSUM
      FLOOPEND
!$OMP END PARALLEL DO
#endif
      KE = 0.5_dl*KE * kinetic_norm() / nvol
      PE = PE / nvol

#ifdef USEMPI
      call MPI_Allreduce(MPI_IN_PLACE,GE,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)
      call MPI_Allreduce(MPI_IN_PLACE,PE,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)
      call MPI_Allreduce(MPI_IN_PLACE,KE,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)
#endif

      kecur = KE
      pecur = PE
      gecur = GE
      rhocur = KE + PE + GE

      scalar_rho = KE + GE + PE
    end function scalar_rho

!
! Wrap the field to impose periodic boundary conditions
!
    subroutine wrap_fields()
      integer :: i,j,k, l, slicesize

      integer :: mpierror, mpistatus(MPI_STATUS_SIZE)

#ifdef USEMPI
#ifdef THREEDIM
      fld(1:nfld,0,:,:) = fld(1:nfld,nx,:,:); fld(1:nfld,nx+1,:,:) = fld(1:nfld,1,:,:)
      fld(1:nfld,:,0,:) = fld(1:nfld,:,ny,:); fld(1:nfld,:,ny+1,:) = fld(1:nfld,:,1,:)

      slicesize = nfld*(nx+2)*(ny+2)
      call MPI_sendrecv(fld(1:nfld,0:nx+1,0:ny+1,z_size), slicesize, MPI_DOUBLE_PRECISION, rightz, 0, &
           fld(1:nfld,0:nx+1,0:ny+1,0), slicesize, MPI_DOUBLE_PRECISION, leftz, 0, &
           MPI_COMM_WORLD, mpistatus, mpierror)
      call MPI_sendrecv(fld(1:nfld,0:nx+1,0:ny+1,1), slicesize, MPI_DOUBLE_PRECISION, leftz, 0, &
           fld(1:nfld,0:nx+1,0:ny+1,z_size+1), slicesize, MPI_DOUBLE_PRECISION, rightz, 0, &
           MPI_COMM_WORLD, mpistatus, mpierror)
#endif
#ifdef TWODIM
      Error, fill this in or compiling will fail
#endif
#ifdef ONEDIM
      Error, fill this in or compiling will fail
#endif

#else
#ifdef THREEDIM
      do l=1,nfld
         fld(l,0,:,:) = fld(l,nx,:,:); fld(l,nx+1,:,:) = fld(l,1,:,:)
         fld(l,:,0,:) = fld(l,:,ny,:); fld(l,:,ny+1,:) = fld(l,:,1,:)
         fld(l,:,:,0) = fld(l,:,:,nz); fld(l,:,:,nz+1) = fld(l,:,:,1)
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
#endif
    end subroutine wrap_fields

end module Hamiltonian
