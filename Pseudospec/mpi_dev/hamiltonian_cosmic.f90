!VPRIME_MACRO
!VPRIME_MACRO_VEC

#define POTENTIAL_MACRO_SUM 0.5*fld(1,IRANGE)**2*(1.+g2*fld(2,IRANGE)**2)
#define POTENTIAL_MACRO_NOSUM 0.5*fld(1,LATIND)**2*(1.+g2*fld(2,LATIND)**2)

module Hamiltonian
#include "macros.h"
  use fftw3
  use params

  implicit none

  real(dl), parameter :: g2=100.**2
  real(dl), parameter :: mpl=2.e5
  integer, parameter :: nfld=2

  real(dl), parameter :: phi0 = 1.0093430384226378929425913902459
  real(dl), parameter :: dphi0 = -0.7137133070120812430962278466136
  real(dl), parameter :: H0 = 0.5046715192113189464712956951230
  real(dl), parameter, dimension(nfld) :: fld0 = (/phi0, 0.0/)
  real(dl), parameter, dimension(nfld) :: dfld0 = (/ dphi0, 0.0/)

  integer, parameter :: n_Hamiltonian_terms = 3

! Evolution variables
  real(dl) :: yscl, ysclp

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
       call Hamiltonian_grav_kinetic(dt)
    case(2)
       call Hamiltonian_fields_kin(dt)
    case(3)
       call Hamiltonian_fields_grad_pot(dt)
    case default
       print*,"Undefined Hamiltonian term"
       stop
    end select
  end subroutine Hamiltonian_Split
  
!
! Now define the individual Hamiltonian pieces
!
  subroutine Hamiltonian_grav_kinetic(dt)
    real(dl) :: dt

    yscl = yscl - dt * ysclp *(3._dl/8._dl)
  end subroutine Hamiltonian_grav_kinetic

  subroutine Hamiltonian_fields_kin(dt)
    real*8 :: dt
    real*8 :: KE2
    integer :: i,j,k
    real(dl) :: yscl_loc

#ifdef VECTORIZE
    fld(:,IRANGE) = fld(:,IRANGE) + dt*fldp(:,IRANGE) / yscl**2
    KE2 = sum(fldp(:,IRANGE)**2)
#endif
#ifdef LOOPEVOLVE
    KE2 = 0.
    yscl_loc = yscl
!$OMP PARALLEL DO FIRSTPRIVATE(dt,yscl_loc) REDUCTION(+:KE2)
    FLOOP
      fld(:,LATIND) = fld(:,LATIND) + dt*fldp(:,LATIND) / yscl_loc**2
      KE2 = KE2 + sum(fldp(:,LATIND)**2)
    FLOOPEND
!$OMP END PARALLEL DO
#endif

#ifdef USEMPI
    call MPI_Allreduce(MPI_IN_PLACE, KE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
#endif

    KE2 = KE2 / nvol
    ysclp = ysclp + dt*KE2/yscl**3
  end subroutine Hamiltonian_fields_kin

  subroutine Hamiltonian_fields_grad_pot(dt)
    real(dl) :: dt
    real(dl) :: GE2, PE
    
    integer :: l,i,j,k
    real(dl) :: elap, lap(nfld), m2(nfld)
    real(dl) :: b0,b1,b2,b3
    real(dl) :: acur

    acur = get_scale_factor()
#ifdef SPECTRAL
    GE2 = 0.
    do l=1,nfld
       laplace(IRANGE) = fld(l,IRANGE)
#ifdef THREEDIM
!       call laplacian_3d(nx, ny, nz, laplace, Fk, dk, planf, planb, z_size, z_start_index)
       call laplacian_spectral(nx,ny,nz,laplace,Fk,dk,planf,planb,k_size,k_start_index)
#endif
#ifdef TWODIM
       call laplacian_2d(nx, ny, laplace, Fk, dk, planf, planb)
#endif
#ifdef ONEDIM
       call laplacian_1d(nx, laplace, Fk, dk, planf, planb)
#endif
       laplace = laplace / acur**2
#ifdef VECTORIZE
       fldp(l,IRANGE) = fldp(l,IRANGE) +  &
            dt * yscl**2*(laplace(IRANGE) -  modeldv(fld(1,IRANGE),fld(2,IRANGE),l) )
       GE2 = GE2 - sum(fld(l,IRANGE)*laplace(IRANGE))
#endif
! Have to fill this in for performance testing still
#ifdef LOOPEVOLVE
!$OMP PARALLEL DO FIRSTPRIVATE(yscl,dt) REDUCTION(+:GE2)
       FLOOP
         fldp(l,LATIND) = fldp(l,LATIND) + &
              dt * yscl**2*(laplace(LATIND) - modeldv(fld(1,LATIND),fld(2,LATIND),l) )
         GE2 = GE2 - fld(l,LATIND)*laplace(LATIND)
       FLOOPEND
!$OMP END PARALLEL DO
#endif
      enddo
#ifdef VECTORIZE
    PE = sum(potential(fld(1,IRANGE),fld(2,IRANGE)))
#endif
#ifdef LOOPEVOLVE
    PE=0.
!$OMP PARALLEL DO REDUCTION(+:PE)
    FLOOP
      PE = PE + potential(fld(1,LATIND),fld(2,LATIND))
    FLOOPEND
!$OMP END PARALLEL DO
#endif

#endif

#ifdef DISCRETE
    elap = 1./(acur*dx)**2/cc
! Scale laplacian coefficients to reduce number of multiplies
    b0=c0*elap; b1=c1*elap; b2=c2*elap; b3=c3*elap
    call wrap_fields()
    GE2=0.
    PE=0.
!$OMP PARALLEL DO PRIVATE(lap,m2,l) FIRSTPRIVATE(b0,b1,b2,b3,yscl,dt) REDUCTION(+:PE,GE2)
    FLOOP
      lap(:) = STENCIL(b,LAPLACIAN)
      do l=1,nfld
         m2(l) = modeldv(fld(1,i,j,k),fld(2,i,j,k),l)
      enddo
      fldp(:,i,j,k) = fldp(:,i,j,k) + yscl**2*dt*(lap(:) - m2(:))
      GE2 = GE2 - sum(fld(:,i,j,k)*lap(:))
      PE = PE + potential(fld(1,i,j,k),fld(2,i,j,k))
    FLOOPEND
!$OMP END PARALLEL DO
#endif

#ifdef USEMPI
    call MPI_Allreduce(MPI_IN_PLACE, GE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
    call MPI_Allreduce(MPI_IN_PLACE, PE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
#endif

    PE = PE / nvol
    GE2 = GE2 / nvol
    ysclp = ysclp - yscl*((1._dl/3._dl)*GE2 + 2._dl*PE)*dt

  end subroutine Hamiltonian_fields_grad_pot

  elemental function potential(f1,f2)
    real(dl) :: potential
    real(dl), intent(in) :: f1,f2

    potential = 0.5*f1**2 + 0.5*g2*f1**2*f2**2
  end function potential

  elemental function modeldv(f1,f2,ind)
    real(dl) :: modeldv
    real(dl), intent(in) :: f1,f2
    integer, intent(in) :: ind

    if (ind==1) then
       modeldv = f1 + g2*f1*f2**2
    elseif (ind==2) then
       modeldv = g2*f1**2*f2
    endif
  end function modeldv

  function get_scale_factor()
    real(dl) :: get_scale_factor
    get_scale_factor = yscl**(2./3.)
  end function get_scale_factor

  function get_hubble()
    real(dl) :: get_hubble
    get_hubble = -ysclp / yscl / 4._dl
  end function get_hubble

  function kinetic_norm()
    real(dl) :: kinetic_norm
    kinetic_norm = 1._dl/yscl**4  !1/a^6
  end function kinetic_norm

! This is generic and can be moved elsewhere
  function grav_energy()
    real(dl) :: grav_energy
    grav_energy = -3.*get_hubble()**2
  end function grav_energy

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

  function scalar_rho()
    real(dl) :: scalar_rho
    real(dl) :: GE, PE, KE, rho
    real(dl) :: elap, lap(1:nfld)
    integer :: i,j,k,l
    real(dl) :: acur

    acur = get_scale_factor()
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
      GE = GE / acur**2

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

  subroutine calc_metric()
    real(dl) :: rho, GE, KE, PE
    integer :: i,j,k, l
    real(dl) :: lap(nfld), elap
    real(dl) :: b0, b1, b2, b3
    real(dl) :: acur

    acur = get_scale_factor()

    KE = sum(fldp(:,IRANGE)**2)
    KE = 0.5*KE*kinetic_norm()/nvol
    PE = sum(potential(fld(1,IRANGE),fld(2,IRANGE)))
    PE = PE / nvol

#ifdef SPECTRAL
    GE = 0._dl
    do l=1,nfld
       laplace(IRANGE) = fld(l,IRANGE)
#ifdef THREEDIM
!       call laplacian_3d(nx,ny,nz,laplace,Fk,dk,planf,planb,z_size,z_start_index)
       call laplacian_spectral(nx,ny,nz,laplace,Fk,dk,planf,planb,k_size,k_start_index)
#endif
#ifdef TWODIM
       call laplacian_2d(nx,ny,laplace,Fk,dk,planf,planb)
#endif
#ifdef ONEDIM
       call laplacian_1d(nx, laplace, Fk, dk, planf, planb)
#endif
       GE = GE - sum(fld(l,IRANGE)*laplace(IRANGE))
    enddo
    GE = 0.5_dl*GE / acur**2 / nvol
#endif

#ifdef DISCRETE
    GE = 0.
    call wrap_fields()
    elap = 0.5/dx**2/cc
    FLOOP
      lap = STENCIL(c,LAPLACIAN)
      GE = GE - sum(fld(:,LATIND)*lap(:))
    FLOOPEND
    GE = elap*GE / nvol / acur**2
#endif

#ifdef USEMPI
    call MPI_Allreduce(MPI_IN_PLACE, GE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
    call MPI_Allreduce(MPI_IN_PLACE, PE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
    call MPI_Allreduce(MPI_IN_PLACE, KE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
#endif

    rho = KE + PE + GE
    ysclp = -acur*sqrt(rho*16._dl/3._dl)
  end subroutine calc_metric

!
! Wrap the field to impose periodic boundary conditions
!                                                                               
    subroutine wrap_fields()
      integer :: i,j,k, l

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
    end subroutine wrap_fields

end module Hamiltonian
