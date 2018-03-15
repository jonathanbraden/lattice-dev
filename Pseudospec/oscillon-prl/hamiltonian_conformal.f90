module Hamiltonian
#include "macros.h"
  use fftw3
  use params

  implicit none
  real(dl), parameter :: V0_mpl = 1.e-19, vev_mpl = 0.01_dl

  integer, parameter :: nfld = 1
  real(dl), parameter :: vev = vev_mpl, vev2 = vev**2, vev_i = 1._dl/vev     ! Planck units
  real(dl), parameter :: V0 = (1._dl/72._dl)*vev2
  real(dl), parameter :: mpl2 = (1._dl/72._dl)*vev2/V0_mpl, mpl = mpl2**0.5
!  real(dl), parameter :: mpl2 = 1.e14/7.2_dl, mpl = mpl2**0.5

  real(dl), parameter :: phi0 = 0.08_dl*vev
  real(dl), parameter :: dphi0 = 2.49e-9*vev**2*mpl  ! This is in cosmic time, not conformal
  real(dl), parameter :: H0 = (V0/3._dl)**0.5  ! Fix this, currently an approximation
  real(dl), parameter, dimension(nfld) :: fld0 = (/phi0/)
  real(dl), parameter, dimension(nfld) :: dfld0 = (/dphi0/)

  integer, parameter :: n_Hamiltonian_terms = 3

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
       call Hamiltonian_grav_kinetic(dt)
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

    yscl = yscl - dt * ysclp / 6._dl
  end subroutine Hamiltonian_grav_kinetic

  subroutine Hamiltonian_fields_kin(dt)
    real*8 :: dt
    real*8 :: KE2
    integer :: i,j,k

#ifdef VECTORIZE
    fld(:,IRANGE) = fld(:,IRANGE) + dt*fldp(:,IRANGE) / yscl**2
    KE2 = sum(fldp(:,IRANGE)**2)
#endif
#ifdef LOOPEVOLVE
    KE2 = 0.
! To do: copy in yscl,dt and use as a local variable (to avoid many reads on it)
!$OMP PARALLEL DO FIRSTPRIVATE(dt,yscl) REDUCTION(+:KE2)
    FLOOP
      fld(:,LATIND) = fld(:,LATIND) + dt*fldp(:,LATIND) / yscl**2
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

#ifdef SPECTRAL
    GE2 = 0.
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
!            dt * (yscl**2*laplace(IRANGE) -  yscl**4*modeldv(fld(1,IRANGE),fld(2,IRANGE),l) ) ! call for a two-field model
            dt * (yscl**2*laplace(IRANGE) - yscl**4*modeldv(fld(1,IRANGE),l)) 
       GE2 = GE2 - sum(fld(l,IRANGE)*laplace(IRANGE))
#endif
#ifdef LOOPEVOLVE
!$OMP PARALLEL DO FIRSTPRIVATE(yscl) REDUCTION(+:GE2)
       FLOOP
         fldp(l,LATIND) = fldp(l,LATIND) + &
!              dt * (yscl**2*laplace(LATIND) - yscl**4*modeldv(fld(1,LATIND),fld(2,LATIND),l) ) ! call for a two-field model
              dt * (yscl**2*laplace(LATIND) - yscl**4*modeldv(fld(1,LATIND),l)
         GE2 = GE2 - fld(l,LATIND)*laplace(LATIND)
       FLOOPEND
!$OMP END PARALLEL DO
#endif
    enddo ! end of loop over fields
#ifdef VECTORIZE
!    PE = sum(potential(fld(1,IRANGE),fld(2,IRANGE))) ! call for two-field model
    PE = sum(potential(fld(1,IRANGE)))
#endif
#ifdef LOOPEVOLVE
    PE = 0.
!$OMP PARALLEL DO REDUCTION(+:PE)
    FLOOP
!      PE = PE + potential(fld(1,LATIND),fld(2,LATIND)) ! call for a two-field model
       PE = PE + potential(fld(1,LATIND))
    FLOOPEND
!$OMP END PARALLEL DO
#endif
    GE2 = GE2 / yscl**2
! end of spectral evolution
#endif

#ifdef DISCRETE
    elap = 1./(yscl*dx)**2/cc
! Scale laplacian coefficients to reduce number of multiplies
    b0=c0*elap; b1=c1*elap; b2=c2*elap; b3=c3*elap
    call wrap_fields()
    GE2=0.
    PE=0.
!$OMP PARALLEL DO PRIVATE(lap,m2,l) FIRSTPRIVATE(b0,b1,b2,b3,yscl,dt) REDUCTION(+:PE,GE2)
    FLOOP
      lap(:) = STENCIL(b,LAPLACIAN)
      do l=1,nfld
!         m2(l) = modeldv(fld(1,LATIND),fld(2,LATIND),l) ! for two field model
         m2(l) = modeldv(fld(1,LATIND),l)
      enddo
      fldp(:,LATIND) = fldp(:,LATIND) + yscl**4*dt*(lap(:) - m2(:))
      GE2 = GE2 - sum(fld(:,LATIND)*lap(:))
!      PE = PE + potential(fld(1,LATIND),fld(2,LATIND)) ! for two field model
      PE = PE + potential(fld(1,LATIND))
    FLOOPEND
!$OMP END PARALLEL DO
#endif

#ifdef USEMPI
    call MPI_Allreduce(MPI_IN_PLACE, GE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
    call MPI_Allreduce(MPI_IN_PLACE, PE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
#endif

    PE = PE / nvol
    GE2 = GE2 / nvol
    ysclp = ysclp - yscl**3*(GE2 + 4._dl*PE)*dt
  end subroutine Hamiltonian_fields_grad_pot

!  elemental function potential(f1,f2)
!    real(dl) :: potential
!    real(dl), intent(in) :: f1,f2

!    potential = 0.25*f1**4 + 0.5*g2*f1**2*f2**2
!  end function potential
  elemental function potential(f1)
    real(dl) :: potential
    real(dl), intent(in) :: f1

    potential = V0*(1._dl-(vev_i*f1)**6)**2
  end function potential

  elemental function modeldv(f1,ind)
    real(dl) :: modeldv
    real(dl), intent(in) :: f1
    integer, intent(in) :: ind

    modeldv = -12._dl*(V0*vev_i)*(1._dl-(vev_i*f1)**6)*(vev_i*f1)**5
  end function modeldv

  function get_scale_factor()
    real(dl) :: get_scale_factor
    get_scale_factor = yscl
  end function get_scale_factor

  function get_hubble()
    real(dl) :: get_hubble
    get_hubble = -ysclp / 6._dl / yscl**2
  end function get_hubble

  function kinetic_norm()
    real(dl) :: kinetic_norm
    kinetic_norm = 1._dl/yscl**6
  end function kinetic_norm

! This is generic and can be moved elsewhere
  function grav_energy()
    real(dl) :: grav_energy
    grav_energy = -3.*get_hubble()**2
  end function grav_energy

  subroutine calc_metric()
    real(dl) :: rho, GE, KE, PE
    integer :: i,j,k, l
    real(dl) :: lap(nfld), elap

    KE = sum(fldp(:,IRANGE)**2)
    KE = 0.5*KE*kinetic_norm()/nvol
!    PE = sum(potential(fld(1,IRANGE),fld(2,IRANGE)))  ! two field call
    PE = sum(potential(fld(1,IRANGE)))
    PE = PE / nvol

#ifdef SPECTRAL
    GE = 0._dl
    do l=1,nfld
       laplace(IRANGE) = fld(l,IRANGE)
#ifdef THREEDIM
       call laplacian_3d(nx,ny,nz,laplace,Fk,dk,planf,planb)
#endif
#ifdef TWODIM
       call laplacian_2d(nx,ny,laplace,Fk,dk,planf,planb)
#endif
#ifdef ONEDIM
       call laplacian_1d(nx, laplace, Fk, dk, planf, planb)
#endif
       GE = GE - sum(fld(l,IRANGE)*laplace(IRANGE))
    enddo
    GE = 0.5_dl*GE / yscl**2 / nvol
#endif

#ifdef DISCRETE
    GE = 0.
    call wrap_fields()
    elap = 0.5/dx**2/cc
    FLOOP
      lap = STENCIL(c,LAPLACIAN)
      GE = GE - sum(fld(:,i,j,k)*lap(:))
    FLOOPEND
    GE = elap*GE / nvol / yscl**2
#endif
    rho = KE + PE + GE
    ysclp = -sqrt(yscl**4*12._dl*rho)
  end subroutine calc_metric

!
! Wrap the field to impose periodic boundary conditions
!                                                                               
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
