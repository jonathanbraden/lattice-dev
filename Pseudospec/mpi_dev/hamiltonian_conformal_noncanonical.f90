!
! Required specifications in this Code are the kinetic metric, gradient metric and their derivatives
!

module Hamiltonian
#include "macros.h"
  use fftw3
  use params

  implicit none
#include "macros_gl.h"

  real(dl), parameter :: g2=2.
  real(dl), parameter :: mpl=1.e7/3.
  real(dl), parameter :: zeta=1.
  real(dl), parameter :: zetafac=zeta*(1.+6.*zeta)
  integer, parameter :: nfld=2

! Parameters for implicit RK solver
  integer, parameter :: nvar=2*nfld+2
  real(dl), dimension(nvar,order) :: force_cur
  integer, parameter :: niter=8

  real(dl), parameter :: phi0 = 2.3393837654714997732962993666073
  real(dl), parameter :: dphi0 = -2.7363582010758065274616992909302
  real(dl), parameter :: chi0 = 3.9e-7
  real(dl), parameter :: H0 = 1.9348974397391251388968698880012
  real(dl), parameter, dimension(nfld) :: fld0 = (/phi0,chi0/)
  real(dl), parameter, dimension(nfld) :: dfld0 = (/dphi0,0./)

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

! Variables to store the energy components as I'm evolving
  real(dl), private :: kecur, gecur, pecur

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

! This can probably be optimized better
! And since I'm already doing a RK step, I could just do everything but the
! gradient at once
  subroutine Hamiltonian_fields_kin(dt)
    real(dl) :: dt
    real(dl) :: KE2
    integer :: i,j,k
    real(dl) :: vars(nvar), fac, ftmp, pia

    KE2 = 0._dl
    pia = 0._dl
!$OMP PARALLEL DO PRIVATE(vars, fac, ftmp) FIRSTPRIVATE(dt, yscl, ysclp) REDUCTION(+:KE2, pia)
    FLOOP
      vars(1:nfld) = fld(:,LATIND)
      vars(nfld+1:2*nfld) = fldp(:,LATIND)
      vars(2*nfld+1) = yscl
      vars(2*nfld+2) = ysclp

      call implicit_rk(vars,dt)

      fld(:,LATIND) = vars(1:nfld)
      fldp(:,LATIND) = vars(nfld+1:2*nfld)
      pia = pia + vars(2*nfld+2)
    FLOOPEND
!$OMP END PARALLEL DO

    ysclp = pia / nvol
  end subroutine Hamiltonian_fields_kin

  subroutine implicit_rk(y, dt)
    real(dl) :: y(nvar), dt
    real(dl) :: giter(nvar,order)  ! not thread-safe if I declare it globally

    integer :: k,i
! Iterate to get an approximation for the derivative
    giter = 0.
    do k=1,niter
       giter = matmul(giter,a)
       do i=1,order
          call evalf(y+giter(:,i)*dt, giter(:,i))
       enddo
    enddo
    y = y + matmul(giter,b)*dt  ! For efficiency, should replace with a BLAS call
  end subroutine implicit_rk
!
! Currently for Starobinski model
!
  subroutine evalf(ycur, yprime)
    real(dl), dimension(1:nvar) :: ycur, yprime
    real(dl) :: acur, fac, ftmp, gp1, gp2, ktmp

    acur = get_scale_factor()
    fac = 1.+zeta*ycur(1)**2
    ftmp = fac / (1.+zetafac*ycur(1)**2)
    gp1 = 4.*zeta*ycur(1)*ftmp - 2.*zetafac*ycur(1)*ftmp**2
    gp2 = 2.*zeta*ycur(1)

    ktmp = fac*(ftmp*ycur(nfld+1)**2 + ycur(2*nfld)**2)

    yprime(1) = fac*ftmp*ycur(3) / acur**2
    yprime(2) = fac*ycur(4) / acur**2
    yprime(nfld+1) = - 0.5*(ycur(3)**2*gp1 + ycur(4)**2*gp2) /acur**2
    yprime(2*nfld) = 0.
    yprime(2*nfld+1) = 0.
    yprime(2*nfld+2) = ktmp / acur**3
  end subroutine evalf

  subroutine Hamiltonian_fields_grad_pot(dt)
    real(dl) :: dt
    real(dl) :: GE2, PE
    
    integer :: l,i,j,k
    real(dl) :: elap, lap(nfld), m2(nfld)
    real(dl) :: b0,b1,b2,b3
    real(dl) :: fac1, fac2

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
       fldp(l,IRANGE) = fldp(l,IRANGE) +  &
            dt * (yscl**2*laplace(IRANGE) -  yscl**4*modeldv(fld(1,IRANGE),fld(2,IRANGE),l) )
       GE2 = GE2 - sum(fld(l,IRANGE)*laplace(IRANGE))
    enddo
    PE = sum(potential(fld(1,IRANGE),fld(2,IRANGE)))
#endif

#ifdef DISCRETE
    elap = 1./dx**2/cc
! Scale laplacian coefficients to reduce number of multiplies
    b0=c0*elap; b1=c1*elap; b2=c2*elap; b3=c3*elap
    call wrap_fields()
    GE2=0.
    PE = 0.
!$OMP PARALLEL DO PRIVATE(fac1,fac2,m2,lap) REDUCTION(+:PE,GE2)
    FLOOP
      fac2 = 1.+zetafac*fld(1,LATIND)**2
      fac1 = 1./(1.+zeta*fld(1,LATIND)**2)
      lap(:) = STENCIL(b,LAPLACIAN)
!      lap = elap*fac1*lap
!      lap(1) = fac2*fac1*lap(1)
      do l=1,nfld
         m2(l) = modeldv(fld(1,LATIND),fld(2,LATIND),l)
      enddo
      fldp(:,i,j,k) = fldp(:,LATIND) + dt*(yscl**2*lap(:) - yscl**4*m2(:))  ! modify this appropriately
      GE2 = GE2 - sum(fld(:,LATIND)*lap(:))  ! modify this appropriately
      PE = PE + potential(fld(1,LATIND),fld(2,LATIND))
    FLOOPEND
!$OMP END PARALLEL DO
#endif

#ifdef USEMPI
    call MPI_Allreduce(MPI_IN_PLACE, GE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
    call MPI_Allreduce(MPI_IN_PLACE, PE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
#endif

    PE = PE / nvol
    GE2 = GE2 / yscl**2 / nvol
    ysclp = ysclp - yscl**3*(GE2 + 4._dl*PE)*dt
  end subroutine Hamiltonian_fields_grad_pot

  elemental function potential(f1,f2)
    real(dl) :: potential
    real(dl), intent(in) :: f1,f2
    real(dl) :: norm

    norm = 1./(1.+zeta*f1**2)**2

    potential = norm*(0.25*f1**4 + 0.5*g2*f1**2*f2**2)
  end function potential

  elemental function modeldv(f1,f2,ind)
    real(dl) :: modeldv
    real(dl), intent(in) :: f1,f2
    integer, intent(in) :: ind

    real(dl) :: norm

    norm = 1./(1.+zeta*f1**2)

    if (ind==1) then
       modeldv = norm**2*(f1**3 + g2*f1*f2**2) &
            - norm**3*zeta*f1*(f1**4 + 2.*g2*f1**2*f2**2)
    elseif (ind==2) then
       modeldv = norm**2*g2*f1**2*f2
    endif
  end function modeldv

  function get_scale_factor()
    real(dl) :: get_scale_factor
    get_scale_factor = yscl
  end function get_scale_factor

  function get_hubble()
    real(dl) :: get_hubble
    get_hubble = -ysclp / 6. / yscl**2
  end function get_hubble

  function kinetic_norm()
    real(dl) :: kinetic_norm
    kinetic_norm = 1._dl/yscl**6
  end function kinetic_norm

  subroutine set_ke(kenew)
    real(dl) :: kenew
    kecur = kenew
  end subroutine set_ke
  function kinetic_energy()
    real(dl) :: kinetic_energy
    kinetic_energy = kecur
  end function kinetic_energy
  subroutine set_ge(genew)
    real(dl) :: genew
    gecur = genew
  end subroutine set_ge
  function grad_energy()
    real(dl) :: grad_energy
    grad_energy = gecur
  end function grad_energy
  subroutine set_pe(penew)
    real(dl) :: penew
    pecur = penew
  end subroutine set_pe

! This is generic and can be moved elsewhere
  function grav_energy()
    real(dl) :: grav_energy
    grav_energy = -3.*get_hubble()**2
  end function grav_energy

  subroutine calc_metric()
    real(dl) :: rho, GE, KE, PE
    integer :: i,j,k, l
    real(dl) :: lap(nfld), elap

    real(dl) :: ftmp, fac

    KE = 0._dl
    FLOOP
      fac = 1.+zeta*fld(1,LATIND)**2
      ftmp = fac / (1.+zetafac*fld(1,LATIND)**2)
      KE = KE + fac*(ftmp*fldp(1,LATIND)**2 + fldp(2,LATIND)**2)
    FLOOPEND
    KE = 0.5*KE*kinetic_norm()/nvol

    PE = sum(potential(fld(1,IRANGE),fld(2,IRANGE)))
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
