!
! Required specifications in this Code are the kinetic metric, gradient metric and their derivatives
!
module Hamiltonian
#include "macros.h"
  use, intrinsic :: iso_c_binding
  use fftw3
  use params

  implicit none
#include "macros_gl.h"

!////////
! Some model parameters and useful constants
!////////
  real(dl), parameter :: g2=0.
  real(dl), parameter :: mpl=1.e7/3.
  real(dl), parameter :: zeta=1.
  real(dl), parameter :: zetafac=zeta*(1.+6.*zeta)
  real(dl), parameter :: root23 = (2./3.)**0.5
  integer, parameter :: nfld=2

! Parameters for implicit RK solver
  integer, parameter :: nvar=2*nfld+2
  real(dl), dimension(nvar,order) :: force_cur  ! why do I have this here?
  integer, parameter :: niter=8

  real(dl), parameter :: phi0 = 2.3393837654714997732962993666073
  real(dl), parameter :: dphi0 = -2.7363582010758065274616992909302
  real(dl), parameter :: chi0 = 0.
  real(dl), parameter :: H0 = 1.9348974397391251388968698880012
  real(dl), parameter :: phiend = phi0
!  real(dl), parameter :: phiend = ( 0.5*(sqrt((1.+24.*zeta)*(1.+8.*zeta))-1.)/(1.+6.*zeta)/abs(zeta) )**0.5
  real(dl), parameter, dimension(nfld) :: fld0 = (/phiend,0./)
  real(dl), parameter, dimension(nfld) :: dfld0 = (/dphi0,0./)

  integer, parameter :: n_Hamiltonian_terms = 3

  real(dl) :: yscl, ysclp

  type(C_PTR) :: planf, planb
#ifdef THREEDIM
  real(C_DOUBLE), pointer :: laplace(:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)
  real(C_DOUBLE), allocatable, dimension(:,:,:,:) :: fld, fldp
  real(C_DOUBLE), allocatable, dimension(:,:,:) :: tmp_lat
#endif
#ifdef TWODIM
  real(C_DOUBLE), pointer :: laplace(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:)
  real(C_DOUBLE), allocatable, dimension(:,:,:) :: fld, fldp
#endif
#ifdef ONEDIM
  real(C_DOUBLE), pointer :: laplace(:)
  complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
  real(C_DOUBLE), allocatable, dimension(:,:) :: fld, fldp
#endif

#ifdef THREEDIM
  real(dl), parameter :: c3=1.0_dl, c2=3.0_dl, c1=14.0_dl, c0=-128.0_dl, cc=30.0_dl
#endif
#ifdef TWODIM
!  real(dl), parameter :: c2=0._dl, c1=1._dl, c0=-4._dl
  real(dl), parameter :: c2=1._dl/6._dl, c1=2._dl/3._dl, c0=-10._dl/3._dl
#endif
#ifdef ONEDIM
  Error, fill in laplacian coefficients for one-dimensional problem
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

! Version to do potential and kinetic evolution together
  subroutine Hamiltonian_fields_kin(dt)
    real(dl) :: dt
    integer :: i,j,k
    real(dl) :: vars(nvar), pia, avar, davar

    pia = 0._dl
    avar = yscl
    davar = ysclp
!$OMP PARALLEL DO PRIVATE(vars) FIRSTPRIVATE(dt, avar, davar) REDUCTION(+:pia)
    FLOOP
      vars(1:nfld) = fld(:,LATIND)
      vars(nfld+1:2*nfld) = fldp(:,LATIND)
      vars(2*nfld+1) = avar
      vars(2*nfld+2) = davar

      call implicit_rk(vars,dt)

      fld(:,LATIND) = vars(1:nfld)
      fldp(:,LATIND) = vars(nfld+1:2*nfld)
      pia = pia + vars(2*nfld+2)
    FLOOPEND
!$OMP END PARALLEL DO
#ifdef USEMPI
    call MPI_allreduce(MPI_IN_PLACE, pia, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
#endif
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
    real(dl) :: acur, kfac, ktmp, gp

    acur = get_scale_factor()
    kfac = 1.+zeta*ycur(1)**2
    gp = 2.*zeta*ycur(1)
    ktmp = ycur(nfld+1)**2 + kfac*ycur(nfld+2)**2

    yprime(1) = ycur(3) / acur**2  ! Treat inflaton as canonical
    yprime(2) = kfac*ycur(4) / acur**2
    yprime(nfld+1) = - 0.5*( ycur(4)**2*gp ) /acur**2 - acur**4*modeldv(ycur(1),ycur(2),1)
    yprime(nfld+2) = -acur**4*modeldv(ycur(1),ycur(2),2)
    yprime(2*nfld+1) = 0.
    yprime(2*nfld+2) = ktmp / acur**3 - 4.*ycur(2*nfld+1)**3*potential(ycur(1),ycur(2))
  end subroutine evalf

  subroutine Hamiltonian_fields_grad_pot(dt)
    real(dl) :: dt
    real(dl) :: GE2, PE
    
    integer :: l,i,j,k
    real(dl) :: elap, lap(nfld), m2(nfld)
    real(dl) :: b0,b1,b2,b3, d0,d1,d2,d3
    real(dl) :: fac1, fac2, gp1, gp2, lfac(nfld), grad2(nfld)

#ifdef SPECTRAL
    Warning, implementation of noncanonical hamiltonian not yet implemented with spectral discretization.  This will abort compilation.
#endif

#ifdef DISCRETE
    elap = 1./dx**2/cc
! Scale laplacian coefficients to reduce number of multiplies
    b0=0.5*c0*elap; b1=0.5*c1*elap; b2=0.5*c2*elap; b3=0.5*c3*elap  ! the 0.5 are in here to do the averaging int the GDOTL stencil
    d0=0.25*c0*elap; d1=0.25*c1*elap; d2=0.25*c2*elap; d3=0.25*c3*elap
    call wrap_fields()

    tmp_lat(SIRANGE) = 1./(1.+zeta*fld(1,SIRANGE)**2)
!$OMP PARALLEL DO FIRSTPRIVATE(b0,b1,b2,b3,yscl,dt)
    FLOOP
      fldp(2,LATIND) = fldp(2,LATIND) + yscl**2*dt*( STENCIL(b,GDOTL2) )
    FLOOPEND
!$OMP END PARALLEL DO
    tmp_lat(SIRANGE) = 1.
!$OMP PARALLEL DO FIRSTPRIVATE(b0,b1,b2,b3,yscl,dt)
    FLOOP
      fldp(1,LATIND) = fldp(1,LATIND) + yscl**2*dt*( STENCIL(b,GDOTL1) )
    FLOOPEND
!$OMP END PARALLEL DO

    GE2 = 0.
!!$OMP PARALLEL DO PRIVATE(gp1,gp2,lfac,lap,grad2) FIRSTPRIVATE(b0,b1,b2,b3,d0,d1,d2,d3,yscl,dt) REDUCTION(+:GE2)
    FLOOP
      lfac(2) = 1./(1.+zeta*fld(1,LATIND)**2)
      lfac(1) = 1.
      gp1 = 1. 
      gp2 = -2.*zeta*fld(1,LATIND)*lfac(2)**2

      grad2(1:nfld) = STENCIL(d,GRAD2)
      fldp(1,i,j,k) = fldp(1,LATIND) - yscl**2*dt*(gp2*grad2(2)) 
      GE2 = GE2 + sum(lfac(1:nfld)*grad2(1:nfld))  ! modify this appropriately
    FLOOPEND
!!$OMP END PARALLEL DO
#endif

#ifdef USEMPI
    call MPI_Allreduce(MPI_IN_PLACE, GE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
#endif

    GE2 = 2. * GE2 / yscl**2 / nvol
    ysclp = ysclp - dt*yscl**3*GE2
  end subroutine Hamiltonian_fields_grad_pot

  elemental function potential(f1,f2)
    real(dl) :: potential
    real(dl), intent(in) :: f1,f2

    potential = 0.25*f1**4
  end function potential

  elemental function modeldv(f1,f2,ind)
    real(dl) :: modeldv
    real(dl), intent(in) :: f1,f2
    integer, intent(in) :: ind

    if (ind==1) then
       modeldv = f1**3
    elseif (ind==2) then
       modeldv = 0.
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
  function potential_energy()
    real(dl) :: potential_energy
    potential_energy = pecur
  end function potential_energy

! This is generic and can be moved elsewhere
  function grav_energy()
    real(dl) :: grav_energy
    grav_energy = -3.*get_hubble()**2
  end function grav_energy

  function scalar_rho()
    real(dl) :: scalar_rho
    real(dl) :: GE, PE, KE, rho
    real(dl) :: egrad, fac(1:nfld), grad2(1:nfld)
    integer :: i,j,k
    real(dl) :: r23, z1, zfac

    r23 = root23
    z1 = zeta
    zfac = zetafac

#ifdef SPECTRAL
    Warning, fill in scalar energy density for Spectral evolver
#endif

#ifdef DISCRETE
    call wrap_fields()
    egrad = 0.25/dx**2/cc
    GE = 0._dl
!!$OMP PARALLEL DO PRIVATE(grad2,fac) FIRSTPRIVATE(z1,zfac,c0,c1,c2,c3) REDUCTION(+:GE)
    FLOOP
      fac(1) = 1.
      fac(2) = 1./(1.+z1*fld(1,LATIND)**2)
      grad2(1:nfld) = STENCIL(c,GRAD2)
      GE = GE + sum(fac(1:nfld)*grad2(1:nfld))
    FLOOPEND
!!$OMP END PARALLEL DO
#endif

#ifdef VECTORIZE
    KE = 
    PE = sum(potential(fld(1,IRANGE),fld(2,IRANGE)))
#endif
#ifdef LOOPEVOLVE
    KE = 0.; PE = 0.
!$OMP PARALLEL DO PRIVATE(fac) FIRSTPRIVATE(z1,zfac) REDUCTION(+:PE,KE)
    FLOOP
      fac(1) = 1.+z1*fld(1,LATIND)**2
      KE = KE + fldp(1,LATIND)**2 + fac(1)*fldp(2,LATIND)**2
      PE = PE + potential(fld(1,LATIND),fld(2,LATIND))
    FLOOPEND
!$OMP END PARALLEL DO
#endif

    KE = 0.5*KE * kinetic_norm() / nvol  ! check the normalization here
    PE = PE / nvol
    GE = GE * egrad / nvol / yscl**2

#ifdef USEMPI
    call MPI_Allreduce(MPI_IN_PLACE,GE,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)
    call MPI_Allreduce(MPI_IN_PLACE,PE,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)
    call MPI_Allreduce(MPI_IN_PLACE,KE,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)
#endif

    call set_ke(KE)
    call set_pe(PE)
    call set_ge(GE)
    rho = KE + PE + GE
    rhocur = rho
    scalar_rho = rho
  end function scalar_rho

  subroutine calc_metric()
    real(dl) :: rho
    rho = scalar_rho()
    ysclp = -sqrt(yscl**4*12._dl*rho)
  end subroutine calc_metric

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
