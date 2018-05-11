!
! Module to store the desired evolution scheme for the fields, etc.
!

!
! To do: Add the Minkowski evolution and fixed background evolution
!   make sure I have the macros in the right place
!

! Flag to set damping or not (currently not tested, only added to Minkowski simulation
!#define DAMPING 1
!#define ABSORBX 1
!#define ABSORBY 1
!#define ABSORBZ 1

#define FLDLOOP do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
#define FLDLOOPEND enddo; enddo; enddo 

#define SCL yscl
#define SCLP ysclp

module evolve

  use mpi
  use params
  use latpar
  use scafld
  use homlat
  
!  use define_latpar
!  use define_lat
!  use define_metric

#include "macros.h"
  implicit none

! This stores the current time within each of the split hamiltonian steps.
! This is used when the background evolution is fixed
  real(dl), dimension(1:3) :: tglobal = (/ 0., 0., 0./)

  contains

#ifdef DAMPING
    real(dl) function damping_function(i,j,k)
      integer :: i,j,k

      integer :: ii,jj,kk
      real(dl) :: temp_damp
      real(dl) :: rad2
! Control parameters for the damping (adjust as needed)
      real(dl), parameter :: eta = 0.005
      real(dl), parameter :: rad_damp = 10.*2.5  ! need initial radius in here
      
      rad2 = (nside(1)/2-i)**2+(nside(2)/2-j)**2+(nside(3)/2-k)**2
      if (rad2 .gt. rad_damp**2) then
         damping_function = eta**2 * (rad2**0.5*dx - rad_damp)**2
      else
         damping_function = 0.
      endif

    end function damping_function
#endif

    subroutine step(nst)
      integer :: nst
      call need_get_energies()
      call symp4(dt, nst)
! Now add in the dissipation step, and noise steps
!      call dissipation()
    end subroutine step

!
! Should I pass in the filter I plan to use?
!
    subroutine dissipation()
      real(dl) :: rho_before, rho_after
      real(dl) :: a
      integer :: i

!      if (.not.do$dissipation) return

      a = get_scale_factor()
      rho_before = scalar_en()
      do i=1,fields
!         call convolve(fld(i,IRANGE), dampk)
!         call convolve(fld(i,IRANGE), dampk)
      enddo
      rho_after = scalar_en()

      rho_rad = rho_rad + (rho_before-rho_after)*a**4  ! check this

    end subroutine dissipation
    

    !
    ! This is the basic building block of the symplectic integrator 
    !
    subroutine symp_o2step(dt, w1, w2)
       real(dl) :: dt, w1, w2

       integer :: i
       integer :: n_Hamiltonian_terms=3
 
       do i=2,n_Hamiltonian_terms-1   ! gerneric for a Hamiltonian with 4 pieces
          call Hamiltonian_Split(w1*dt/2._dl,i)
       enddo
       call Hamiltonian_split(w1*dt, n_Hamiltonian_terms)
       do i=n_Hamiltonian_terms-1,2,-1
          call Hamiltonian_Split(w1*dt/2._dl,i)
       enddo
       call Hamiltonian_Split((w1+w2)*dt/2._dl,1)
       return
    end subroutine symp_o2step

    !
    ! Here is a very simple simplectic integrator for the homogeneous background
    ! It needs to be modified if I change my dynamical variables
    subroutine symp2(dt, nsteps)
       real :: dt
       integer :: nsteps

       integer :: j

! Doe the actual integrations.
       call Hamiltonian_Split(dt/2._dl, 1)
       do j=1,nsteps-1
          call symp_o2step(dt,1._dl,1._dl)
       enddo
       call symp_o2step(dt,1._dl,0._dl)

    end subroutine symp2

    subroutine symp4(dt, nsteps)
      real(dl) :: dt
      integer :: nsteps

      integer :: i
      
      real(dl), parameter :: w1 = 1._dl/(2._dl-2._dl**(1._dl/3._dl))
      real(dl), parameter :: w0 = 1._dl - 2._dl*w1

      call Hamiltonian_Split(w1*dt/2._dl,1)
      do i=1,nsteps
         call symp_o2step(dt, w1, w0)
         call symp_o2step(dt, w0, w1)
         if (i.eq.nsteps) then
            call symp_o2step(dt, w1, 0._dl)
         else
            call symp_o2step(dt, w1, w1)
         endif
      enddo

    end subroutine symp4

    subroutine symp6(dt, nsteps)
      real(dl) :: dt
      integer :: nsteps

      integer :: i

      real(dl), parameter :: w1 = -1.17767998417887_dl
      real(dl), parameter :: w2 = 0.235573213359357_dl
      real(dl), parameter :: w3 = 0.784513610477560_dl
      real(dl), parameter :: w0 = 1._dl-2._dl*(w1+w2+w3)

      call Hamiltonian_Split(w3*dt/2._dl, 1)
      do i=1,nsteps
         call symp_o2step(dt, w3, w2)
         call symp_o2step(dt, w2, w1)
         call symp_o2step(dt, w1, w0)
         call symp_o2step(dt, w0, w1)
         call symp_o2step(dt, w1, w2)
         call symp_o2step(dt, w2, w3)
         if (i.eq.nsteps) then
            call symp_o2step(dt, w3, 0._dl)
         else
            call symp_o2step(dt, w3, w3)
         endif
      enddo
    end subroutine symp6

!
! To do: fill in the coefficients (see PyCOOL paper)
!
    subroutine symp8(dt, nsteps)
      real :: dt
      integer :: nsteps

      integer :: j
      real(dl), parameter :: w1 = 0.
      real(dl), parameter :: w2 = 0.
      real(dl), parameter :: w3 = 0.
      real(dl), parameter :: w4 = 0.
      real(dl), parameter :: w5 = 0.
      real(dl), parameter :: w6 = 0.
      real(dl), parameter :: w7 = 0.
      real(dl), parameter :: w0 = 1._dl - 2._dl*(w1+w2+w3+w4+w5+w6+w7)

      call Hamiltonian_Split(w7*dt/2._dl,1)

      do j=1,nsteps
         call symp_o2step(dt, w7, w6)
         call symp_o2step(dt, w6, w5)
         call symp_o2step(dt, w5, w4)
         call symp_o2step(dt, w4, w3)
         call symp_o2step(dt, w3, w2)
         call symp_o2step(dt, w2, w1)
         call symp_o2step(dt, w1, w0)
         call symp_o2step(dt, w0, w1)
         call symp_o2step(dt, w1, w2)
         call symp_o2step(dt, w2, w3)
         call symp_o2step(dt, w3, w4)
         call symp_o2step(dt, w4, w5)
         call symp_o2step(dt, w5, w6)
         call symp_o2step(dt, w6, w7)
         if (j.eq.nsteps) then
            call symp_o2step(dt, w7, 0._dl)
         else
            call symp_o2step(dt, w7, w7)
         endif
      enddo

    end subroutine symp8

    subroutine Hamiltonian_Split(dt, hamind)
       real :: dt
       integer :: hamind

       select case(hamind)
       case(1)
         call Hamiltonian_metric_kin(dt)
       case(2)
         call Hamiltonian_fields_kin(dt)
       case(3)
         call Hamiltonian_fields_grad_pot(dt)  ! this is called the fewest times
       case default
         if (mpirank.eq.0) print*,"Undefined Hamiltonian term"
         stop
       end select  
    end subroutine Hamiltonian_Split 

! The stuff that's in Conformal will only work in Minkowski space!!!!
#ifdef MINK
    subroutine Hamiltonian_fields_kin(dt)
      real(dl) :: dt

      real(dl) :: a
      integer :: i,j,k

      FLDLOOP
        fld(:,i,j,k) = fld(:,i,j,k) + dt*fldp(:,i,j,k)
      FLDLOOPEND
    end subroutine Hamiltonian_fields_kin

    subroutine Hamiltonian_fields_grad_pot(dt)
      real(dl) :: dt

      real(dl), dimension(1:fields) :: lap
      real(dl) :: elap
      real(dl) :: a
      real(dl) :: gamma
      integer :: i,j,k

! This is setting the boundary conditions to be periodic
! What I should really do is set the outer boundaries on the whole cube to be absorbing here
! This actually requires a bit of thought, since the phi update is local, so if I just update the momenta outside the boundary, nothing will change
      call wrap_field(fld)

      elap = 1./(dx)**2/cc

      FLDLOOP
        lap(:) = STENCIL(c,LAPLACE)
        lap(:) = lap(:) * elap

#ifdef DAMPING
        gamma = dt*0.5*damping_function(i,j,k) ! Put this in a model file or something, so it doesn't have to be inlined and it's easy to change the way damping works
        fldp(:,i,j,k) = fldp(:,i,j,k)*(1.-gamma) + dt*(lap - modeldv(fld,i,j,k))
        fldp(:,i,j,k) = fldp(:,i,j,k)/(1.+gamma)   
#endif
#ifndef DAMPING
        fldp(:,i,j,k) = fldp(:,i,j,k) + dt*(lap - modeldv(fld,i,j,k))

#ifdef ABSORBX
! This is the second order one-sided estimate for the first derivative
        if (i.eq.1) then
           fldp(:,i,j,k) = ( -0.5*fld(:,i+2,j,k) + 2.*fld(:,i+1,j,k)-1.5*fld(:,i,j,k) ) / dx
        else if (i.eq.nside(1)) then
           fldp(:,i,j,k) = ( -0.5*fld(:,i-2,j,k) + 2.*fld(:,i-1,j,k)-1.5*fld(:,i,j,k) ) / dx
        endif
#endif
#ifdef ABSORBY
        if (j.eq.1) then
           fldp(:,i,j,k) = ( -0.5*fld(:,i,j+2,k) + 2.*fld(:,i,j+1,k)-1.5*fld(:,i,j,k) ) / dx
        else if (j.eq.nside(2)) then
           fldp(:,i,j,k) = ( -0.5*fld(:,i,j-2,k) + 2.*fld(:,i,j-1,k)-1.5*fld(:,i,j,k) ) / dx
        endif
#endif
#ifdef ABSORBZ
        if (k.eq.1) then
           fldp(:,i,j,k) = ( -0.5*fld(:,i,j,k+2) + 2.*fld(:,i,j,k+1)-1.5*fld(:,i,j,k) ) / dx
        else if (k.eq.nside(3)) then
           fldp(:,i,j,k) = ( -0.5*fld(:,i,j,k-2) + 2.*fld(:,i,j,k-1)-1.5*fld(:,i,j,k) ) / dx
        endif
#endif
#endif

! And here is the equivalent thing with the first-order estimate for the derivative
! check plus and minus signs here
!        if (i.eq.1) then
!           fldp(:,i,j,k) = ( fld(:,i+1,j,k) - fld(:,i,j,k) ) / dx
!        else if (i.eq.nside(1)) then
!           fldp(:,i,j,k) = ( fld(:,i-1,j,k) - fld(:,i,j,k) ) / dx
!        endif
!        if (j.eq.1) then
!           fldp(:,i,j,k) = ( fld(:,i,j+1,k) - fld(:,i,j,k) ) / dx
!        else if (j.eq.nside(2)) then
!           fldp(:,i,j,k) = ( fld(:,i,j-1,k) - fld(:,i,j,k) ) / dx
!        endif
!        if (k.eq.1) then
!           fldp(:,i,j,k) = ( fld(:,i,j,k+1) - fld(:,i,j,k) ) / dx
!        else if (k.eq.nside(3)) then
!           fldp(:,i,j,k) = ( fld(:,i,j,k-1) - fld(:,i,j,k) ) / dx
!        endif
!#endif
      FLDLOOPEND

    end subroutine Hamiltonian_fields_grad_pot

    subroutine Hamiltonian_metric_kin(dt)
      real(dl) :: dt
    end subroutine Hamiltonian_metric_kin

#else
#ifdef FIXEDBG
#ifdef CONFORMAL
    subroutine Hamiltonian_fields_kin(dt)
      real(dl) :: dt

      real(dl) :: afac
      integer :: i,j,k

      afac = integral_aneg2(tglobal(1),dt)

      FLDLOOP
        fld(:,i,j,k) = fld(:,i,j,k) + afac*fldp(:,i,j,k)
      FLDLOOPEND

      tglobal(1) = tglobal(1) + dt
    end subroutine Hamiltonian_fields_kin

    subroutine Hamiltonian_fields_grad_pot(dt)
      real(dl) :: dt

      real(dl), dimension(1:fields) :: lap
      real(dl) :: elap
      real(dl) :: a2, a4
      integer :: i,j,k

      call wrap_field(fld)

      a4 = integral_a4(tglobal(2), dt)
      a2 = integral_a2(tglobal(2), dt)

      elap = a2 / dx**2 / cc 

      FLDLOOP
        lap(:) = STENCIL(c,LAPLACE)
        fldp(:,i,j,k) = fldp(:,i,j,k) + (elap*lap(:) - a4*modeldv(fld,i,j,k))
      FLDLOOPEND

      tglobal(2) = tglobal(2) + dt
    end subroutine Hamiltonian_fields_grad_pot

#else  ! this is for the case of cosmic time
    subroutine Hamiltonian_fields_kin(dt)
      real(dl) :: dt

      real(dl) :: afac
      integer :: i,j,k

      afac = integral_aneg3(tglobal(1),dt)

      FLDLOOP
        fld(:,i,j,k) = fld(:,i,j,k) + afac*fldp(:,i,j,k)
      FLDLOOPEND
      tglobal(1) = tglobal(1) + dt
    end subroutine Hamiltonian_fields_kin

    subroutine Hamiltonian_fields_grad_pot(dt)
      real(dl) :: dt

      real(dl), dimension(1:fields) :: lap
      real(dl) :: elap
      real(dl) :: a1, a3
      integer :: i,j,k

      call wrap_field(fld)

      a1 = integral_a1(tglobal(2), dt)
      a3 = integral_a3(tglobal(2), dt)
      elap = a1/dx**2/cc

      FLDLOOP
        lap(:) = STENCIL(c,LAPLACE)
        fldp(:,i,j,k) = fldp(:,i,j,k) + (elap*lap(:) - a3*modeldv(fld,i,j,k))
      FLDLOOPEND  

      tglobal(2) = tglobal(2) + dt
    end subroutine Hamiltonian_fields_grad_pot

#endif
    ! This gravitational evolution is empty in either cosmic or conformal time
    subroutine Hamiltonian_metric_kin(dt)
      real(dl) :: dt
      tglobal(3) = tglobal(3) + dt
    end subroutine Hamiltonian_metric_kin

!
! In order to allow for nontrivial metrics on fields space, etc, this should really be transferred to some sort of model file
!
! These are the subroutines for conformal time
!
#else  ! This following is the code for a dynamical homogeneous background 
#ifdef CONFORMAL
    subroutine Hamiltonian_fields_kin(dt)
       real :: dt

       integer :: i,j,k
       real(dl) :: KE2

       KE2 = 0._dl

       FLDLOOP 
         fld(:,i,j,k) = fld(:,i,j,k) + dt*fldp(:,i,j,k) / SCL**2
         KE2 = KE2 + sum(fldp(:,i,j,k)**2)
       FLDLOOPEND   

       call MPI_Allreduce(MPI_IN_PLACE, KE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
       KE2 = KE2/dble(n**3)

       SCLP = SCLP + dt * KE2 / SCL**3
    end subroutine Hamiltonian_fields_kin



!
! Fix all the factors of a lying around here
!
! This is the only subroutine where I need to wrap the fields, so do it
!
    subroutine Hamiltonian_fields_grad_pot(dt)
       real(dl) :: dt

       integer :: i,j,k
       real(dl) :: GE2, PE2, V2  ! accumulator for gradient energy
       real(dl), dimension(1:fields) :: lap
       real(dl) :: a, elap

       call wrap_field(fld)

       a = get_scale_factor()
       elap = 1./(a*dx)**2/cc  ! check this normalization

       GE2 = 0._dl
       PE2 = 0._dl

!
! The laplacian is given by D/(a*dx)**2 where D has the various c coefficients
! The Gradient energy

       FLDLOOP
         lap(:) = STENCIL(c,LAPLACE)
         lap(:) = lap * elap
         V2 = modelv(fld,i,j,k)  ! change call in model.f90
         fldp(:,i,j,k) = fldp(:,i,j,k) + dt*SCL**4*(lap - modeldv(fld,i,j,k))
         GE2 = GE2 - sum(lap(:)*fld(:,i,j,k))
         PE2 = PE2 + V2
! Progressively add the values of the Laplacian (or alternatively gradient) here
       FLDLOOPEND

       call MPI_Allreduce(MPI_IN_PLACE, GE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
       call MPI_Allreduce(MPI_IN_PLACE, PE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)

       GE2 = GE2/dble(n**3)
       PE2 = PE2/dble(n**3)

       SCLP = SCLP - SCL**3*(GE2 + 2._dl*PE2)*dt

    end subroutine Hamiltonian_fields_grad_pot



    subroutine Hamiltonian_metric_kin(dt)
       real(dl) :: dt

       SCL = SCL - dt * SCLP / 6._dl
    end subroutine Hamiltonian_metric_kin

#else   ! evolution scheme in cosmic time  (this part is working, testing below)
    subroutine Hamiltonian_fields_kin(dt)
      real(dl) :: dt

      integer :: i,j,k
      real(dl) :: KE2

      KE2 = 0._dl
      FLDLOOP
        fld(:,i,j,k) = fld(:,i,j,k) + dt*fldp(:,i,j,k)/SCL**2
        KE2 = KE2 + sum(fldp(:,i,j,k)**2) 
      FLDLOOPEND

      call MPI_Allreduce(MPI_IN_PLACE, KE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      KE2 = KE2 / dble(n**3)

      SCLP = SCLP + dt * KE2 / SCL**3 
    end subroutine Hamiltonian_fields_kin

! This one is definitely still wrong
    subroutine Hamiltonian_fields_grad_pot(dt)
      real(dl) :: dt

      integer :: i,j,k
      real(dl), dimension(1:fields) :: lap
      real(dl) :: elap, a      
      real(dl) :: GE2, PE2

!
! I'm computing gradients here, so wrap the fields
!
      call wrap_field(fld)

      GE2 = 0.
      PE2 = 0. 
      a = get_scale_factor()
      elap = 1./(a*dx)**2/cc

      FLDLOOP
        lap(:) = STENCIL(c,LAPLACE)
        lap(:) = lap*elap
        V = modelv(fld,i,j,k)
        fldp(:,i,j,k) = fldp(:,i,j,k) + dt*SCL**2*( lap(:) - modeldv(fld,i,j,k))
        GE2 = GE2 - sum(lap(:)*fld(:,i,j,k))
        PE2 = PE2 + V
      FLDLOOPEND

      call MPI_Allreduce(MPI_IN_PLACE, GE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, PE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      GE2 = GE2/dble(n**3)
      PE2 = PE2/dble(n**3)
      SCLP = SCLP - SCL*( PE2 + (1._dl/3._dl)*GE2)*dt
    end subroutine Hamiltonian_fields_grad_pot

    subroutine Hamiltonian_metric_kin(dt)
      real ::dt
! I have a factor of n**3 that Zhiqi doesn't have (which also finds it's way into above thisgs.  I'm sotring a different SCLP
      SCL = SCL - dt*(3._dl/8._dl)*SCLP
    end subroutine Hamiltonian_metric_kin

#endif
#endif
#endif

#ifdef SCALEDFIELDS
!
! These will need to be commented
!

!
! These are the evolutions of rescaled variables (so I'm not comparing large and small numbers)
!
    subroutine Hamiltonian_fields_kin(dt)
      real(dl) :: dt

      integer :: i,j,k
      real(dl) :: KE2

      KE2 = 0._dl
      FLDLOOP
        fld(:,i,j,k) = fld(:,i,j,k) + dt*fldp(:,i,j,k)
        KE2 = KE2 + sum(fldp(:,i,j,k)**2) 
      FLDLOOPEND

      call MPI_Allreduce(MPI_IN_PLACE, KE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      KE2 = KE2 / dble(n**3)

      SCLP = SCLP + dt * KE2 / SCL
    end subroutine Hamiltonian_fields_kin

! This one is definitely still wrong, I need to fix the normalization of V and GE
    subroutine Hamiltonian_fields_grad_pot(dt)
      real(dl) :: dt

      integer :: i,j,k
      real(dl), dimension(1:fields) :: lap
      real(dl) :: elap, a      
      real(dl) :: GE2, PE2

!
! I'm computing gradients here, so wrap the fields
!
      call wrap_field(fld)

      GE2 = 0.
      PE2 = 0. 
      a = get_scale_factor()
      elap = 1./(a*dx)**2/cc

      FLDLOOP
        lap(:) = STENCIL(c,LAPLACE)
        lap(:) = lap*elap 
        V = modelv(fld,i,j,k)
! potential part is normalized wrong
        fldp(:,i,j,k) = fldp(:,i,j,k) + dt*SCL*( lap(:) - modeldv(fld,i,j,k))
        GE2 = GE2 - sum(lap(:)*fld(:,i,j,k))
        PE2 = PE2 + V
      FLDLOOPEND

      call MPI_Allreduce(MPI_IN_PLACE, GE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, PE2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      GE2 = GE2/dble(n**3)
      PE2 = PE2/dble(n**3)
      SCLP = SCLP - ( PE2 + (1._dl/3._dl)*GE2)*dt/SCL
    end subroutine Hamiltonian_fields_grad_pot

    subroutine Hamiltonian_metric_kin(dt)
      real ::dt
      SCL = SCL - dt*(3._dl/8._dl)*SCLP
    end subroutine Hamiltonian_metric_kin
#endif

end module evolve
