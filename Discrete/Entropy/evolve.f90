!
! Module to store the desired evolution scheme for the fields, etc.
!

!
! To do: Add the Minkowski evolution and fixed background evolution
!   make sure I have the macros in the right place
!

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

  contains

    subroutine step(nst)
      integer :: nst
      call need_get_energies()
      call symp2(dt, nst)
! Now add in the dissipation step, and noise steps
      call dissipation()
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

!
! In order to allow for nontrivial metrics on fields space, etc, this should really be transferred to some sort of model file
!
! These are the subroutines for conformal time
!

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

#else   
! evolution scheme in cosmic time  (this part is working, testing below)
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
