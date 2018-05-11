!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!MODULE: Integrator Yoshida
!  Implementation of Yoshida symplectic integration scheme for lattice
!
!>@author
!> Jonathan Braden, University College London
!
!!TODO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Yoshida
#ifdef USE_MPI
  use mpi  ! Do I need this?  Use in a safer way
  !use mpi_f08
#endif
  use params, only : dl
  implicit none
contains

  !>@brief
  !> Implements an O(dt^2) accurate step of the symplectic integrator, 
  !> including operator fusion
  subroutine symp_o2step(dt,w1,w2)
    real(dl), intent(in) :: dt, w1, w2
    integer :: i
    integer :: n_Ham_terms = 3  ! Adjust this for Minkowski

    do i=2,n_Ham_terms-1
       call Hamiltonian_Split(0.5_dl*w1*dt,i)
    enddo
    call Hamiltonian_Split(w1*dt,n_Ham_terms)
    do i=n_Ham_terms,2,-1
       call Hamiltonian_Split(0.5_dl*w1*dt,i)
    enddo
    call Hamiltonian_Split(0.5_dl*(w1+w2)*dt,1)
  end subroutine symp_o2step

  !>@Brief
  !> Implements symplectic integrator logic for arbitrary choices of
  !> stepping parameters w
  !
  !>@todo Write this subroutine
  subroutine symp(dt,ns,w)
    real(dl), intent(in) :: dt
    integer, intent(in) :: ns
    real(dl), dimension(:), intent(in) :: w
  end subroutine symp

  !>@brief
  !> Second order accurate symplectic integration with time step dt for ns steps
  subroutine symp2(dt,ns)
    real(dl), intent(in) :: dt
    integer, intent(in) :: ns

    integer :: i
    call Hamiltonian_Split(0.5_dl*dt,1)
    do i=1,ns-1
       call symp_o2step(dt,1._dl,1._dl)
    enddo
    call symp_o2step(dt,1._dl,0._dl)
  end subroutine symp2

  subroutine symp4(dt,ns)
    real(dl), intent(in) :: dt
    integer, intent(in) :: ns

    ! Parameters controlling the integrator
    real(dl), parameter :: w1 = 1._dl/(2._dl-2._dl**(1._dl/3._dl))
    real(dl), parameter :: w0 = 1._dl - 2._dl*w1
    ! real(dl), parameter, dimension(1:2) :: w = (/0._dl,0._dl/)
    integer :: i

    !!!!!!!!!! Remove the nested if statement in here !!!!!!!!!!!!!!!!!!
    call Hamiltonian_Split(0.5*w1*dt,1)
    do i=1,ns
       call symp_o2step(dt,1._dl,1._dl)
       if (i == ns) then
          call symp_o2step(dt,w1,0._dl)
       else
          call symp_o2step(dt,w1,w1)
    enddo
    ! call symp_o2step(dt,1._dl,1._dl)
    ! call symp_o2step(dt,w1,0._dl)
 end do

 !>@brief
 !> Implementation of 6th order Yoshida integrator
 subroutine symp6(dt,ns)
   real(dl), intent(in) :: dt
   integer, intent(in) :: ns
   
   real(dl), parameter :: w1 = -1.17767998417887_dl
   real(dl), parameter :: w2 = 0.235573213359357_dl
   real(dl), parameter :: w3 = 0.784513610477560_dl
   real(dl), parameter :: w0 = 1._dl-2._dl*(w1+w2+w3)
   integer :: i

   call Hamiltonian_Split(0.5_dl*w3*dt,1)
   do i=1,ns
      call symp_o2step(dt,w3,w2)
      call symp_o2step(dt,w2,w1)
      call symp_o2step(dt,w1,w0)
      call symp_o2step(dt,w0,w1)
      call symp_o2step(dt,w1,w2)
      call symp_o2step(dt,w2,w3)
      if (i == ns) then
         call symp_o2step(dt,w3,0._dl)
      else
         call symp_o2step(dt,w3,w3)
      endif
   enddo
 end subroutine symp6

 subroutine symp8(dt,ns)
   real(dl), intent(in) :: dt
   integer, intent(in) :: ns

   integer :: i
 end subroutine symp8

end module Yoshida
