!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!MODULE: Hamiltonian
!  Split evolution operator for canonical scalar fields in conformal time
!
!>@author
!> Jonathan Braden, University College London
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Hamiltonian
#ifdef USE_MPI
  use mpi
  ! use mpi_f08
#endif
  use params, only : dl
  implicit none

  integer, parameter :: n_Ham_terms = 3

contains

  subroutine Hamiltonian_Split(dt, ind)
    real(dl), intent(in) :: dt
    integer, intent(in) :: ind

    select case(hamind)
    case(1)
       call Hamiltonian_metric_kin(dt)
    case(2)
       call Hamiltonian_fields_kin(dt)
    case(3)  ! Called the fewest times by symplectic integrator
       call Hamiltonian_fields_grad_pot(dt)
    case default
#ifdef USE_MPI
       if (mpirank == 0) print*,"Undefined Hamiltonian term"
#else
       print*,"Undefined Hamiltonian term"
#endif
       stop
    end select
  end subroutine Hamiltonian_Split

  subroutine Hamiltonian_fields_kin(dt)
    real(dl), intent(in) :: dt
  end subroutine Hamiltonian_fields_kin

  subroutine Hamiltonian_fields_grad_pot(dt)
    real(dl), intent(in) :: dt
  end subroutine Hamiltonian_fields_grad_pot

  subroutine Hamiltonian_metric_kin(dt)
    real(dl), intent(in) :: dt
  end subroutine Hamiltonian_metric_kin

end module Hamiltonian
