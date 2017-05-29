!
! Module to call MPI subroutines that automatically turns them
! into non-MPI calls for the case when not running in mpi mode
!
module mpiutils
  implicit none
  include "mpif.h"

  INTERFACE 
     module procedure Allreduce_dble, Allreduce_int
  end INTERFACE

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do sum reductions of variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine Allreduce_dble(vars)
    real(dl) :: vars

#ifdef USEMPI
    integer :: nvar, mpierror
    
    nvar = size(vars)
    call MPI_Allreduce(MPI_IN_PLACE, vars, nvar, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
#endif
  end subroutine ALLREDUCE

  subroutine Allreduce_int(vars)
    integer :: vars
  end subroutine ALLREDUCE_INT

end module mpiutils
