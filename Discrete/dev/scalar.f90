!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!MODULE: Scalar
!
! Stores information about the scalar field and it momentum
!
!>@author
!> Jonathan Braden, University College London
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!ISSUES!!!!!!!!!!
! Require knowledge about the lattice to wrap the fields
!
! Since it currently requires the number of fields, it also requires
! knowing information specifically about the scalars.
! Would be nice to make this more general to allow for easy implementation
! of gauge fields, etc.

module Scalar
#ifdef USE_MPI
  use mpi
  !use mpi_f08
#endif
  use constants, only : dl
  implicit none

  real(dl), allocatable, dimension(:,:,:,:) :: fld, fldp

contains

  subroutine wrap_field_periodic_mpi(f)
    real(dl), dimension(nfld,SIRANGE), intent(inout) :: f
  end subroutine wrap_field_periodic_mpi

  subroutine wrap_field_periodic(f)
    real(dl), dimension(nfld,SIRANGE), intent(inout) :: f
  end subroutine wrap_field_periodic

end module Scalar
