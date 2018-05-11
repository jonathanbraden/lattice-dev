!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!MODULE : Lattice
!
! Stores information about the lattice dimensions and its
! MPI distribution over processes
!
!>@author
!> Jonathan Braden, University College London
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Lattice
#ifdef USE_MPI
  use mpi
  !use mpi_f08
#endif
  
  implicit none

  integer, dimension(1:3) :: nsize
  integer, dimension(1:3) :: istart, iend, padst, paden
  integer, dimension(1:3) :: fstart, fend
  integer, dimension(1:3) :: isize, fsize, sisize

#ifdef USE_MPI
  integer :: mpirank, mpisize
  integer :: right_z_rank, left_z_rank, right_y_rank, left_y_rank
  integer :: ydivs, zdivs
  integer :: lzdiv, lydiv
#endif

contains

  !>@brief
  !> Compute the properties of the sliced lattice if using MPI.
  !> If not using MPI will simply return
  subroutine slicelat(pad,mom)
    integer, intent(in) :: pad
    logical, intent(in) :: mom  ! switch to pad momentum (e.g. for noncanonical kin terms
#ifdef USE_MPI
    
#endif
  end subroutine slicelat

end module Lattice
