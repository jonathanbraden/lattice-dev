!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!MODULE : Evolve
!  Evolves the scalar field equations
!
!>@author
!> Jonathan Braden, University College London
!
!!TODO
!>@todo Add Minkowski evolution equations
!>@todo Factor out the time integrator
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module evolve
  use mpi  
  use params
  use latpar
  use scafld
  use homlat

  implicit none

contains

  subroutine step(nst)
    integer, intent(in) :: nst
    call need_get_energies()
    call symp2(dt,nst)
!  Add dissipation
    call dissipation()
  end subroutine step
