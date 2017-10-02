!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! A Module for Performing an SVT decomposition of a 3x3 matrix
!
!>@author
!> Jonathan Braden, University College London
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module SVT_Decomposition

contains

  subroutine scalar_potential(s,hij,n)
    real(dl), dimension(:n,:n,:n), intent(out) :: a
    real(dl), dimension(2,3,:n,:n,:n), intent(in) :: hij
    integer, intent(in) :: n
    
    complex(dl), dimension(:nn/2+1,:n,:n) :: Fk
    
  end subroutine scalar_potential

  subroutine vector_potential(v,hij,n)
    real(dl), dimension(3,:n,:n,:n), intent(out) :: v
    real(dl), dimension(2,3,:n,:n,:n), intent(in) :: hij
    integer, intent(in) :: n
  end subroutine vector_potential

  
  
end module SVT_Decomposition
