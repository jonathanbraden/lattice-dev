!#define PSPEC
!#define FD
!#define FAST_GW

module GW_data
#ifdef PSPEC
  use fftw3
#endif
!  use constants, only dl
  implicit none
  integer, parameter :: dl = kind(1.d0)
  
  integer :: n
  real(dl), dimension(:,:,:,:,:), allocatable :: hij, pi_hij
#ifdef PSPEC
  type(transformPair3D) :: tForm
#endif
  
  type grav_waves
     real(dl), dimension(:,:,:,:,:), allocatable :: hij, pi_hij
  end type grav_waves
  
contains

  subroutine create_gw(gw,n)
    type(grav_waves), intent(inout) :: gw
    integer, intent(in) :: n
    allocate( gw%hij(2,3,1:n,1:n,1:n), gw%pi_hij(2,3,1:n,1:n,1:n) )
  end subroutine create_gw

  subroutine initialise_gw(n)
    integer, intent(in) :: n
    allocate( hij(2,3,1:n,1:n,1:n),pi_hij(2,3,1:n,1:n,1:n) )
  end subroutine initialise_gw
  
end module GW_data
