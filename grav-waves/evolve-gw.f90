!#define PSPEC
!#define FD
!#define FAST_GW

module GW_Evolution
!  use constants, only dl
  implicit none
  integer, parameter :: dl = kind(1.d0)
  
contains
  ! TO DO: Add scale factor piece in here
  subroutine gw_kinetic_cosmic(h,pi_h,a,dt,n)
    real(dl), dimension(1:2,1:3,1:n,1:n,1:n), intent(out) :: h
    real(dl), dimension(1:2,1:3,1:n,1:n,1:n), intent(in) :: pi_h
    real(dl), intent(in) :: a,dt
    integer, intent(in) :: n

    h = h + pi_h*dt/a**3
  end subroutine gw_kinetic_cosmic

  subroutine gw_source_cosmic(pi_h,h,phi,a,dt,n)
    real(dl), dimension(1:2,1:3,1:n,1:n,1:n), intent(out) :: pi_h
    real(dl), dimension(1:2,1:3,1:n,1:n,1:n), intent(in) :: h
    real(dl), dimension(1:n,1:n,1:n), intent(in) :: phi
    real(dl), intent(in) :: a,dt
    integer, intent(in) :: n

    ! First get laplacian piece
#ifdef PSPEC
    
#else
    
#endif
    
    ! Then compute the source (might be smarter to just pass in the source?  But then won't vectorise properly)
  end subroutine gw_source_cosmic

  subroutine gw_kinetic_conformal(h,pi_h,a,dt)
    real(dl), dimension(:,:,:,:,:), intent(out) :: h
    real(dl), dimension(:,:,:,:,:), intent(in) :: pi_h
    real(dl), intent(in) :: a,dt

    h = h + pi_h*dt
  end subroutine gw_kinetic_conformal

  subroutine gw_source_conformal(pi_h,h,phi,a,dt)
    real(dl), dimension(:,:,:,:,:), intent(out) :: pi_h
    real(dl), dimension(:,:,:,:,:), intent(in) :: h
    real(dl), dimension(:,:,:) :: phi
    real(dl), intent(in) :: a,dt
  end subroutine gw_source_conformal
  
    ! These are space saving subroutines where we only evolve one diagonal and one off diagonal piece on the grounds that for homogeneous and isotropic evolution, the statistical averages are the same for any choice of them

#ifdef FAST_GW
    subroutine gw_kinetic_fast(h,pi_h,dt)
      real(dl), dimension(2,:,:,:), intent(out) :: h
      real(dl), dimension(2,:,:,:), intent(in) :: pi_h
      real(dl), intent(in) :: dt
      
      h = h + pi_h*dt
    end subroutine gw_kinetic_fast

    subroutine gw_source_fast(pi_h,h,phi,dt)
      real(dl), dimension(2,:,:,:), intent(out) :: pi_h
      real(dl), dimension(2,:,:,:), intent(in) :: h
      real(dl), dimension(:,:,:), intent(in) :: phi
      real(dl), intent(in) :: dt
    end subroutine gw_source_fast
#endif
    
end module GW_Evolution
