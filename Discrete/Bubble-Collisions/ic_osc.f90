!
! Initialize an oscillon profile according to the Hindmarsh paper
!

module init_cond

  use mpi       ! to transfer the fluctuations between processors
  use modparam  ! I probably don't need this
  use params
  use model
  use latpar

#include "macros.h"

  implicit none

contains

  subroutine init_fields(f,fp)
    real(dl), dimension(1:fields,SIRANGE) :: f, fp

    integer :: i, j, k, ii, jj, kk
    real(dl) :: amp
    real(dl) :: r0
    
    amp = -1.
    r0 = 2.**0.5*3.0

    call make_oscillon(f, fp, (/ n/2, n/2, n/2 /), r0, amp)
    call make_oblong(f,fp, (/ n/2, n/2, n/2 /), (/ r0, 1.3*r0, 1.7*r0 /), amp)
    
  end subroutine init_fields

  subroutine make_oblong(f, fp, center, rinit, amp)
    real(dl), dimension(1:fields, SIRANGE) :: f, fp 
    integer, dimension(1:3) :: center
    real(dl) :: amp

    integer :: i,j,k, ii,jj,kk
    real(dl), dimension(1:3) :: rinit
    real(dl) :: rad

    do k=istart(3), iend(3); kk=center(3)-k
       do j=istart(2), iend(2); jj=center(2)-j
          do i=istart(1), iend(1); ii = center(1)-i
             rad = dble(ii**2/rinit(1)**2 + jj**2/rinit(2)**2 + kk**2/rinit(3)**2)

             f(phi,i,j,k) = (1.-amp*exp(-rad))
             fp(phi,i,j,k) = 0.
          enddo
       enddo
    enddo

  end subroutine make_oblong

  subroutine make_oscillon(f, fp, center, rinit, amp)
    real(dl), dimension(1:fields, SIRANGE) :: f, fp
    integer, dimension(1:3) :: center
    real(dl) :: rinit, amp

    integer :: i,j,k, ii,jj,kk
    real(dl) :: rad

    do k=istart(3),iend(3); kk=center(3)-k
       do j=istart(2),iend(2); jj=center(2)-j
          do i=istart(1),iend(1); ii=center(1)-i
             rad = dble(ii**2+jj**2+kk**2)*dx**2

             f(phi,i,j,k) = (1-amp*exp(-rad/rinit**2))
             fp(phi,i,j,k) = 0.
          enddo
       enddo
    enddo

  end subroutine make_oscillon

end module init_cond
