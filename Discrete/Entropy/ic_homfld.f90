module init_cond

  use modparam
  use params
  use model
  use latpar

#include "macros.h"

  implicit none

  integer, parameter :: numbubble = 1
  integer, dimension(3,1:numbubble) :: bubpos

contains

  subroutine init_fields(f, fp)
    real(dl), dimension(fields, SIRANGE) :: f,fp
    real(dl), dimension(fields) :: m2eff

    integer :: i

    m2eff = modelm2eff()
    m2eff = m2eff - 2.25*H0**2

    do i=1,fields
       call sample(tmp, -0.25, m2eff(i))
       f(i,IRANGE) = tmp + fld0(i)
    enddo

    do i=1,fields
       call sample(tmp, 0.25, m2eff(i))
       fp(i,IRANGE) = tmp + dfld0(i)
    enddo

!    call wrap_field(f)

  end subroutine init_fields

  !
! Randomly sample a gaussian random field with the appropriate spectrum
!
  subroutine sample(f, gamma, m2eff)
    real, dimension(IRANGE) :: f
    real :: gamma
    real :: m2eff

    integer*8 :: plan

#ifndef KCUT_FAC
#define KCUT_FAC 2.0 !2.0
#endif

    integer, parameter :: os = 16, nos = n * os**2
    real, parameter :: dxos = dx/os, dkos = dk/(2*os), kcut = KCUT_FAC*nn*dk/2.0
    complex, parameter :: w = (0.0, twopi)

    real ker(nos), a(nn), p(nn)
    integer i, j, k, l; real kk, norm

    ! calculate (oversampled) radial profile of convolution kernel
    do k = 1,nos; kk = (k-0.5)*dkos
       ker(k) = kk*(kk**2 + m2eff)**gamma * exp(-(kk/kcut)**2)
    end do

    call dfftw_plan_r2r_1d(plan,nos,ker,ker,FFTW_RODFT10,FFTW_ESTIMATE)
    call dfftw_execute(plan); call dfftw_destroy_plan(plan)

    norm = 0.5/(n**3 * sqrt(twopi*dk**3) * mpl) * (dkos/dxos)
    do k = 1,nos; ker(k) = norm * ker(k)/k; end do

       ! initialize 3D convolution kernel (using linear interpolation of radial profile)
       !$omp parallel do
       do k = istart(3),iend(3); do j = istart(2),iend(2); do i = istart(1),iend(1)
          kk = sqrt(real(i-nn)**2 + real(j-nn)**2 + real(k-nn)**2) * os; l = floor(kk)

          if (l > 0) then
             f(i,j,k) = ker(l) + (kk-l)*(ker(l+1)-ker(l))
          else
             f(i,j,k) = (4.0*ker(1)-ker(2))/3.0
          end if
       end do; end do; end do

       ! convolve kernel with delta-correlated Gaussian noise
       call p3dfft_ftran_r2c(f, Fk)

       !$omp parallel do
       do k = fstart(3),fend(3); do j = fstart(2),fend(2)
          call random_number(a); call random_number(p)
          Fk(:,j,k) = sqrt(-2.0*log(a)) * exp(w*p) * Fk(:,j,k)
       end do; end do

       call p3dfft_btran_c2r(Fk, f)  

     end subroutine sample

end module init_cond
