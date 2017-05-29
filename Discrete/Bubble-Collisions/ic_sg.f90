!
! Initialize two walls parallel to the ( , ) - plane (decide on this)
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

    integer :: i
    integer :: numwalls
    logical :: supos

    integer :: nx

!    real(dl), dimension(1:nside(1),1:nside(2)) :: posfluc

    nx = nside(3)

!    call initwall(f,fp,7*nside(1)/8,phitrue, phifalse, 2., 0.9, .false.)
    call initwall(f, fp, nx/2-nx/4, phitrue, phifalse, 2.**0.5, 0.05, .false.)
    call initwall(f, fp, nx/2+nx/4, phifalse, phitrue, 2.**0.5, -0.05, .true.)

    f(:,:,:,:) = f(:,:,:,:) - phifalse

  end subroutine init_fields

  subroutine get_wall_pos

  end subroutine get_wall_pos

! Start the breather at it's maximum amplitude, so that it's field momenta are zero
! Alternatively, could just initialize it so that the field is zero and we input the velocity
  subroutine init_breather(f, fp, pos, vparam)
    real(dl), dimension(fields, SIRANGE) :: f, fp
    integer :: pos
    real(dl) :: vparam

    real(dl) :: gamma
    real(dl) :: xcur

    integer :: i,j,k,kk

    gamma = 1./sqrt(1. + vparam**2)

    do k=istart(3),iend(3)
       do j=istart(2),iend(2)
          do i=istart(1),iend(1)
             kk = k-pos
             xcur = kk*dx
             f(phi,i,j,k) = 4.*atan(1./(vparam*cosh(gamma*xcur)))
          enddo
       enddo
    enddo
    fp(phi,:,:,:) = 0.

  end subroutine init_breather

  subroutine init_sgkink(f, fp, pos, phileft, phiright, width, vwall, superpose)
    real(dl), dimension(fields, SIRANGE) :: f, fp
    integer :: pos
    real(dl) :: phileft, phiright
    real(dl) :: width
    logical :: superpose
    real(dl) :: vwall

    real(dl),dimension(1:nside(1),1:nside(2)) :: posflu
    integer :: i,j,k, kk
    real(dl) :: phistep, ftmp, fptmp, tmpexp
    real(dl) :: xpos, gamma

    if (mpirank.eq.0) then
       call wall_fluc(posfluc, 1.)
    endif
    call MPI_Bcast(posfluc, nside(1)*nside(2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)

    gamma = 1./sqrt(1. - vwall**2)

    do k=istart(3), iend(3)
       do j=istart(2), iend(2)
          do i=istart(1), iend(1)
             kk = k-pos
             xpos = dble(kk)*dx + posfluc(i,j)
             tmpexp = exp(gamma*xpos)


             ftmp = 4.*atan(tmpexp)
             fptmp = -4.*tmpexp*gamma*vwall/(1.+tmpexp**2)

             if (superpose) then
                f(phi,i,j,k) = f(phi,i,j,k) + ftmp
                fp(phi,i,j,k) = fp(phi,i,j,k) + fptmp
             else
                f(phi,i,j,k) = ftmp
                fp(phi,i,j,k) = fptmp
             endif
          enddo
       enddo
    enddo

  end subroutine init_sgkink

  !
  ! Initialization subroutine a plane wall in the simulation, using a tanh profile for now
  ! 
  ! To do: I can actually distribute the 2-d position arrays more efficiently
  ! (such as with p3dfft) so that each process only gets the part of the array it needs
  ! e.g.  If I'm slicing the y-z plane and the walls move along x, then each mpiprocess
  ! only needs a small portion of the 2-d array.
  ! This can probably be accomplished easily

  subroutine initwall( f, fp, pos, phileft, phiright, width, vwall, superpose )
    real(dl), dimension(fields, SIRANGE) :: f, fp
    integer :: pos
    real(dl) :: phileft, phiright
    real(dl) :: width
    logical :: superpose
    real(dl) :: vwall

    real(dl), dimension(1:nside(1),1:nside(2)) :: posfluc
    integer :: i,j,k, ii
    real(dl) :: phistep, ftmp, fptmp, phiave
    real(dl) :: xpos, gamma

    phistep = (phiright - phileft ) /2.
    phiave = (phiright + phileft) / 2.

! Add in the call to get the fluctuations in the wall here
    if (mpirank.eq.0) then
       call wall_fluc(posfluc , 1.)
    endif
    call MPI_Bcast(posfluc, nside(1)*nside(2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)

    gamma = 1._dl/sqrt(1-vwall**2)

! Now actually initialize the bubble wall
    do k=istart(3), iend(3)
       do j=istart(2),iend(2)
          do i=istart(1),iend(1)
             ii = k - pos  ! align wall perp to z direction
             xpos = dble(ii)*dx + posfluc(i,j)

             ftmp = phistep*tanh(xpos*gamma/width) + phiave
             fptmp = -(phistep*vwall*gamma/width) / cosh(gamma*xpos/width)**2

             if (superpose) then
                f(phi,i,j,k) = f(phi,i,j,k) + ftmp
                fp(phi,i,j,k) = fp(phi,i,j,k) + fptmp
             else
                f(phi,i,j,k) = ftmp
                fp(phi,i,j,k) = fptmp
             endif
          enddo
       enddo
    enddo

  end subroutine initwall

!
! Add a feature (such as a Gaussian bump) to the wall
!
  subroutine wall_shape()

  end subroutine wall_shape

  !
  ! Add fluctuations to the location of the wall (essentially add fluctuations
  ! to pos given in the above subroutine
  ! for the 1-D wall this is accomplished by just initializing a spectrum of 
  ! fourier amplitudes
  !
  ! Input : amp - the amplitude of the fluctuations (should be fixed to allow a spectrum)
  !
  ! To do : change amp (which is assuming a scale inv. spectrum) to actually be a spectrum
  subroutine wall_fluc(fluc, amp)
!    real(dl), dimension(1: :: spec
    real(dl), dimension(nside(1),nside(2)) :: fluc
    real(dl) :: amp 

    complex, dimension(nside(1)/2+1,nside(2)) :: fluc_fourier
    complex, parameter :: w = (0.0, twopi)
    integer*8 :: plan
    integer :: i,j, l, ii,jj
    real(dl) :: rad
    real(dl) :: norm
    real(dl) :: a(nnside(1)), p(nnside(1))
    real(dl) :: dk_effective

    real(dl) :: flucamp


!    norm = amp/(dk)   ! Normalize to make independent of 2-d box size
!    print*,'norm is ',norm
    norm = lside(1)*amp/twopi/2.**0.5
    print*,"norm = ",norm

    dk_effective = twopi/lside(1)

! Start by creating the spectrum of fluctuations I wish to initialize
!
! Currently, I'm just assuming a 1/sqrt(2k) spectrum
    do j=1,nside(2); if (j<=nnside(2)) then; jj=j-1; else; jj=nside(2)+1-j; endif
       do i=1,nnside(1); ii=i-1
          rad = sqrt(dble(ii**2+jj**2))
          l = floor(rad)
          if (l > 0) then
             fluc_fourier(i,j) = norm/sqrt(2.*rad*dk_effective)*exp(-(2.*rad/dble(nnside(1)))**2)
          else
             fluc_fourier(i,j) = 0.   ! don't add a zero-mode (more generally I could)
          endif
       enddo
    enddo

    do j=1,nside(2)
       call random_number(a)
       call random_number(p)

       fluc_fourier(:,j) = sqrt(-2.0*log(a))*exp(w*p) * fluc_fourier(:,j) / dble(nside(1)*nside(2))
    enddo

    ! Now perform the inverse FT to generate the fluctuations in x
    call dfftw_plan_dft_c2r_2d(plan, nside(1), nside(2), fluc_fourier, fluc, FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(plan, fluc_fourier, fluc)
    call dfftw_destroy_plan(plan)

    flucamp = 0.
    do j=1,nside(2)
       do i=1,nside(1)
          flucamp = flucamp + fluc(i,j)**2
       enddo
    enddo
    print*,"flucamp is ",flucamp
    flucamp = flucamp / dble(nside(1)*nside(2))

    print*,"Total fluctuation amplitude is delx^2 = ",flucamp

  end subroutine wall_fluc

  subroutine sample(f, gamma, m2eff)
    real, dimension(IRANGE) :: f
    real :: gamma
    real :: m2eff

    integer*8 :: plan

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
