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

!    call get_wall_pos

!    do i=1,numwalls
!       supos = .true.
!       if (i.eq.1) supos = .false.

!       xloc = 
!       call initwall(f, fp, n/2-n/4, , , supos)
!    enddo

    call initwall(f, fp, n/2-n/4, phitrue, phifalse, 2., 0.25, .false.)
    call initwall(f, fp, n/2+n/4, phifalse, phitrue, 2., 0.25, .true.)

    f(:,:,:,:) = f(:,:,:,:) - phifalse

  end subroutine init_fields

  subroutine get_wall_pos

  end subroutine get_wall_pos

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

    real(dl), dimension(1:n,1:n) :: posfluc
    integer :: i,j,k, ii
    real(dl) :: phistep, ftmp, fptmp, phiave
    real(dl) :: xpos, gamma

    phistep = (phiright - phileft ) /2.
    phiave = (phiright + phileft) / 2.

! Add in the call to get the fluctuations in the wall here
    if (mpirank.eq.0) then
       call wall_fluc(posfluc , 500.)
    endif
    call MPI_Bcast(posfluc, n*n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)

    gamma = 1._dl/sqrt(1-vwall**2)

! Now actually initialize the bubble wall
    do k=istart(3), iend(3)
       do j=istart(2),iend(2)
          do i=istart(1),iend(1)
             ii = i - pos
             xpos = dble(ii)*dx + posfluc(j,k)

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
    real(dl), dimension(n,n) :: fluc
    real(dl) :: amp 

    complex, dimension(n/2+1,n) :: fluc_fourier
    complex, parameter :: w = (0.0, twopi)
    integer*8 :: plan
    integer :: i,j, l
    real(dl) :: rad
    real(dl) :: norm
    real(dl) :: a(nn), p(nn)

    norm = amp/(dk**0.5) 

! Start by creating the spectrum of fluctuations I wish to initialize
!
! Currently, I'm just assuming a 1/sqrt(2k) spectrum
    do j=1,n
       do i=1,nn
          rad = sqrt(dble((i-nn)**2+(j-nn)**2))
          l = floor(rad)
          if (l > 0) then
!             fluc_fourier(i,j) = spec(l) + (rad-dble(l))*(spec(l+1)-spec(l))
             fluc_fourier(i,j) = amp/sqrt(2.*rad*dk)
          else
             fluc_fourier(i,j) = 0.   ! don't add a zero-mode (more generally I could)
          endif
       enddo
    enddo

    do j=1,n
       call random_number(a)
       call random_number(p)

       fluc_fourier(:,j) = sqrt(-2.0*log(a))*exp(w*p) * fluc_fourier(:,j) / dble(n**2)
    enddo

    ! Now perform the inverse FT to generate the fluctuations in x
    call dfftw_plan_dft_c2r_2d(plan, n, n, fluc_fourier, fluc, FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(plan, fluc_fourier, fluc)
    call dfftw_destroy_plan(plan)

  end subroutine wall_fluc

  subroutine wall_spec

  end subroutine wall_spec

end module init_cond
