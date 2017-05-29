module init_cond

  use modparam
  use params
  use model
  use latpar
  use SHTOOLS  ! needed for the spherical harmonics (probably can change to use a library

#include "macros.h"

  implicit none

  integer, parameter :: numbubble = 1
  integer, dimension(3,1:numbubble) :: bubpos

contains

  subroutine init_fields(f, fp)
    real(dl), dimension(fields, SIRANGE) :: f,fp

    integer, dimension(1:3) :: bcent
    integer :: i
    logical :: dosup

    call getbubpos()
    do i=1, numbubble
       bcent(:) = bubpos(:,i)
       dosup = .true.
       if (i.eq.1) dosup = .false.
       call initbubble(f, fp, bcent, rinit, 0., dosup)
    enddo

    f(:,:,:,:) = f(:,:,:,:) + phifalse
  end subroutine init_fields

  subroutine getbubpos()
!    bubpos(1,1) = n/2 + n/5
!    bubpos(1,2) = n/2 - n/5
    bubpos(1,:) = n/2
    bubpos(2,:) = n/2
    bubpos(3,:) = n/2
  end subroutine getbubpos

  subroutine initbubble(f, fp, center, radinit, tinit, superpose )
    real(dl), dimension(1:fields, SIRANGE) :: f, fp
    integer, dimension(3) :: center
    real(dl) :: radinit, tinit
    logical :: superpose

    real*8, allocatable, dimension(:,:) :: fluc
    real(dl) :: rval
    integer :: r
    integer :: i,j,k,ii,jj,kk, lat, long
    real(dl) :: tmpf, tmpd
    real(dl) :: ang1, ang2, radfluc
    real(dl) :: random
    real(dl) :: integralcheck

    integer :: nsampsize

!    print*,atan2(0.,1.)

    nsampsize = n
    allocate(fluc(1:nsampsize,1:nsampsize))

    if (mpirank.eq.0) then
       call sample_sphericalharm(fluc, nsampsize)
    endif
    call MPI_Bcast(fluc, nsampsize*nsampsize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
    print*,maxval(fluc(:,:))

!    integralcheck=0.
!    do i=1,nsampsize
!       integralcheck = integralcheck + sum(fluc(i,:))
!    enddo

    do k=istart(3),iend(3); kk=k-center(3)
       do j=istart(2),iend(2); jj = j-center(2)
          do i=istart(1),iend(1); ii=i-center(1)

             rval = dble((ii**2+jj**2+kk**2))*dx**2
             rval = rval**0.5

             ang2= atan2(dble(jj),dble(ii))
             ang1 = atan2(-dble(kk), sqrt(dble(ii**2+jj**2)))  + twopi/4._dl
! Interpolate from the grid computed using spherepack
! To do: check the range returned by atan2 and adjust as needed
!  Also, I should really perform some sort of interpolation here
             lat = floor(ang1*dble(nsampsize)*2._dl/twopi)
             if (lat.le.0) lat=1
             if (lat.gt.nsampsize) print*,"Error, lat = ",lat, ang1, i, j ,k

             long = floor(ang2*dble(nsampsize)/twopi) + 1
             if (long.le.0) long = nsampsize + long

             if (long.gt.nsampsize) print*,"error, long = ",long, ang2, i,j,k

             radfluc = fluc(lat,long)

             if (long.gt.nsampsize .or.lat.gt.nsampsize) print*,"too big",i,j,k
             if (long.le.0 .or. lat.le.0) print*,"too small",i,j,k

             tmpf = initprofile(rval+radfluc, radinit, tinit)
             tmpd = initderiv(rval+radfluc, radinit, tinit)

             if (superpose) then
                f(phi,i,j,k) = f(phi,i,j,k) + tmpf
                fp(phi,i,j,k) = fp(phi,i,j,k) + tmpd
             else
                f(phi,i,j,k) = tmpf
                fp(phi,i,j,k) = tmpd
             endif
          enddo
       enddo
    enddo

    deallocate(fluc)

  end subroutine initbubble

  !
  ! The initial profile of the bubble
  ! Eventually, I probably just want to evaluate the instanton eqn.
  ! numerically and then evolve it
  !
  function initprofile(rad, r0, t0)
    real(dl) :: initprofile
    real(dl) :: r0, t0, rad

    real(dl) :: s, temp

    real(dl) :: phistep, phiave

    ! Here I'll put in the initial bubble profile

    phistep = (phifalse-phitrue)/2._dl
    phiave = (-phifalse+phitrue)/2._dl

    if (rad**2.gt.t0**2) then
       s = (rad**2-t0**2)**0.5 - r0
       ! Kosowsky profile
       initprofile = phistep*tanh(0.5*s) + phiave
    else
       initprofile = phitrue-phifalse ! I add the false vacuum on later
    endif

    ! Giblin's profile

  end function initprofile

  function initderiv(rad, r0, t0)
    real :: initderiv
    real :: r0, t0, rad

    real :: s, temp

    if (rad**2.gt.t0**2) then
       s = (rad**2-t0**2)**0.5
       initderiv = -0.5*t0/s
       s = (rad**2-t0**2)**0.5 - r0
       initderiv = initderiv / cosh(0.5*s)**2
    else
       initderiv = 0.  ! check this later
    endif

  end function initderiv

!
! This subroutine requires spherepack.  I can probably actually do this myself but why bother
!
  subroutine sample_sphericalharm(fluc, nsamp)
    integer :: nsamp
    real*8, dimension(1:nsamp,1:nsamp) :: fluc

    integer :: gridsize, nlat, nlong 
    real*8, dimension(:,:,:), allocatable :: alm

    real*8, dimension(1:2) :: a
    complex*8, parameter :: w = (0.0, twopi)
    complex*8 :: gaussrand

    integer :: maxdeg
    integer :: i,j

    maxdeg = (nsamp-2)/2
    allocate(alm(1:2,1:maxdeg+1,1:maxdeg+1))

    gridsize = nsamp

    alm = 0.

    ! Start by getting the a_lm's (need to determine the spectrum)
    ! Exclude the L=0 and 1 modes from being sampled
    do i=3,maxdeg
       do j=1,i
          call random_number(a)
          gaussrand = sqrt(-2.*log(a(1)))*exp(w*a(2))
          a(1) = real(gaussrand)
          a(2) = aimag(gaussrand)
          alm(:,i,j) = a / dble(maxdeg)
       enddo
    enddo
    
    ! Now transform to a grid (I can specify more options here, specifically the normalization of the spherical harmonics)
    call MakeGridDH(fluc, gridsize, alm, maxdeg)
    deallocate(alm)

  end subroutine sample_sphericalharm

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

end module init_cond
