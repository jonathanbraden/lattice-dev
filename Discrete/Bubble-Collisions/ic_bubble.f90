module init_cond

  use modparam
  use params
  use model
  use latpar

#include "macros.h"

  implicit none

  integer, parameter :: numbubble = 2
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
       call initbubble(f, fp, bcent, rinit, 1.3*rinit, dosup)
    enddo

    print*, phifalse, phitrue
    print*, initprofile(100.,0.,0.)

    f(:,:,:,:) = f(:,:,:,:) + phifalse
  end subroutine init_fields

  subroutine getbubpos()
    bubpos(1,1) = n/2 + n/5
    bubpos(1,2) = n/2 - n/5
    bubpos(2,:) = n/2
    bubpos(3,:) = n/2
  end subroutine getbubpos

  subroutine initbubble(f, fp, center, radinit, tinit, superpose )
    real(dl), dimension(1:fields, SIRANGE) :: f, fp
    integer, dimension(3) :: center
    real(dl) :: radinit, tinit
    logical :: superpose

    real(dl) :: rval
    integer :: r
    integer :: i,j,k,ii,jj,kk
    real(dl) :: tmpf, tmpd

    do k=istart(3),iend(3); kk=k-center(3)
       do j=istart(2),iend(2); jj = j-center(2)
          do i=istart(1),iend(1); ii=i-center(1)

             rval = dble((ii**2+jj**2+kk**2))*dx**2
             rval = rval**0.5

             tmpf = initprofile(rval, radinit, tinit)
             tmpd = initderiv(rval, radinit, tinit)

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

    phistep = (phifalse - phitrue)/2._dl
    phiave = (-phifalse + phitrue)/2._dl

    ! Here I'll put in the initial bubble profile

    if (rad**2.gt.t0**2) then
       s = (rad**2-t0**2)**0.5 - r0
       ! Kosowsky profile
       initprofile = phistep*tanh(0.5*s)+phiave
    else
       initprofile = phitrue - phifalse ! I add the false vacuum on later
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



end module init_cond
