!
! To do : add proper initialization in two-dimensions and one-dimension
!

module init_fields

  implicit none

  integer, parameter :: nbubble = 2
  integer, dimension(3,nbubble) :: center

contains

! Make this independent of the number of dimension
  subroutine init_fields(fld, fldp)
    real(C_DOUBLE), dimension(nfld, SIRANGE) :: fld
    real(C_DOUBLE), dimension(nfld, IRANGE) :: fldp

    integer :: m

    do m=1,nbubble
       if (m == 1) then
          call init_bubble(fld, fldp, center(:,m), .false.)
       else
          call init_bubble(center(:,m), .true.)
       endif
    enddo

    call add_bulk_fluctuations(fld,fldp)
  end subroutine init_fields

  subroutine init_bubble(f, fp, cent, , superpose)
    real(C_DOUBLE), dimension(nfld,SIRANGE), intent(inout) :: f
    real(C_DOUBLE), dimension(nfld,IRANGE), intent(inout) :: fp
    integer, dimension(3), intent(in) :: cent
    logical, intent(in) :: superpose

    integer :: i,j,k, ii,jj,kk
    real*8 :: rad

    do k= ; kk =  - cent(3)
       do j= ; jj = j-cent(2)
          do i= ; ii = i-cent(1)
             rad = dx*dble(ii**2+jj**2+kk**2)
          enddo
       enddo
    enddo

  end subroutine init_bubble
 
  subroutine init_bubble_fromfile(f,fp,cent, ,superpose)
    real(C_DOUBLE), dimension(nfld,SIRANGE), intent(inout) :: f
    real(C_DOUBLE), dimension(nfld,IRANGE), intent(inout) :: fp
    integer, dimension(3), intent(in) :: cent
    logical, intent(in) :: superpose

    integer :: i,j,k,ii,jj,kk
    real*8 :: tmpf, tmpd
    integer, int_rval2

    do k= ; kk= - cent(3)
       do j= ; jj=j-cent(2)
          do i= ;; ii=i-cent(1)
             int_rval2 = ii**2+jj**2+kk**2
             int_rval2 = (1024/n)**2*int_rval2+1 ! fix this, nsample?

             tmpf = profile(int_rval2)
             tmpd = 0.

             if (superpose) then
                f(1,LATIND) = f(1,LATIND) + tmp
                fp(1,LATIND) = fp(1,LATIND) + tmpd
             else
                f(1,LATIND) = tmp
                fp(1,LATIND) = tmpd
             endif
          enddo
       enddo
    enddo
  end subroutine init_bubble_fromfile

  subroutine add_bulk_fluctuations(f,fp)
    real(C_DOUBLE), dimension(nfld,SIRANGE), intent(inout) :: f
    real(C_DOUBLE), dimension(nfld,IRANGE), intent(inout) :: fp

    integer :: m
    
    do m=1,nfld
       call sample(0.25,1.)
       f(m,IRANGE) = f(m,IRANGE) + laplace(IRANGE)
       call sample(-0.25,1.)
       fp(m,IRANGE) = fp(m,IRANGE) + laplace(IRANGE)
    enddo
  end subroutine add_bulk_fluctuations

end module init_fields
