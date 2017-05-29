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

  subroutine generate_fluctuations_powerlaw(gamma, m2eff, fluc)
    real(C_DOUBLE), intent(in) :: gamma, m2eff
    real(C_DOUBLE), dimension(IRANGE), intent(out) :: fluc

    type(C_PTR) :: plan_sin, plan1
    integer, parameter :: os = 16, nos = max(nx,ny,nz)*os**2
    real(C_DOUBLE), parameter :: dxos = dx/os, dkos=dk/(2.*os), kcut = KCUT_FAC*min(nnx,nny,nnz)*dk/2.0
    complex, parameter :: w = (0.0, twopi)

    real(C_DOUBLE) :: ker(nos), a(nnx), p(nnx)
    

  end subroutine generate_fluctuations_powerlaw

  subroutine generate_fluctuations(fluc, spectrum)
    real(C_DOUBLE), dimension(IRANGE), intent(out) :: fluc
    real(C_DOUBLE), dimension(), intent(in) :: spectrum

    
  end subroutine generate_fluctuations

!!!!!!!!!!!!!!
! Use linear perturbation theory to compute initial spectrum of fluctuations
!!!!!!!!!!!!!!
  subroutine generate_linear_spectrum()

  end subroutine generate_linear_spectrum

end module init_fields
