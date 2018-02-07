module fftw3
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

! Define some generic interfaces
  interface allocate_fftw_array
     module procedure allocate_1d_array, allocate_2d_array, allocate_3d_array
  end interface allocate_fftw_array

  contains
   
    subroutine allocate_1d_array(L, arr, Fk)
      integer :: L

      real(C_DOUBLE), pointer :: arr(:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)

      type(C_PTR) :: fptr, fkptr
      integer :: LL

      LL = L/2+1

      fptr = fftw_alloc_real(int(L, C_SIZE_T))
      call c_f_pointer(fptr, arr, [L])
      fkptr = fftw_alloc_complex(int(L, C_SIZE_T))
      call c_f_pointer(fkptr, Fk, [LL])

    end subroutine allocate_1d_array
    
    subroutine allocate_2d_array(L,M, arr, Fk)
      integer :: L,M

      real(C_DOUBLE), pointer :: arr(:,:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:)

      type(C_PTR) :: fptr, fkptr
      integer :: LL

      LL = L/2+1
! Check orderings in here      
      fptr = fftw_alloc_real(int(L*M, C_SIZE_T))
      call c_f_pointer(fptr, arr, [L,M])
      fkptr = fftw_alloc_complex(int(LL*M, C_SIZE_T))
      call c_f_pointer(fkptr, Fk, [LL,M])

    end subroutine allocate_2d_array

    subroutine allocate_3d_array(L,M,N, arr, Fk)
      integer :: L,M,N
      
      real(C_DOUBLE), pointer :: arr(:,:,:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)

      type(C_PTR) :: fptr, fkptr
      integer :: LL

      LL = L/2+1

      fptr = fftw_alloc_real(int(L*M*N, C_SIZE_T))
      call c_f_pointer(fptr, arr, [L,M,N])
      fkptr = fftw_alloc_complex(int(LL*M*N, C_SIZE_T))
      call c_f_pointer(fkptr, Fk, [LL,M,N])

    end subroutine allocate_3d_array

    subroutine allocate_3d_real_array(L,M,N,arr)
      integer :: L,M,N
      
      real(C_DOUBLE), pointer :: arr(:,:,:)
      type(C_PTR) :: fptr
      
      fptr = fftw_alloc_real(int(L*M*N, C_SIZE_T))
      call c_f_pointer(fptr,arr,[L,M,N])
    end subroutine allocate_3d_real_array

    subroutine allocate_3d_fourier_array(L,M,N,Fk)
      integer :: L,M,N

      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)
      type(C_PTR) :: fkptr
      integer :: LL

      LL = L/2+1
      fkptr = fftw_alloc_complex(int(LL*M*N, C_SIZE_T))
      call c_f_pointer(fkptr, Fk, [LL,M,N])
    end subroutine allocate_3d_fourier_array

! To add : allocate and initialize inplace/outofplace FFTW plans
! To add : extend this to work with MPI

!!!!!!!!!!!!!!
! Mathematical Subroutines, compute various derivatives
!!!!!!!!!!!!!!
    subroutine derivative_1d(n,f,Fk,dk)
      integer :: n
      real*8 :: dk
      real(C_DOUBLE), pointer :: f(:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
      integer :: i,nn

      type(C_PTR) :: planf, planb

      nn = n/2+1
      planf = fftw_plan_dft_r2c_1d(n, f, Fk, FFTW_ESTIMATE)
      planb = fftw_plan_dft_c2r_1d(n, Fk, f, FFTW_ESTIMATE)
      call fftw_execute_dft_r2c(planf, f, Fk)
      do i=1,nn
         Fk(i) = (i-1)*dk*(0,1.)*Fk(i)
      enddo

      call fftw_execute_dft_c2r(planb, Fk, f)
      call fftw_destroy_plan(planf)
      call fftw_destroy_plan(planb)
      f = f / dble(n)

    end subroutine derivative_1d

    subroutine laplacian_1d(n, f, Fk, dk, planf, planb)
      integer :: n
      real(C_DOUBLE), pointer :: f(:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
      real*8 :: dk

      integer :: i,nn

      type(C_PTR) :: planf, planb

      nn = n/2+1
!      planf = fftw_plan_dft_r2c_1d(n, f, Fk, FFTW_ESTIMATE)
!      planb = fftw_plan_dft_c2r_1d(n, Fk, f, FFTW_ESTIMATE)
      call fftw_execute_dft_r2c(planf, f, Fk)
      do i=1,nn
         Fk(i) = -((i-1)*dk)**2*Fk(i)
      enddo

      call fftw_execute_dft_c2r(planb, Fk, f)
!      call fftw_destroy_plan(planf)
!      call fftw_destroy_plan(planb)
      f = f / dble(n)
    end subroutine laplacian_1d

    subroutine laplacian_2d(n1, n2, f, Fk, dk, planf, planb)
      integer :: n1, n2
      real(C_DOUBLE), pointer :: f(:,:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:)
      real*8 :: dk
      type(C_PTR) :: planf, planb

      real*8 :: rad2
      integer :: nn1, nn2
      integer :: i,j,ii,jj

      nn1 = n1/2+1; nn2=n2/2+1
!      planf = fftw_plan_dft_r2c_2d(n2, n1, f, Fk, FFTW_ESTIMATE)
!      planb = fftw_plan_dft_c2r_2d(n2, n1, Fk, f, FFTW_ESTIMATE)
      call fftw_execute_dft_r2c(planf, f, Fk)
      do j=1,n2; if (j>nn2) then; jj = n2+1-j; else; jj=j-1; endif
         do i=1,nn1
            rad2 = dble((i-1)**2 + jj**2)
            Fk(i,j) = -rad2*dk**2*Fk(i,j)
         enddo
      enddo

      call fftw_execute_dft_c2r(planb, Fk, f)
!      call fftw_destroy_plan(planf)
!      call fftw_destroy_plan(planb)
      f = f / dble(n1) / dble(n2)
    end subroutine laplacian_2d

    subroutine laplacian_3d(n1, n2, n3, f, Fk, dk, planf, planb)
      integer :: n1, n2, n3
      real(C_DOUBLE), pointer :: f(:,:,:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)
      real*8 :: dk
      type(C_PTR) :: planf, planb

      real*8 :: rad2
      integer :: nn1, nn2, nn3
      integer :: i,j,k,ii,jj,kk

      nn1 = n1/2+1; nn2=n2/2+1; nn3=n3/2+1
!      planf = fftw_plan_dft_r2c_3d(n3, n2, n1, f, Fk, FFTW_ESTIMATE)
!      planb = fftw_plan_dft_c2r_3d(n3, n2, n1, Fk, f, FFTW_ESTIMATE)
      call fftw_execute_dft_r2c(planf, f, Fk)
!$OMP PARALLEL DO PRIVATE(ii,jj,kk,rad2) FIRSTPRIVATE(n1,n2,n3,nn3,nn2,nn1,dk)
      do k=1,n3; if (k>nn3) then; kk = n3+1-k; else; kk=k-1; endif
      do j=1,n2; if (j>nn2) then; jj = n2+1-j; else; jj=j-1; endif
         do i=1,nn1
            rad2 = dble((i-1)**2) + dble(jj**2) + dble(kk**2)
            Fk(i,j,k) = -rad2*dk**2*Fk(i,j,k)
         enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      call fftw_execute_dft_c2r(planb, Fk, f)
!      call fftw_destroy_plan(planf)
!      call fftw_destroy_plan(planb)
      f = f / dble(n1) / dble(n2) / dble(n3)
    end subroutine laplacian_3d

    real*8 function grad_energy_1d(n, f, Fk, dk)
      integer :: n
      real(C_DOUBLE), pointer :: f(:)
      complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
      real*8 :: dk

      integer :: i,ii,nn
      real*8 :: GE

      type(C_PTR) :: planf, planb

      nn = n/2+1
      planf = fftw_plan_dft_r2c_1d(n, f, Fk, FFTW_ESTIMATE)
      planb = fftw_plan_dft_c2r_1d(n, Fk, f, FFTW_ESTIMATE)
      call fftw_execute_dft_r2c(planf, f, Fk)

      GE = 0.
      do i=1,n
         if (i<=nn) then; ii = i-1; else; ii=n+1-i; endif
         GE = GE + Fk(ii+1)*conjg(Fk(ii+1)) * (ii*dk)**2
      enddo
      GE = GE / n / n  ! add normalizations to average and correct unnormalized inverse DFT

      grad_energy_1d = 0.5*GE
    end function grad_energy_1d

  end module fftw3
