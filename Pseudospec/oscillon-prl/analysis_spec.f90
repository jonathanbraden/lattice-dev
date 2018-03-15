!
! Module with various useful analysis routines
!

module analysis
  use fftw3
  use params

!  use params (or store these params in the FFTW module?)

  implicit none

!  interface init_spectrum
!     module procedure init_spectrum_1d, init_spectrum_2d, init_spectrum_3d
!  end interface init_spectrum

  integer, private :: nkx, nky, nkz
  integer :: ns, ns1, ns2

!  interface spectrum
!     module procedure spectrum_1d, spectrum_2d, spectrum_3d
!  end interface spectrum

!  interface cross_spectrum
!     module procedure crossspec_1d, crossspec_2d, crossspec_3d
!  end interface cross_spectrum

contains

  subroutine init_spectrum_1d(nx)
    integer :: nx

    nkx = nx/2+1

    ns1 = nkx
  end subroutine init_spectrum_1d

  subroutine init_spectrum_2d(nx, ny)
    integer :: nx, ny

    nkx = nx/2+1
    nky = ny

    ns2 = (nx/2+1)**2+(ny/2+1)**2
    ns2 = floor(sqrt(dble(ns)+1.))
  end subroutine init_spectrum_2d

  subroutine init_spectrum_3d(nx, ny, nz)
    integer :: nx, ny, nz

    nkx = nx/2+1
    nky = ny
    nkz = nz

    ns = (nx/2+1)**2+(ny/2+1)**2+(nz/2+1)**2
    ns = floor(sqrt(dble(ns)+1.))    
  end subroutine init_spectrum_3d

  subroutine spectrum_3d(S, f, Fk, plan)
    real(dl), dimension(1:ns) :: S
!    real(C_DOUBLE), pointer :: f(:,:,:)
!    complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)
    real(C_DOUBLE) :: f(:,:,:)
    complex(C_DOUBLE_COMPLEX) :: Fk(:,:,:)
    type(C_PTR), optional :: plan

    type(C_PTR) :: plan2
    integer :: i,j,k,ii,jj,kk
    real(dl) :: p, c(2), W(ns)
    integer :: l
    integer :: n1,n2,n3,nn1,nn2,nn3 ! find a better solution to get OMP working

    n1=nx;n2=ny;n3=nz;nn1=nnx;nn2=nny;nn3=nnz
! Do the FFT here
!    plan = fftw_plan_dft_r2c_3d(nz,ny,nx,f,Fk,FFTW_ESTIMATE)
    if (present(plan)) then
       call fftw_execute_dft_r2c(plan, f, Fk)
    else
       plan2 = fftw_plan_dft_r2c_3d(nz,ny,nx,f,Fk,FFTW_ESTIMATE)
       call dfftw_execute_dft_r2c(plan2, f, Fk)
       call dfftw_destroy_plan(plan2)
    endif

    W = 0.
    S = 0.
!$OMP PARALLEL DO PRIVATE(ii,jj,kk,p,l,c) FIRSTPRIVATE(n1,n2,n3,nn1,nn2,nn3) REDUCTION(+:W,S)
    do k=1,n3; if (k<=nn3) then; kk=k-1; else; kk=n3+1-k; endif
       do j=1,n2; if (j<=nn2) then; jj=j-1; else; jj=n2+1-j; endif
          do i=1,n1; if (i<=nn1) then; ii=i-1; else; ii=n1+1-i; endif

             p = sqrt(dble(ii**2 + jj**2 + kk**2)); l=floor(p)
             c = (1.0 - (/l-p,l+1-p/)**2)**2

             S(l+1:l+2) = S(l+1:l+2) + c*Fk(ii+1,j,k)*conjg(Fk(ii+1,j,k))
             W(l+1:l+2) = W(l+1:l+2) + c
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO
    
    where (W /= 0.0) S = S/W!/dble(nx)/dble(ny)/dble(nz)
  end subroutine spectrum_3d

  subroutine spectrum_2d(S, f, Fk, plan)
    real(dl), dimension(1:ns) :: S
    real(C_DOUBLE), pointer :: f(:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:)
    type(C_PTR), optional :: plan

    type(C_PTR) :: plan2
    integer :: i,j,ii,jj
    real(dl) :: p, c(2), W(ns)
    integer :: l

    if (present(plan)) then
       call fftw_execute_dft_r2c(plan, f, Fk)
    else
       plan2 = fftw_plan_dft_r2c_2d(ny,nx,f,Fk,FFTW_ESTIMATE)
       call dfftw_execute_dft_r2c(plan2, f, Fk)
       call dfftw_destroy_plan(plan2)
    endif

    W = 0.
    S = 0.
    do j=1,ny; if (j<=nny) then; jj=j-1; else; jj=ny+1-j; endif
       do i=1,nx; if (i<=nnx) then; ii=i-1; else; ii=nx+1-i; endif

          p = sqrt( dble(ii)**2 + dble(jj)**2 ) ; l=floor(p)
          c = (1.0 - (/l-p,l+1-p/)**2)**2

          S(l+1:l+2) = S(l+1:l+2) + c*Fk(ii+1,j)*conjg(Fk(ii+1,j))
          W(l+1:l+2) = W(l+1:l+2) + c
       enddo
    enddo
    
    where (W /= 0.0) S = S/W/dble(nx)/dble(ny)
  end subroutine spectrum_2d

  subroutine spectrum_1d(S)
    real, dimension(1:ns), intent(out) :: S
    S=0.

  end subroutine spectrum_1d

  subroutine crossspec_3d( Fk1, Fk2, S_real, S_imag)
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:) :: Fk1, Fk2
    real(dl) :: S_real(ns), S_imag(ns)

    real(dl) :: W(ns), p, c(2)
    complex(C_DOUBLE_COMPLEX) :: atmp
    integer :: i,j,k,ii,jj,kk, l
    
! This is a hack to get OMP working
    integer :: n1,n2,n3,nn1,nn2,nn3 ! find a better solution to get OMP working
    n1=nx;n2=ny;n3=nz;nn1=nnx;nn2=nny;nn3=nnz

    W = 0.
    S_real = 0.
    S_imag = 0.
!$OMP PARALLEL DO PRIVATE(ii,jj,kk,p,l,c) FIRSTPRIVATE(n1,n2,n3,nn1,nn2,nn3) REDUCTION(+:W,S_real,S_imag)
    do k=1,n3; if (k<=nn3) then; kk=k-1; else; kk=nz+1-k; endif
       do j=1,n2; if (j<=nn2) then; jj=j-1; else; jj=ny+1-j; endif
          do i=1,n1; if (i<=nn1) then; ii=i-1; else; ii=nx+1-i; endif
             p = sqrt(dble(ii)**2 + dble(jj)**2 + dble(kk)**2); l=floor(p)
             c = (1.0 - (/l-p,l+1-p/)**2)**2
             atmp = Fk1(ii+1,j,k)*conjg(Fk2(ii+1,j,k))

             S_real(l+1:l+2) = S_real(l+1:l+2) + c*real(atmp)
             S_imag(l+1:l+2) = S_imag(l+1:l+2) + c*aimag(atmp)
             W(l+1:l+2) = W(l+1:l+2) + c
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO

    where (W/=0.) S_real = S_real / W !/ (dble(nx)*dble(ny)*dble(nz))
    where (W/=0.) S_imag = S_imag / W !/ (dble(nx)*dble(ny)*dble(nz))
  end subroutine crossspec_3d

  subroutine crossspec_2d()

  end subroutine crossspec_2d

  subroutine crossspec_1d()

  end subroutine crossspec_1d

! More subroutines to add

! Combined Fourier moment, CDF, spectrum code
! Sorting subroutine
! CDF subroutine
! Wavelet subroutine

end module analysis
