module output
  use fftw3
  use Hamiltonian
  implicit none

#ifdef THREEDIM
  real*8, allocatable :: emtens(:,:,:,:), grad2(:,:,:)
  complex(C_DOUBLE_PRECISION), pointer :: Fk2(:,:,:)
#endif

#ifdef TWODIM
  real*8, allocatable :: emtens(:,:,:)
  complex(C_DOUBLE_PRECISION), pointer :: Fk2(:,:)
#endif

#ifdef ONEDIM
  real*8, allocatable :: emtens(:,:)
  complex(C_DOUBLE_PRECISION), pointer :: Fk2(:)
#endif

contains

  subroutine initialize_output_arrays_3d(nx,ny,nz)
    integer, intent(in) :: nx,ny,nz

    allocate(emtens(1:2,1:nx,1:ny,1:nz))
    allocate(grad2(1:nx,1:ny,1:nz))
    call allocate_3d_fourier_array(nx,ny,nz,Fk2)
  end subroutine initialize_output_arrays_3d

  subroutine initialize_output_arrays_2d()

  end subroutine initialize_output_arrays_2d

  subroutine initialize_output_arrays_1d()

  end subroutine initialize_output_arrays_1d

!!!!!!!
! User supplied subroutine to compute any desired variables
!!!!!!!
  subroutine compute_variables()
    integer :: LATIND
    integer :: l

    emtens(1,IRANGE) = 0.5*sum(fldp(1:nfld,IRANGE)**2,DIM=1)
    emtens(2,IRANGE) = 0.5*sum(fldp(1:nfld,IRANGE)**2,DIM=2)

    do l=1,nfld
       laplace(IRANGE) = fld(l,IRANGE)
       call gradient_squared_3d_spectral( (/nx,ny,nz/) , laplace, Fk, Fk2, grad2, , planf, planb)
       emtens(1,IRANGE) = emtens(1,IRANGE) + 0.5*grad2(IRANGE)
       emtens(2,IRANGE) = emtens(2,IRANGE) - (0.5/3.)*grad2(IRANGE)
    enddo

!    FLOOP
!       emtens(1,LATIND) = emtens(1,LATIND) + potential()
!       emtens(2,LATIND) = emtens(2,LATIND) + potential()
!    FLOOPEND

  end subroutine compute_variables

  subroutine dump_rho(time)
    real(dl) :: time
    
    integer :: l
    real(dl) :: GE, PE, KE, rho, mom
    real(dl) :: elap, lap(nfld)
    integer :: LATIND
    real(dl) :: acur, fac1, fac2

#ifdef SPECTRAL
    GE = 0._dl
    do l=1,nfld
       laplace = fld(l,IRANGE)
       
    enddo
#endif

#ifdef DISCRETE
    call wrap_fields()
#endif

end module output
