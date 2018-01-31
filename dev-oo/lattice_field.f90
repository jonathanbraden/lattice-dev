#define SCALAR_TYPE_3D Scalar_3d
#define SCALAR_TYPE_2D Scalar_2d
#define SCALAR_TYPE_1D Scalar_1d

#define LATTICE_TYPE   Lattice_Params_Sym

! TO DO : Add padding ghost cells to the lattice storage

module Lattice_Fields
  use Lattice_Spec
  implicit none

  type SCALAR_TYPE_3D
     ! Pointer to the lattice specs vs. carrying around a copy
!     type(LATTICE_TYPE), pointer  
     type(LATTICE_TYPE) :: lat
     integer :: nfld
     real(dl), dimension(:,:,:,:), allocatable :: fld, fldp
  end type SCALAR_TYPE_3D

  type SCALAR_TYPE_2D
     integer :: nfld
     real(dl), dimension(:,:,:), allocatable :: fld, fldp
  end type SCALAR_TYPE_2D

  type SCALAR_TYPE_1D
     integer :: nfld
     real(dl), dimension(:,:), allocatable :: fld, fldp
  end type SCALAR_TYPE_1D

contains

  subroutine create_scalar_3d(sc,nf,lat)
    type(SCALAR_TYPE_3D), intent(out) :: sc
    integer, intent(in) :: nf
    type(LATTICE_TYPE), intent(in) :: lat

    sc%lat = copy_lattice_params_s(lat)
    sc%nfld = nf
    allocate(sc%fld(1:nf,1:lat%n,1:lat%n,1:lat%n))
    allocate(sc%fldp(1:nf,1:lat%n,1:lat%n,1:lat%n))
  end subroutine create_scalar_3d

  subroutine create_scalar_2d(sc,nf,lat)
    type(SCALAR_TYPE_2D), intent(out) :: sc
    integer, intent(in) :: nf
    type(LATTICE_TYPE), intent(in) :: lat
    
    sc%lat = copy_lattice_params_s(lat)
    sc%nfld = nf
    allocate(sc%fld(1:nf,1:lat%n,1:lat%n))
    allocate(sc%fldp(1:nf,1:lat%n,1:lat%n))
  end subroutine create_scalar_2d

  subroutine create_scalar_1d(sc,nf,lat)
    type(SCALAR_TYPE_1D), intent(out) :: sc
    integer, intent(in) :: nf
    type(LATTICE_TYPE), intent(in) :: lat

    sc%lat = copy_lattice_params_s(lat)
    sc%nfld = nf
    allocate(sc%fld(1:nf,1:lat%n))
    allocate(sc%fldp(1:nf,1:lat%n))
  end subroutine create_scalar_1d()
  
end module Lattice_Fields
