module Lattice_Spec
  use constants

  implicit none

  !>@brief
  !> Stores information about a rectangular (but anisotropic) lattice
  type Lattice_Params_Asym
     integer :: dim
     integer, dimension(:), allocatable :: n, nn
     real(dl), dimension(:), allocatable :: dx, len, dk
  end type Lattice_Params_Asym

  !>@brief
  !> Stores information about a cubic (and isotropic) lattice
  type Lattice_Params_Sym
     integer :: dim
     integer :: n, nn
     real(dl) :: dx, len, dk
  end type Lattice_Params_Sym

contains
  
  !>@brief
  !> Create a new lattice definition for an asymmetric rectangular lattice
  subroutine create_lattice_params_a(latpar,nd,n,dx,len)
    type(Lattice_Params_Asym), intent(out) :: latpar
    integer, intent(in) :: nd
    integer, dimension(1:nd), intent(in) :: n
    real(dl), dimension(1:nd), intent(in) :: dx,len
  
    latpar%dim = nd
    allocate( latpar%n(1:nd), latpar%nn(1:nd), latpar%dx(1:nd), latpar%len(1:nd), latpar%dk(1:nd) )
    latpar%n = n; latpar%dx = dx; latpar%len = len
    latpar%dk = twopi / len; latpar%nn = n/2+1
  end subroutine create_lattice_params_a

  !>@brief
  !> Create a new lattice definition of a symmetric cubic lattice
  subroutine create_lattice_params_s(latpar,nd,n,dx,len)
    type(Lattice_Params_Sym), intent(out) :: latpar
    integer, intent(in) :: nd, n
    real(dl), intent(in) :: dx, len

    latpar%dim = nd; latpar%n = n; latpar%dx = dx; latpar%len = len
    latpar%dk = twopi / l; latpar%nn = n/2+1
  end subroutine create_lattice_params_s

  !>@brief
  !> Copy the lattice definition of a symmetric cubic lattice into a new lattice definition
  function copy_lattice_params_s(old) return(lat)
    type(Lattice_Params_Sym) :: lat
    type(Lattice_Params_Sym), intent(in) :: old

    call create_lattice_params_s(lat,old%dim,old%n,old%dx,old%len)
  end function copy_lattice_params_s

  !>@brief
  !> Copy the lattice definition of an asymmetric cubic lattice into a new lattice definition
  function copy_lattice_params_a(old) return(lat)
    type(Lattice_Params_Asym) :: lat
    type(Lattice_Params_Asym), intent(in) :: old

    call create_lattice_params_a(lat,old%dim,old%n,old%dx,old%len)
  end function copy_lattice_params_a

  !Move this somewhere else
  subroutine wrap_field(fld)
    real(dl), dimension(:,:,:), intent(inout) :: fld
  end subroutine wrap_field

end module Lattice_Spec
