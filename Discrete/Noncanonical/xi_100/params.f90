module params
  use, intrinsic :: iso_c_binding
  integer, parameter :: dl = kind(1.d0)
  real(dl), parameter :: twopi = 6.2831853071795864769252867665590

  integer(C_INTPTR_T), parameter :: nx=256, ny=256, nz=256
  integer, parameter :: nnx=nx/2+1, nny=ny/2+1, nnz=nz/2+1
  integer, parameter :: nn = min(nnx,nny,nnz)

  real(dl), parameter :: nvol = dble(nx)*dble(ny)*dble(nz)
  real(dl), parameter :: len = 200.
  real(dl), parameter :: dx = len/dble(nx)
  real(dl), parameter :: dk = twopi / len

  integer(C_INTPTR_T) :: z_size, z_start_index, k_size, k_start_index
  integer :: leftz,rightz, mpisize, mpirank
  integer :: mpierror

end module params
