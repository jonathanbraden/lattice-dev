module params
  integer, parameter :: dl = kind(1.d0)
  real(dl), parameter :: twopi = 6.2831853071795864769252867665590_dl

  integer, parameter :: nx=64, ny=64, nz=64
  integer, parameter :: nnx=nx/2+1, nny=ny/2+1, nnz=nz/2+1
  integer, parameter :: nn = min(nnx,nny,nnz)

  real(dl), parameter :: nvol = dble(nx)*dble(ny)*dble(nz)
  real(dl), parameter :: len = 0.01_dl / (6.8041381743977173e-4)  ! Fix this to be an actual computation
  real(dl), parameter :: dx = len/dble(nx)
  real(dl), parameter :: dk = twopi / len

end module params
