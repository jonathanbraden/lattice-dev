!
! This is a temporary module until I get the full OO working
!
module latpar
  use mpi
  use p3dfft
  use params

  implicit none

  integer, dimension(1:3) :: nsize
  integer, dimension(1:3) :: istart, iend, padst, paden
  integer, dimension(1:3) :: fstart, fend
  integer, dimension(1:3) :: isize, fsize, sisize
  integer :: mpirank, mpisize
  integer :: rightzrank, leftzrank, rightyrank, leftyrank
  integer :: ydivs, zdivs
  integer :: lzdiv, lydiv  ! these probably aren't needed globally

  integer :: tstep
  real :: time
  real :: conftime

  contains

    subroutine slicelat()

      integer, dimension(1:2) :: dims
      integer :: ny, nz

      nsize(:) = nside(:)
      nz = nside(3)
      ny = nside(2)

      ! determine how many slices in y and z directions
      ydivs = (mpisize/nz) + 1
      zdivs = (mpisize/ydivs)
      
      ! Find indices in y and z direction of lattice on this process
      lydiv = mpirank/zdivs
      lzdiv = mpirank - (lydiv*zdivs)

      ! Find ranks of neighbouring lattice sites
      if (lydiv.eq. (ydivs-1)) then
         rightyrank = lzdiv
      else
         rightyrank = mpirank + zdivs
      endif

      if (lydiv .eq. 0) then
         leftyrank = (ydivs - 1)*zdivs + lzdiv
      else
         leftyrank = mpirank - zdivs
      endif

      if (lzdiv .eq. (zdivs - 1)) then
         rightzrank = mpirank - (zdivs - 1)
      else
         rightzrank = mpirank + 1
      endif

      if (lzdiv .eq. 0) then
         leftzrank = mpirank + (zdivs - 1)
      else
         leftzrank = mpirank - 1
      endif

      ! Now set up P3DFFT
      dims(1) = ydivs
      dims(2) = zdivs

      call p3dfft_setup(dims, nsize(1), nsize(2), nsize(3))
      call get_dims(istart, iend, isize, 1)
      call get_dims(fstart, fend, fsize, 2)

      ! The extra layer is to do double wrapping with the emt
      ! Add in some options here
      padst(:) = istart(:) - 1 - 1
      paden(:) = iend(:) + (p-n) + 1 
      sisize(:) = paden(:) - padst(:) + 1

!      print*, 'rank is ',mpirank, istart, iend

    end subroutine slicelat

end module latpar
