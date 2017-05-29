module initialize
#include "macros.h"
  use mpi
  use fftw3
  use params
  use hamiltonian  ! needed for laplace and Fk, ugly solution
#ifdef OMP
  use omp_lib
#endif

  implicit none
!
! To Do : Fix up the boot_MPI and start_threads so that I can do hybrid

contains

  subroutine boot_MPI(rank, size)
    integer, intent(out) :: rank, size
    integer :: mpierror, error_code, thread_level

#ifdef OMP
    call MPI_Init_Thread(MPI_THREAD_FUNNELED, thread_level,mpierror)
    if(mpierror /= MPI_SUCCESS) then
       print*,"Error starting MPI, quitting"
       call MPI_Abort(MPI_COMM_WORLD, error_code, mpierror)
    endif
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,mpierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,mpierror)

    if (rank == 0) print*,"Available level of MPI threading support is ", thread_level," required level for FFTW is ",MPI_THREAD_FUNNELED

    if (thread_level >= MPI_THREAD_FUNNELED) then
       call start_threads_omp(rank)
    else
       if (rank == 0) print*,"Insufficient Threading Support, Not Using"
    endif
#else
    call MPI_Init(mpierror)

    if(mpierror /= MPI_SUCCESS) then
       print*,"Error starting MPI, quitting"
       call MPI_Abort(MPI_COMM_WORLD, error_code, mpierror)
    endif
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,mpierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,mpierror)
#endif
 
    if (rank == 0) then
       print*,"Total number of processes is ", size, "master is ",rank
    endif
    call fftw_mpi_init()
  end subroutine boot_MPI

  subroutine start_threads_omp(rank)
    integer, intent(in) :: rank
    integer :: terror
    terror = fftw_init_threads()
    if (rank == 0) print*,"Threading error code is ",terror
    terror = omp_get_max_threads()
    if (rank == 0) print*,"Number of threads is ",terror
    call fftw_plan_with_nthreads(omp_get_max_threads())
  end subroutine start_threads_omp

  subroutine slice_lattice(L,M,N,zsize, zstart)
    integer(C_INTPTR_T), intent(in) :: L,M,N
    integer(C_INTPTR_T), intent(out) :: zsize, zstart
    integer(C_INTPTR_T) :: loc_size, LL

#ifdef USEMPI
    LL = L/2+1
#ifdef THREEDIM
    loc_size = fftw_mpi_local_size_3d(N,M,LL,MPI_COMM_WORLD, zsize,zstart)
#endif
#ifdef TWODIM
    print*,"Error, 2D MPI slicing not yet implemented"
#endif
#ifdef ONEDIM
    print*,"Error, 1D MPI slicing not yet implemented"
#endif
! Set up the MPI grid
    rightz = mpirank+1
    leftz = mpirank-1
    if ( mpirank == (mpisize-1) ) rightz = 0
    if ( mpirank == 0 ) leftz = mpisize - 1
#else
    zsize = N
    zstart = 0
#endif

    print*,"size ",z_size

! Add some ghost cells if we're using finite-difference derivatives
#ifdef SPECTRAL
    allocate(fld(1:nfld,IRANGE))
    allocate(fldp(1:nfld,IRANGE))
#endif
#ifdef DISCRETE
    allocate(fld(1:nfld,SIRANGE))
    allocate(fldp(1:nfld,IRANGE))
#endif
  end subroutine slice_lattice

  subroutine set_mpi_neighbours()
    
  end subroutine set_mpi_neighbours

  subroutine allocate_arrays_3d(n1,n2,n3,z_len, z_off, k_len, k_off)
    integer(C_INTPTR_T), intent(in) :: n1,n2,n3
    integer(C_INTPTR_T), intent(out) :: z_len, z_off, k_len, k_off

    integer(C_INTPTR_T) :: local_size, LL

    call allocate_3d_array(nx,ny,nz,laplace,Fk)
#ifdef USEMPI
    planf = fftw_mpi_plan_dft_r2c_3d(nz,ny,nx,laplace,Fk,MPI_COMM_WORLD,FFTW_MEASURE+FFTW_DESTROY_INPUT)
    planb = fftw_mpi_plan_dft_c2r_3d(nz,ny,nx,Fk,laplace,MPI_COMM_WORLD,FFTW_MEASURE+FFTW_DESTROY_INPUT)
#else
    planf = fftw_plan_dft_r2c_3d(nz, ny, nx, laplace, Fk, FFTW_MEASURE+FFTW_DESTROY_INPUT) 
    planb = fftw_plan_dft_c2r_3d(nz, ny, nx, Fk, laplace, FFTW_MEASURE+FFTW_DESTROY_INPUT) 
#endif

    LL = n1/2+1
    local_size = fftw_mpi_local_size_3d(n3,n2,LL,MPI_COMM_WORLD,z_len,z_off)
    k_len = z_len
    k_off = z_off
  end subroutine allocate_arrays_3d


  subroutine allocate_arrays_3d_transpose(n1,n2,n3,z_len, z_off, k_len, k_off)
    integer(C_INTPTR_T), intent(in) :: n1,n2,n3
    integer(C_INTPTR_T), intent(out) :: z_len, z_off, k_len, k_off

    integer(C_INTPTR_T) :: local_size, LL

    call allocate_3d_array_mpi_transpose(nx,ny,nz,laplace,Fk)
    planf = fftw_mpi_plan_dft_r2c_3d(nz,ny,nx,laplace,Fk,MPI_COMM_WORLD,FFTW_MEASURE+FFTW_DESTROY_INPUT+FFTW_MPI_TRANSPOSED_OUT)
    planb = fftw_mpi_plan_dft_c2r_3d(nz,ny,nx,Fk,laplace,MPI_COMM_WORLD,FFTW_MEASURE+FFTW_DESTROY_INPUT+FFTW_MPI_TRANSPOSED_IN)

    LL = n1/2+1
    local_size = fftw_mpi_local_size_3d_transposed(n3,n2,LL,MPI_COMM_WORLD,z_len,z_off,k_len,k_off)
  end subroutine allocate_arrays_3d_transpose

  subroutine allocate_arrays_2d()

  end subroutine allocate_arrays_2d

  subroutine allocate_arrays_1d()

  end subroutine allocate_arrays_1d

end module initialize
