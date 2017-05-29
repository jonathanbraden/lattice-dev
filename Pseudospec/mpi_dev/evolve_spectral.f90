#define POTENTIAL_MACRO_SUM 0.25*fld(1,IRANGE)**4
#define POTENTIAL_MACRO_NOSUM 0.25*fld(1,LATIND)**4

program lattice
#include "macros.h"
#ifdef USEMPI
  use mpi
#endif
  use fftw3
  use params
  use initialize
  use integrator
  use hamiltonian
  use analysis
#ifdef OMP
  use omp_lib
#endif

  implicit none

! Time Stepping Properties
  integer :: j, jj
  integer, parameter :: nstep = 2**10
  integer, parameter :: stepsize = 2**3
  real(dl), parameter :: tstep = dx/10.

  integer :: terror

  integer(kind=8) :: ti1,ti2,clock_rate
  real*8 :: tr1,tr2

  complex(C_DOUBLE_COMPLEX), pointer :: Fk2(:,:,:)

  integer :: error_code, thread_level
!  integer :: mpirank, mpisize
  integer(C_INTPTR_T) :: loc_size

#ifdef USEMPI
  call boot_MPI(mpirank, mpisize)
#else
  mpirank = 0
  mpisize = 1
#endif
  call slice_lattice(nx,ny,nz,z_size, z_start_index)
  if (mpirank == 0) then
     print*,"Running with ", nstep*stepsize," steps and ",nstep," output steps"
     print*,"Lattice size is ",nx,ny,nz
  endif
! these two lines are in initialize.f90 now
  allocate(tmp_lat(SIRANGE))  ! uncomment this for noncanonical case

  call cpu_time(tr1);call system_clock(ti1)
! Initialize arrays for doing spectral laplacian
#ifdef THREEDIM
!  call allocate_arrays_3d(nx,ny,ny,z_size,z_start_index)
  call allocate_arrays_3d_transpose(nx,ny,nz,z_size, z_start_index, k_size, k_start_index)
  call allocate_3d_fourier_array(nx,ny,nz,Fk2)  ! why isn't this transposed?
#endif
#ifdef TWODIM
!  call allocate_2d_array(nx, ny, laplace, Fk)
  call allocate_fftw_array(nx,ny,laplace,Fk)
  planf = fftw_plan_dft_r2c_2d(ny,nx,laplace,Fk,FFTW_PATIENT+FFTW_DESTROY_INPUT)
  planb = fftw_plan_dft_c2r_2d(ny,nx,Fk,laplace,FFTW_PATIENT+FFTW_DESTROY_INPUT)
#endif
#ifdef ONEDIM
!  call allocate_1d_array(nx, laplace, Fk)
  call allocate_fftw_array(nx,laplace,Fk)
  planf = fftw_plan_dft_r2c_1d(nx,laplace,Fk,FFTW_PATIENT+FFTW_DESTROY_INPUT)
  planb = fftw_plan_dft_c2r_1d(nx,Fk,laplace,FFTW_PATIENT+FFTW_DESTROY_INPUT)
#endif
  call cpu_time(tr2);call system_clock(ti2,clock_rate)
  if (mpirank == 0) print*,"FFTW setup :", dble(ti2-ti1)/clock_rate, tr2-tr1

!!!!!!!!!!!!!!!!!!!!!!
! Set up initial conditions, including field fluctuations
!!!!!!!!!!!!!!!!!!!!!!
  call cpu_time(tr1);call system_clock(ti1)
  call init_output()
  call init_fields(1.)
  call make_output(0.)
  call cpu_time(tr2);call system_clock(ti2,clock_rate)
  if (mpirank == 0) print*,"Initialization time :",dble(ti2-ti1)/clock_rate, tr2-tr1

  call cpu_time(tr1);call system_clock(ti1)
  do j=1,nstep
     if (mpirank == 0) print*,"step ", j
     call symp4(tstep, stepsize)
     call make_output(j*stepsize*tstep)
  enddo
  call cpu_time(tr2);call system_clock(ti2,clock_rate)
  if (mpirank == 0) print*,"evolution time :",dble(ti2-ti1)/clock_rate, tr2-tr1
  if (mpirank == 0) print*,"time per step :",dble(ti2-ti1)/clock_rate/dble(nstep*stepsize)

#ifdef USEMPI
  call MPI_Finalize(mpierror)
  print*,"Error in MPI_Finalize",mpierror
#endif

  contains

    subroutine tick(ti, tr, rate)
      integer, intent(out) :: ti, rate
      real*8, intent(out) :: tr

      call cpu_time(tr)
      call system_clock(ti,rate)
    end subroutine tick

    subroutine write_fields(time)
      integer :: i
      real(dl) :: time

      integer :: j,k
      real(dl) :: lap(1:nfld)
      j=ny/2; k=nz/2

      laplace = fld(2,IRANGE)
!      call laplacian_3d(nx, ny, nz,laplace,Fk,dk, planf, planb, z_size, z_start_index)
      call laplacian_spectral(nx,ny,nz,laplace,Fk,dk,planf,planb,k_size,k_start_index)

!      do j=2, ny-1
      do i=2,nx-1
         k=nz/2; j=ny/2
         lap = 1./dx**2/cc*STENCIL(c,LAPLACIAN)
#ifdef THREEDIM
         write(99,*) i*dx, j*dx, fld(:,i,j,k), fldp(:,i,j,k), laplace(i,j,k), lap(:)
#endif
#ifdef TWODIM
         write(99,*) time, i*dx, fld(:,i,j), fldp(:,i,j)
#endif
#ifdef ONEDIM
         write(99,*) time, i*dx, fld(:,i), fldp(:,i)
#endif
      enddo; !write(99,*); enddo
      write(99,*)

    end subroutine write_fields

    subroutine init_fields(vparam)
      real(dl) :: vparam

      integer :: i, j
      real(dl) :: x, x0, y, y0
      real(dl) :: ptemp, gamma
      integer :: nseed
      integer, allocatable :: seed(:)

#ifdef SMOOTHPROFILE
      gamma = 1. / (1.+vparam**2)**0.5
      x0 = len/2.
      y0 = len/2.

      do i=1,nx
         x=i*dx-x0
         ptemp = vparam*cosh(gamma*x)
#ifdef THREEDIM
         fld(:,i,:,:) = 4.*atan(1./ptemp)
#endif
#ifdef TWODIM
         fld(:,i,:) = 4.*atan(1./ptemp)
#endif
#ifdef ONEDIM
         fld(:,i) = 4.*atan(1./ptemp)
#endif
      enddo
      fldp = 0.
#endif

!#ifdef RANDOMIC
      call random_seed(SIZE=nseed)
      allocate(seed(nseed))
      seed = 37*(/ (i-1, i=1,nseed) /)
      call random_seed(PUT=seed)
      deallocate(seed)

      do j=1,nfld
         call sample(-0.25, 3.*fld0(1)**2)
         fld(j,IRANGE) = fld0(j) + laplace(IRANGE)
         call sample(0.25, 3.*fld0(1)**2)
         fldp(j,IRANGE) = dfld0(j) + laplace(IRANGE)
      enddo
!#endif
      yscl = 1.
      call calc_metric()
    end subroutine init_fields

    subroutine dump_rho(time)
      real(dl) :: time
      integer :: l

      real(dl) :: GE, PE, KE, rho, mom
      real(dl) :: elap, lap(nfld)
      integer :: i,j,k
      real(dl) :: acur, fac1,fac2
      real(dl) :: fld_ave, fldp_ave
      
      acur = get_scale_factor()

#ifdef TEMP

#ifdef SPECTRAL
      GE = 0._dl
      do l=1,nfld
         laplace(IRANGE) = fld(l,IRANGE)
#ifdef THREEDIM
!         call laplacian_3d(nx, ny, nz,laplace,Fk,dk, planf, planb, z_size, z_start_index)
         call laplacian_spectral(nx,ny,nz,laplace,Fk,dk,planf,planb,k_size,k_start_index)
#endif
#ifdef TWODIM
         call laplacian_2d(nx, ny, laplace, Fk, dk, planf, planb)
#endif
#ifdef ONEDIM
         call laplacian_1d(nx, laplace, Fk, dk, planf, planb)
#endif

#ifdef VECTORIZE
         GE =  GE - sum(fld(l,IRANGE)*laplace(IRANGE))
#endif
#ifdef LOOPEVOLVE
!$OMP PARALLEL DO REDUCTION(+:GE)
         FLOOP
           GE = GE - fld(l,LATIND)*laplace(LATIND)
         FLOOPEND
!$OMP END PARALLEL DO
#endif
      enddo
      GE = 0.5_dl*GE / nvol
#endif

#ifdef DISCRETE
      call wrap_fields()
      elap = 0.5/dx**2/cc
      GE = 0._dl
!$OMP PARALLEL DO PRIVATE(lap) REDUCTION(+:GE)
      FLOOP
        lap = STENCIL(c,LAPLACIAN)
        GE = GE - sum(fld(:,LATIND)*lap(:))
      FLOOPEND
!$OMP END PARALLEL DO
      GE = elap*GE / nvol
#endif

#ifdef VECTORIZE
!      PE = sum(potential(fld(1,IRANGE)))
!      PE = sum(potential(fld(1,IRANGE),fld(2,IRANGE)))
      PE = sum(POTENTIAL_MACRO_SUM)
      KE = sum(fldp(:,IRANGE)**2)
#endif
#ifdef LOOPEVOLVE
      PE=0.
!$OMP PARALLEL DO REDUCTION(+:PE)
      FLOOP
!        PE = PE + potential(fld(1,LATIND)) 
!        PE = PE + potential(fld(1,LATIND),fld(2,LATIND))
        PE = PE + POTENTIAL_MACRO_NOSUM
      FLOOPEND  
!$OMP END PARALLEL DO
#endif 

#ifdef LOOPEVOLVE
      KE=0.
!$OMP PARALLEL DO REDUCTION(+:KE)      
      FLOOP
        KE = KE + sum(fldp(:,LATIND)**2)
      FLOOPEND
!$OMP END PARALLEL DO
#endif

#ifdef DO_NONCANONICAL
      KE = 0._dl
!$OMP PARALLEL DO REDUCTION(+:KE)
      FLOOP
        fac1 = 1.+zeta*fld(1,LATIND)**2
        fac2 = fac1 / (1.+zetafac*fld(1,LATIND)**2)
        KE = KE + fac1*(fac2*fldp(1,LATIND)**2 + fldp(2,LATIND)**2)
      FLOOPEND
!$OMP END PARALLEL DO
      GE = 0.
      PE = 0.
!$OMP PARALLEL DO REUCTION(+:GE,PE)
      FLOOP
        fac1 = 1.+zeta*fld(1,LATIND)**2
        fac2 = fac1 / (1.+zetafac*fld(1,LATIND)**2)
        GE = GE + blargh * STENCIL(c,GRAD_SQUARED)
        PE = PE + POTENTIAL_MACRO_NOSUM
      FLOOPEND
!$OMP END PARALLEL DO

#endif
      KE = 0.5_dl*KE * kinetic_norm() / nvol
      PE = PE / nvol
      GE = GE / acur**2
      
#ifdef USEMPI
      call MPI_Allreduce(MPI_IN_PLACE,GE,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)
      call MPI_Allreduce(MPI_IN_PLACE,PE,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)
      call MPI_Allreduce(MPI_IN_PLACE,KE,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)
#endif

      rho = KE + PE + GE

#endif  ! This is the endif that finished the commented part above out.  I've moved this into the Hamiltonian

      rho = scalar_rho()
      KE = kinetic_energy()
      GE = grad_energy()
      PE = potential_energy()

      if (mpirank == 0) then
#ifdef THREEDIM
      write(98,'(30(ES22.15,2X))') time, acur, rho, KE, PE, GE, grav_energy(), &
           (rho+grav_energy())/rho, -ysclp**2/12._dl/acur**4, &
           sum(fld(1,IRANGE))/nvol, sum(fldp(1,IRANGE))/nvol
#endif
#ifdef TWODIM
      write(98,*) time, rho, KE, PE, GE, grav_energy(), sum(fld(1,:,:))/nvol, sum(fld(2,:,:))/nvol
#endif
#ifdef ONEDIM
      write(98,*) time, rho, KE, PE, GE, grav_energy(), &
           (rho+grav_energy())/rho, sum(fld(1,:))/nvol, sum(fld(2,:))/nvol
#endif
   endif

    end subroutine dump_rho

!
! Randomly sample a gaussian random field with the appropriate spectrum
!
#define KCUT_FAC 1.
    subroutine sample(gamma, m2eff, spec)
!      real(C_DOUBLE), pointer :: f(:,:,:)
      real(dl) :: gamma
      real(dl) :: m2eff
      real(dl), optional :: spec
      
      type(C_PTR) :: plan_sin
      type(C_PTR) :: plan1
      
      integer, parameter :: os = 16, nos = max(nx,ny,nz)*os**2
      real, parameter :: dxos = dx/os, dkos = dk/(2*os), kcut = KCUT_FAC*min(nnx,nny,nnz)*dk/2.0
      complex, parameter :: w = (0.0, twopi)

      real(dl) :: ker(nos), a(nnx), p(nnx)
      integer, allocatable :: seed(:)
      integer :: nseed

      integer :: i, j, k, l; 
      real(dl) :: kk
#ifdef THREEDIM
      real(dl), parameter :: norm = 0.5/(nvol*(twopi*dk**3)**0.5*mpl)*(dkos/dxos)
#endif
#ifdef TWODIM
      real(dl), parameter :: norm = 0.5/nvol/mpl*(dkos/dxos)  ! Work this one out
#endif
#ifdef ONEDIM
      real(dl), parameter :: norm = 0.5/nvol/mpl*(dkos/dxos)  ! Work this one out
#endif
      
    ! calculate (oversampled) radial profile of convolution kernel
      do k = 1,nos; kk = (k-0.5)*dkos
         ker(k) = kk*(kk**2 + m2eff)**gamma * exp(-(kk/kcut)**2)
      end do
      
!      Assign the kernel here
!      plan = fftw_plan_r2r_1d(nos, ker, ker, FFTW_RODFT10, FFTW_ESTIMATE)
!      call fftw_execute
      
      plan_sin = fftw_plan_r2r_1d(nos,ker,ker,FFTW_RODFT10,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
      call fftw_execute_r2r(plan_sin, ker, ker)
      call fftw_destroy_plan(plan_sin)

      do k = 1,nos; ker(k) = norm * ker(k)/k; end do
       ! initialize 3D convolution kernel (using linear interpolation of radial profile)
         FLOOP
#ifdef THREEDIM
#ifdef USEMPI
            kk = sqrt(dble(i-nn)**2 + dble(j-nn)**2 + dble(k+z_start_index-nn)**2) * os
#else
            kk = sqrt(real(i-nn)**2 + real(j-nn)**2 + real(k-nn)**2) * os
#endif
#endif
#ifdef TWODIM
            kk = sqrt(dble(i-nn)**2 + dble(j-nn)**2)
#endif
#ifdef ONEDIM
            kk = sqrt(dble(i-nn)**2)
#endif
            l = floor(kk)
            
            if (l > 0) then
               laplace(LATIND) = ker(l) + (kk-l)*(ker(l+1)-ker(l))
            else
#ifdef THREEDIM
               laplace(LATIND) = (4.0*ker(1)-ker(2))/3.0
#endif
#ifdef TWODIM
               laplace(LATIND) = (4.0*ker(1)-ker(2))/3.0  ! check this one
#endif
#ifdef ONEDIM
               laplace(LATIND) =
#endif
            end if
      FLOOPEND

       ! convolve kernel with delta-correlated Gaussian noise
!         plan1 = fftw_plan_dft_r2c_3d(nz,ny,nx,f,Fk, FFTW_ESTIMATE)
#ifdef USEMPI
      call fftw_mpi_execute_dft_r2c(planf, laplace, Fk)
#else
      call fftw_execute_dft_r2c(planf, laplace, Fk)
#endif
!         call fftw_destroy_plan(plan1)

#ifdef USEMPI
      do k=1,z_size; do j=1,ny
#else
      do k=1,nz; do j=1,ny
#endif
         call random_number(a); call random_number(p)
         Fk(:,j,k) = sqrt(-2.0*log(a)) * exp(w*p) * Fk(:,j,k)
      enddo; enddo

!         plan1 = fftw_plan_dft_c2r_3d(nz, ny, nx, Fk, f, FFTW_ESTIMATE)
#ifdef USEMPI
      call fftw_mpi_execute_dft_c2r(planb, Fk, laplace)
#else
      call fftw_execute_dft_c2r(planb, Fk, laplace)
#endif
!         call fftw_destroy_plan(plan1)

    end subroutine sample

    subroutine init_output()
      if (mpirank == 0) then
         open(unit=99,file="field_values_spec.out")
         open(unit=98,file="energy_spec.out")
         open(unit=97,file="spectrum.out")
      endif
#ifdef THREEDIM
      call init_spectrum_3d(nx,ny,nz, k_size, k_start_index)
#endif
#ifdef TWODIM
      call init_spectrum_2d(nx,ny)
#endif
#ifdef ONEDIM
      call init_spectrum_1d(nx)
#endif
    end subroutine init_output

    subroutine make_output(time)
      real(dl) :: time
      integer :: i
      real(dl) :: spec(ns,2*nfld+2)
      
!      call write_fields(time)
      call dump_rho(time) 

#ifdef THREEDIM
      laplace(IRANGE) = fld(1,IRANGE)
      call spectrum_3d(spec(:,1),laplace, Fk, planf)
      Fk2=Fk
      laplace(IRANGE) = fldp(1,IRANGE)
      call spectrum_3d(spec(:,2), laplace, Fk, planf)
      call crossspec_3d(Fk, Fk2, spec(:,3),spec(:,4))
#endif
#ifdef TWODIM
      call spectrum_2d(spec, laplace, Fk, planf)
#endif
#ifdef ONEDIM
      call spectrum_1d(spec, laplace, Fk, planf)
#endif

      if (mpirank == 0) then
         do i=1,ns
            write(97,'(30(ES22.15,2x))') (i-1)*dk, spec(i,:)
         enddo
         write(97,*)
      endif
      
    end subroutine make_output
  end program lattice
