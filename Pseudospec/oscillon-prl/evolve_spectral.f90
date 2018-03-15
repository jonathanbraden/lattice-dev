!
! Module to store the desired evolution scheme for the fields, etc.
!

!
! To do : Add packing and unpacking of phase space variables
!

!
! To do: Add the Minkowski evolution and fixed background evolution
!   make sure I have the macros in the right place
!

! Add some macros to allow more seamless switching between numbers of dimensions

!#define THREADS 1
#define OMP 1
!#define USEMPI

!#define VECTOR(A,B) (/A,B/)
! Define the model parameters
!#define DVDPHI sin(fld(l,IRANGE))
!#define POTENTIAL -cos(fld(IRANGE))
!#define DVDPHI(PHI,PSI) VECTOR((PHI)**3 + g2*(PHI)*(PSI)**2, g2*(PHI)**2)
!#define POTENTIAL 0.25*fld(1,IRANGE)**4 + 0.5*g2*fld(1,IRANGE)**2*fld(2,IRANGE)**2
!#define M2I(PHI,PSI) VECTOR((PHI)**2+g2*(PSI)**2, g2*(PHI)**2)
!#define NUMFIELD 2

#define EXACT_HUB 1

program lattice
#include "macros.h"
  use fftw3
  use params
  use hamiltonian
  use analysis
#ifdef OMP
  use omp_lib
#endif

  implicit none

! Time Stepping Properties
  integer :: j, jj
  integer, parameter :: nstep=2**0, stepsize = 2**0 !nstep = 2**14, stepsize = 2**5
  real(dl), parameter :: tstep = dx/10._dl

  integer :: terror

  integer(kind=8) :: ti1,ti2,clock_rate
  real*8 :: tr1,tr2

  complex(C_DOUBLE_COMPLEX), pointer :: Fk2(:,:,:)

! This is currently just for debugging
  real(C_DOUBLE), dimension(nx,ny,nz) :: grad_squared_s, grad_squared_d

  print*,"Running with ", nstep*stepsize," steps and ",nstep," output steps, symp6"
  print*,"Lattice size is ", nx,ny,nz

  print*,"mpl is ",mpl
  print*,"vev is ",vev
  print*,"V0 is ", V0, " and in planck units ", V0/mpl2
  print*,"mphi in Planck units ", 1._dl/mpl
  print*,"dt is ",tstep

! Begin by initializing FFTW
#ifdef THREADS
  terror = fftw_init_threads()
  call fftw_plan_with_nthreads(4)
#endif
#ifdef OMP
  terror = fftw_init_threads()
  print*,"Error code is ",terror
  terror = omp_get_max_threads()
  print*,"Num Threads = ",terror
  !  call fftw_plan_with_nthreads(omp_get_max_threads())
  call fftw_plan_with_nthreads(1)
#endif

! initialize arrays for doing FFT
  call cpu_time(tr1); call system_clock(ti1)
  call initialize_arrays()
  call cpu_time(tr2); call system_clock(ti2,clock_rate)
  print*,"FFTW setup :", dble(ti2-ti1)/clock_rate, tr2-tr1

  call cpu_time(tr1); call system_clock(ti1)
  call init_output()
!  call init_fields(phi0,dphi0,mpl)
  call init_fields_w_fluc()
  call cpu_time(tr2); call system_clock(ti2,clock_rate)
  print*,"Field Initialization time :",dble(ti2-ti1)/clock_rate, tr2-tr1

  call cpu_time(tr1); call system_clock(ti1)
  call make_output(0.)
  call cpu_time(tr2); call system_clock(ti2,clock_rate)
  print*,"Initial Output Time :",dble(ti2-ti1)/clock_rate, tr2-tr1

  call cpu_time(tr1); call system_clock(ti1)
  do j=1,nstep
     print*,"step ", j
     call symp6(tstep, stepsize)
     call make_output(j*stepsize*tstep)
  enddo
  call cpu_time(tr2); call system_clock(ti2,clock_rate)
  print*,"Total Evolution Time :",dble(ti2-ti1)/clock_rate, tr2-tr1

  contains

    subroutine initialize_threads()

    end subroutine initialize_threads

    subroutine initialize_arrays()
#ifdef THREEDIM
      call allocate_fftw_array(nx,ny,nz,laplace,Fk)
      planf = fftw_plan_dft_r2c_3d(nz, ny, nx, laplace, Fk, FFTW_MEASURE) 
      planb = fftw_plan_dft_c2r_3d(nz, ny, nx, Fk, laplace, FFTW_MEASURE) 

      call allocate_3d_fourier_array(nx,ny,nz,Fk2)
#endif

#ifdef TWODIM
      call allocate_fftw_array(nx,ny,laplace,Fk)
      planf = fftw_plan_dft_r2c_2d(ny,nx,laplace,Fk,FFTW_PATIENT+FFTW_DESTROY_INPUT)
      planb = fftw_plan_dft_c2r_2d(ny,nx,Fk,laplace,FFTW_PATIENT+FFTW_DESTROY_INPUT)
#endif

#ifdef ONEDIM
      call allocate_fftw_array(nx,laplace,Fk)
      planf = fftw_plan_dft_r2c_1d(nx,laplace,Fk,FFTW_PATIENT+FFTW_DESTROY_INPUT)
      planb = fftw_plan_dft_c2r_1d(nx,Fk,laplace,FFTW_PATIENT+FFTW_DESTROY_INPUT)
#endif
    end subroutine initialize_arrays

    subroutine symp8(dt, nsteps)
      real*8 :: dt
      integer :: nsteps

      integer :: j
      real*8, parameter :: w1 = 0.74167036435061295344822780
      real*8, parameter :: w2 = -0.40910082580003159399730010
      real*8, parameter :: w3 = 0.19075471029623837995387626
      real*8, parameter :: w4 = -0.57386247111608226665638733
      real*8, parameter :: w5 = 0.29906418130365592384446354
      real*8, parameter :: w6 = 0.33462491824529818378495798
      real*8, parameter :: w7 = 0.31529309239676659663205666
      real*8, parameter :: w0 = 1._dl - 2._dl*(w1+w2+w3+w4+w5+w6+w7)

      call Hamiltonian_Split(w7*dt/2._dl,1)

      do j=1,nsteps
         call symp_o2step(dt, w7, w6)
         call symp_o2step(dt, w6, w5)
         call symp_o2step(dt, w5, w4)
         call symp_o2step(dt, w4, w3)
         call symp_o2step(dt, w3, w2)
         call symp_o2step(dt, w2, w1)
         call symp_o2step(dt, w1, w0)
         call symp_o2step(dt, w0, w1)
         call symp_o2step(dt, w1, w2)
         call symp_o2step(dt, w2, w3)
         call symp_o2step(dt, w3, w4)
         call symp_o2step(dt, w4, w5)
         call symp_o2step(dt, w5, w6)
         call symp_o2step(dt, w6, w7)
         if (j.eq.nsteps) then
            call symp_o2step(dt, w7, 0._dl)
         else
            call symp_o2step(dt, w7, w7)
         endif
      enddo

    end subroutine symp8

    subroutine symp6(dt, nsteps)
      real(dl) :: dt
      integer :: nsteps

      real(dl), parameter :: w3=0.784513610477560_dl
      real(dl), parameter :: w2=0.235573213359357_dl
      real(dl), parameter :: w1=-1.177679984177887_dl
      real(dl), parameter :: w0=1._dl - 2._dl*(w1+w2+w3)

      integer :: j

      call Hamiltonian_Split(w3*dt/2._dl,1)
      do j=1,nsteps
         call symp_o2step(dt,w3,w2)
         call symp_o2step(dt,w2,w1)
         call symp_o2step(dt,w1,w0)
         call symp_o2step(dt,w0,w1)
         call symp_o2step(dt,w1,w2)
         call symp_o2step(dt,w2,w3)
         if (j==nsteps) then
            call symp_o2step(dt,w3,0._dl)
         else
            call symp_o2step(dt,w3,w3)
         endif
      enddo
    end subroutine symp6

    subroutine symp4(dt, nsteps)
      real(dl) :: dt
      integer :: nsteps

      real(dl), parameter :: w1 = 1._dl/(2._dl-2._dl**(1._dl/3._dl))
      real(dl), parameter :: w0 = 1._dl - 2._dl*w1

      integer :: j
      
      call Hamiltonian_Split(w1*dt/2._dl,1)
      do j=1,nsteps
         call symp_o2step(dt, w1,w0)
         call symp_o2step(dt,w0,w1)
         if (j==nsteps) then
            call symp_o2step(dt,w1,0._dl)
         else
            call symp_o2step(dt,w1,w1)
         endif
      enddo
    end subroutine symp4

    subroutine symp2(dt,nsteps)
      real(dl) :: dt
      integer :: nsteps

      integer :: j
      call Hamiltonian_Split(dt/2._dl, 1)
      do j=1,nsteps-1
         call symp_o2step(dt,1._dl,1._dl)
      enddo
      call symp_o2step(dt, 1._dl,0._dl)
    end subroutine symp2

    subroutine write_fields(time)
      integer :: i
      real(dl) :: time

      integer :: j,k
      real(dl) :: lap(1:nfld), grad_sq(1:nfld)
      real(dl), dimension(1:3) :: GE
      j=ny/2; k=nz/2

      laplace = fld(1,IRANGE)
      GE(1) = 0.5_dl*gradient_energy_3d(nx,ny,nz,laplace,Fk,dk,planf)
      laplace = fld(1,IRANGE)
      call gradient_squared_3d_spectral([nx,ny,nz],laplace,Fk,Fk2,grad_squared_s,dk,planf,planb)
      laplace = fld(1,IRANGE)
      call laplacian_3d(nx, ny, nz,laplace,Fk,dk, planf, planb)
      GE(2) = -0.5_dl*sum(fld(1,IRANGE)*laplace(IRANGE))/nvol
      GE(3) = 0.5_dl*sum(grad_squared_s(IRANGE))/nvol
      print*,"Mean gradient energies are ",GE
      
      do j=2,ny-1
         k=nz/2; i=nx/2
         lap = (1._dl/dx**2/cc)*(STENCIL(c,LAPLACIAN))
         grad_sq = (0.5_dl/dx**2/cc)*(STENCIL(c,GRAD2)) ! check factors here
#ifdef THREEDIM
!         write(99,*) i*dx, j*dx, fld(:,i,j,k), fldp(:,i,j,k), laplace(i,j,k), lap(:)
         write(99,'(30(ES22.15,2x))') i*dx, j*dx, fld(:,i,j,k), fldp(:,i,j,k), grad_squared_s(i,j,k), grad_sq(1), laplace(i,j,k), lap(1)
#endif
#ifdef TWODIM
         write(99,*) time, i*dx, fld(:,i,j), fldp(:,i,j)
#endif
#ifdef ONEDIM
         write(99,*) time, i*dx, fld(:,i), fldp(:,i)
#endif
      enddo;
      write(99,*)

    end subroutine write_fields

    subroutine init_fields(phi0,dphi0,mpl)
      real(dl), intent(in) :: phi0, dphi0,mpl
      fld = phi0; fldp = dphi0; yscl = 1._dl
#ifdef EXACT_HUB
      call calc_metric()
#else
      ysclp = -6._dl*H0  ! approximate Hubble (used for testing)
#endif
    end subroutine init_fields

    subroutine init_mean_fields(phi0,dphi0)
      real(dl), intent(in) :: phi0, dphi0
      fld = phi0; fldp = dphi0; yscl = 1._dl
      call calc_metric()
    end subroutine init_mean_fields

    ! This thing is hiding horrible nonlocality, fix it
    subroutine init_fluctuations(mpl)
      real(dl), intent(in) :: mpl
      integer :: j

      do j=1,nfld
         call sample(-0.25, 3.*phi0**2)
         fld(j,IRANGE) = fld(j,IRANGE) + laplace
         call sample(0.25, 3.*phi0**2)
         fldp(j,IRANGE) = fldp(j,IRANGE) + laplace
      enddo
      yscl = 1._dl
      call calc_metric()
    end subroutine init_fluctuations

    subroutine init_fields_w_fluc()
      integer :: i, j
      real(dl) :: x, x0, y, y0
      real(dl) :: ptemp, gamma
      integer :: nseed
      integer, allocatable :: seed(:)

      call random_seed(SIZE=nseed)
      allocate(seed(nseed))
      seed = 37*(/ (i-1, i=1,nseed) /)
      call random_seed(PUT=seed)
      deallocate(seed)

      do j=1,nfld
         call sample(-0.25, 3.*fld0(1)**2)
         fld(j,IRANGE) = fld0(j) + laplace
         call sample(0.25, 3.*fld0(1)**2)
         fldp(j,IRANGE) = dfld0(j) + laplace
      enddo
      yscl = 1.
      call calc_metric()
    end subroutine init_fields_w_fluc
    
!
! Randomly sample a gaussian random field with the appropriate spectrum
!
! Check this doesn't screw up when the k_x mode is homogeneous, since in this case we aren't explicitly forcing the correct symmetry of the matrix to ensure a real field.  Probably doesn't matter
#define KCUT_FAC 0.5
    subroutine sample(gamma, m2eff, spec)
!      real(C_DOUBLE), pointer :: f(:,:,:)
      real(dl) :: gamma
      real(dl) :: m2eff
      real(dl), optional :: spec
      
      type(C_PTR) :: plan_sin
      type(C_PTR) :: plan1
      
      integer, parameter :: os = 16, nos = max(nx,ny,nz)*os**2
      real, parameter :: dxos = dx/os, dkos = dk/(2*os), kcut = KCUT_FAC*min(nnx,nny,nnz)*dk/2._dl
      complex, parameter :: w = (0._dl, twopi)

      real(dl) :: ker(nos), a(nnx), p(nnx)
      integer, allocatable :: seed(:)
      integer :: nseed

      integer :: i, j, k, l; 
      real(dl) :: kk

      real(dl) :: rad2

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
      
      plan_sin = fftw_plan_r2r_1d(nos,ker,ker,FFTW_RODFT10,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
      call fftw_execute_r2r(plan_sin, ker, ker)
      call fftw_destroy_plan(plan_sin)

      do k = 1,nos; ker(k) = norm * ker(k)/k; end do
       ! initialize 3D convolution kernel (using linear interpolation of radial profile)
         FLOOP
#ifdef THREEDIM
            kk = sqrt(real(i-nn)**2 + real(j-nn)**2 + real(k-nn)**2) * os
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
               laplace(LATIND) = (4._dl*ker(1)-ker(2))/3._dl
#endif
#ifdef TWODIM
               laplace(LATIND) = (4._dl*ker(1)-ker(2))/3._dl  ! check this one
#endif
#ifdef ONEDIM
               laplace(LATIND) = (4._dl*ker(1)-ker(2))/3._dl
#endif
            end if
      FLOOPEND

      ! convolve kernel with delta-correlated Gaussian noise
      call fftw_execute_dft_r2c(planf, laplace, Fk)

      do k=1,nz; do j=1,ny
         call random_number(a); call random_number(p)
         Fk(:,j,k) = sqrt(-2.0*log(a)) * exp(w*p) * Fk(:,j,k)
      enddo; enddo

      call fftw_execute_dft_c2r(planb, Fk, laplace)  
!      call fftw_execute_dft_r2c(planf, laplace, Fk)
    end subroutine sample

#ifdef OSCILLON_START
    subroutine init_fields(r0,eps,phi0)
      real(dl), intent(in) :: r0, eps, phi0
      integer :: i,j,k
      real(dl) :: r2

      do k=1,nz; do j=1,ny; do i=1,nx
         r2 = ( (i-nx/2)**2 + (j-ny/2)**2/(1+eps)**2 + (k-nz/2)**2) * (dx/r0)**2
         fld(1,i,j,k) = 1. + phi0*exp(-0.5*r2)
      enddo; enddo; enddo
      fldp = 0._dl

      yscl = 1.
      call calc_metric()
    end subroutine init_fields
#endif

    subroutine dump_rho(time)
      real(dl) :: time
      integer :: l

      real(dl) :: GE, PE, KE, rho, mom
      real(dl) :: elap, lap(nfld)
      integer :: i,j,k
      real(dl) :: acur, fac1,fac2
      
      acur = get_scale_factor()

#ifdef SPECTRAL
      GE = 0._dl
      do l=1,nfld
         laplace = fld(l,IRANGE)
#ifdef THREEDIM
         call laplacian_spectral(nx, ny, nz,laplace,Fk,dk, planf, planb)
#endif
#ifdef TWODIM
         call laplacian_spectral(nx, ny, laplace, Fk, dk, planf, planb)
#endif
#ifdef ONEDIM
         call laplacian_spectral(nx, laplace, Fk, dk, planf, planb)
#endif

#ifdef VECTORIZE
         GE =  GE - sum(fld(l,IRANGE)*laplace(IRANGE))
#endif
#ifdef LOOPEVOLVE
         FLOOP
           GE = GE - fld(l,LATIND)*laplace(LATIND)
         FLOOPEND  
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
        GE = GE - sum(fld(:,i,j,k)*lap(:))
      FLOOPEND
!$OMP END PARALLEL DO
      GE = elap*GE / nvol
#endif

#ifdef VECTORIZE
      PE = sum(potential(fld(1,IRANGE))) !sum(potential(fld(1,IRANGE),fld(2,IRANGE)))
      KE = sum(fldp(:,IRANGE)**2)
#endif
#ifdef LOOPEVOLVE
      PE=0.
!$OMP PARALLEL DO REDUCTION(+:PE)
      FLOOP
        PE = PE + potential(fld(1,LATIND)) !potential(fld(1,LATIND),fld(2,LATIND))
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
      KE = 0.5_dl*KE * kinetic_norm() / nvol
      PE = PE / nvol

      GE = GE / acur**2
      
      rho = KE + PE + GE

#ifdef THREEDIM
      write(98,'(30(ES22.15,2X))') time, acur, rho, KE, PE, GE, grav_energy(), &
           (rho+grav_energy())/rho, get_hubble(), &
           sum(fld(1,IRANGE))/nvol, sum(fldp(1,IRANGE))/nvol
#endif
#ifdef TWODIM
      write(98,*) time, rho, KE, PE, GE, grav_energy(), sum(fld(1,:,:))/nvol, sum(fld(2,:,:))/nvol
#endif
#ifdef ONEDIM
      write(98,*) time, rho, KE, PE, GE, grav_energy(), &
           (rho+grav_energy())/rho, sum(fld(1,:))/nvol, sum(fld(2,:))/nvol
#endif

    end subroutine dump_rho

    subroutine init_output()
      open(unit=99,file="field_values_spec.out")
      open(unit=98,file="energy_spec.out")
      open(unit=97,file="spectrum.out")
#ifdef THREEDIM
      call init_spectrum_3d(nx,ny,nz)
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
      
      call write_fields(time)
      call dump_rho(time)
      
!      laplace(IRANGE) = fld(2,IRANGE)
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
      do i=1,ns
         write(97,'(30(ES22.15,2x))') (i-1)*dk, spec(i,:)
      enddo
      write(97,*)
      
    end subroutine make_output
  end program lattice
