program oolattice

  use mpi
  use p3dfft
  use params
  use modparam
  use model
  use init_cond
  use latpar
  use scafld
  use evolve
  use iofile

#include "macros.h"

  implicit none

  integer :: i
  integer :: outsteps, stepsize
!  real :: dt
!  integer :: mpierror
!  integer :: mpisize, mpirank

  integer :: rssize
  integer, allocatable :: rs(:)

  outsteps = tt/nt  ! only for debugging purposes
  stepsize = nt

! Start MPI
  call MPI_Init(mpierror)
  call MPI_Comm_size(MPI_COMM_WORLD, mpisize, mpierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpierror)

  call slicelat()   ! modify this subroutine

!  print*,"lattice in ",padst, paden

#ifndef SEED
#define SEED 2
#endif

!
! Check that this procedure is thread safe!!!!
!
  call random_seed(SIZE=rssize)
  allocate(rs(rssize))
  rs = SEED + rssize*mpirank +(/ (i-1, i=1,rssize) /)
  seednum = SEED
  call random_seed(PUT=rs)


! Now initialize the fields
  call allocate_fields()
  allocate(tmp(IRANGE))
  allocate(tmp2(IRANGE))
  allocate(Fk(FRANGE))

! As a check, find out what the ranges are for the Fourier transforms
!  print*,"rank is ",mpirank," with fourier ranges ", fstart,fend

  call init_output()

!
! To do: add option to resume an old run and put in a different module (so that we can adjust as needed, specifically the read_metric piece
!
!
  if (resume$run) then
     do i=1,fields
        call read_check_pt("fld_"//int2str(i), fld(i, IRANGE), istart, iend)
        call read_check_pt("fldp_"//int2str(i), fldp(i, IRANGE), istart, iend)
     enddo
     if (mpirank.eq.0) then
        open(unit=99,file="checkpt.dat")
        read(99,*) yscl, ysclp
        read(99,*) time
        read(99,*) outframenum
        close(99)
     endif
     
     call MPI_Bcast(yscl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
     call MPI_Bcast(ysclp, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
     call MPI_Bcast(time, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
     call MPI_Bcast(outframenum, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror)
     
  else
     outframenum = 0
     call init_metric()
     call init_fields(fld(:,:,:,:), fldp(:,:,:,:))

     call need_get_energies()
     call correct_metric()
     call make_output()
  endif

! Do the actual time-evolution
  do i=1,outsteps
     call step(stepsize)
     time = time + stepsize*dt
     outframenum = i
     call make_output()
     if (mpirank.eq.0) print*, "step ",outframenum
  enddo

! Do any additional output that is needed
  if (write$checkpt) then
     do i=1,fields
        call make_check_pt("fld_"//int2str(i), fld(i,IRANGE), istart, iend)
        call make_check_pt("fldp_"//int2str(i), fldp(i,IRANGE), istart, iend)
     enddo
     if (mpirank.eq.0) then
        open(unit=99, file="checkpt.dat")
        write(99,*) yscl, ysclp
        write(99,*) time
        write(99,*) outframenum
        close(99)
     endif
  endif


! Stop MPI
  call p3dfft_clean()
  call MPI_Finalize(mpierror)

end program oolattice
