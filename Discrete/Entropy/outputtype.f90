!
! A class designed specifically for output of the lattice code
!
! I define and outputbuffer type and then the corresponding operations to create the data for the output buffer and to flush it to file
!

! Define some preprocessors in case I decide to define a type (instead of just declaring variables) to store all the output arrays
! If I do, justt append name% in front of everything and include passing the object into the subroutines below
#define IDX idx
#define IDXF idxf
#define IDXFM idxfm
#define VNAMES VNAMES

module OutputType

  use mpi
  use params
  use p3dfft
  use model
  use modparam
  use analysis

  implicit none

  character(10) :: fpos, fstat

!  Type OutputBuffer
  logical :: get$PSD
  logical :: get$CDF
  logical :: get$fouriermoms
  logical :: get$fourierdist
  
  integer, parameter :: numkbins = ns
  integer :: numbins
  integer :: nummoms
  real(dl), parameter :: pdfsig = 3.
  
  integer, parameter :: maxstats = 2*fields + 15
  integer, parameter :: maxfstats = fields + 1
  
  integer :: idx, idxf, idxfm
  
  !# Add some preprocessors here to avoid using extra storage space
  real(dl), dimension(1:numkbins, 1:nummoms, 1:2, 1:maxfstats) :: PSD_MOMS
  real(dl), dimension(1:numkbins, 1:numbins, 1:maxfstats) :: PSD_DIST
  real(dl), dimension(1:n+1, 1:maxstats) :: CDF
  real(dl), dimension(1:ns, 1:maxstats) :: PSD
  character(12), dimension(1:maxstats) :: DVAR
!  end Type OutputBuffer

!
! Store the current file numbers to be used for output
!
  logical :: validfiles  ! checks to make sure output files have been set
  integer :: ampfilen, phasefilen, refilen, imfilen
  integer :: psdfilen, cdffilen

  contains

    subroutine dump_buffer(v, t, f, st, en, dolog, dofstats
      character(*) :: v
      real(dl) :: t
      integer, dimension(1:3) :: st, en
      real(dl), dimension(st(1):en(1),st(2):en(2),st(3):en(3)) :: f
      logical :: dolog
      logical :: dofstats
      
      IDX = IDX + 1
      VNAME(IDX) = v

      if (dofstats) then
         IDXF = IDXF + 1
         IDXFM = IDXFM + 1
         call get_fourier_stats(f, st, en
      endif

      call get_realspace_stats(f, st, en

    end subroutine dump_buffer

!
! Input : file numbers for the various output files
!
    subroutine flush_buffer(t, scl, hub, 
      real(dl) :: t, scl, hub

      logical :: o

      if (mpirank.ne.0) return

!# Add a preprocessor
      if (output$psd) then
         inquire(psdfilen, opened=o) 
         if (.not.o) then
            open(unit=psdfilen, file=trime(psdname), position=fpos, status=fstat)
            call head(psdfilen, (/"t   ","a   ","H    ","k    ",VNAMES(1:IDX)/))
         endif

         do k=1,ns
            write(psdfilen,'(32g25.16e2)') t, scl, hub, (k-1)*dk, PSD(k,1:IDX)
         enddo
         write(psdfilen,'(g)') "",""
      endif


!# Add a preprocessor
      if (output$cdf) then
         inquire(cdffilen, opened=o)
         if (.not.o) then
            open(unit=cdffilen, file=trim(cdfname), position=fpos, status=fstat)
            call head(cdffilen, (/"t     ","percentile ",VNAMES(1:IDX)/))
         endif

         do k=1,n+1
            write(cdffilen,'(32g25.16e2)') t, real(k-1)/n, CDF(k,1:IDX)           
         enddo

         write(cdffilen,'(g)') "",""
         flush(cdffilen)
      endif

!# Add a preprocessor call
      if (IDXF .gt. 0) then   ! check if we've stored anything to output
         if (output$fourierprob) then
            do i=1,n
               write(ampfilen,'(32g)') t, (real(k-1)/ ), (k,1:IDXF)
               write(phasefilen,'(32g)') t, real(k-1)*twopi/dble(n) - twopi/2._dl, (k,1:IDXF)
               write(refilen,'(32g)') t, real(k-n/2)/, (k,1:IDXF)
               write(imfilen,'(32g)') t, real(k-n/2)/, (k,1:IDXF)
            enddo
            write(ampfilen,'(g)') "",""
            write(phasefilen,'(g)') "",""
            write(refilen,'(g)') "",""
            write(imfilen,'(g)') "",""
         endif
      endif

      if (IDXFM .gt. 0) then
         if (output$fouriermoms) then
            do k=1,ns
               write(fourmomfn,'(32g)') t, k*dk, , , ,
            enddo
            write(fourmomfn,'(g)') "", ""
         endif

      endif

      IDX = 0
      IDXF = 0
      IDXFM = 0

    end subroutine flush_buffer

!
! This isn't finished
!
    subroutine open_outfile(filenum, filename)
      integer :: filenum
      character(*) :: filename

      logical :: o

      inquire(filenum, opened=o)
      if (.not.o) then
         open(unit=filenum, file=trim(filename), position=fpow, status=fstat)
      endif

    end subroutine open_outfile

!
! Sets the output files for the next data dump
!
    subroutine set_output_files

    end subroutine set_output_files

    subroutine get_fourier_stats

    end subroutine get_fourier_stats

    subroutine get_realspace_stats

    end subroutine get_realspace_stats

end module OutputType
