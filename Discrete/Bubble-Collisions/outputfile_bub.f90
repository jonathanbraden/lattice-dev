module OutputFile

  use mpi
  use params
  use p3dfft
  use model
  use modparam
  use analysis  ! to do: split this up into analysis for the output and analysis in general

  implicit none

  logical, parameter :: output$cdf = .false.
  logical, parameter :: output$psd = .false.
  logical, parameter :: output$pdf = .false.
  logical, parameter :: output$fourierprob = .false.
  logical, parameter :: output$fouriermoms = .false.

  character(10) :: fpos, fstat
  real(dl), allocatable, dimension(:,:,:) :: bov

  integer, parameter :: numkbins = ns
!  integer, parameter :: maxstats = 2*fields+15
  integer, parameter :: maxfstats = fields + 1

!
! This type provides a common area in which to store various statistics that we wish to output
! This has been implemented to ease the readability of user output routines (and also to permit easier integration
! of new postprocessing into the code) 
!
  logical :: get$PSD
  logical :: get$fouriermoms
  logical :: get$fourierdist
     
  integer :: numbins
  integer :: nummoms

  integer :: idx, idxf
!# Add some preprocessors here so I'm not creating extra arrays unless I need them
! real(dl), dimension(1:numkbins, 1:nummoms, 1:2, 1:maxfstats) :: PSD_MOMS
! real(dl), dimension(1:numkbins, 1:numbins, 1:maxfstats) :: PSD_DIST
!  real(dl), dimension(1:ns, 1:maxstats) :: PSD
!  real(dl), dimension(1:n+1, 1:maxstats) :: CDF
!  character(12), dimension(1:maxstats) :: DVAR

  contains

!
! To do's: Combine a lot of the individual subroutines into a single subroutine that only iterates through a single loop (c.f computing spectra, moments and distributions)
!


!
! spectrumangle  ! get variance in power spectrum angles
!

!
! The tricky thing here is that I'm actually going to be doing analysis with the data, so I need the arrays to analyze
! Thus, putting all of the decision structure in here will be tricky.
!
    subroutine outputdata
      

    end subroutine outputdata


!
! Required subroutines
!   - open files
!   - decide what to output (currently in step)
!   - combine all of my stupid different dump subroutines into 1


!
! It will be better to put this all in a structure designed for output and then do it in a more object oriented way
!

!
! Modify this so various things go in different files
!
! To do : make PSP, CDFN, VNAME optional arguments
!
! Input for subroutine
!    v : Variable name (a string)
!    t : time (put in whatever you want to store, such as cosmic, conformal, or ln(a)
!    f : the array with the field values to store data from
!    st : starting indices for f
!    en : ending indices for array f
!    VNAME : array to store variable names in
!    PSP : array to store computed power spectrum in
!    CDFN : array to store computed CDF in
!    idx : integer labelling the various variables being stored (different one for each final final to output) - stored dynamically to reduce user input
!    dolog : (boolean) store PSD as log10 or not
!
    subroutine dumpbuffer(v, t, f, st, en, VNAME, PSP, CDFN, count, dolog, AMPBIN, PHASEBIN, REBIN, IMBIN, FMBINS)
      character(*) :: v
      real(dl) :: t
      integer, dimension(1:3) :: st, en 
      real(dl), dimension(st(1):en(1),st(2):en(2),st(3):en(3)) :: f
      character(*), intent(inout), dimension(:) :: VNAME
      real(dl), dimension(:,:), intent(inout) :: PSP
      real(dl), dimension(:,:), intent(inout) :: CDFN
!      real, dimension(1:ns,:), intent(inout) :: PSP
!      real, dimension(1:n+1,:), intent(inout) :: CDFN
      integer,optional, intent(inout) :: count
      logical :: dolog
      integer, optional, dimension(:,:) :: AMPBIN, PHASEBIN, REBIN, IMBIN
      real(dl), optional, dimension(:,:,:) :: FMBINS

      real(dl), dimension(ns) :: S
      real(dl), dimension(n+1) :: C

      integer k

      if (present(count)) then
         count = count+1
         VNAME(count) = v
      endif

! get the power spectrum
      if (output$psd) then
!         call spectrum(f, S)
         call spectrum_bin(f,S)
         if (dolog) PSP(:,count) = log10(S(:))
      endif

      if (output$fourierprob) then
         if (present(AMPBIN) .and. present(PHASEBIN)) then
            call bin_fourier_prewhiten(f, st, en, AMPBIN(:,count),PHASEBIN(:,count),REBIN(:,count), IMBIN(:,count),n)
         endif
      endif

      if (output$fouriermoms) then
         if (present(FMBINS) ) then
            call fourier_moments(f, st, en, FMBINS(1:ns, 1:4, 1:2), 4)
         endif
      endif

      if (output$cdf) then
         call getcdf(f, tmp2, C, size(f,1)*size(f,2)*size(f,3))  ! get rid of reference to tmp2, find a better way to store temp arrays
         CDFN(:,count) = C(:)             
      endif

    end subroutine dumpbuffer

!
! From the point of view of my code, many of the input parameters above are pretty redundant, so here it might be better
! to simply include a simplified dump subroutine, so the calls are shorter
!
! In particular, if several sets of files are being written (say N), simply write N output subroutines, one for each set of files.
! It's probably actually best to write this in a user specified-module, since then it separates out what needs to be user-supplied.
!


!
! Flush the data obtained from dump into a file
!
! This needs some work in order to also be able to dump the RG'd stuff
!
    subroutine flushbuffer(t, scl, hub, VNAMES, PSD, CDF, nstat, psdfilen, psdname, cdffilen, cdfname, AMPS, ampfilen, PHASES, phasefilen, RPART, refilen, IPART, imfilen, FOURMOMS, fourmomfn) 
      real :: t, scl, hub
      character(*), dimension(1:nstat), intent(in) :: VNAMES
      real, dimension(1:ns ,1:nstat), intent(in)  :: PSD
      real, dimension(1:n+1,1:nstat), intent(in)  :: CDF
      integer, intent(in) :: nstat
      integer :: psdfilen, cdffilen
      character(*) :: psdname, cdfname
      integer, optional, dimension(1:n,1:nstat) :: AMPS, PHASES, RPART, IPART
      integer, optional :: ampfilen, phasefilen, refilen, imfilen
      real, optional, dimension(1:ns,4,2) :: FOURMOMS
      integer, optional :: fourmomfn 

      integer :: k,i
      logical :: o

      if (mpirank.ne.0) return

      if (output$psd) then
         inquire(psdfilen, opened=o)
         if (.not.o) then
            open(unit=psdfilen, file=trim(psdname), position=fpos, status=fstat)
            call head(psdfilen, (/"t   ","a   ","H    ","k    ", VNAMES(1:nstat)/))  ! fix this line
         endif

         do k=1,ns
            write(psdfilen,'(32g25.16e2)') t, scl, hub, (k-1)*dk, PSD(k,1:nstat)
         enddo
         write(psdfilen,'(g)') "",""  ! change this as needed
         flush(psdfilen)
      endif

      if (output$cdf) then
         inquire(cdffilen, opened=o)
         if (.not.o) then
            open(unit=cdffilen, file=trim(cdfname), position=fpos, status=fstat)
            call head(cdffilen, (/"t    ","percentile ",VNAMES(1:nstat)/))
         endif

         do k=1,n+1
            write(cdffilen,'(32g25.16e2)') t, real(k-1)/n, CDF(k,1:nstat)
         enddo
         write(cdffilen,'(g)') "",""
         flush(cdffilen)
      endif

      if (present(AMPS)) then
      if (output$fourierprob) then
         inquire(ampfilen ,opened=o)
         if (.not.o) then
            open(unit=ampfilen, file="AMPS.dat", position=fpos, status=fstat)
            call head(ampfilen, (/"t     ","amp/power",VNAMES(1:nstat)/))
         endif

         inquire(phasefilen ,opened=o)
         if (.not.o) then
            open(unit=phasefilen, file='PHASES.dat', position=fpos, status=fstat)
            call head(phasefilen, (/"t    ","Angle",VNAMES(1:nstat)/))
         endif

         inquire(refilen, opened=o)
         if (.not.o) then
            open(unit=refilen, file='REALAMP.dat', position=fpos, status=fstat)
            call head(refilen, (/"t     ","amp/pow**0.5",VNAMES(1:nstat)/))
         endif

         inquire(imfilen, opened=o)
         if (.not.o) then
            open(unit=imfilen, file='IMAMP.dat', position=fpos, status=fstat)
            call head(imfilen, (/"t     ","amp/pow**0.5",VNAMES(1:nstat)/))
         endif

         do k=1,n
            write(ampfilen,'(32g)') t, (real(k-1)/3.), AMPS(k,1:nstat)
            write(phasefilen,'(32g)') t, real(k-1)*twopi/dble(n)-twopi/2._dl, PHASES(k,1:nstat)
            write(refilen,'(32g)') t, real(k-n/2)/1.5, RPART(k,1:nstat)
            write(imfilen,'(32g)') t, real(k-n/2)/1.5, IPART(k,1:nstat)
         enddo
         write(ampfilen,'(g)') "",""
         write(phasefilen,'(g)') "",""
         write(refilen,'(g)') "",""
         write(imfilen,'(g)') "",""
      endif
      endif

      if (present(FOURMOMS) .and. present(fourmomfn) ) then
      if (output$fouriermoms) then
         inquire(fourmomfn, opened=o)
         if (.not.o) then
            open(unit=fourmomfn, file="FOURIER_MOMS", position=fpos, status=fstat)
            call head(fourmomfn, (/"t    ","k    ", "<Re f> ","<IM  f>","<|f|^2> ","Re f^3   ","Im f^3   ","|f|^4  "/) )
         endif
         if (mpirank.eq.0) then
            do k=1,ns
               write(fourmomfn, '(32g)') t, (k-0.5)*dk, FOURMOMS(k, 1,:), FOURMOMS(k, 2, 1) &
                    , FOURMOMS(k, 3, :), FOURMOMS(k, 4, 1)
            enddo
            write(fourmomfn, '(g)') "", ""
         endif

      endif   
      endif

    end subroutine flushbuffer

!
! Flush buffer for the RG flow analysis
!
! Parameter :: fld - fld to perform smoothing to
!              sfld - array to store smoothed field
!              rfld - ????
!
    subroutine dump_rgbuffer(fld, sfld, rfld, st, en, count, t, vname, dolog)
      character(12) :: vname
      integer, dimension(1:3) :: st, en
      real, dimension(st(1):en(1),st(2):en(2),st(3):en(3)) :: fld, sfld, rfld
      integer, intent(inout) :: count
      real :: t
      logical :: dolog

      count = count + 1
!      VNAMES(count) = trim(vname)

      ! Start by performing the smoothing
      call blockrg(fld(st(1):en(1),st(2):en(2),st(3):en(3)), sfld(st(1):en(1),st(2):en(2),st(3):en(3)), st, en)

      ! Now get the desired statistics


    end subroutine dump_rgbuffer

    subroutine flush_rgbbuffer

    end subroutine flush_rgbbuffer

!
! To reduce the number of FFT required and also reduce the number of times I'm going through loops, gather all the operations that require working in Fourier space here
!
! This can probably be simplified considerably by defining a new type to store all of the arrays I have, then just passing it in and out
! I can even include logical parameters in the type to decide which things I need to compute when I call this subroutine
! This is probably the best thing to do from a usability standpoint (but maybe not speed)
!
    subroutine get_fourier_stats(f, st, en, famp, fst, fen, S_weight, four_moms, nummoms, amps, phases, recomp, imcomp, nbins)
      integer, dimension(1:3) :: st, en, fst, fen
      real*8, dimension(st(1):en(1), st(2):en(2), st(3):en(3)) :: f
      complex, dimension(fst(1):fen(1),fst(2):fen(2),fst(3):fen(3)) :: famp
      real*8, dimension(1:ns) :: S_weight
      integer :: nummoms
      real, dimension(1:ns,1:nummoms,1:2), optional :: four_moms
      real, dimension(1:ns+1, 1:nbins), optional :: amps, phases, recomp, imcomp

      integer :: i,j,k, ii,jj,kk
      integer :: kint, l, m
      real*8 :: keff

      real*8, dimension(1:ns) :: S_bin, W
      real*8, dimension(1:2) :: c
      integer, dimension(1:ns) :: nkcount
      integer :: numk, nbins

      real*8 :: ik, rk, amp2
      real*8, parameter :: sigsize = 3.

      numk = ns
! Start by getting the Fourier spectrum
      call p3dfft_ftran_r2c(f, famp)

! Since I need the spectrum for normalizing stuff, compute it here (unnecessary if I only bin the Fourier modes)
      call spectrum_bin_nofft(famp, S_bin)   ! write this subroutine 

! Initialize the various statistics we want to store
      S_weight = 0.
      W = 0.
      four_moms = 0.
      nkcount = 0
      amps = 0
      phases = 0
      recomp = 0
      imcomp = 0

! Now successively compute the desired statistics.  Eventually I'll want to write some _in_loop subroutines to make it look nicer
      do k=fst(3),fen(3); if (k <= nn) then; kk = k-1; else; kk = n + 1 - k; endif
         do j=fst(2),fen(2); if (j <= nn) then; jj=j-1; else; jj=n+1-j; endif
            do i=fst(1),fen(1); if (i <= nn) then; ii=i-1; else; ii=n+1-i; endif

               keff = sqrt(dble(ii**2+jj**2+kk**2))
               kint = floor(keff) + 1
               l = kint - 1
               nkcount(kint) = nkcount(kint) + 1
               amp2 = famp(i,j,k)*conjg(famp(i,j,k))

!# preprocessor for getting PSD
               ! Start by getting volume weighted Fourier spectrum (put in a preprocessor somewhere)
               c = (1.0 - (/ dble(l)-keff , dble(l+1)-keff/)**2)**2
               W(kint:kint+1) = W(kint:kint+1) + c
               S_weight(kint:kint+1) = S_weight(kint:kint+1) + c*amp2

!# preprocessor for getting spectral moments
               ! Now get the moments of the spectrum
               rk = real(famp(i,j,k))
               ik = aimag(famp(i,j,k))
               do m=1,nummoms, 2
                  four_moms(kint,m ,1) = four_moms(kint,m,1) + rk**m
                  four_moms(kint,m,2) = four_moms(kint,m,2) + ik**m
                  four_moms(kint,m+1, 1) = four_moms(kint,m+1,1) + amp2**((m+1)/2)
               enddo

!# preprocessor for getting binned distributions
               rk = real(famp(i,j,k))
               ik = aimag(famp(i,j,k))
               ! Normalize the amplitudes to the power
               amp2 = amp2 / S_bin(kint)
               l = floor(amp2*dble(nbins)/sigsize) + 1
               if (l.gt.nbins) l=nbins   ! make sure we don't overshoot the bins
               amps(kint,l) = amps(kint,l) + 1

               ! Now study the phases
               amp2 = atan2(ik, rk)
               l = floor(amp2*dble(nbins)/twopi) + nbins/2 + 1
               if (l.gt.nbins) print*, "Error in assigning Fourier phases"; l=nbins
               phases(kint,l) = phases(kint,l) + 1

               ! Finally get the rescaled real and imaginary parts
               ik = ik / sqrt(amp2)
               rk = rk / sqrt(amp2)

               l = floor(rk*dble(nbins)/(2.*sigsize)) + 1 + nbins/2
               if (l.gt.nbins) l=nbins
               if (l.lt.1) l = 1
               recomp(kint,l) = recomp(kint,l) + 1

               l = floor(ik*dble(nbins)/(2.*sigsize)) + 1 + nbins/2
               if (l.gt.nbins) l=nbins
               if (l.lt.1) l = 1
               imcomp(kint,l) = imcomp(kint,l) + 1

            enddo
         enddo
      enddo

! Now sum up the totals over all processors
! # Add preprocessors corresponding to the above (as required)
      call MPI_reduce(MPI_IN_PLACE, nkcount, numk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror)

!# Add spectrum preprocessor)
      call MPI_reduce(MPI_IN_PLACE, W, numk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
      call MPI_reduce(MPI_IN_PLACE, S_weight, numk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierror)
      if (mpirank.eq.0) then
         where( W /= 0. ) S_weight = S_weight / W / dble(n)**6
      endif

!# preprocessor for spectral moments
      call MPI_reduce(MPI_IN_PLACE, four_moms, numk*nummoms*2, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror)
      if (mpirank.eq.0) then
         do i=1,numk
            if ( nkcount(i) .ne. 0) then
               four_moms(i,:,:) = four_moms(i,:,:) / dble(nkcount(i))
            endif
         enddo
      endif

!# preprocessor for spectral distributions (binned)
      call MPI_reduce(MPI_IN_PLACE, amps, numk*nbins, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror)
      call MPI_reduce(MPI_IN_PLACE, phases, numk*nbins, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror)
      call MPI_reduce(MPI_IN_PLACE, recomp, numk*nbins, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror)
      call MPI_reduce(MPI_IN_PLACE, imcomp, numk*nbins, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierror)

      if (mpirank.eq.0) then
         
      endif

    end subroutine get_fourier_stats

!
! Similar to the Fourier thing, but to prevent doing multiple loops over real space
!
! Ideally, I could simply do a single big loop over the entire cube and as I'm going along compute most of the desired statistics (I could for example bin, but I can't get the CDF)
!
    subroutine get_realspace_stats

    end subroutine get_realspace_stats


    subroutine bin_fourier_prewhiten_in_loop

    end subroutine bin_fourier_prewhiten_in_loop

!
! Get binned pdf of Fourier amplitudes and phases averaged over all Fourier modes.  Here, we prewhiten, so the amplitudes are always in (-1,1) and phases in (0,2pi), so no need to do anything fancy like find mins and maxes
!
    subroutine bin_fourier_prewhiten(f, st, en, amps, phases, recomp, imcomp, nbins)
      integer, dimension(1:3) :: st, en
      real, dimension(st(1):en(1),st(2):en(2),st(3):en(3)) :: f
!      integer, dimension(1:3) :: fs, fe
!      complex, dimension(fs(1):fe(1),fs(2):fe(2),fs(3):fe(3)) :: famp
      real*8, dimension(1:ns) :: spec
      integer :: nbins
      integer, dimension(1:nbins), intent(out) :: amps, phases, recomp, imcomp
 
      integer :: ii,jj,kk, i, j, k
      real*8 :: keff, tmpamp, tmpphase
      integer :: kint, l

      amps=0
      phases=0
      recomp=0
      imcomp=0

      call spectrum_bin(f,spec)   ! check how much this changes if I don't do a binning of the spectrum
      call p3dfft_ftran_r2c(f,Fk)

      do k=fstart(3), fend(3); if (k <= nn) then; kk=k-1; else; kk=n+1-k; endif
         do j=fstart(2), fend(2); if (j <= nn) then; jj=j-1; else; jj=n+1-j; endif
            do i=1, n; if (i <= nn) then; ii=i-1; else; ii=n+1-i;  endif
               if ((ii+1) .lt. fstart(1) .or. (ii+1) .gt. fend(1)) cycle

!               keff = get_keff(ii,jj,kk)
!               kinteger = get_kinteger(ii,jj,kk)
               keff = sqrt(dble(ii**2+jj**2+kk**2))
               kint = floor(keff) + 1

               tmpamp = Fk(ii+1,j,k)*conjg(Fk(ii+1,j,k))
               tmpamp = tmpamp / spec(kint)
!               print*,"tmpamp",tmpamp
               l = floor(tmpamp*dble(nbins)/3.) +1  ! I need to tune this parameter (currently largest bin is 3sigma (probably too big))
!               print*,"l = ",l, kint
               if (l.gt.nbins) l=nbins
               amps(l) = amps(l) + 1

               tmpphase =  atan2(aimag(Fk(ii+1,j,k)), real(Fk(ii+1,j,k))) 
               l = floor(tmpphase*dble(nbins)/twopi) + nbins/2 + 1
!               if (l.lt.1 .or. l.gt.nbins) print*,"error",l, tmpphase, kint, i,j,k, Fk(ii+1,j,k), fend
               if (l.gt.nbins) l = nbins
               phases(l) = phases(l) + 1

               tmpamp=real(Fk(ii+1,j,k))/spec(kint)**0.5
               l = floor(tmpamp*dble(nbins)/6.) + 1 + nbins/2
               if (l.gt.nbins) l=nbins
               if (l.lt.1) l=1
               recomp(l)=recomp(l)+1

               tmpamp = aimag(Fk(ii+1,j,k))/spec(kint)**0.5
               l = floor(tmpamp*dble(nbins)/6.) + 1 + nbins/2
               if (l.gt.nbins) l=nbins
               if (l.lt.1) l=1
               imcomp(l) = imcomp(l)+1

            enddo
         enddo
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, amps, nbins, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, phases, nbins, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, recomp, nbins, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, imcomp, nbins, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierror)

    end subroutine bin_fourier_prewhiten

!
! Same as the above subroutine, except we also bin into slices of |k|
!
! Finish this!!!!
!
    subroutine bin_fourier_slice_prewhiten(famp, fs, fe, spectrum, amps, phases, kbins, nbins)
      integer, dimension(1:3) :: fs, fe
      complex, dimension(fs(1):fe(1),fs(2):fe(2),fs(3):fe(3)) :: famp
      real*8, dimension(1:ns) :: spectrum
      integer :: kbins, nbins
      integer, dimension(kbins, nbins) :: phases, amps

      integer :: ii, jj, kk, i, j, k
      integer :: l, kint
      real*8 :: keff, tmpamp, tmpphase

      phases = 0
      amps = 0
      
      do kk=fstart(3),fend(3)
         do jj=fstart(2),fend(2)
            do ii=fstart(1),fend(1)
               keff = sqrt(dble(ii**2+jj**2+kk**2))
               kint = floor(keff) + 1

               tmpamp = famp(ii+1,j,k)*conjg(famp(ii+1,j,k))/spectrum(kint)
               l = floor(tmpamp*dble(nbins)/6.)
               amps(kint,l) = amps(kint, l)+1

               phases(kint, l) = phases(kint, l) + 1
            enddo
         enddo
      enddo
      
      
      call MPI_Allreduce(MPI_IN_PLACE, amps, nbins*kbins, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, phases, nbins*kbins, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierror)

      
    end subroutine bin_fourier_slice_prewhiten

!
! Get the kurtosis, etc. for the Fourier coeffs (the spectrum is just the variance, although good to check)
!

    subroutine fourier_moments(f, st, en, spec_moms, mombins, famp)
      integer, dimension(1:3) :: st, en
      real*8, dimension(st(1):en(1),st(2):en(2),st(3):en(3)) :: f
      integer :: mombins
      real*8, dimension(1:ns, 1:mombins, 1:2) :: spec_moms
      complex, optional, dimension(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) :: famp

      real*8 :: keff
      integer :: kint, l

      integer, dimension(1:ns) :: ncount
      integer :: kbins

      real*8 :: rk, ik
      integer :: i,j,k,ii,jj,kk

! Fix this, it won't work if I pass in Fk
      if (.not. present(famp)) then
         call p3dfft_ftran_r2c(f, Fk)
      endif

      kbins = ns

      ncount = 0
      spec_moms = 0.

      do k=fstart(3), fend(3); if (k <= nn) then; kk=k-1; else; kk=n+1-k; endif
         do j=fstart(2), fend(2); if (j <= nn) then; jj=j-1; else; jj=n+1-j; endif
            do i=1, n; if (i <= nn) then; ii=i-1; else; ii=n+1-i;  endif
               if ((ii+1) .lt. fstart(1) .or. (ii+1) .gt. fend(1)) cycle
! Get the desired Fourier mode
               keff = sqrt(dble(ii**2+jj**2+kk**2))
               kint = floor(keff) + 1
               ncount(kint) = ncount(kint) + 1

               if (kint.gt.ns) then
                  if (mpirank.eq.0) print*,"kint out of bounds"
               endif

               rk = real(Fk(ii+1,j,k))
               ik = aimag(Fk(ii+1,j,k))

               do l=1,mombins, 2
                  spec_moms(kint, l, 1) = spec_moms(kint, l, 1) + rk**l
                  spec_moms(kint, l, 2) = spec_moms(kint, l, 2) + ik**l
               enddo
               do l=2,mombins, 2
                  spec_moms(kint, l, 1) = spec_moms(kint, l, 1) + Fk(ii+1,j,k)*conjg(Fk(ii+1,j,k))**(l/2)
               enddo
            enddo
         enddo
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, ncount, kbins, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, spec_moms, kbins*mombins*2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      ! This won't work since spec_moms and count have different sizes !!!

      do i =1,ns
         if (ncount(i) .ne. 0) then
            spec_moms(i,:,:) = spec_moms(i,:,:)/dble(ncount(i))
         endif
      enddo

    end subroutine fourier_moments
    
!
! Compute the generating function for moments of the 1-point distribution E[e^{tX}]
! How do I tell a function to return an array of unknown size? 
!
    subroutine generating_function(f, st, en, slims, ssteps)
      integer, dimension(1:3) :: st, en
      real*8, dimension(st(1):en(1),st(2):en(2), st(3):en(3)) :: f
      integer :: ssteps
      real, dimension(1:2) :: slims

      real*8 :: s, ds, tmpf
      integer :: i, j, k
      integer :: l
      real*8, dimension(1:ssteps) :: genfunc

      ds = (slims(2)-slims(1)) / dble(ssteps)
      genfunc = 0.

      do k=st(3),en(3)
         do j=st(2),en(2)
            do i=st(1),en(1)
               tmpf = f(i,j,k)
               do l=1,ssteps
                  s = ds*dble(l)
                  genfunc(l) = genfunc(l) + exp(s*tmpf)
               enddo
            enddo
         enddo
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, genfunc, ssteps, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      genfunc = genfunc / dble(n**3)   ! can make this more independent of knowledge of the cube (at expense of being slower?) 

    end subroutine generating_function


!
! Write out the data cube into a BOV file to read using VisIt
!
! Input: v - name of variable
!        frame - frame number to label data with
!        t - time to use as an index
!        fld - lattice of field values to output
!        st - starting indices of lattice
!        en - ending indices of lattice
!        bovfn file number of use for the BOV file (will change to SILO eventually) (removed)
!
    subroutine output_cube( v, frame, t, lat, st, en)  !, bovfn)
      character(*) :: v
      integer :: frame   ! labels identifying the time,etc. for the file
      real :: t
      integer, dimension(1:3) :: st, en
      real(dl), dimension(st(1):en(1),st(2):en(2),st(3):en(3)), intent(in) :: lat

      integer :: bovfn


      character(256) :: buffer
      integer :: db

      integer, dimension(1:3) :: latsize

      integer :: rawfh   ! stores handle for MPI file used to output the box
      integer :: datsize  ! number of data elements to store on each processor

      integer :: i,j, k,ii, jj, kk  ! loop variables

      integer(kind=MPI_OFFSET_KIND) :: disp
      real(dl) :: tsize  ! used to set file view by giving the size of a real number

      bovfn = 99  ! open on unit 99

      latsize(:) = en(:)-st(:)+1
      datsize = latsize(1)*latsize(2)*latsize(3)/ds**3
      disp = mpirank*sizeof(lat(istart(1),istart(2),istart(3)))*datsize

      write(buffer,'(a,a,a,i4.4,a)') 'output/',trim(v),'-',frame,'.raw'
      call MPI_File_open(MPI_COMM_WORLD, buffer, MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, rawfh, mpierror)

      call MPI_file_set_view(rawfh, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, mpierror)
      if (ds .eq. 1) then
         call MPI_file_write(rawfh, lat(st(1):en(1),st(2):en(2),st(3):en(3)), datsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, mpierror)
      else
         do k = st(3),en(3), ds
            kk = (k - st(3))/ds + 1   ! where am I computing k?
            do j = st(2),en(2), ds
               jj = (j - st(2))/ds + 1
               do i = st(1),en(1), ds
                  ii = (i-st(1))/ds + 1
                  bov(ii,jj,kk) = lat(i,j,k)
               enddo
            enddo
         enddo
         call MPI_File_write(rawfh, bov, datsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, mpierror)
      endif

      call MPI_file_close(rawfh, mpierror)


      ! Open BOV header file (see below for writing)
      if (mpirank .eq. 0) then
         write (buffer,'(a,a,a,i4.4,a)') 'output/',trim(v),'-',frame,'.bov'
         open(bovfn,file=buffer)
         write(buffer,'(a,a,i4.4,a)') trim(v),'-',frame,'.raw'
!      endif
      ! Finally, write the BOV header file
!      if (mpirank.eq.0) then
         write(bovfn,'(a,g25.16e2)')  'TIME:        ', t
         write(bovfn,'(2a)')          'VARIABLE:    ', v
         write(bovfn,'(a,3g25.16e2)') 'BRICK_ORIGIN:', 0.0, 0.0, 0.0
         write(bovfn,'(a,3g25.16e2)') 'BRICK_SIZE:  ', n*dx, n*dx, n*dx
         write(bovfn,'(2a)')          'DATA_FILE:   ', trim(buffer)
         write(bovfn,'(a,3i5)')       'DATA_SIZE:   ', n/ds, n/ds, n/ds
         write(bovfn,'(2a)')          'DATA_FORMAT: ', "DOUBLE"
         write(bovfn,'(2a)')          'DATA_ENDIAN: ', "LITTLE"   ! bug check this, it depends on the processor
         write(bovfn,'(2a)')          'CENTERING:   ', "zonal"    ! change if needed
!         write(bovfn,'(a,i5)')        'BYTE_OFFSET: ', 4    ! check if this is needed (since it's fortran) (don't need or visit screws up
         close(bovfn)
      endif


    end subroutine output_cube

    subroutine output_slice( v, frame, t, slice, st, en)  !, bovfn)
      character(*) :: v
      integer :: frame   ! labels identifying the time,etc. for the file
      real :: t
      integer, dimension(1:2) :: st, en
      real(dl), dimension(st(1):en(1),st(2):en(2)), intent(in) :: slice

      integer :: bovfn


      character(256) :: buffer
      integer :: db

      integer, dimension(1:2) :: latsize

      integer :: rawfh   ! stores handle for MPI file used to output the box
      integer :: datsize  ! number of data elements to store on each processor

      integer :: i,j, k,ii, jj, kk  ! loop variables

      integer(kind=MPI_OFFSET_KIND) :: disp
      real(dl) :: tsize  ! used to set file view by giving the size of a real number

      bovfn = 99  ! open on unit 99

      latsize(:) = en(:)-st(:)+1
      datsize = latsize(1)*latsize(2)
      disp = mpirank*sizeof(slice(st(1),st(2)))*datsize

      write(buffer,'(a,a,a,i4.4,a)') 'output/',trim(v),'-',frame,'.raw'
      call MPI_File_open(MPI_COMM_WORLD, buffer, MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, rawfh, mpierror)

      call MPI_file_set_view(rawfh, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, mpierror)
      call MPI_file_write(rawfh, slice(st(1):en(1),st(2):en(2)), datsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, mpierror)

      call MPI_file_close(rawfh, mpierror)


      ! Open BOV header file (see below for writing)
      if (mpirank .eq. 0) then
         write (buffer,'(a,a,a,i4.4,a)') 'output/',trim(v),'-',frame,'.bov'
         open(bovfn,file=buffer)
         write(buffer,'(a,a,i4.4,a)') trim(v),'-',frame,'.raw'
!      endif
      ! Finally, write the BOV header file
!      if (mpirank.eq.0) then
         write(bovfn,'(a,g25.16e2)')  'TIME:        ', t
         write(bovfn,'(2a)')          'VARIABLE:    ', v
         write(bovfn,'(a,3g25.16e2)') 'BRICK_ORIGIN:', 0.0, 0.0
         write(bovfn,'(a,3g25.16e2)') 'BRICK_SIZE:  ', n*dx, n*dx
         write(bovfn,'(2a)')          'DATA_FILE:   ', trim(buffer)
         write(bovfn,'(a,3i5)')       'DATA_SIZE:   ', latsize(1), latsize(2)
         write(bovfn,'(2a)')          'DATA_FORMAT: ', "DOUBLE"
         write(bovfn,'(2a)')          'DATA_ENDIAN: ', "LITTLE"   ! bug check this, it depends on the processor
         write(bovfn,'(2a)')          'CENTERING:   ', "zonal"    ! change if needed
!         write(bovfn,'(a,i5)')        'BYTE_OFFSET: ', 4    ! check if this is needed (since it's fortran) (don't need or visit screws up
         close(bovfn)
      endif


    end subroutine output_slice


    subroutine output_rg_cube()

    end subroutine output_rg_cube

    subroutine output_smooth_cube(v,t,frame)
      character(*) ::v
      real*8 :: t
      integer :: frame 

      integer :: bovfn
      character(256) :: buffer

      bovfn=99

      if (mpirank.eq.0) then
         write (buffer,'(a,a,a,i4.4,a)') 'output/',trim(v),'-',frame,'.bov'
         open(bovfn,file=buffer)
         write(buffer,'(a,a,i4.4,a)') trim(v),'-',frame,'.raw'
                                   
                                                                                                                
                                                                                                                                                   
         write(bovfn,'(a,g25.16e2)')  'TIME:        ', t
         write(bovfn,'(2a)')          'VARIABLE:    ', v
         write(bovfn,'(a,3g25.16e2)') 'BRICK_ORIGIN:', 0.0, 0.0, 0.0
         write(bovfn,'(a,3g25.16e2)') 'BRICK_SIZE:  ', n*dx, n*dx, n*dx
         write(bovfn,'(2a)')          'DATA_FILE:   ', trim(buffer)
         write(bovfn,'(a,3i5)')       'DATA_SIZE:   ', n/ds, n/ds, n/ds
         write(bovfn,'(2a)')          'DATA_FORMAT: ', "DOUBLE"
         write(bovfn,'(2a)')          'DATA_ENDIAN: ', "LITTLE"   ! bug check this, it depends on the processor                                                                  
         write(bovfn,'(2a)')          'CENTERING:   ', "zonal"    ! change if needed                                                                                             
!         write(bovfn,'(a,i5)')        'BYTE_OFFSET: ', 4    ! check if this is needed (since it's fortran) (don't need or visit screws up                                       
         close(bovfn)
      endif

    end subroutine output_smooth_cube


!!!!!
!
! Here store some subroutines temporarily.  I really want these in the analysis module
!
! These are auxilliary routines for sorting, making cdfs, etc.
!
!!!!!

    subroutine mymerge(f, ft, fs)
      integer :: fs
      real, dimension(1:fs), intent(inout) :: f
      real, dimension(1:fs) :: ft
      
      call mergesort(f,ft)
      
    end subroutine mymerge
    
! local merge sort
    subroutine mergesort(f, ft)
      real, dimension(:), intent(inout) :: f
      real, dimension(1:size(f)), intent(inout) :: ft
      
      !        print*,"In mergesort, rank = ",mpirank
      
      call mergesortinner(f, ft)
        !        print*,"Done mergesort, rank = ",mpirank
    end subroutine mergesort
      
    recursive subroutine mergesortinner(f, ft)
      real, dimension(:), intent(inout) :: f
      real, dimension(1:size(f)), intent(inout) :: ft
      
      integer fs, mid, cur, left, right 
      
      !        print*,"In innter mergesort, rank = ",mpirank
      
      fs = size(f)
      if (fs .lt. 2) return
      
      mid = fs/2
      call mergesortinner(f(1:mid), ft(1:mid))
      call mergesortinner(f((mid+1):fs), ft((mid+1):fs))
      
      cur = 1; left = 1
      right = mid + 1
      
      do while (cur .le. fs)
         if (left .gt. mid) then
            ft(cur) = f(right)
            right = right + 1
         else if (right .gt. fs) then
            ft(cur) = f(left)
            left = left + 1
         else
            if (f(left) .le. f(right)) then
               ft(cur) = f(left)
               left = left + 1
            else
               ft(cur) = f(right)
               right = right + 1
            end if
         end if
         
         cur = cur + 1
      end do
      
      f = ft
    end subroutine mergesortinner


    !  Simple subroutine to estimate the CDF.  Needs to be improved
!  Currently, I sort each cube, then find partitions for each one and average them
!  At some point, I should really do another sort or something
!
    subroutine getcdf(f, ft, pivots, asize)
      integer :: asize
      real, dimension(1:asize), intent(inout) :: f
      real, dimension(1:size(f)), intent(inout) :: ft
      real, dimension(n+1), intent(inout) :: pivots

      integer :: dsize, i

!      if (mpirank.eq.0) print*,"Before sort ",f
      ! start by sorting the cube on each process
      call mergesort(f, ft)
!      if (mpirank.eq.0) print*,"After sort ",f

      ! Now get the median values for the CDF
      dsize = size(f) / n
      do i=1,n
         pivots(i) = f((i-1)*dsize+1)
      enddo
      pivots(n+1) = f(asize)

      ! Now average over all the processes, and get min/max for the 0% and 100% ones
      call MPI_Allreduce(MPI_IN_PLACE, pivots(2:n), n-1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, pivots(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE,pivots(n+1), 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierror)
      ! As a test, I should compute say std. dev. of the pivots, works as a test of ergodicity as well

      do i=2,n
         pivots(i) = pivots(i) / mpisize
      enddo

!      if (mpirank.eq.0) print*,"pivots are ",pivots

    end subroutine getcdf

    subroutine make_bov_header(v, frame)
      character(*) :: v
      integer :: frame

      if (mpirank.eq.0) then

      endif

    end subroutine make_bov_header


  ! identify yourself
  subroutine head(fd, vars)
        integer(4) fd; character(*) :: vars(:)
        character(512) :: buffer; integer a, b, c, l
        character(*), parameter :: rev = "$Revision: 2.0 $"

        if (mpirank .ne. 0) return
        if (resume$run) return
        
        ! behold the horror that is Fortran string parsing
        a = index(rev, ": ") + 2
        b = index(rev, ".")
        c = index(rev, " $", .true.) - 1
        buffer = rev(a:b-1); read (buffer, '(i12)') a
        buffer = rev(b+1:c); read (buffer, '(i12)') b
        
        ! ID string
        write (fd,'(a,4(i0,a),4(g,a))') "# This is DEFROST revision ", a-1, ".", b, " (", fields, " fields, ", n, "^3 grid) L=",len," dt = ",dt, " mpl = ",mpl," seed = ",seednum

#ifdef _OPENMP
        write (fd,'(a,i0)') "# Using OpenMP with max threads: ", OMP_GET_MAX_THREADS()
#endif

        write (fd,'(a,i0,a,i0)') "# MPI rank: ", mpirank, "; total nodes: ", mpisize
        write (fd,'(a,i0,a,i0,a)') "# Using ", ydivs, " y division(s) and ", zdivs, " z division(s)"
        write (fd,'(a,i0,a,i0)') "# Local y division: ", lydiv, "; local z division: ", lzdiv
        write (fd,'(a,i0,a,i0)') "# neighbor ranks: ", leftyrank, " <- y -> ", rightyrank
        write (fd,'(a,i0,a,i0)') "# neighbor ranks: ", leftzrank, " <- z -> ", rightzrank

        ! model summary
        call modelsummary(fd)
        
        ! variable list
        write (buffer,'(a,g12.12",",32(g24.12","))') "OUTPUT:", adjustr(vars); l = index(buffer, ',', .true.);
        write (fd,'(2a)') "# ", repeat('=', l)
        write (fd,'(2a)') "# ", buffer(1:l-1)//';'
        write (fd,'(2a)') "# ", repeat('=', l)
      end subroutine head

end module OutputFile
