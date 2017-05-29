module output_analysis
#ifdef USEMPI
  use mpi
#endif
  use analysis

  implicit none

  logical, private :: output$cdf
  logical, private :: output$psd
  logical, private :: output$fourierdist
  logical, private :: output$fouriermoms

  integer, private :: psd_filenum, cdf_filenum, fourdist_filenum, fourmom_filenum
  character(80), private :: psd_filename, cdf_filename, fourdist_filename, fourmom_filename
contains

! Edit this to control the output that is generated
  subroutine dump_outputbuffers( v, t, f, st, en)
    character(*) :: v
    real(dl) :: t
    integer, dimension(1:3) :: st, en
#ifdef THREEDIM
    real(C_DOUBLE), pointer
#endif
#ifdef TWODIM

#endif
#ifdef ONEDIM

#endif
  end subroutine dump_outputbuffers

  subroutine set_dump_options(needpsd, needcdf, needfmoms, needfdist)
    logical, optional :: needpsd, needcdf, needfmoms, needfdist

    if (present(needpsd)) output$psd = needpsd
    if (present(needcdf)) output$cdf = needcdf
    if (present(needfmoms)) output$fouriermoms = needfmoms
    if (present(needfdist)) output$fourierdist = needfdist
  end subroutine set_dump_options

  subroutine reset_dump_options()
    output$psd = .false.
    output$cdf = .false.
    output$fouriermoms = .false.
    output$fourierdist = .false.
  end subroutine reset_dump_options
    subroutine set_output_files(psdfile, cdffile, fourdistfile, fourmomfile)
      integer, optional :: psdfile, cdffile
      integer, optional, dimension(1:4) :: fourdistfile
      integer, optional :: fourmomfile
      
      if (present(psdfile)) psdfilenum = psdfile
      if (present(cdffile)) cdffilenum = cdffile
      if (present(fourdistfile)) fourdistfilenum = fourdistfile
      if (present(fourmomfile)) fourmomfilenum = fourmomfile
    end subroutine set_output_files

    subroutine set_output_file_names(psdfile, cdffile, fdistfile, fmomfile)
      character(*),optional :: psdfile, cdffile, fdistfile(1:4), fmomfile

      if (present(psdfile)) psdfilename = psdfile
      if (present(cdffile)) cdffilename = cdffile
      if (present(fdistfile)) fourdistfilename = fdistfile
      if (present(fmomfile)) fourmomfilename = fmomfile
    end subroutine set_output_file_names  

!
! Subroutines that need to be either debugged or written
!
    subroutine get_fourier_stats_3d(f, st, en, famp, fst, fen)
      integer, dimension(1:3) :: st, end, fst, fen
      real, dimension(st(1):en(1),st(2):en(2),st(3):en(3)), intent(in) :: f
      complex, dimension(fst(1):fen(1),fst(2):fen(2),fst(3):fen(3)), intent(in) :: famp

      integer :: i,j,k,ii,jj,kk, k_loc
      real :: p
      integer :: l

      do k=fst(3),fen(3); k_loc=k+z_start_index; if (k_loc > nn3) then; kk=-(n3+1-k_loc); else; kk=k_loc-1; endif
         do j=fst(2),fen(2); if (j<=nn2) then; jj=j-1; else; jj=n2+1-j; endif
            do i=fst(1),fen(1); if (i<=nn1) then; ii=i-1; else; ii=n1+1-i; endif

               p = sqrt(dble(ii**2+jj**2+kk**2)); l=floor(p)
               c = (1.0 - (/l-p,l+1-p/)**2)**2

#ifdef COMPUTE_PSD
               S(l+1:l+2) = S(l+1:l+2) + c*Fk(ii+1,j,k)*conjg(Fk(ii+1,j,k))
               W(l+1:l+2) = W(l+1:l+2) + c
#endif

#ifdef COMPUTE_FOURMOMS
               
#endif

#ifdef COMPUTE_FOURDIST

#endif

            enddo
         enddo
      enddo

    end subroutine get_fourier_stats_3d

    subroutine get_fourier_stats_2d(f, st, en, famp, fst, fen)
      integer, dimension(1:2) :: st, end, fst, fen
      real, dimension(st(1):en(1),st(2):en(2)) :: f
      complex, dimension(fst(1):fen(1),fst(2):fen(2)) :: famp
    end subroutine get_fourier_stats_2d

    subroutine get_fourier_stats_1d(f, st, en, famp, fst, fen)
      integer :: st, en, fst, fen
      real, dimension(st:en) :: f
      complex, dimension(fst:fen) :: famp
    end subroutine get_fourier_stats_1d

    subroutine output_cube(fld)
      real(), dimension(IRANGE) :: fld


    end subroutine output_cube

end module output_analysis
