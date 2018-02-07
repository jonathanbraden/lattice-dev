module output_analysis
#ifdef USEMPI
  use mpi
#endif
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

end module output_analysis
