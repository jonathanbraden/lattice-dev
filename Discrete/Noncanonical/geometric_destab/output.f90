module output
  use analysis

  contains
    subroutine init_output()
      call open_files()

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

    subroutine open_files()
      if (mpirank == 0) then
         open(unit=98,file="energy_spec.out")
         open(unit=97,file="spectrum.out")
         open(unit=99,file="field_values_spec.out")
      endif
    end subroutine open_files

! This will replace make_output
    subroutine dump_output(time)
      real(dl), intent(in) :: time
      
      integer :: i
      real(dl) :: spec(ns,2*nfld+2)

      call dump_rho(time)
      
#ifdef THREEDIM
      laplace(IRANGE) = fld(1,IRANGE)
      call spectrum_3d(spec(:,1),laplace, Fk, planf)
      Fk2 = Fk
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
    end subroutine dump_output

    subroutine write_output_to_file()
      if (mpirank == 0) then
         do i=1,ns
            write(97,'(30(ES22.15,2x))') (i-1)*dk, spec(i,:)
         enddo
         write(97,*)
      endif
    end subroutine write_output_to_file

end module output
