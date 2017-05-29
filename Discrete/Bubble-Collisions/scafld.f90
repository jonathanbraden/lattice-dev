module scafld

  use mpi
  use p3dfft
  use modparam
  use params
  use model
  use latpar

#include "macros.h"

  implicit none

  real, allocatable, dimension(:,:,:,:) :: fld, fldp
  real :: yscl, ysclp  ! FRW metric  Turn this into a type  (Coformal time, a da/dtau*factors, cosmic time a^1.5 da^(1.5)/dt*factors
  real :: rho_rad, rho_mat

  logical, dimension(fields) :: xboundcond, yboundcond, zboundcond

  contains

    subroutine allocate_fields()
      allocate(fld(1:fields,SIRANGE))
      allocate(fldp(1:fields,SIRANGE))
    end subroutine allocate_fields


    subroutine init_metric()
#ifdef CONFORMAL
      yscl = 1.
      ysclp = -6.*H0*yscl!*mpl2
#else
      yscl = 1.
      ysclp = -4.*H0*yscl!*mpl2
#endif
    end subroutine init_metric

    subroutine init_radiation()
      rho_rad = 0.
    end subroutine init_radiation

    subroutine init_matter()
      rho_mat = 0.
    end subroutine init_matter

#ifdef MINK
    real(dl) function get_scale_factor()
      get_scale_factor = 1.
    end function get_scale_factor

    real(dl) function get_hubble()
      get_hubble = 0.
    end function get_hubble
#else
#ifndef FIXEDBG
    real(dl) function get_scale_factor()
      
#ifdef CONFORMAL
      get_scale_factor = yscl
#else
      get_scale_factor = yscl**(2./3.)
#endif
      return
    end function get_scale_factor

! Returns the value of the physical Hubble parameter (adot/a) (not a'/a for Conformal case)
    real(dl) function get_hubble()
#ifdef CONFORMAL
      get_hubble = -ysclp / yscl / 6._dl / yscl! last piece is to get H not a' !/mpl2
#else
      get_hubble = -ysclp / yscl / 4._dl !/ mpl2
#endif
      return
    end function get_hubble
#endif
#endif


!    subroutine wrap_field(hr)
    subroutine wrap_field_periodic(hr)
      real, dimension(fields,SIRANGE) :: hr

      ! the x dimension is non-distributed        
      hr(:,0,:,:) = hr(:,nside(1),:,:); hr(:,nside(1)+1,:,:) = hr(:,1,:,:)
      
      ! the y dimension may be distributed
      if (ydivs .eq. 1) then
         hr(:,:,0,:) = hr(:,:,nsize(2),:); hr(:,:,nside(2)+1,:) = hr(:,:,1,:)
      else
         call MPI_Sendrecv(hr(:,:,iend(2),:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, rightyrank, 0, &
              hr(:,:,istart(2)-1,:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, leftyrank, 0, &
              MPI_COMM_WORLD, mpistatus, mpierror)
         call MPI_Sendrecv(hr(:,:,istart(2),:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, leftyrank, 0, &
              hr(:,:,iend(2)+1,:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, rightyrank, 0, &
              MPI_COMM_WORLD, mpistatus, mpierror)
      endif
      
      ! the z dimension is distributed
      call MPI_Sendrecv(hr(:,:,:,iend(3)), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, rightzrank, 0, &
           hr(:,:,:,istart(3)-1), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, leftzrank, 0, &
           MPI_COMM_WORLD, mpistatus, mpierror)
      call MPI_Sendrecv(hr(:,:,:,istart(3)), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, leftzrank, 0, &
           hr(:,:,:,iend(3)+1), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, rightzrank, 0, &
           MPI_COMM_WORLD, mpistatus, mpierror)
      
    end subroutine wrap_field_periodic

    subroutine set_boundary_conditions(xbc, ybc, zbc)
      logical, dimension(fields) :: xbc, ybc, zbc

      xboundcond = xbc
      yboundcond = ybc
      zboundcond = zbc
    end subroutine set_boundary_conditions

!
! This subroutine won't work properly if the y-direction is distributed.  This needs to be fixed
!
!    subroutine wrap_field_antiperiodic(hr)
    subroutine wrap_field(hr)
      real, dimension(fields, SIRANGE) :: hr

      integer :: i

! Start with the xboundaries (which aren't distributed
      do i=1,fields
         if (xboundcond(i)) then
            hr(i,0,:,:) = 0. - hr(i,nside(1),:,:)
            hr(i,nside(1)+1,:,:) = 0. - hr(i,1,:,:)
         else
            hr(i,0,:,:) = hr(i,nside(1),:,:)
            hr(i,nside(1)+1,:,:) = hr(i,1,:,:)
         endif
      enddo


      do i=1,fields
         if (ydivs .eq. 1) then
            if (yboundcond(i)) then
               hr(i,:,0,:) = 0. - hr(i,:,nside(2),:)
               hr(i,:,nside(2)+1,:) = 0. - hr(i,:,1,:)
            else
               hr(i,:,0,:) = hr(i,:,nside(2),:)
               hr(i,:,nside(2)+1,:) = hr(i,:,1,:)
            endif

         else
            call MPI_Sendrecv(hr(:,:,iend(2),:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, rightyrank, 0, &
                 hr(:,:,istart(2)-1,:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, leftyrank, 0, &
                 MPI_COMM_WORLD, mpistatus, mpierror)
            call MPI_Sendrecv(hr(:,:,istart(2),:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, leftyrank, 0, &
                 hr(:,:,iend(2)+1,:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, rightyrank, 0, &
                 MPI_COMM_WORLD, mpistatus, mpierror)

            if (yboundcond(i)) then
               if (istart(2).eq.1) then
                  hr(i,:,(istart(2)-1),:) = -hr(i,:,(istart(2)-1),:)
               elseif (iend(2).eq.nside(2)) then
                  hr(i,:,(iend(2)+1),:) = -hr(i,:,(iend(2)+1),:)
               endif
            endif
         endif
      enddo

! Start by sending the data assuming it's periodic.  Then, afterwards multiply the approprate stuff by -1 if it's antiperiodic

      ! the z dimension is distributed
      call MPI_Sendrecv(hr(:,:,:,iend(3)), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, rightzrank, 0, &
           hr(:,:,:,istart(3)-1), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, leftzrank, 0, &
           MPI_COMM_WORLD, mpistatus, mpierror)
      call MPI_Sendrecv(hr(:,:,:,istart(3)), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, leftzrank, 0, &
           hr(:,:,:,iend(3)+1), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, rightzrank, 0, &
           MPI_COMM_WORLD, mpistatus, mpierror)

      ! Now, if we're on the edges of the lattice, apply the antiperiodic b.c.s
      if (istart(3).eq.1) then
         do i=1,fields
            if (zboundcond(i)) then
               hr(i,:,:,(istart(3)-1)) = 0. - hr(i,:,:,(istart(3)-1))
            endif
         enddo
      endif

      if (iend(3).eq.nside(3)) then
         do i=1,fields
            if (zboundcond(i)) then
               hr(i,:,:,(iend(3)+1)) = 0. - hr(i,:,:,(iend(3)+1))
            endif
         enddo
      endif
      
    end subroutine wrap_field


    subroutine wrap2_field(up)
      real, dimension(fields,SIRANGE) :: up

      ! x-direction
      up(:,-1,:,:) = up(:,n-1,:,:)
      up(:,n+2,:,:) = up(:,2,:,:)

      ! ydirection
      if (ydivs .eq. 1) then
         up(:,:,-1,:) = up(:,:,n-1,:)
         up(:,:,n+2,:) = up(:,:,2,:)
      else
         call MPI_Sendrecv(up(:,:,iend(2)-1,:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, rightyrank, 0, &
              up(:,:,istart(2)-2,:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, leftyrank, 0, &
              MPI_COMM_WORLD, mpistatus, mpierror)
         call MPI_Sendrecv(up(:,:,istart(2)+1,:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, leftyrank, 0, &
              up(:,:,iend(2)+2,:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, rightyrank, 0, &
              MPI_COMM_WORLD, mpistatus, mpierror)
      end if

      call MPI_Sendrecv(up(:,:,:,iend(3)-1), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, rightzrank, 0, &
           up(:,:,:,istart(3)-2), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, leftzrank, 0, &
           MPI_COMM_WORLD, mpistatus, mpierror)
      call MPI_Sendrecv(up(:,:,:,istart(3)+1), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, leftzrank, 0, &
           up(:,:,:,iend(3)+2), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, rightzrank, 0, &
           MPI_COMM_WORLD, mpistatus, mpierror)

    end subroutine wrap2_field

!
! Read in and store checkpoint files for the fields themselves
!
    subroutine read_check_pt(v, lat, st, end)
      character(*) :: v
      integer, dimension(1:3) :: st, end
      real :: lat(st(1):end(1),st(2):end(2),st(3):end(3))

      integer, dimension(1:3) :: latsize
      integer :: rawfh
      character(256) :: buffer
      integer :: mpierror, mpistat(MPI_STATUS_SIZE)
      integer :: datsize
      integer(kind=MPI_OFFSET_KIND) :: disp
      real :: tsize


      write(buffer,'(a,a)') trim(v),'.raw'
      call MPI_File_open(MPI_COMM_WORLD, buffer, MPI_MODE_RDONLY, MPI_INFO_NULL, rawfh, mpierror)

      latsize = end - st + 1
      datsize = latsize(1)*latsize(2)*latsize(3)
      disp = mpirank*sizeof(tsize)*datsize

      call MPI_file_set_view(rawfh, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, mpierror)

      call MPI_FILE_READ(rawfh, lat, datsize, MPI_DOUBLE_PRECISION, mpistat, mpierror)

      call MPI_FILE_CLOSE(rawfh, mpierror)

    end subroutine read_check_pt
    
    
    subroutine make_check_pt(v, lat, st, end)
      character(*) :: v
      integer, dimension(1:3) :: st, end
      real :: lat(st(1):end(1),st(2):end(2),st(3):end(3))

      integer, dimension(1:3) :: latsize
      character(256) :: buffer
      integer :: rawfh, mpierror

      integer :: datsize
      integer(kind=MPI_OFFSET_KIND) :: disp

      real :: tsize

      ! make the file name
      write(buffer,'(a,a)') trim(v),'.raw'
      call MPI_File_open(MPI_COMM_WORLD, buffer, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, rawfh, mpierror)

      ! Now set the view on the file before writing
      latsize = end - st + 1
      datsize = latsize(1)*latsize(2)*latsize(3)
      disp = mpirank*sizeof(tsize)*datsize
      call MPI_file_set_view(rawfh, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, mpierror)

      ! Finally write the file
      call MPI_FILE_WRITE(rawfh, lat, datsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, mpierror)

      ! Now close the file
      call MPI_FILE_CLOSE(rawfh, mpierror)

    end subroutine make_check_pt

end module scafld
