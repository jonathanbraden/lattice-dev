module homlat

! I need to make the ci's accessible in here !!!

  use mpi
  use modparam
  use model
  use latpar
  use scafld  ! can I get out of doing this?

#include "macros.h"

  implicit none

  real(dl) :: graden, poten, kinen, gradenlap

  real :: G, T, V
  real, dimension(1:fields) :: KINE1, GRADE1

  real, dimension(1:fields) :: FDOTF, FLAPF, FVPRI, FDFD
  real, dimension(1:fields) :: KINE, GRADE, POTE
  real, dimension(1:fields) :: fdot, flap
  real :: ekin, egrad, epot
!  real :: rhoave

  logical :: calc$ge, calc$pe, calc$ke, calc$gelap
  logical :: calc$gecomp, calc$pecomp, calc$kecomp

  contains

!
! The coefficient getting subroutines might be better somewhere else?
!
    subroutine get_energy_coeffs()
      real :: a
      real :: norm

      a = get_scale_factor()
!
! Careful with all of my definitions here, it's easy to screw something up
!
      ekin = 0.5_dl/a**6
      egrad = 0.25_dl/(a*dx)**2/cc
      epot = 0.5_dl

    end subroutine get_energy_coeffs

    subroutine get_encomp_coeffs

    end subroutine get_encomp_coeffs

    subroutine need_get_energies()
      calc$ge = .true.
      calc$pe = .true.
      calc$ke = .true.
      calc$gelap = .true.
    end subroutine need_get_energies

    subroutine need_get_energy_comp()
      calc$gecomp = .true.
      calc$pecomp = .true.
      calc$kecomp = .true.
    end subroutine need_get_energy_comp

    real*8 function grad_en()
      integer :: i,j,k
      real(dl) :: gtmp, GE
      real(dl) :: a, enorm

      if (calc$ge) then

         call wrap_field(fld)
         
         a = get_scale_factor()
         enorm = 0.25_dl/(a*dx)**2/cc
         GE = 0.
         ILOOP
         gtmp = STENCIL(c,GRAD2)
         GE = GE + gtmp
         ENDLOOP
         call MPI_Allreduce(MPI_IN_PLACE, GE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
         GE = enorm*GE/rnvol

         graden = GE
         calc$ge = .false.
      endif

      grad_en = graden
      return
    end function grad_en

!
! Computes gradient energy over the lattice using the Laplacian definition
!
    real*8 function grad_en_wlap()
      integer :: i,j,k
      real(dl) :: gtmp, GE
      real(dl), dimension(1:fields) :: lap
      real(dl) :: a, enorm

      if (calc$gelap) then

         call wrap_field(fld)

         a = get_scale_factor()
         enorm = -0.5_dl/(a*dx)**2/cc
         GE = 0._dl
         ILOOP
         lap(:) = STENCIL(c,LAPLACE)
         lap(:) = fld(:,i,j,k)*lap(:) 
         gtmp = sum(lap(:))
         GE = GE + gtmp
         ENDLOOP
         call MPI_Allreduce(MPI_IN_PLACE, GE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
         GE = enorm*GE/rnvol

         gradenlap = GE
         calc$gelap = .false.
      endif

      grad_en_wlap = gradenlap
    end function grad_en_wlap

    real*8 function pot_en()
      integer :: i,j,k
      real(dl) :: vtmp, PE
      real(dl) :: enorm, a


      if (calc$pe) then
         PE = 0._dl

         a = get_scale_factor()
         enorm = 0.5_dl

         ILOOP
         vtmp = modelv(fld,i,j,k)
         PE = PE + vtmp
         ENDLOOP
         call MPI_Allreduce(MPI_IN_PLACE, PE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
         PE = enorm*PE/rnvol

         poten = PE
         calc$pe = .false.
      endif

         pot_en = poten
         return
    end function pot_en

    real*8 function kinetic_en()
      integer :: i,j,k
      real(dl) :: ttmp, KE
      real(dl) :: enorm, a

      if (calc$ke) then
         KE = 0._dl
         a = get_scale_factor()
         enorm = 0.5_dl/a**6
         
         ILOOP
         ttmp = sum(fldp(:,i,j,k)**2)
         KE = KE + ttmp
         ENDLOOP
         call MPI_Allreduce(MPI_IN_PLACE, KE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)

         KE = enorm*KE/rnvol
         kinen = KE
      endif

      kinetic_en = kinen
      calc$ke = .false.
      return
    end function kinetic_en

    real*8 function grav_en()
      real :: a, H

      grav_en = -3._dl*get_hubble()**2
    end function grav_en

    real*8 function scalar_en
      scalar_en = kinetic_en() + grad_en() + pot_en()
    end function scalar_en

    real*8 function hamiltonian()
      hamiltonian = grad_en() + kinetic_en() + pot_en() + grav_en()
    end function hamiltonian

    subroutine get_gradients_local(i,j,k)
      integer :: i,j,k
      integer :: ifld

      do ifld = 1, fields
         GRADE1(:) = STENCIL(c,GR2)*egrad
      enddo
      G = sum(GRADE1(:))
    end subroutine get_gradients_local

    subroutine get_kinetics_local(ii,jj,kk)
      integer :: ii,jj,kk

      KINE1(:) = fldp(:,ii,jj,kk)**2*ekin
      T = sum(KINE1(:))
    end subroutine get_kinetics_local

    subroutine get_potential_local(ii,jj,kk)
      integer :: ii,jj,kk

      V = modelv(fld,ii,jj,kk)*epot
    end subroutine get_potential_local

    subroutine get_gradients_dt_local(i,j,k)
      integer :: i,j,k

    end subroutine get_gradients_dt_local

    real*8 function get_latave(latfld, st, en)
      integer, dimension(1:3) :: st, en
      real(dl), dimension(st(1):en(1),st(2):en(2),st(3):en(3)) :: latfld

      integer :: i,j,k
      real :: tmpave

      tmpave = 0._dl
      do k=st(3),en(3)
         do j=st(2),en(2)
            do i=st(1),en(1)
               tmpave = tmpave + latfld(i,j,k)
            enddo
         enddo
      enddo

      tmpave = tmpave / dble(nvol)

      call MPI_Allreduce(MPI_IN_PLACE, tmpave, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)

      get_latave = tmpave
      return

    end function get_latave

    real*8 function get_latmom(latfld, st, en, fldave, order)
      integer, dimension(1:3) :: st, en
      real(dl), dimension(st(1):en(1),st(2):en(2),st(3):en(3)) :: latfld
      real(dl) :: fldave  ! average of the field over the lattice
      integer :: order

      integer :: i,j,k
      real :: tmpave

      tmpave = 0._dl
      do k=st(3),en(3)
         do j=st(2),en(2)
            do i=st(1),en(1)
               tmpave = tmpave + (latfld(i,j,k)-fldave)**order
            enddo
         enddo
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, tmpave, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      get_latmom = tmpave/dble(nvol)
      return
    end function get_latmom

end module homlat
