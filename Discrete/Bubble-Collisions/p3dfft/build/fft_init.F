! This file is part of P3DFFT library
!
! Version 2.2 beta
!
! Copyright (C) 2006-2008 Dmitry Pekurovsky
!
!    P3DFFT is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

!----------------------------------------------------------------------------

c This file contains routines intended for initializing 1D FFT oprations
c as well as one-time initialization and clean-up

c Initialize backward complex-to-complex FFT 

      subroutine plan_b_c1(X,stride_x1,stride_x2,Y,stride_y1,
     &    stride_y2,N,m)

      use p3dfft
      use fft_spec
      implicit none

      integer stride_x1,stride_x2,stride_y1,stride_y2,N,m,fflag
      complex(mytype) X(N*m),Y(N*m)

#ifdef FFTW

#ifdef DOUBLE_PREC
      call dfftw_plan_many_dft(plan1_bc,1,N,m,X,NULL,stride_x1,stride_x2,
     &     Y,NULL,stride_y1,stride_y2, FFTW_BACKWARD,fftw_flag)
#else
      call sfftw_plan_many_dft(plan1_bc,1,N,m,X,NULL,stride_x1,stride_x2,
     &     Y,NULL,stride_y1,stride_y2, FFTW_BACKWARD,fftw_flag)
#endif
#endif
      return
      end

      subroutine plan_b_c2(X,stride_x1,stride_x2,Y,stride_y1,
     &    stride_y2,N,m)

      use p3dfft
      use fft_spec
      implicit none

      integer stride_x1,stride_x2,stride_y1,stride_y2,N,m,fflag
      complex(mytype) X(N*m),Y(N*m)

#ifdef FFTW

#ifdef DOUBLE_PREC
      call dfftw_plan_many_dft(plan2_bc,1,N,m,X,NULL,stride_x1,stride_x2,
     &     Y,NULL,stride_y1,stride_y2, FFTW_BACKWARD,fftw_flag) 
#else
      call sfftw_plan_many_dft(plan2_bc,1,N,m,X,NULL,stride_x1,stride_x2,
     &     Y,NULL,stride_y1,stride_y2, FFTW_BACKWARD,fftw_flag )
#endif
#endif
      return
      end



      subroutine init_b_c(X,stride_x1,stride_x2,Y,stride_y1,
     &    stride_y2,N,m)

      use p3dfft
      use fft_spec
      implicit none


      integer stride_x1,stride_x2,stride_y1,stride_y2,N,m,fflag
      complex(mytype) X(N*m),Y(N*m)


#ifdef ESSL

#ifdef DOUBLE_PREC
      call dcft (1,X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m,-1,1.0d0,
     1        caux1,cnaux,caux2,cnaux)  
#else
      call scft (1,X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m,-1,1.0,
     1        caux1,cnaux,caux2,cnaux)  
#endif

#endif
      return
      end


!!--------------------------------------------------------------
c Initialize backward complex-to-real FFT 

      subroutine plan_b_c2r(X,dimx,Y,dimy,N,m)

      use p3dfft
      use fft_spec
      implicit none

      integer dimx,dimy,N,m,fflag
      complex(mytype) X((N/2+1)*m)
      real(mytype) Y(N*m)

#ifdef FFTW

#ifdef DOUBLE_PREC
      call dfftw_plan_many_dft_c2r(plan1_bcr,1,N,m,
     &        X,NULL,1,dimx, Y,NULL,1,dimy,fftw_flag)
#else
      call sfftw_plan_many_dft_c2r(plan1_bcr,1,N,m,
     &        X,NULL,1,dimx, Y,NULL,1,dimy,fftw_flag)
#endif
#endif

      return
      end

      subroutine init_b_c2r(X,dimx,Y,dimy,N,m)

      use p3dfft
      use fft_spec
      implicit none

      integer dimx,dimy,N,m,fflag
      complex(mytype) X((N/2+1)*m)
      real(mytype) Y(N*m)

#ifdef ESSL

#ifdef DOUBLE_PREC
      call dcrft(1,X,dimx,Y,dimy,N,m,-1,1.0d0,
     &      raux1, rnaux1,raux2,rnaux2,raux3,1)  
#else
      call scrft(1,X,dimx,Y,dimy,N,m,-1,1.0,
     &      raux1, rnaux1,raux2,rnaux2,raux3,1)  
#endif
#endif

      return
      end

!!--------------------------------------------------------------
c Initialize forward complex-to-complex FFT 

      subroutine plan_f_c1(X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m) 
      
      use p3dfft
      use fft_spec
      implicit none

      integer N,m,stride_x1,stride_x2,stride_y1,stride_y2,fflag
      complex(mytype) X(N*m),Y(N*m)

#ifdef FFTW

      fflag = fftw_flag 

#ifdef DOUBLE_PREC
      call dfftw_plan_many_dft(plan1_fc,1,N,m, X,NULL,stride_x1,stride_x2,
     &  Y,NULL,stride_y1,stride_y2,FFTW_FORWARD,fflag)
#else
      call sfftw_plan_many_dft(plan1_fc,1,N,m, X,NULL,stride_x1,stride_x2,
     &  Y,NULL,stride_y1,stride_y2,FFTW_FORWARD,fflag)
#endif
#endif
      return
      end


      subroutine plan_f_c2(X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m) 
      
      use p3dfft
      use fft_spec
      implicit none

      integer N,m,stride_x1,stride_x2,stride_y1,stride_y2,fflag
      complex(mytype) X(N*m),Y(N*m)

#ifdef FFTW

      fflag = fftw_flag

#ifdef DOUBLE_PREC
      call dfftw_plan_many_dft(plan2_fc,1,N,m, X,NULL,stride_x1,stride_x2,
     &  Y,NULL,stride_y1,stride_y2,FFTW_FORWARD,fflag)
#else
      call sfftw_plan_many_dft(plan2_fc,1,N,m, X,NULL,stride_x1,stride_x2,
     &  Y,NULL,stride_y1,stride_y2,FFTW_FORWARD,fflag)
#endif
#endif
      return
      end


      subroutine init_f_c(X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m) 
      
      use p3dfft
      use fft_spec
      implicit none

      integer N,m,stride_x1,stride_x2,stride_y1,stride_y2,fflag
      complex(mytype) X(N*m),Y(N*m)

#ifdef ESSL

#ifdef DOUBLE_PREC
      call dcft(1,X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m,1,1.0d0, 
     1        caux1,cnaux,caux2,cnaux)
#else
      call scft(1,X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m,1,1.0, 
     1        caux1,cnaux,caux2,cnaux)
#endif

#endif

      return
      end


!!--------------------------------------------------------------
c Initialize forward real-to-complex FFT 

      subroutine plan_f_r2c(X,dimx,Y,dimy,N,m)

      use p3dfft
      use fft_spec

      integer dimx,dimy,N,m,fflag
      real(mytype) X(N*m)
      complex(mytype) Y((N/2+1)*m)

#ifdef FFTW

#ifdef DOUBLE_PREC
      call dfftw_plan_many_dft_r2c(plan1_frc,1,N,m,
     &     X,NULL,1,dimx, Y,NULL,1,dimy,fftw_flag)
#else
      call sfftw_plan_many_dft_r2c(plan1_frc,1,N,m,
     &     X,NULL,1,dimx, Y,NULL,1,dimy,fftw_flag)
#endif
#endif

      return
      end

      subroutine init_f_r2c(X,dimx,Y,dimy,N,m)

      use p3dfft
      use fft_spec

      integer dimx,dimy,N,m,fflag
      real(mytype) X(N*m)
      complex(mytype) Y((N/2+1)*m)

#ifdef ESSL

#ifdef DOUBLE_PREC
      call drcft(1,X,dimx,Y,dimy,N,m,1,1.0d0,
     &      raux1, rnaux1,raux2,rnaux2,raux3,1)           
#else
      call srcft(1,X,dimx,Y,dimy,N,m,1,1.0,
     &      raux1, rnaux1,raux2,rnaux2,raux3,1)           
#endif

#endif

      return
      end

!!--------------------------------------------------------------
c One-time initialization
C Initialize work arrays for ESSL

      subroutine init_work(nx,ny,nz)
      
      use fft_spec
      integer nx,ny,nz,err

#ifdef ESSL

#ifdef DOUBLE_PREC
      integer nyz
      nyz = max(ny,nz)
      if(nyz .le. 2048) then
         cnaux = 20000
      else 
         cnaux = 20000+2.28*nyz
      endif
      if(nyz .ge. 252) then
         cnaux = cnaux+(2*nyz+256)*64
      endif
      if(nx .le. 4096) then
         rnaux1 = 22000
         rnaux2 = 20000
      else
         rnaux1 = 20000+1.64*nx
         rnaux2=20000+1.14*nx
      endif
#else
      integer nyz
      nyz = max(ny,nz)
      if(nyz .le. 8192) then
         cnaux = 20000
      else 
         cnaux = 20000+1.14*nyz
      endif
      if(nyz .ge. 252) then
         cnaux = cnaux +(nyz+256)*64
      endif

      if(nx .le. 16384) then
         rnaux1 = 25000
         rnaux2 = 20000
      else
         rnaux1 = 20000+0.82*nx
         rnaux2=20000+0.57*nx
      endif
#endif

      allocate(caux1(cnaux),stat=err)
      if(err .ne. 0) then
         print *,'Error allocating caux1'
      endif
      allocate(caux2(cnaux),stat=err)
      if(err .ne. 0) then
         print *,'Error allocating caux2'
      endif
      allocate(raux1(rnaux1),stat=err)
      if(err .ne. 0) then
         print *,'Error allocating raux1'
      endif
      allocate(raux2(rnaux2),stat=err)
      if(err .ne. 0) then
         print *,'Error allocating raux2'
      endif

#endif

      return
      end

!!--------------------------------------------------------------
c Clean-up routines for FFTW

      subroutine clean_x1

      use fft_spec
#ifdef FFTW
#ifdef DOUBLE_PREC
      call dfftw_destroy_plan(plan1)      
#else
      call sfftw_destroy_plan(plan1)      
#endif
#endif
      
      return
      end

!!--------------------------------------------------------------

      subroutine clean_x2

      use fft_spec

#ifdef FFTW
#ifdef DOUBLE_PREC
      call dfftw_destroy_plan(plan2)      
#else
      call sfftw_destroy_plan(plan2)      
#endif
#endif
      
      return
      end

!!--------------------------------------------------------------

      subroutine clean_x3

#ifdef FFTW
      use fft_spec
#ifdef DOUBLE_PREC
      call dfftw_destroy_plan(plan3)      
#else
      call sfftw_destroy_plan(plan3)      
#endif
#endif
      
      return
      end

!!--------------------------------------------------------------
c Release work arrays for ESSL

      subroutine free_work
      
      use fft_spec

#ifdef ESSL
      deallocate(caux1)      
      deallocate(caux2)
      deallocate(raux1)
      deallocate(raux2)
#endif

      return
      end
