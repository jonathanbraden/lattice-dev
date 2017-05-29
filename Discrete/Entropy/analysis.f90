module analysis

  use modparam
  use params
  use model
  use mpi; use p3dfft

  implicit none

  contains


!!!!!!!!!!!!!!!!!!!!
!
! Here we'll include several definitions of the entropy to compare
!
!!!!!!!!!!!!!!!!!!!!

!
! Power spectrum obtained by binning
!
    subroutine spectrum_bin(f, S)
      real :: f(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)), S(ns) 
      integer :: W(ns)

      integer :: i,j,k,ii,jj,kk,l
      real*8 :: p

      call p3dfft_ftran_r2c(f, Fk)
      S = 0._dl
      W = 0

      do k=fstart(3),fend(3); if (k<=nn) then; kk=k-1; else; kk=n+1-k; endif
         do j=fstart(2),fend(2); if (j<=nn) then; jj=j-1; else; jj=n+1-j; endif
            do i=1,n; if (i<=nn) then; ii=i-1; else; ii=n+1-i; endif
               if ((ii+1) .lt. fstart(1) .or. (ii+1) .gt. fend(1)) cycle

               p = sqrt(dble(ii**2+jj**2+kk**2));
               l = floor(p) + 1
               
               S(l) = S(l) + Fk(ii+1,j,k)*conjg(Fk(ii+1,j,k))
               W(l) = W(l) + 1
            enddo
         enddo
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, S, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, W, ns, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierror)  

      where (W /= 0) S = S/dble(W)

    end subroutine spectrum_bin

    subroutine spectrum_bin_nofft(famp, S)
      complex, dimension(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) :: famp
      real(dl), dimension(1:ns) :: S
      integer, dimension(1:ns) :: kcount

      integer :: i,j,k,ii,jj,kk
      real(dl) :: p
      integer :: l

      S = 0._dl
      kcount = 0

      do k=fstart(3),fend(3); if (k<=nn) then; kk=k-1; else; kk=n+1-k; endif
         do j=fstart(2),fend(2); if (j<=nn) then; jj=j-1; else; jj=n+1-j; endif
            do i=fstart(1),fend(1); if (i<=nn) then; ii=i-1; else; ii=n+1-i; endif
               if (i .ge. nn) cycle

               p = sqrt(dble(ii**2+jj**2+kk**2));
               l = floor(p) + 1
               
               S(l) = S(l) + Fk(i,j,k)*conjg(Fk(i,j,k))
               kcount(l) = kcount(l) + 1
            enddo
         enddo
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, S, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, kcount, ns, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierror)  

      where (kcount /= 0) S = S/dble(kcount)

    end subroutine spectrum_bin_nofft
    

! one-sided power spectrum density estimator (output valid only on rank-0 node)
    subroutine spectrum(f, S)
        real :: f(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)), S(ns), W(ns)
        integer :: i, j, k, ii, jj, kk, l; real p, c(2)

        call p3dfft_ftran_r2c(f, Fk)
        
        S = 0.0; W = 0.0

        do k = fstart(3),fend(3); if (k <= nn) then; kk = k-1; else; kk = n+1-k; end if
        do j = fstart(2),fend(2); if (j <= nn) then; jj = j-1; else; jj = n+1-j; end if
        do i = 1,n;               if (i <= nn) then; ii = i-1; else; ii = n+1-i; end if
                if ((ii+1) .lt. fstart(1) .or. (ii+1) .gt. fend(1)) cycle

                p = sqrt(real(ii**2 + jj**2 + kk**2)); l = floor(p)
                
                c = (1.0 - (/l-p,l+1-p/)**2)**2
                
                S(l+1:l+2) = S(l+1:l+2) + c * Fk(ii+1,j,k)*conjg(Fk(ii+1,j,k))
                W(l+1:l+2) = W(l+1:l+2) + c
        end do; end do; end do

        call MPI_Allreduce(MPI_IN_PLACE, S, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
        call MPI_Allreduce(MPI_IN_PLACE, W, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
        
        where (W /= 0.0) S = S/W/dble(n)**6  
    end subroutine spectrum

    subroutine spectrum_nofft(famp, S)
      complex, dimension(fstart(1):fend(1), fstart(2):fend(2), fstart(3):fend(3)) :: famp
      real(dl), dimension(1:ns) :: S

      real(dl), dimension(1:ns) :: W
      integer :: l
      integer :: i,j,k,ii,jj,kk
      real(dl) :: p, c(2)


      S = 0.
      W = 0.

      do k = fstart(3),fend(3); if (k<=nn) then; kk=k-1; else; kk=n+1-k;endif
         do j=fstart(2),fend(2); if (j<=nn) then; jj=j-1; else; jj=n+1-j; endif
            do i=fstart(1),fend(1); if (i<=nn) then; ii=i-1; else; ii=n+1-i; endif
               p = sqrt(real(ii**2+jj**2+kk**2))
               l=floor(p)

               c = (1. - (/l-p,l+1-p/)**2)**2

               S(l+1:l+2) = S(l+1:l+2) + c*famp(i,j,k)*conjg(famp(i,j,k))
               W(l+1:l+2) = W(l+1:l+2) + c
            enddo
         enddo
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, S, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, W, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)

      where (W /= 0.) S = S/W/dble(n)**6

    end subroutine spectrum_nofft

    subroutine crossspec_nofft(famp1, famp2, S1, S2)
      complex, dimension(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)) :: famp1, famp2
      real :: W(ns), S1(ns), S2(ns)

      real(dl) :: p
      integer :: l
      complex(dl) :: atmp
      real(dl) :: c(2)
      integer :: i,j,k, ii,jj,kk

      W = 0.
      S1 = 0.
      S2 = 0.

      do k=fstart(3),fend(3); if (k<=nn) then; kk=k-1; else; kk=n+1-k; endif
         do j=fstart(2),fend(2); if (j<=nn) then; jj=j-1; else; jj=n+1-j; endif
            do i=fstart(1),fend(1); if (i<=nn) then; ii=i-1; else; ii = n+1-i ; endif
               if (i >= nn ) cycle

               p = sqrt(dble(ii**2+jj**2+kk**2)); l=floor(p)
               c = (1.0 - (/l-p,l+1-p/)**2)**2
               atmp = famp1(i,j,k)*conjg(famp2(i,j,k))

               S1(l+1:l+2) = S1(l+1:l+2) + c*real(atmp)
               S2(l+1:l+2) = S2(l+1:l+2) + c*aimag(atmp)
               W(l+1:l+2) = W(l+1:l+2) + c

            enddo
         enddo
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, W, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, S1, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, S2, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)

      where (W /= 0.) S1 = S1/W!/dble(n)**6  
      where (W /= 0.) S2 = S2/W!/dble(n)**6  

    end subroutine crossspec_nofft

!
! Compute the cross-spectrum of two variables
! (In computing time, this will be more efficient if I FT and store first instead of many times
! Easiest way to implement this is to include flags for whether or not to FT each variable.
!
!
! July 8: To do: make S a complex variable !!!!
!
    subroutine crossspec(f1, f2, S1, S2)
      real, dimension(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)) :: f1, f2
      real :: W(ns), S1(ns), S2(ns)
      complex*8 :: Stmp(ns)

      integer :: i,j,k,ii,jj,kk,l
      real :: p, c(2)

      call p3dfft_ftran_r2c(f1,Fk)
      call p3dfft_ftran_r2c(f2,Fk2)

      Stmp = 0.0
      W = 0.0

      do k = fstart(3),fend(3); if (k <= nn) then; kk = k-1; else; kk = n+1-k; end if
         do j = fstart(2),fend(2); if (j <= nn) then; jj = j-1; else; jj = n+1-j; end if
            do i = 1,n;               if (i <= nn) then; ii = i-1; else; ii = n+1-i; end if
               if ((ii+1) .lt. fstart(1) .or. (ii+1) .gt. fend(1)) cycle

               p = sqrt(real(ii**2 + jj**2 + kk**2)); l = floor(p)
                
               c = (1.0 - (/l-p,l+1-p/)**2)**2
                
               Stmp(l+1:l+2) = Stmp(l+1:l+2) + 0.5 * c * ( Fk(ii+1,j,k)*conjg(Fk2(ii+1,j,k)) + conjg(Fk(ii+1,j,k))*Fk2(ii+1,j,k) )
               W(l+1:l+2) = W(l+1:l+2) + c
            end do
         end do
      end do

      S1 = real(Stmp)
      S2 = aimag(Stmp)

      call MPI_Allreduce(MPI_IN_PLACE, S1, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, S2, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, W, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)

      where (W/=0.0) S1 = S1/W!/dble(n)**6
      where (W/=0.0) S2 = S2 / W !/dble(n)**6
    end subroutine crossspec
      

!
! Get spectrum and phase spectra
!
    subroutine fullspec(f, S, ang)
      real :: f(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)), S(ns), W(ns)
      complex :: tmpcmp
      real :: ang(ns)  ! this shouldn't be complex
      integer :: i,j,k, ii,jj,kk, l
      real :: p, c(2), powmod

      call p3dfft_ftran_r2c(f, Fk)

      S=0.0; W=0.0; ang=0.0

      ! Now get out the desired power spectrum
      do k=fstart(3),fend(3); if (k<=nn) then; kk=k-1; else; kk=n+1-k; endif
      do j=fstart(2),fend(2); if (j<=nn) then; jj=j-1; else; jj=n+1-j; endif
      do i=1,n;               if (i<=nn) then; ii=i-1; else; ii=n+1-i; endif
         if ((ii+1) .lt. fstart(1) .or. (ii+1) .gt. fend(1)) cycle
         
         p=sqrt(real(ii**2+jj**2+kk**2)); l=floor(p)

         c = (1.0 - (/l-p,l+1-p/)**2)**2

         powmod = Fk(ii+1,j,k)*conjg(Fk(ii+1,j,k))

         S(l+1:l+2) = S(l+1:l+2) + c*powmod
 
         ang(l+1:l+2) = ang(l+1:l+2) + c*Fk(ii+1,j,k)/sqrt(powmod)

         W(l+1:l+2) = W(l+1:l+2) + c
      enddo; enddo;enddo

      call MPI_Allreduce(MPI_IN_PLACE, S, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, W, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, ang, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)

      where (W /= 0.0) S = S/W/real(n)**6
      where (W /= 0.0) ang = ang/W    ! check this

      ! We also want probability distributions, so get these here.
      ! probably best to find dist of abs(f_k) and (theta_k)
      ! The hardest part here is the fact that I'm going to have different wavenumbers on different processors,
      ! making this a pain in the ass

      

    end subroutine fullspec

    subroutine specstats
! Compute vaious correlators for the spectra (such as 4 pt and 3 pt functions, etc
      
      real*8, dimension(1:ns) :: S1, S2

      call MPI_Allreduce(MPI_IN_PLACE, S1, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
      call MPI_Allreduce(MPI_IN_PLACE, S2, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
    end subroutine specstats


    subroutine latave(lat, st, en, ave)
      real, intent(inout) :: ave
      integer, dimension(1:3) :: st, en
      real, dimension(st(1):en(1),st(2):en(2),st(3):en(3)) :: lat
      real :: tmpave
      integer :: i,j,k

      tmpave = 0.
      do k=st(3),en(3); do j=st(2),en(2); do i=st(1),en(1)
         tmpave = tmpave + lat(i,j,k)
      enddo; enddo; enddo
      tmpave = tmpave / real(n)**3

      call MPI_Allreduce(MPI_IN_PLACE, tmpave, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)

      ave = tmpave

    end subroutine latave
      

!
! This will output the trajectories of several different grid points to a file
! This will help determine how much of the behaviour is just a tracking of initial conditions
!
! This subroutine takes a 3-d lattice of variables stored in fld
! and outputs the smoothed version (stored in sfld) with a particular
! choice of storage
! This will be easy to upgrade to a full Haar decomposition and will be done
! in the future
!
! IMPORTANT : This subroutine assumes only a single direction has been distributed by MPI
! This needs to be fixed in the more general case
!
!
subroutine blockrg(fld, sfld, st, en)

!  real, dimension(:,:,:), intent(inout) :: fld, sfld
! fix this to use assumed shape arrays somehow
  integer, dimension(1:3) :: st, en
  real, dimension(st(1):en(1), st(2):en(2), st(3):en(3)) :: fld, sfld

  integer :: i,j,k
  integer :: ii,jj,kk
  integer :: nb, stepb

  real :: tmpval

  integer :: nstart, nend, nzstart, zstep


!
! We want to minimize the computations, so this will be done recursively
! although the subroutine isn't explicitly recursive

!
! Start with the first RG smoothing and store it in sfld
! The temp variable here is used to allow the smoothed field to be stored in the original array in the event we pass fld and sfld as the same array (actually check that this is ok)
!
  do k=istart(3), istart(3)-1+isize(3)/2     !!! use the array indices that have been passed in!!!!!!
     do j=1,n/2
        do i=1,n/2
           tmpval = 0.0
           do kk=0,1; do jj=0,1; do ii=0,1
              tmpval =  tmpval + fld( 2*(i-1)+ii+1, 2*(j-1)+jj+1, istart(3)+2*(k-istart(3))+kk)
           enddo; enddo; enddo
           sfld(i,j,k) = tmpval / 8.
        enddo
     enddo
  enddo

!  print*,"done first loop on ",mpirank

!
! Now do the remaining RG smoothings
! At this moment it requires some externally computed "control parameters" to decide for example how many smoothings, etc.
! are needed before we quit.
!

  nstart = 1
  nend = n/2
  nzstart = 1

  xstart(1) = 1
  xend(1) = n/2
  kend(1) = isize(3)/2

  do nb = 2,numblock

     stepb = 2**nb
     if (stepb.gt.isize(3)) then
        zstep=0
        kend(nb) = 1
     else
        zstep = 1
        kend(nb) = isize(3)/stepb
     endif

     do k=istart(3), istart(3)-1+kend(nb)
        do j=1, n/stepb
           do i=1, n/stepb
              tmpval = 0.0
              do kk=0,zstep; do jj=0,1; do ii=0,1
                 tmpval = tmpval + sfld(nstart+(i-1)*2+ii, nstart+(j-1)*2+jj, istart(3)+(k-istart(3))*2+kk)
              enddo; enddo; enddo
              sfld(i+nend, j+nend, k) = tmpval / 8.
           enddo
        enddo
     enddo

     xstart(nb) = nend+1
     xend(nb) = nend + n/stepb
     nstart = nend + 1
     nend = nend + n/stepb

  enddo


!  if (mpirank.eq.0) print*,"done second loop"

!
! Now, the sums in the x and y directions are complete (since these aren't distributed)
! We now need to finish the sums over the distributed z-directions
!
  do nb = rgstart, numblock
     call MPI_Allreduce(MPI_IN_PLACE, sfld(xstart(nb):xend(nb), xstart(nb):xend(nb), istart(3)), ncount(nb), &
          MPI_DOUBLE_PRECISION, MPI_SUM, rgcomm(nb), mpierror)
  enddo

!  if (mpirank.eq.0) print*,"done allreduce"


end subroutine blockrg

!
! Input: lat : full cube of desired value
!        slat : RG smoothed cube (ie. averaged values, changed on Aug. 13 to use averaged values)
!        st, end : start and end indices for lat and slat
!        mean :: average value (ie. normalization constant) of lattice

subroutine getrgmoms(lat, slat, st, end, mean, mom, getlog)
  integer, dimension(1:3) :: st, end
  real, dimension(st(1):end(1),st(2):end(2),st(3):end(3)) :: lat, slat
  real :: mean
  real, dimension(1:maxmom,0:numblock, 2), intent(inout) :: mom
  logical :: getlog

  integer :: ii, jj, i, j, k, nb

! Do the full cube separately
      mom = 0.0   ! pass in this variable in the future
      do i=istart(1),iend(1); do j=istart(2),iend(2); do k=istart(3),iend(3)
         do ii=1,maxmom
            mom(ii,0,1) = mom(ii,0,1) + lat(i,j,k)**ii
            if (getlog) mom(ii,0,2) = mom(ii,0,2) + ( log(lat(i,j,k)/mean) )**ii
         enddo
      enddo;enddo;enddo
      do ii=1,maxmom
         mom(ii,0,1) = mom(ii,0,1) / mean**ii
      enddo

! Do the RG smoothed cubes now
      do nb = 1, numblock
         do i= xstart(nb),xend(nb); do j=xstart(nb),xend(nb); do k=istart(3),istart(3)-1+kend(nb)
            do ii=1,maxmom
               mom(ii,nb,1) = mom(ii,nb,1) + slat(i,j,k)**ii
!               if (getlog) mom(ii,nb,2) = mom(ii,nb,2) + ( log(tmp(i,j,k)/(8.**nb*mean)) )**ii
               if (getlog) mom(ii,nb,2) = mom(ii,nb,2) + ( log(slat(i,j,k)/mean) )**ii
            enddo
         enddo;enddo;enddo
         do ii=1,maxmom
!            mom(ii,nb,1) = mom(ii,nb,1) / 8.**(ii*nb) / mean**ii
            mom(ii,nb,1) = mom(ii,nb,1) * 8.**nb  ! do this so everything is normalized properly later when I divide by n**3?)
            mom(ii,nb,2) = mom(ii,nb,2) * 8.**nb
            mom(ii,nb,1) = mom(ii,nb,1) / mean**ii

         enddo
         if (nb .ge. rgstart) then   ! does this need to be fixed? Aug. 13? isize(3)/2**nb = # of processors that need to be averaged? I think it stays
            do ii=1,maxmom; do jj=1,2
               mom(ii,nb,jj) = mom(ii,nb,jj) * isize(3) / 2.**nb
            enddo; enddo
         endif

      enddo

      ! Now combine the results off of each processor
      call MPI_Allreduce(MPI_IN_PLACE, mom, maxmom*(numblock+1)*2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)

      mom = mom / n**3

end subroutine getrgmoms


subroutine getlatmoms(lat, st, end, norm, nmom, mom)
  integer, dimension(1:3) :: st, end
  real, dimension(st(1):end(1),st(2):end(2),st(3):end(3)) :: lat
  real :: norm   ! I could just compute this in here
  integer :: nmom
  real, dimension(1:nmom,2), intent(inout) :: mom

  integer :: ii,jj,i,j,k, nb
  real :: latave

  mom=0.0
  
  ! first get the lattice average
  do k=st(3),end(3)
     do j=st(2),end(2)
        do i=st(1),end(1)
           mom(1,1) = mom(1,1) + lat(i,j,k)
        enddo
     enddo
  enddo
  
  call MPI_Allreduce(MPI_IN_PLACE, mom(1,1), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
  mom(1,1) = mom(1,1) / n**3
  latave = mom(1,1)

  do nb=2,nmom
     do k=st(3),end(3)
        do j=st(2),end(2)
           do i=st(1),end(1)
              mom(nb,1) = mom(nb,1) + lat(i,j,k)**nb
              mom(nb,2) = mom(nb,2) + (lat(i,j,k)-latave)**nb
           enddo
        enddo
     enddo
     mom(nb,1) = mom(nb,1) / norm**nb
     mom(nb,2) = mom(nb,2) / norm**nb
  enddo

  call MPI_Allreduce(MPI_IN_PLACE, mom(2:nmom,1:2), (nmom-1)*2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)

  mom = mom / n**3


end subroutine getlatmoms

! Maps between the indices for the smoothed field (with 1,2,3,... in new coords) and the original
! box (stored in tmp), for a given level of the RG flow
!
    subroutine rgindex(ir,jr,kr,rg1, rg2)
      integer, intent(in) :: rg1, rg2
      integer, intent(inout) :: ir, jr, kr
      integer :: stepdiff

      if (rg2.lt.rg1) return

      stepdiff = 2**(rg2-rg1)
      ir = ir/stepdiff + xstart(rg2)
      jr = jr/stepdiff + xstart(rg2)

      if (rg2.ge.rgstart) then
         kr = istart(3)
      else
         kr = kr/stepdiff + istart(3)  ! check this one
      endif  

    end subroutine rgindex

    subroutine win_init
      
      integer :: i, j, k, l, ii, jj, kk
      real :: norm, scl
      integer :: size
      
      integer :: nwin

      nwin = 1

      ! start by allocating space for all of the FT window functions
!      allocate(wind(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3),1:nwin))
      
      
!      do l = 1,nwin
         ! First compute the real space version of the window
         ! In here we put the desired real space window function
!         do k=istart(3),iend(3); if (k <= nn) then kk=k-1; else; kk=n+1-k; endif
!            do j=istart(2), iend(2); if (j<=nn) then jj=j-1; else; jj=n+1-j; endif
!               do i=istart(1), iend(1); if (i<=nn) then ii=i-1; else; ii=n+1-i; endif
!                  scl = 5.
!                  tmp(i,j,k) = exp(-(ii**2+jj**2+kk**2)/(2*scl))
!               enddo
!            enddo
!         enddo
         
!         norm = sum(tmp)
!         call MPI_Allreduce(MPI_IN_PLACE, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
         
!         tmp = tmp / norm
         
         ! Now compute the Fourier space version by FFTing
!         call p3dfft_ftran_r2c(tmp, wind(:,:,:,l))
         
!      enddo
      

    end subroutine win_init
    
    subroutine win_smooth(fld, st, end, wnum)
      
      integer, dimension(3) :: st, end
      real, dimension(st(1):end(1), st(2):end(2), st(3):end(3)), intent(inout) :: fld
      integer :: wnum
      
      integer :: i,j,k
      
      ! begin by transforming desired field
      call p3dfft_ftran_r2c(fld, Fk)
      
      ! Now multiply by the desired window
      ! Here, we are using the particular storage structure of FFTW
      !
      ! Make sure this is correct
      do k=fstart(3), fend(3)
         do j=fstart(2), fend(2)
            do i=fstart(1), fend(1)           
!               Fk(i,j,k) = Fk(i,j,k) * wind(i,j,k, wnum)
            enddo
         enddo
      enddo
  
      ! Now FT backwards to get the desired field (perhaps don't store it in original field?
!      fld = p3dfft_btran_c2r(Fk, fld) 


    end subroutine win_smooth

!
! Convolves the array fld with the window function win.
! The st,en give start/end locations of the arrays.
! fs,fe are the start end wavenumbers for the FFT'd arrays
!
    subroutine convolve(fld, win, st, en, fs, fe)
      integer, dimension(1:3), intent(in) :: st, en, fs, fe
      real, dimension(st(1):en(1),st(2):en(2),st(3):en(3)), intent(inout):: fld
      real, dimension(st(1):en(1),st(2):en(2),st(3):en(3)), intent(inout) :: win
!      complex*8, dimension(fs(1):fe(1),fs(2):fe(2),fs(3):fe(3)), intent(out) :: fld_fft, win_fft    

      real :: sp(ns)

      integer :: i, j, k

      call p3dfft_ftran_r2c(fld, Fk)
      call p3dfft_ftran_r2c(win, Fk2)
    
      ! Perform the Convolution
      do i=fs(1),fe(1)
         do j=fs(2),fe(2)
            do k=fs(3),fe(3)
               Fk(i,j,k) = Fk(i,j,k)*Fk2(i,j,k)
            enddo
         enddo
      enddo


      ! Now invert the transformation
      call p3dfft_btran_c2r(Fk, fld)
      call p3dfft_btran_c2r(Fk2, win)

      fld = fld/real(n**3)

  end subroutine convolve


! This subroutine assumes that the window has already been fft's and thus doesn't need to perform the extra fft
  subroutine convolve_fft(fld, fld_fft, win_fft,st, en, fs, fe)
    integer, dimension(1:3) :: st, en, fs, fe
    real, dimension(st(1):en(1),st(2):en(2),st(3):en(3)),intent(inout) :: fld
    complex, dimension(fs(1):fe(1),fs(2):fe(2),fs(3):fe(3)) :: fld_fft
    complex, dimension(fs(1):fe(1),fs(2):fe(2),fs(3):fe(3)), intent(in) :: win_fft
    
    ! transform the field
    call p3dfft_ftran_r2c(fld, fld_fft)
!    call p3dfft_ftran_r2c(fld, Fk)
    ! Now do the multiplication
    fld_fft = fld_fft*conjg(win_fft)
    call p3dfft_btran_c2r(fld_fft, fld)
    fld = fld / dble(n)**6
    
  end subroutine convolve_fft

!
! This subroutine initializes the window kernel that I want
! Currently, the only option is a gaussian.  This can be extended later
!
! I'm currently assuming all lattices are the same size in sim (since I use nn in here)
!
  subroutine make_kernel_fft(sclfac, win, win_fft, st, en, fs, fe)
    integer, dimension(1:3), intent(in) :: st, en, fs, fe
    real :: sclfac
!    real, dimension(st(1):en(1),st(2):en(2),st(3):en(3)) :: win
    real, dimension(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)) :: win
    complex, dimension(fs(1):fe(1),fs(2):fe(2),fs(3):fe(3)), intent(inout) :: win_fft

    integer :: i,j,k
    real :: rad, wnorm
    
    do k=st(3),en(3); do j=st(2),en(2); do i=st(1),en(1)
        rad = sqrt( real(i-nn)**2 + real(j-nn)**2 + real(k-nn)**2 )   ! change if window size changes
        win(i,j,k) = exp(-rad**2/sclfac**2)
     enddo; enddo; enddo
     
     wnorm = sum(win)
     call MPI_ALLREDUCE(MPI_IN_PLACE, wnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
     
     win = win / wnorm
     
     call p3dfft_ftran_r2c(win, win_fft)     
   end subroutine make_kernel_fft

   subroutine make_kernel_real(sclfac, win, st, en)
     real :: sclfac
     integer, dimension(1:3) :: st, en
     real, dimension(st(1):en(1),st(2):en(2),st(3):en(3)), intent(out) :: win

     integer :: i,j,k, ii,jj,kk
     real :: rad, wnorm
    
     wnorm =0.
     do k=st(3),en(3)
        kk=k-1
        do j=st(2),en(2)
           jj=j-1
           do i=st(1),en(1)
              ii=i-1
              
              rad = sqrt( real(ii)**2 + real(jj)**2 + real(kk)**2 )*dx   ! change if window size changes
              !             rad = sqrt( real(i-1)**2 + real(j-1)**2 + real(k-1)**2 )*dx
              win(i,j,k) = exp(-rad**2/sclfac**2)
              wnorm = wnorm + win(i,j,k)
              
           enddo
        enddo
     enddo
    
     call MPI_ALLREDUCE(MPI_IN_PLACE, wnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
     
     win = win / wnorm   

   end subroutine make_kernel_real

  function int2str(I)
    integer, intent(in)::I
    character(LEN=16) INT2STR
    write(INT2STR,*) I
    INT2STR = Trim(ADJUSTL(INT2STR))
  end function int2str


end module analysis
