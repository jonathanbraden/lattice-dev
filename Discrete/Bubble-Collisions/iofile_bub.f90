module iofile

  use modparam
  use params
  use model
  use latpar
  use scafld
  use homlat
  use analysis
  use OutputFile

#include "macros.h"

  implicit none

!!!!!
  !
  ! Declare some variables that are only used for the output of various user defined quantities in here instead of parameters.f90
  !
!!!!!

!!!!!!!!
  ! USER DEFINED parameters
!!!!!!!

  integer, parameter :: maxstats=2*fields+15    ! user defined, gives the maximum number of statistics, ideally compute this in some way, maybe through a user-supplied function

  !
  ! Flags determining what is output and what isn't.  This is simply for user convenience
  !
  logical, parameter :: output$bov = .true.

  !  logical, parameter :: output$psd = .true.
  !  logical, parameter :: output$cdf = .true.
  !  logical, parameter :: output$pdf = .false.  ! Don't change this, it's broken

  logical, parameter :: output$fld = .true.
  logical, parameter :: output$set = .true.
  logical, parameter :: output$pot = .false.
  logical, parameter :: get$npart = .false.
  logical, parameter :: get$fluc = .false.

  logical, parameter :: output$encomp = .false.
  logical, parameter :: get$emtensor = .false.

  logical, parameter :: get$crosscorr = .false.
  logical, parameter :: datslice = .false.

  logical, parameter :: get$trajectories = .false.
  logical, parameter :: get$entropy = .false.
  logical, parameter :: get$rgtraj = .false.
  logical, parameter :: get$rgflow = .false. .or. get$rgtraj
  logical, parameter :: get$tseries = .false.

  ! Make these better later
  logical, parameter :: output$any = output$fld .or. output$set .or. output$pot .or. output$encomp
  logical, parameter :: output$crv = (output$psd.or.output$cdf) .and. (output$gnu .or. output$vis)

!!!!!!!!
  ! Temporary storage arrays
!!!!!!!!

  real, allocatable, dimension(:,:,:,:) :: rhocomp   ! array to store various decompositions of the energy
  !  real, allocatable, dimension(:,:,:)   :: tmp, tmp2         ! temporary array to store various statistics
  real, allocatable, dimension(:,:,:,:) :: emtens    ! stores components of the emtensor 
  real, allocatable, dimension(:,:,:,:) :: pp        ! homogeneous parts of stress-energy
  integer, parameter :: rho=1,prs=2

!  complex, allocatable :: Fk(:,:,:), Fk2(:,:,:)

!!!!!!!!
  ! Storage arrays for spectra,etc to output
!!!!!!!!
  integer, parameter :: numspecmoms = 4
  character(12), dimension(1:maxstats) :: DVAR
  real, dimension(1:n+1,1:maxstats)    :: CDF
  real, dimension(1:ns,1:maxstats)     :: PSD
  integer, dimension(1:n, 1:maxstats)  :: PAMP
  integer, dimension(1:n, 1:maxstats)  :: PPHASE
  integer, dimension(1:n, 1:maxstats)  :: PREAL, PIMAG
  real, dimension(1:ns, 1:numspecmoms, 1:2) :: SPECMOMS

!!!!!!!!
  ! File numbers for various output files
  ! To do: decide which ones to put in here
!!!!!!!!

  integer :: rgfn, rgdistfn
  integer :: emtpfn
  integer :: trajfn(1:numtstat)
  integer :: logfn, partfn, entfn, emtfn, logefn, logvfn
  integer :: timefn
  integer :: psdfn, cdffn, pdffn
  integer :: rgdifffn, rgratfn
  integer :: rgtrajfn
  integer :: encompfn(3)
  integer :: rgcdffn(numblock), rgrafn(numblock), rgdiffn(numblock)
  integer :: rgmomfn(numrgstat)
  integer :: ampfn, phasefn, refn, imfn
  integer :: fouriermomfn

  integer :: curfilenum

!  character(10) :: fpos, fstat

  contains


    subroutine correct_metric()
      real(dl) :: rho

      rho = grad_en() + kinetic_en() + pot_en()
#ifdef CONFORMAL
      ysclp = -yscl*6.*rho**0.5/sqrt3
#else      
      ysclp = -yscl**(2./3.)*4.*rho**0.5/sqrt3
#endif

    end subroutine correct_metric

  !
  ! Subroutine to open files, etc. for later output
  !
  subroutine init_output()  
    call allocate_output_arrays()
    call open_output_files()
  end subroutine init_output

  subroutine allocate_output_arrays()
    if (output$encomp) then
       allocate(rhocomp(2*fields+2,IRANGE))
    endif

    if (output$set) then
       allocate(pp(1:2,IRANGE))
    endif

    if (get$emtensor) then
       allocate(emtens(IRANGE,1:10))
    endif

    if (output$bov .and. ds .ne. 1) then
       allocate(bov(isize(1)/ds, isize(2)/ds,isize(3)/ds))
    endif

!    allocate(Fk(FRANGE))
    allocate(Fk2(FRANGE))

  end subroutine allocate_output_arrays

  subroutine open_output_files()
    character(100), parameter :: scratchdir='/mnt/scratch-3week/jbraden/'
    integer i
    ! only the master process should be opening stuff
    if (mpirank.ne.0) return

    curfilenum = 14

    ! Determine whether or not new files should be opened or old ones appended to
    if (resume$run) then
       fstat = "UNKNOWN"
       fstat = trim(fstat)
       fpos = "APPEND"
       fpos = trim(fpos)
    else
       fstat = "UNKNOWN"
       fstat = trim(fstat)
       fpos = "REWIND"
       fpos = trim(fpos)
    endif


    ! In order to not have so many file numbers to set in the parameters file, start a loop here, only the master node needs the file numbers anyway
    logfn = curfilenum
    open(unit=logfn, file="LOG.out", position=fpos, status=fstat)
    call head(logfn, (/"t    ", "a    ", "H    ", "<rho>", "<P>  ","delrho^n", "<delphi^n>","<delpsi^n>","<phidot^n>","<psidot^n>","KE","PE","GELAP","GE","H"/))
    curfilenum = curfilenum+1

    logefn = curfilenum
    open(unit=logefn, file="LOG_energy", position=fpos, status=fstat)
    call head(logefn, (/"t    ","a    ","H     ","<rho>","<E_kin>","<E_pot>","<E_grad>","KE_phi","KE_psi","GE_phi","GE_psi","0.5m^2phi^2>"/))
    curfilenum = curfilenum+1

    logvfn = curfilenum
    open(unit=logvfn, file="LOG_vir", position=fpos, status=fstat)
    call head(logvfn, (/"t    ","a    ","H     ","<rho>","P"/))
    curfilenum = curfilenum+1

    if (output$psd) then
       psdfn = curfilenum
       curfilenum = curfilenum + 1
    endif

    if (output$cdf) then
       cdffn = curfilenum
       curfilenum = curfilenum + 1
    endif

    if (output$fourierprob) then
       ampfn = curfilenum
       curfilenum = curfilenum + 1
       phasefn = curfilenum
       curfilenum = curfilenum + 1
       refn = curfilenum
       curfilenum = curfilenum+1
       imfn = curfilenum
       curfilenum = curfilenum + 1
    endif

    if (output$fouriermoms) then
       fouriermomfn = curfilenum
       curfilenum = curfilenum + 1
    endif

!    if (output$pdf) then
!       pdffn = curfilenum
!       do i=1,numstat   ! fix this !!!!  (I'd rather not have to hardcode the # of statistics)
!          curfilenum = curfilenum + 1
!       enddo
!    endif

    if (output$encomp) then
       encompfn(1) = curfilenum
       curfilenum = curfilenum + 1
       encompfn(2) = curfilenum
       curfilenum = curfilenum + 1
       encompfn(3) = curfilenum
       curfilenum = curfilenum + 1
    endif

! Also to do, after I sample the trajectories, print the local indices of the trajectory being followed
    if (get$trajectories) then
       do i=1,numtstat
          trajfn(i) = curfilenum
          open(unit=trajfn(i), file="LOG_"//int2str(i), position=fpos, status=fstat)
          write(trajfn(i),*) "# tracking ", numtraj," lattice points"
          call head(trajfn(i), (/"t    ", "a    ", "H    ", "x", "y","z","statistic"/))
          curfilenum = curfilenum+1
       enddo
    endif

    if (get$rgtraj) then
       rgtrajfn = curfilenum
       open(unit=rgtrajfn, file="RG_TRAJ", position=fpos, status=fstat)
!       write(rgtrajfn,*) "# tracking", numtraj," lattice points"
       call head(rgtrajfn, (/"t    ", "a    ","H     ","RG_level  ", "rho_tot  "/))
       curfilenum = curfilenum+1
    endif

! Now open a file to output particle numbers into
    if (get$npart) then
       partfn = curfilenum
       open(unit=partfn, file="LOG_npart", position=fpos, status=fstat)
       call head(partfn, (/"t     ", "a       ", "H     ", "n_phi", "n_psi", "Sc_phi", "Sq_phi", "Sc_psi", "Sq_psi","m2_phi", "m2_psi"/))
       curfilenum = curfilenum + 1
    endif

! Open the entropy file if required
    if (get$entropy) then
       entfn = curfilenum
       open(unit=entfn, file="LOG_S", position=fpos, status=fstat)
       call head(entfn, (/"t    ","a     ","H    ", "<rho>  ", "S_etot","S_etotcheck","S_nohom","S_nohomcheck",&
            "S_nohom-hom","S_nohom-fluc"/))
       curfilenum = curfilenum + 1
    endif

    if (get$emtensor) then
       emtfn = curfilenum
       open(unit=emtfn, file="PSD_EM", position=fpos, status=fstat)
       call head(emtfn, (/"t    ","a     ","H     ","k     ", DVAR_EM(1:10)/))
       
       curfilenum = curfilenum + 1

       emtpfn = curfilenum
       open(unit=emtpfn, file="CDF_EM", position=fpos, status=fstat)
       call head(emtpfn, (/"t      ","a    ","H     ","Percentile  ", DVAR_EM(1:10)/))
       curfilenum = curfilenum+1
      
    endif

    if (get$rgflow) then
        rgfn = curfilenum
	open(unit=rgfn, file="RG_entropy", position=fpos, status=fstat)
	call head(rgfn,(/"t      ","a    ","H     ","entropies"/))
        curfilenum = curfilenum + 1

        open(unit=rgfn+1,file="RG_ent", position=fpos, status=fstat)
        call head(rgfn+1,(/"RG Level ","t      ","a    ","rhoave","rhohom ","S_etot    ","sumxi_tot","S_e-hom","sumxi_-hom"/))
        curfilenum = curfilenum+1

        ! Open the rg CDF and PDF files, move this elsewhere in the future
        do i=1,numblock
           rgcdffn(i) = curfilenum
           curfilenum = curfilenum + 1
           
           rgrafn(i) = curfilenum
           curfilenum = curfilenum + 1

           rgdiffn(i) = curfilenum
           print*,rgdiffn(i)
           curfilenum = curfilenum + 1
        enddo

        do i=1,numrgstat
           rgmomfn(i) = curfilenum
           curfilenum = curfilenum + 1
        enddo
     endif

     if (get$tseries) then
        timefn = curfilenum

        open(unit=timefn, file="T_SERIES", position=fpos, status=fstat)
        call head(timefn,(/"t        ","statistic "/))
        curfilenum = curfilenum + 1
     endif

  end subroutine open_output_files

  !
  ! User supplied subroutine to create all of the desired output
  !
  subroutine make_output()
    call wrap_field(fld)
    call get_variables()
    call output_reg()

!    call get_emt()
!    call output_emt()

!    call get_entot()
!    call output_entot()

!    call get_homogeneous_fields()
    call output_log_file()

    call output_bov()  ! to speed up output, put this somewhere else

  end subroutine make_output

  subroutine get_emt()

  end subroutine get_emt


!
! To fix, take care of tout, a, H in subroutine
!   Deal with FLDSCL (and appropriate norm on derivatives)
!

#ifdef MINK
#define FLDSCL 1.
#define FLDPSCL 1.
#ifdef CONFORMAL
#undef CONFORMAL
#endif
#else
#ifdef CONFORMAL
#define FLDSCL a
#define FLDPSCL 1/a
#else
#define FLDSCL a**1.5
#define FLDPSCL 1/a**1.5
#endif
#endif

  subroutine output_bov()
    real :: Q
    real :: a, H
    real :: tout

    integer :: i

    if (.not.output$bov) return

    tout = time
    a = get_scale_factor()
    H = get_hubble()

    Q = FLDSCL
    do i=1,fields
       tmp = Q*fld(1,IRANGE)
       call output_cube( "phi"//trim(int2str(i)), outframenum, tout, tmp, istart, iend)
    enddo
    Q = 1._dl/rhoave
!    tmp = log(Q*pp(rho,IRANGE))
    tmp = Q*pp(rho,IRANGE)
    call output_cube("rho_norm",outframenum, tout, tmp, istart, iend)
  end subroutine output_bov

  subroutine output_reg()
    real :: tout, a, H

    real :: Q
    integer :: i
    integer :: idx

!
! Fix these
!
    tout = time  ! from latpar
    a = get_scale_factor()
    H = get_hubble()
 
    idx = 0
    ! Start by outputting the fields
    if (output$fld) then

       Q = 1.0; if (oscale) Q = FLDSCL
       do i=1,fields
          tmp = Q*fld(i,IRANGE)
          call p3dfft_ftran_r2c(tmp, Fk)
          call dumpbuffer("phi_"//trim(int2str(i)), tout, tmp, istart, iend, DVAR, PSD, CDF, idx, .true., PAMP, PPHASE, PREAL, PIMAG)
       enddo

       Q = 1.0; if (oscale) Q = FLDPSCL
       do i=1,fields
          tmp = Q*fldp(i,IRANGE)
          call dumpbuffer("dphi_"//trim(int2str(i)), tout, tmp, istart, iend, DVAR, PSD, CDF, idx, .true.)
       enddo
    endif

    if (get$npart) then
        ! Finishing copying the particle number code here
       if (mpirank.eq.0) then
          Q=1.0; if (oscale) Q= 0. ! check this !!!
          do i=1,fields
             fldspec(:,i) = 10.**PSD(:,i)
             flddspec(:,i) = 10.**PSD(:,fields+i)
          enddo
          call getnumpart(Q,a)   ! do I have a in here?
          
          do i=1,fields
             idx = idx+1
             DVAR(idx) = "npartfix_"//int2str(i)
             PSD(:,idx) = npart(:,i)
          enddo
       endif

    endif

    if (output$set) then
       Q = 1.0; if (oscale) Q = 1./rhoave
       tmp = Q*pp(rho,IRANGE)
       call dumpbuffer("rho_ave",tout,tmp, istart, iend, DVAR, PSD, CDF, idx, .true.)

       tmp = log(pp(rho,IRANGE)) - log(rhoave)
       call dumpbuffer("logrho_ave",tout,tmp,istart,iend,DVAR,PSD,CDF,idx,.true., PAMP, PPHASE, PREAL, PIMAG, SPECMOMS)

       if (get$fluc) then
          tmp = (Q*pp(rho,IRANGE))**2
          call dumpbuffer("rho_fluc", tout, tmp, istart, iend, DVAR, PSD, CDF, idx, .true.)

          tmp = (log(pp(rho,IRANGE)/rhoave))**2
          call dumpbuffer("logrho_fluc", tout, tmp, istart, iend, DVAR, PSD, CDF, idx, .true.)
       endif

       tmp = Q*pp(prs,IRANGE)
       call dumpbuffer("prs", tout, tmp, istart, iend, DVAR, PSD, CDF, idx, .true.)

       tmp = Q*(pp(prs,IRANGE)+pp(rho,IRANGE))
       call dumpbuffer("prs+rho",tout,tmp,istart,iend,DVAR,PSD, CDF,idx,.true.)

       tmp = Q*(3.*pp(prs,IRANGE)-pp(rho,IRANGE))
       call dumpbuffer("3P-rho",tout,tmp,istart,iend,DVAR,PSD,CDF,idx,.true.)

       tmp = pp(prs,IRANGE)/pp(rho,IRANGE)
       call dumpbuffer("EOS",tout,tmp,istart,iend,DVAR,PSD,CDF,idx,.true.)       
    endif

! Fix this value of Q, it's wrong
    if (output$pot) then
#ifdef CONFORMAL
       Q = 0.5
#else
       Q = a**2/2.0   ! check this (it's different in conformal time)
#endif
       tmp = Q*pp(rho,IRANGE)
       call laplace(tmp,tmp)
       call dumpbuffer("PSI",tout,tmp,istart,iend,DVAR,PSD,CDF,idx,.true.)
    endif

    call flushbuffer(tout, a, H, DVAR, PSD, CDF, idx, psdfn, "PSD", cdffn, "CDF", PAMP(1:n,1:idx), ampfn, PPHASE(1:n,1:idx), phasefn, PREAL(1:n,1:idx), refn,PIMAG(1:n,1:idx), imfn, SPECMOMS(1:ns, 1:4, 1:2), fouriermomfn)

  end subroutine output_reg

  subroutine get_variables()
    integer :: i,j,k
    real(dl) :: T, G, V

    call need_get_energies()
    call get_energy_coeffs()

    ILOOP
    !
    ! Start by getting the various components of energy that we'll need
    ! These are saved by homlat
      call get_gradients(i,j,k)
      call get_kinetics(i,j,k)
      call get_potential(i,j,k)
      call get_gradients_dt(i,j,k)

      if (output$encomp) then
         call get_entot_inloop(i,j,k)
      endif
      if (output$set) then
!         call get_pp_inloop(i,j,k)
         T = sum(fldp(:,i,j,k))**2
         G = STENCIL(c,GRAD2)
         V = modelv(fld,i,j,k)

         pp(rho,i,j,k) = ekin*T + egrad*G + epot*V
         pp(prs,i,j,k) = ekin*t - (1./3.)*egrad*G - epot*V
      endif
    ENDLOOP
    
    rhoave = grad_en() + pot_en() + kinetic_en()
    prsave = kinetic_en() - pot_en() - (1./3.)*grad_en()  
  end subroutine get_variables

  subroutine get_pp_inloop(ii,jj,kk)
    integer :: ii,jj,kk
    real(dl) :: T,G,V
    integer :: i,j,k

    i=ii
    j=jj
    k=kk

    call get_energy_coeffs()
    T = ekin*sum(fldp(:,ii,jj,kk)**2)
    G = egrad*STENCIL(c,GRAD2)
    V = epot*modelv(fld,ii,jj,kk)

    pp(rho,ii,jj,kk) = T + G + V
    pp(prs,ii,jj,kk) = T - (1._dl/3._dl)*G - V
  end subroutine get_pp_inloop

  subroutine get_entot_inloop(ii,jj,kk)
    integer :: ii,jj,kk

    rhocomp(1:fields,ii,jj,kk) = KINE1(1:fields)
    rhocomp(fields+1:2*fields,ii,jj,kk) = GRADE1(1:fields)
    rhocomp(2*fields+1,ii,jj,kk) = V
    rhocomp(2*fields+2,ii,jj,kk) = 0.    ! check this one (3 parts in defrost code)
  end subroutine get_entot_inloop

  subroutine output_entot()

  end subroutine output_entot

  subroutine output_emt()

  end subroutine output_emt

  subroutine output_log_file()
    integer, parameter :: nummom = 1
    real(dl), dimension(1:fields, 1:nummom) :: fldmom, fldpmom
    real(dl), dimension(1:nummom) :: rhomom

    integer :: i, j
    real(dl) :: current_avefld
    real(dl) :: a

    a = get_scale_factor()
    ! Start by getting some field variances, etc.
    ! To do: write a separate subroutine to do this
    do i=1,fields
       fldmom(i,1) = FLDSCL*get_latave(fld(i,IRANGE),istart,iend)
       fldpmom(i,1) = FLDPSCL*get_latave(fldp(i,IRANGE), istart, iend)
    enddo
    tmp = log(pp(rho,IRANGE)/scalar_en())
    rhomom(1) = get_latave(tmp,istart,iend)
    current_avefld = rhomom(1)
    do j=2,nummom
       rhomom(j) = get_latmom(tmp, istart, iend, current_avefld, j)
    enddo

    do i=1,fields
       tmp = fld(i,IRANGE) * FLDSCL
       do j=2,nummom
          current_avefld = fldmom(i,1)
          fldmom(i,j) = get_latmom(tmp(IRANGE) ,istart,iend, current_avefld ,j)
       enddo
       tmp = fldp(i,IRANGE) * FLDPSCL
       do j=2,nummom
          current_avefld = fldpmom(i,1)
          fldpmom(i,j) = get_latmom(tmp(IRANGE),istart,iend, current_avefld,j)
       enddo
    enddo

! I need to precompute stuff here or else it get's stuck in the following loop
!
! To avoid reading through the data multiple times, it's easiest to just compute all of these simultaneously (by calling a single function and doing only 1 data loop)
! (Actually I only need to comput the potential and gradient terms together in this way
    tmp=kinetic_en()
    tmp=grad_en()
    tmp=pot_en()
    tmp=grad_en_wlap()

    if (mpirank.eq.0) then
       write(logfn,'(32g)') time, get_scale_factor(), get_hubble(), scalar_en(), prsave,   &
            rhomom(:), fldmom(1,:), fldpmom(1,:)
       write(logefn,'(32g)') time, get_scale_factor(), get_hubble(), scalar_en(), kinetic_en(), pot_en(), grad_en(), grad_en_wlap(), grav_en(), hamiltonian()
    endif

  end subroutine output_log_file

  subroutine output_enlog_file()
!    write(logefn,'(32g)') time, get_scale_factor(), get_hubble(), get_rho(), kinetic_en(), pot_en(), grad_en(), grad_en_wlap(), grav_en(), hamiltonian()
  end subroutine output_enlog_file

  subroutine output_vir_log_file()
!    write(logvfn,'(32g)') time, get_scale_factor(), get_hubble()
  end subroutine output_vir_log_file

  subroutine make_slice()

  end subroutine make_slice

end module iofile
