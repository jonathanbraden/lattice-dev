!#define GETFIELDMOMS 1

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

!  integer, parameter :: maxstats=2*fields+15    ! user defined, gives the maximum number of statistics, ideally compute this in some way, maybe through a user-supplied function

  !
  ! Flags determining what is output and what isn't.  This is simply for user convenience
  !
  logical, parameter :: output$bov = .false.

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

#ifdef GETFIELDMOMS
  integer,parameter :: nummomlog = 4
  real(dl) :: fldmom(fields,nummomlog), fldpmom(fields, nummomlog), rhomom(nummomlog)
#endif

  real(dl), dimension(fields) :: phiave, phidotave
  real(dl) :: lnrhoave, eosave, rhoave, prsave

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
       allocate(pp(1:6,IRANGE))
    endif

    if (get$emtensor) then
       allocate(emtens(1:10,IRANGE))
    endif

    if (output$bov .and. ds .ne. 1) then
       allocate(bov(isize(1)/ds, isize(2)/ds,isize(3)/ds))
    endif

    allocate(Fk2(FRANGE))
!    allocate(Fk3(FRANGE))
!    allocate(Fk4(FRANGE))

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
    call head(logfn, [character(len=namelength) :: "t    ", "a    ", "H    ", "<rho>", "<P>  ","<lnrho>","<EOS>","<E_KIN>","E_POT>","<E_grad>","GRAVEN","Ham"])
    curfilenum = curfilenum+1

    logefn = curfilenum
    open(unit=logefn, file="LOG_energy", position=fpos, status=fstat)
    call head(logefn, [character(len=namelength) :: "t    ","a    ","H     ","<rho>","dellnrho","delphi","delphidot"])
    curfilenum = curfilenum+1

    logvfn = curfilenum
    open(unit=logvfn, file="LOG_vir", position=fpos, status=fstat)
    call head(logvfn, [character(len=namelength) :: "t    ","a    ","H     ","<rho>","P"])
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
          call head(trajfn(i), [character(len=namelength) :: "t    ", "a    ", "H    ", "x", "y","z","statistic"])
          curfilenum = curfilenum+1
       enddo
    endif

    if (get$rgtraj) then
       rgtrajfn = curfilenum
       open(unit=rgtrajfn, file="RG_TRAJ", position=fpos, status=fstat)
!       write(rgtrajfn,*) "# tracking", numtraj," lattice points"
       call head(rgtrajfn, [character(len=namelength) :: "t    ", "a    ","H     ","RG_level  ", "rho_tot  "])
       curfilenum = curfilenum+1
    endif

! Now open a file to output particle numbers into
    if (get$npart) then
       partfn = curfilenum
       open(unit=partfn, file="LOG_npart", position=fpos, status=fstat)
       call head(partfn, [character(len=namelength) :: "t     ", "a       ", "H     ", "n_phi", "n_psi", "Sc_phi", "Sq_phi", "Sc_psi", "Sq_psi","m2_phi", "m2_psi"])
       curfilenum = curfilenum + 1
    endif

! Open the entropy file if required
    if (get$entropy) then
       entfn = curfilenum
       open(unit=entfn, file="LOG_S", position=fpos, status=fstat)
       call head(entfn, [character(len=namelength) :: "t    ","a     ","H    ", "<rho>  ", "S_etot","S_etotcheck","S_nohom","S_nohomcheck",&
            "S_nohom-hom","S_nohom-fluc"])
       curfilenum = curfilenum + 1
    endif

    if (get$emtensor) then
       emtfn = curfilenum
       open(unit=emtfn, file="PSD_EM", position=fpos, status=fstat)
       call head(emtfn, [character(len=namelength) :: "t    ","a     ","H     ","k     ", DVAR_EM(1:10)])
       
       curfilenum = curfilenum + 1

       emtpfn = curfilenum
       open(unit=emtpfn, file="CDF_EM", position=fpos, status=fstat)
       call head(emtpfn, [character(len=namelength) :: "t      ","a    ","H     ","Percentile  ", DVAR_EM(1:10)])
       curfilenum = curfilenum+1
      
    endif

    if (get$rgflow) then
        rgfn = curfilenum
	open(unit=rgfn, file="RG_entropy", position=fpos, status=fstat)
	call head(rgfn,[character(len=namelength) :: "t      ","a    ","H     ","entropies"])
        curfilenum = curfilenum + 1

        open(unit=rgfn+1,file="RG_ent", position=fpos, status=fstat)
        call head(rgfn+1,[character(len=namelength) :: "RG Level ","t      ","a    ","rhoave","rhohom ","S_etot    ","sumxi_tot","S_e-hom","sumxi_-hom"])
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
        call head(timefn,[character(len=namelength) :: "t        ","statistic "])
        curfilenum = curfilenum + 1
     endif

  end subroutine open_output_files

  !
  ! User supplied subroutine to create all of the desired output
  !
  subroutine make_output()
    call wrap_field(fld)
    call wrap_field(fldp)
    call get_variables()

    call output_reg()
!    call get_emt()  ! this needs to be done in the get_variables subroutine
!    call output_emt()

!    call get_entot()
!    call output_entot()

!    call get_homogeneous_fields()
    call output_log_file()

!    call output_bov()  ! to speed up output, put this somewhere else

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
    tmp = log(Q*pp(rho,IRANGE))
    call output_cube("lnrho",outframenum, tout, tmp, istart, iend)
  end subroutine output_bov

  subroutine output_reg()
    real :: tout, a, H

    real :: Q
    integer :: i
!
! Initialize the output stuff
!
    call set_output_files(psdfn, cdffn, (/ ampfn, phasefn, refn, imfn /), fouriermomfn )
    call set_output_file_names( "PSD","CDF", [character(len=6) :: "AMPS","PHASES","REDIST","IMDIST"], "FOURMOMS" )
    call set_flush_options(.true.,.true.,.true.,.true.)
    call reset_dump_options()
 
!
! Fix these
!
    tout = time  ! from latpar
    a = get_scale_factor()
    H = get_hubble()
 
! Move these into outputfile.f90 (this is very nonlocal having them here)
    idx = 0
    idxf = 0
    idxfd = 0
    idxcdf = 0

    call set_dump_options(needfmoms=.false., needfdist=.false.)

    tmp = log(pp(rho,IRANGE)) - log(rhoave) - lnrhoave
    call set_dump_options(needfmoms=.false., needfdist=.false.)
    call dump_outputbuffers("logrho_ave",tout,tmp,istart,iend)
    call reset_dump_options()
    Fk2=Fk

    tmp = -3.*H*( pp(prs,IRANGE)/pp(rho,IRANGE) - prsave/rhoave) +   &
         ( fldp(1,IRANGE)*pp(3,IRANGE) + fldp(2,IRANGE)*pp(4,IRANGE) + pp(5,IRANGE) + pp(6,IRANGE) ) / pp(rho,IRANGE) / a**3
    call set_dump_options(needfmoms=.false., needfdist=.false.)
    call dump_outputbuffers("dlnrho",tout,tmp,istart,iend)

    call dump_crossspec("lnr_dlnr",tout,Fk,Fk2)

    call flush_outputbuffers(tout, a, H)
  end subroutine output_reg

  subroutine get_variables()
!
! To make a future version with OpenMP more efficient, define
! some local storage variables
!
    integer :: i,j,k
    real(dl) :: T, G, V
    real(dl) :: KE, GE, PE
    real(dl) :: a
    real(dl) :: kefac, pefac, gefac
    real(dl) :: lap_coeff, piphi_coeff
! Coefficients for for em tensor
    real(dl), dimension(fields) :: par_t, par_x, par_y, par_z
    real(dl) :: emtnorm

    real(dl) :: lnrtmp, eostmp, rhotmp, prstmp

    integer :: l

    call need_get_energies()
    call get_energy_coeffs()

    a = get_scale_factor()

    KE = 0.
    PE = 0.
    GE = 0.

    phiave = 0.
    phidotave = 0.

    lnrtmp = 0.
    eostmp = 0. 
    rhotmp = 0.
    prstmp = 0.

    rhoave = 0.
    prsave = 0.

    gefac = egrad
    pefac = epot
    kefac = ekin

    lap_coeff = 1./(a*dx)**2/cc
    piphi_coeff = 0.5/(a*dx)**2/cc
    
    if (mpirank.eq.0) print*,"ratios in variables laplacian ", lap_coeff/gefac," cross grad ", piphi_coeff/gefac

!OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,T,G,V,l) FIRSTPRIVATE(a,gfac,pfac,kfac,KE,PE,GE) SHARED(pp)
!OMP DO
    ILOOP
      phiave = phiave + fld(:,i,j,k)
      phidotave = phidotave + fldp(:,i,j,k)
    !
    ! Start by getting the various components of energy that we'll need
    ! These are saved by homlat
      call get_gradients_local(i,j,k)
      call get_kinetics_local(i,j,k)
      call get_potential_local(i,j,k)
      call get_gradients_dt_local(i,j,k)

      if (output$encomp) then
         call get_entot_inloop(i,j,k)
      endif
      if (output$set) then

         T = sum(fldp(:,i,j,k)**2)
         G = STENCIL(c,GRAD2)
         V = modelv(fld,i,j,k)

         pp(rho,i,j,k) = kefac*T + gefac*G + pefac*V
         pp(prs,i,j,k) = kefac*T - (1./3.)*gefac*G - pefac*V
         rhoave = rhotmp + pp(rho,i,j,k)
         prsave = prstmp + pp(prs,i,j,k)

         pp(3:4,i,j,k) = lap_coeff* (STENCIL(c,LAPLACE)) 
         pp(5:6,i,j,k) = piphi_coeff* (STENCIL(c,PI_PHI))

         GE = GE + G
         PE = PE + V
         KE = KE + T
         lnrtmp = lnrtmp + log(pp(rho,i,j,k))
         eostmp = eostmp + pp(prs,i,j,k)/pp(rho,i,j,k)
         rhotmp = rhotmp + pp(rho,i,j,k)
         prstmp = prstmp + pp(prs,i,j,k)
      endif
      if (get$emtensor) then
         ! Start by getting all the relevant pieces
         ! Currently, I'm just using a low-order spatial derivative
         par_t(:) = fldp(:,i,j,k)/a**3   ! phi_dot
         par_x(:) = fld(:,i+1,j,k) - fld(:,i-1,j,k)  !grad_x phi * normalization
         par_y(:) = fld(:,i,j+1,k) - fld(:,i,j-1,k)
         par_z(:) = fld(:,i,j,k+1) - fld(:,i,j,k-1)

         emtnorm = pp(rho,i,j,k) + pp(prs,i,j,k)
         emtens(2,i,j,k) = sum(par_t(:)*par_x(:))
         emtens(3,i,j,k) = sum(par_t(:)*par_y(:))
         emtens(4,i,j,k) = sum(par_t(:)*par_z(:))
         emtens(1,i,j,k) = sum(emtens(:,i,j,k))

         emtens(1:4,i,j,k) = emtens(1:4,i,j,k) / emtnorm        

      endif

#ifdef GETFIELDMOMENTS      
      do l=1,nummomlog
         fldmom(:,l) = fldmom(:,l) + fld(:,i,j,k)**l
         fldpmom(:,l) = fldpmom(:,l) + fldp(:,i,j,k)**l
         rhomom(l) = rhomom(l) + log(pp(rho,i,j,k))**2
      enddo
#endif

    ENDLOOP
!OMP END DO
!OMP END PARALLEL
    
    if (output$set) then
       call MPI_allreduce(MPI_IN_PLACE, GE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
       call MPI_allreduce(MPI_IN_PLACE, PE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
       call MPI_allreduce(MPI_IN_PLACE, KE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
       call MPI_allreduce(MPI_IN_PLACE, eostmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
       call MPI_allreduce(MPI_IN_PLACE, lnrtmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
       call MPI_allreduce(MPI_IN_PLACE, rhotmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
       call MPI_allreduce(MPI_IN_PLACe, prstmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)

       GE = egrad * GE / dble(n)**3
       PE = epot * PE / dble(n)**3
       KE = ekin * KE / dble(n)**3
       eostmp = eostmp / dble(n)**3
       lnrtmp = lnrtmp / dble(n)**3
       rhotmp = rhotmp / dble(n)**3
       prstmp = prstmp / dble(n)**3

       call set_energies(kinetic=KE, gradient=GE, potential=PE)
       eosave = eostmp
       lnrhoave = lnrtmp
       rhoave = rhotmp
       prsave = prstmp
    endif

    call MPI_allreduce(MPI_IN_PLACE, phiave, fields, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
    call MPI_allreduce(MPI_IN_PLACE, phidotave, fields, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
    phiave = phiave / dble(n)**3
    phidotave = phidotave / dble(n)**3

! Again, these should be computed in the above loop, rather than after the fact
    rhoave = grad_en() + pot_en() + kinetic_en()
    prsave = kinetic_en() - pot_en() - (1./3.)*grad_en()  
    lnrhoave = lnrhoave - log(rhoave)
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
    integer :: i, j
    real(dl) :: current_avefld
    real(dl) :: a

    a = get_scale_factor()

! I need to precompute stuff here or else it get's stuck in the following loop
!
! To avoid reading through the data multiple times, it's easiest to just compute all of these simultaneously (by calling a single function and doing only 1 data loop)
! (Actually I only need to comput the potential and gradient terms together in this way
    tmp=kinetic_en()
    tmp=grad_en()
    tmp=pot_en()
    tmp=grad_en_wlap()

    if (mpirank.eq.0) then
#ifdef GETFIELDMOMS
       rhomom(1) = rhomom(1) - log(rhoave)
       rhomom(2) = rhomom(2) 
       rhomom(3) = 
       rhomom(4) = 
       write(logefn,'(32g25.16e2)') time, get_scale_factor(), get_hubble(), &
            rhomom(:),  &
            fldmom(1,:), fldmom(2,:),  fldpmom(1,:), fldpmom(2,:) 
#else
       write(logfn,'(32g25.16e2)') time, get_scale_factor(), get_hubble(), scalar_en(), prsave, lnrhoave, eosave, kinetic_en(), pot_en(), grad_en(), grav_en(), hamiltonian(), phiave, phidotave/get_scale_factor()**3
#endif
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
