!#define MINK 1
!#define FIXEDBG
!#define CONFORMAL 1

#define IRANGE istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)
#define SIRANGE padst(1):paden(1),padst(2):paden(2),padst(3):paden(3)
#define FRANGE fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)

#define ILOOP do k=istart(3),iend(3); do j=istart(2),iend(2); do i=istart(1),iend(1)
#define SILOOP do k=padst(3),paden(3); do j=padst(2),paden(2); do i=padst(1),paden(1)
#define ENDLOOP enddo; enddo; enddo

! I need to add in the adjustments to the values of the effective wavenumbers
! and I need to test/fix this preprocessor
#define FLOOP (F1LOOP(k,3); F1LOOP(j,2); F1LOOP(i,1))
#define F1LOOP(a,i) do (a)=fstart((i)),fend((i)); if ((a) <= nn) then; (a)(a)=(a)-1; else; (a)(a)=n+1-(a); endif

! stencil operators (implemented as preprocessor macros)
#define HR(x,y,z) fld(:,i+(x),j+(y),k+(z))
#define LAPLACE(x,y,z) fld(:,i+(x),j+(y),k+(z))
#define GRAD2(x,y,z) sum((fld(:,i+(x),j+(y),k+(z))-fld(:,i,j,k))**2)
#define GR2(x,y,z) (fld(:,i+(x),j+(y),k+(z))-fld(:,i,j,k))**2

#define PI_PHI(x,y,z) ( (fldp(:,i+(x),j+(y),k+(z))-fldp(:,i,j,k))*(fld(:,i+(x),j+(y),k+(z))-fld(:,i,j,k)) )

#define RANK0(O) (O(0,0,0))
#define RANK1(O) (O(-1,0,0) + O(1,0,0) + O(0,-1,0) + O(0,1,0) + O(0,0,-1) + O(0,0,1))
#define RANK2(O) (O(-1,-1,0) + O(1,-1,0) + O(-1,1,0) + O(1,1,0) + O(-1,0,-1) + O(1,0,-1) + O(-1,0,1) + O(1,0,1) + O(0,-1,-1) + O(0,1,-1) + O(0,-1,1) + O(0,1,1))
#define RANK3(O) (O(-1,-1,-1) + O(1,-1,-1) + O(-1,1,-1) + O(1,1,-1) + O(-1,-1,1) + O(1,-1,1) + O(-1,1,1) + O(1,1,1))

! pgf90 has a line length limit of 2048 which is too short for the expanded stencil
#ifdef __PGI
#define STENCIL(C,O) ((C/**/0) * RANK0(O) + (C/**/1) * RANK1(O) + (C/**/2) * RANK2(O) + (C/**/3) * RANK3(O))
#else
#ifdef __GFORTRAN__
#define STENCIL(C,O) ((C/**/0) * RANK0(O) + (C/**/1) * RANK1(O) + (C/**/2) * RANK2(O) + (C/**/3) * RANK3(O))
#else
#define STENCIL(C,O) ((C ## 0) * RANK0(O) + (C ## 1) * RANK1(O) + (C ## 2) * RANK2(O) + (C ## 3) * RANK3(O))
#endif
#endif
