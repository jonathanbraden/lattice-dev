!#define SPECTRAL 1
#define DISCRETE 1

#define THREEDIM 1
!#define TWODIM 1
!#define ONEDIM 1

#define LOOPEVOLVE 1
!#define VECTORIZE 1

!#define CONFORMAL 1
!#define COSMIC 1
!#define MINKOWSKI 1

#ifdef THREEDIM
#define IRANGE 1:nx,1:ny,1:nz
#define SIRANGE 0:nx+1,0:ny+1,0:nz+1
#define LATIND i,j,k
#define NPTS dble(nx)*dble(ny)*dble(nz)
#define FLOOP do k=1,nz; do j=1,ny; do i=1,nx
#define FLOOPEND enddo; enddo; enddo
!#define FKLOOP do k=1,nz; if (k<=nnz) then; kk=k-1; else; kk= ;endif	\
!  do j=1,ny; if (j<nny) then; jj=j-1; else; jj= ;endif			\
!  do i=1,nnx; ii=i-1
!#define FKLOOPEND enddo; enddo; enddo

#define LAPLACIAN(x,y,z) fld(:,i+(x),j+(y),k+(z))
#define GRAD2(x,y,z) ( fld(:,i+(x),j+(y),k+(z))-fld(:,i,j,k) )**2
#define GRAD12(x,y,z) (fld(1,i+(x),j+(y),k+(z))-fld(1,i,j,k))*(fld(2,i+(x),j+(y),k+(z))-fld(2,i,j,k))
#define GDOTL1(x,y,z) ( tmp_lat(i+(x),j+(y),k+(z))+tmp_lat(i,j,k) )*(fld(1,i+(x),j+(y),k+(z))-fld(1,i,j,k))
#define GDOTL2(x,y,z) ( tmp_lat(i+(x),j+(y),k+(z))+tmp_lat(i,j,k) )*(fld(2,i+(x),j+(y),k+(z))-fld(2,i,j,k))
#define GDOTL3(x,y,z) ( tmp_lat(i+(x),j+(y),i+(z))+tmp_lat(i,j,k) )*(fld(3,i+(x),j+(y),k+(z))-fld(3,i,j,k)) 

#define RANK0(O) (O(0,0,0))
#define RANK1(O) (O(-1,0,0) + O(1,0,0) + O(0,-1,0) + O(0,1,0) + O(0,0,-1) + O(0,0,1))
#define RANK2(O) (O(-1,-1,0) + O(1,-1,0) + O(-1,1,0) + O(1,1,0) + O(-1,0,-1) + O(1,0,-1) + O(-1,0,1) + O(1,0,1) + O(0,-1,-1) + O(0,1,-1) + O(0,-1,1) + O(0,1,1))
#define RANK3(O) (O(-1,-1,-1) + O(1,-1,-1) + O(-1,1,-1) + O(1,1,-1) + O(-1,-1,1) + O(1,-1,1) + O(-1,1,1) + O(1,1,1))
!#define STENCIL(C,O) ((C/**/0)*RANK0(O) + (C/**/1)*RANK1(O) + (C/**/2)*RANK2(O) + (C/**/3)*RANK3(O))

#define STENCIL(C,O) ((C ## 0)*RANK0(O) + (C ## 1)*RANK1(O) + (C ## 2)*RANK2(O) + (C ## 3)*RANK3(O))

#define RANKFLD0(O,fld)
#define RANKFLD1(O,fld)
#define RANKFLD2(O,fld)
#define RANKFLD3(O,fld)
#define STENCIL_FLD(C,O) ((C ## 0)*RANKFLD0(O) + (C ## 1)*RANKFLD1(O) + (C ## 2)*RANKFLD2(O) + (C ## 3)*RANKFLD3(O))  
#endif

#ifdef TWODIM
#define IRANGE 1:nx,1:ny
#define SIRANGE 0:nx+1,0:ny+1
#define LATIND i,j
#define NPTS dble(nx)*dble(ny)
#define FLOOP do j=1,ny; do i=1,nx
#define FLOOPEND enddo;enddo
#define LAPLACIAN(x,y)
#define RANK0(O) (O(0,0))
#define RANK1(O) (O(1,0) + O(-1,0) + O(0,1) + O(0,-1))
#define RANK2(O) (O(1,1) + O(1,-1) + O(-1,1) + O(-1,-1))
#define STENCIL(C,O) ((C/**/0)*RANK0(O) + (C/**/1)*RANK1(O) + (C/**/2)*RANK2(O))
#endif

#ifdef ONEDIM
#define IRANGE 1:nx
#define SIRANGE 0:nx+1
#define LATIND i
#define NPTS dble(nx)
#define FLOOP do i=1,nx
#define FLOOPEND enddo
#endif
