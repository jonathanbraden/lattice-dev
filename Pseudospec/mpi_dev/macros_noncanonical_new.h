#define FLD_DIFF(a,x,y,z) (fld(a,i+(x),j+(y),k+(z))-fld(a,i,j,k))
#define VAR_AVE(variable,x,y,z)  (variable(i+(x),j+(y),k+(z))+variable(i,j,k))

#define LAPLACE_A(x,y,z,a) fld(a,i+(x),j+(y),k+(z))
#define GRAD_AB(x,y,z,a,b) (FLD_DIFF(a,x,y,z))*(FLD_DIFF(b,x,y,z))
#define GRAD2_A(x,y,z,a) (FLD_DIFF(a,x,y,z))**2
#define MET_DOT_LAP(x,y,z,a,v) (VAR_AVE(v,x,y,z))*(FLD_DIFF(a,x,y,z))
#define LAPLACIAN(x,y,z) fld(:,i+(x),j+(y),k+(z))

#define RK0(O,flds...) ( O(0,0,0,##flds) )
#define RK1(O,flds...) ( O(-1,0,0,##flds) + O(1,0,0,##flds) + O(0,-1,0,##flds) + O(0,1,0,##flds) + O(0,0,-1,##flds) + O(0,0,1,##flds) )
#define RK2(O,flds...) ( O(-1,-1,0,##flds) + O(1,-1,0,##flds) + O(-1,1,0,##flds) + O(1,1,0,##flds) + O(-1,0,-1,##flds) + O(1,0,-1,##flds) + O(-1,0,1,##flds) + O(1,0,1,##flds) + O(0,-1,-1,##flds) O(0,1,-1,##flds) + O(0,-1,1,##flds) + O(0,1,1,##flds) )
#define RK3(O,flds...) ( O(fld,-1,-1,-1) + O(fld,1,-1,-1) + O(fld,-1,1,-1) + O(fld,1,1,-1) + O(fld,-1,-1,1) + O(fld,1,-1,1) + O(fld,-1,1,1) + O(fld,1,1,1))
#define STENCIL_F(C,O,flds...) ((C##0)*RK0(O,##flds) + (C##1)*RK1(O,##flds) + (C##2)*RK2(O,##flds) + (C##3)*RK2(O,##flds))
