!!!!!!!!!!!!!!!!!!!
!
! Lattice-based Gauss-Legendre Integrator
!
!!!!!!!!!!!!!!!!!!!

#if ORDER_GL==2
#define STEPS_GL integer, paraemter :: s=3
#define BUTCHER_TABLEAU real, parameter :: a(s,s) = 
#elif ORDER_GL==4
#define BUTCHER_TABLEAU 
#elif ORDER_GL==6
#define BUTCHER_TABLEAU
#elif ORDER_GL==8
#define BUTCHER_TABLEAU
#elif ORDER_GL==10
#define BUTCHER_TABLEAU
#endif

module integrator_gl

  real, dimension(IRANGE, ) :: derivs

contains

  subroutine gl_step(ycur, dt)
    real, dimension(), intent(inout) :: y
    real, intent(in) :: dt
    integer, parameter :: niter=8
    STEPS_GL
    BUTCHER_TABLEAU

    g = 0.0
    do k=1,niter
       g = DGEMM  !matmul(g,a)
       do i=1,s
          call eval_derivs()
       enddo
    enddo
    ycur = ycur + DGEMM   !matmul(g,b)*dt
  end subroutine gl_step

  subroutine derivs()
    
#if DISCRETIZATION==1
    
#else  ! default to finite-difference
    call wrap_fields(fld)
    call wrap_fields(fldp)
#endif
    

  end subroutine derivs

end module integrator_gl
