!
! Module storing the time-integration subroutines.
!

module integrator

  use hamiltonian

contains

  subroutine symp2(dt, nsteps)
    real*8 :: dt
    integer :: nsteps

    integer :: j
   
    call Hamiltonian_Split(dt/2._dl, 1)
    do j=1,nsteps-1
       call symp_o2step(dt,1._dl,1._dl)
    enddo
    call symp_o2step(dt,1._dl,0._dl)
  end subroutine symp2

  subroutine symp4(dt, nsteps)
    real*8 :: dt
    integer :: nsteps

    integer :: i
    
    real*8, parameter :: w1 = 1._dl/(2._dl-2._dl**(1._dl/3._dl))
    real*8, parameter :: w0 = 1._dl - 2._dl*w1

    call Hamiltonian_Split(w1*dt/2._dl,1)
    do i=1,nsteps
       call symp_o2step(dt, w1, w0)
       call symp_o2step(dt, w0, w1)
       if (i.eq.nsteps) then
          call symp_o2step(dt, w1, 0._dl)
       else
          call symp_o2step(dt, w1, w1)
       endif
    enddo

  end subroutine symp4

  subroutine symp6(dt, nsteps)
    real*8 :: dt
    integer :: nsteps

    integer :: i

    real*8, parameter :: w1 = -1.17767998417887_dl
    real*8, parameter :: w2 = 0.235573213359357_dl
    real*8, parameter :: w3 = 0.784513610477560_dl
    real*8, parameter :: w0 = 1._dl - 2._dl*(w1+w2+w3)

    call Hamiltonian_Split(w3*dt/2._dl, 1)
    do i=1,nsteps
       
    enddo

  end subroutine symp6

  subroutine symp8(dt, nsteps)
    real*8 :: dt
    integer :: nsteps
    
    integer :: j
    real*8, parameter :: w1 = 0.74167036435061295344822780
    real*8, parameter :: w2 = -0.40910082580003159399730010
    real*8, parameter :: w3 = 0.19075471029623837995387626
    real*8, parameter :: w4 = -0.57386247111608226665638733
    real*8, parameter :: w5 = 0.29906418130365592384446354
    real*8, parameter :: w6 = 0.33462491824529818378495798
    real*8, parameter :: w7 = 0.31529309239676659663205666
    real*8, parameter :: w0 = 1._dl - 2._dl*(w1+w2+w3+w4+w5+w6+w7)
    
    call Hamiltonian_Split(w7*dt/2._dl,1)
    
    do j=1,nsteps
       call symp_o2step(dt, w7, w6)
       call symp_o2step(dt, w6, w5)
       call symp_o2step(dt, w5, w4)
       call symp_o2step(dt, w4, w3)
       call symp_o2step(dt, w3, w2)
       call symp_o2step(dt, w2, w1)
       call symp_o2step(dt, w1, w0)
       call symp_o2step(dt, w0, w1)
       call symp_o2step(dt, w1, w2)
       call symp_o2step(dt, w2, w3)
       call symp_o2step(dt, w3, w4)
       call symp_o2step(dt, w4, w5)
       call symp_o2step(dt, w5, w6)
       call symp_o2step(dt, w6, w7)
       if (j.eq.nsteps) then
          call symp_o2step(dt, w7, 0._dl)
       else
          call symp_o2step(dt, w7, w7)
       endif
    enddo
  end subroutine symp8

!
! Gauss-Legendre Integrators of various orders
!
! 4th order implicit Gauss-Legendre integrator
  subroutine gl4(y, dt)
    integer, parameter :: s = 2
    real :: y(nvar), g(nvar,s), dt; integer i, k
    
    ! Butcher tableau for 4th order Gauss-Legendre method
    real, parameter :: a(s,s) = reshape( (/ 0.25, 0.25 - 0.5/sqrt(3.0), 0.25 + 0.5/sqrt(3.0), 0.25 /), [s,s] )
    real, parameter ::   b(s) = (/ 0.5, 0.5 /)
        
    ! iterate trial steps
    g = 0.0
    do k = 1,16
       g = matmul(g,a)
       do i = 1,s
          call evalf(y + g(:,i)*dt, g(:,i))
       end do
    end do
        
    ! update the solution
    y = y + matmul(g,b)*dt
  end subroutine gl4

  ! 6th order implicit Gauss-Legendre integrator
  subroutine gl6(y, dt)
    integer, parameter :: s = 3
    real :: y(nvar), g(nvar,s), dt
    integer :: i, k
        
    ! Butcher tableau for 6th order Gauss-Legendre method
    real, parameter :: a(s,s) = reshape( (/ &
         5.0/36.0, 2.0/9.0 - 1.0/sqrt(15.0), 5.0/36.0 - 0.5/sqrt(15.0), &
         5.0/36.0 + sqrt(15.0)/24.0, 2.0/9.0, 5.0/36.0 - sqrt(15.0)/24.0, &
         5.0/36.0 + 0.5/sqrt(15.0), 2.0/9.0 + 1.0/sqrt(15.0), 5.0/36.0 /) &
         , [s,s])
    real, parameter ::   b(s) = (/ 5.0/18.0, 4.0/9.0, 5.0/18.0/)
        
    ! iterate trial steps
    g = 0.0
    do k = 1,16
       g = matmul(g,a)
       do i = 1,s
          call evalf(y + g(:,i)*dt, g(:,i))
       end do
    end do
        
    ! update the solution
    y = y + matmul(g,b)*dt
  end subroutine gl6

  ! 8th order implicit Gauss-Legendre integrator
  subroutine gl8(y, dt)
    integer, parameter :: s = 4
    real y(nvar), g(nvar,s), dt; integer i, k
        
    ! Butcher tableau for 8th order Gauss-Legendre method
    real, parameter :: a(s,s) = reshape( (/ &
         0.869637112843634643432659873054998518D-1, -0.266041800849987933133851304769531093D-1, &
         0.126274626894047245150568805746180936D-1, -0.355514968579568315691098184956958860D-2, &
         0.188118117499868071650685545087171160D0,   0.163036288715636535656734012694500148D0,  &
         -0.278804286024708952241511064189974107D-1,  0.673550059453815551539866908570375889D-2, &
         0.167191921974188773171133305525295945D0,   0.353953006033743966537619131807997707D0,  &
         0.163036288715636535656734012694500148D0,  -0.141906949311411429641535704761714564D-1, &
         0.177482572254522611843442956460569292D0,   0.313445114741868346798411144814382203D0,  &
         0.352676757516271864626853155865953406D0,   0.869637112843634643432659873054998518D-1 /), [s,s] )
    real, parameter ::   b(s) = (/ &
         0.173927422568726928686531974610999704D0,   0.326072577431273071313468025389000296D0,  &
         0.326072577431273071313468025389000296D0,   0.173927422568726928686531974610999704D0  /)
        
    ! iterate trial steps
    g = 0.0
    do k = 1,16
       g = matmul(g,a)
       do i = 1,s
          call evalf(y + g(:,i)*dt, g(:,i))
       end do
    end do
        
    ! update the solution
    y = y + matmul(g,b)*dt
  end subroutine gl8

  subroutine gl10( y, dt )
    real*8 :: y(nvar)
    real*8 :: dt

    integer, parameter :: s = 5
    real*8 :: g(nvar,s)

    ! Butcher tableau for 8th order Gauss-Legendre method
    real*8, parameter :: a(s,s) = reshape( (/ &
         0.5923172126404727187856601017997934066D-1, -1.9570364359076037492643214050884060018D-2, &
         1.1254400818642955552716244215090748773D-2, -0.5593793660812184876817721964475928216D-2, &
         1.5881129678659985393652424705934162371D-3,  1.2815100567004528349616684832951382219D-1, &
         1.1965716762484161701032287870890954823D-1, -2.4592114619642200389318251686004016630D-2, &
         1.0318280670683357408953945056355839486D-2, -2.7689943987696030442826307588795957613D-3, &
         1.1377628800422460252874127381536557686D-1,  2.6000465168064151859240589518757397939D-1, &
         1.4222222222222222222222222222222222222D-1, -2.0690316430958284571760137769754882933D-2, &
         4.6871545238699412283907465445931044619D-3,  1.2123243692686414680141465111883827708D-1, &
         2.2899605457899987661169181236146325697D-1,  3.0903655906408664483376269613044846112D-1, &
         1.1965716762484161701032287870890954823D-1, -0.9687563141950739739034827969555140871D-2, &
         1.1687532956022854521776677788936526508D-1,  2.4490812891049541889746347938229502468D-1, &
         2.7319004362580148889172820022935369566D-1,  2.5888469960875927151328897146870315648D-1, &
         0.5923172126404727187856601017997934066D-1 /) , [s,s])
    real, parameter :: b(s) = (/ &
         1.1846344252809454375713202035995868132D-1,  2.3931433524968323402064575741781909646D-1, &
         2.8444444444444444444444444444444444444D-1,  2.3931433524968323402064575741781909646D-1, &
         1.1846344252809454375713202035995868132D-1 /)

      integer :: i,k

      g = 0.
      do k=1,16
         g = matmul(g,a)
         do i=1,s
            call evalf( y+g(:,i)*dt, g(:,i) )
         enddo
      enddo
      y = y + matmul(g,b)*dt

    end subroutine gl10

end module integrator
