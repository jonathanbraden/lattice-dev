!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!MODULE : Metric
!
! Stores information on the smooth FRW metric
!
!>@author
!> Jonathan Braden, University College London
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!ISSUES!!!!
! Now I need to have H0 defined somewhere to initialize the metric

module Metric
  use constants, only : dl
  implicit none

  real(dl) :: yscl, ysclp

contains

  subroutine init_metric_conformal()
    yscl = 1._dl
    ysclp = -6.*H0*yscl
  end subroutine init_metric_conformal

  subroutine init_metric_cosmic()
    yscl = 1._dl
    ysclp = -4.*H0*yscl
  end subroutine init_metric_cosmic

end module Metric
