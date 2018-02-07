module Metric_Model_Conformal

  use params
  implicit none

  type Grav_Metric
     real(dl) :: yscl, ysclp
  end type Grav_Metric

  real(dl) :: yscl, ysclp
  
contains
  
  function get_scale_factor() return(scl)
    real(dl) :: scl
    get_scale_factor = yscl
  end function get_scale_factor

  function get_hubble() return(hub)
    real(dl) :: hub
    hub = -ysclp / 6._dl / yscl**2
  end function get_hubble

  function kinetic_norm()
    real(dl) :: kinetic_norm
    kinetic_norm = 1._dl / yscl**6
  end function kinetic_norm

  subroutine calc_metric()
    real(dl) :: rho
    rho = scalar_rho()
    ysclp = -sqrt(yscl**4*12._dl*rho)
  end subroutine calc_metric
  
end module Metric_Model_Conformal
