module mod_physics
  use mod_io, only: nx, ny, evap, precip
  implicit none

  real, parameter :: rho_water  = 1000.0  ! density of water [kg/m3]
  real, parameter :: accum_time = 3600.0  ! accumulation time of data [s] (Typically 1 hour = 3600 s for hourly ERA5)
contains 

  subroutine calculate_physics(dt, q, w_model, source_mask)
    real, intent(in)    :: dt 
    real, intent(inout) :: q(:,:)           ! colored moisture (tracer)
    real, intent(inout) :: w_model(:,:)     ! total moisture in the model
    real, intent(in)    :: source_mask(:,:)

    real                :: r, e_flux, p_flux
    integer             :: i, j 

    do j = 1, ny
      do i = 1, nx 
        ! unit conversion (m/h -> kg/m2/s)
        e_flux = ( -evap(i,j) * rho_water ) / accum_time
        p_flux = ( precip(i,j) * rho_water ) / accum_time

        if (e_flux < 0.0) e_flux = 0.0
        if (p_flux < 0.0) p_flux = 0.0

        ! calc the current r 
        r = q(i,j) / max(w_model(i,j), 1.0e-12)
        ! water vapor budget
        w_model(i,j) = w_model(i,j) + dt * (e_flux - p_flux)
        ! tracer budget
        q(i,j) = q(i,j) + dt * (e_flux * source_mask(i,j) - p_flux * r)

        if (w_model(i,j) < 0.0) w_model(i,j) = 0.0
        if (q(i,j) < 0.0) q(i,j) = 0.0
        
      end do ! i
    end do ! j

  end subroutine calculate_physics

end module mod_physics
