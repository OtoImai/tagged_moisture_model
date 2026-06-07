module mod_physics
  use mod_io, only: nx, ny, evap, precip
  implicit none

contains 

  subroutine calculate_physics(dt, q, w_model, source_mask)
    real, intent(in)    :: dt 
    real, intent(inout) :: q(:,:)
    real, intent(inout) :: w_model(:,:)
    real, intent(in)    :: source_mask(:,:)

    real                :: r
    integer             :: i, j 

    do j = 1, ny
      do i = 1, nx 
        e_rate = ( -evap(i,j) * 1000.0 ) / 3600.0
        p_rate = ( precip(i,j) * 1000.0 ) / 3600.0

        if (e_rate < 0.0) e_rate = 0.0
        if (p_rate < 0.0) p_rate = 0.0

        if (tcwv(i,j) > 1.0e-3) then
          tracer_ratio = tracer(i,j) / tcwv(i,j)
            if (tracer_ratio > 1.0) tracer_ratio = 1.0
        else
          tracer_ratio = 0.0
        endif ! (tcwv(i,j) > 1.0e-3)
        
        tracer(i,j) = tracer(i,j) + dt * (e_rate * source_mask(i,j) - p_rate * tracer_ratio)

        if (tracer(i,j) < 0.0) tracer(i,j) = 0.0
        if (tracer(i,j) > tcwv(i,j)) tracer(i,j) = tcwv(i,j)    

      end do ! i
    end do ! j

  end subroutine calculate_physics

end module mod_physics
