module mod_advect
  use mod_io, only: nx, ny, lat, u_eff, v_eff
  implicit none

  real, parameter :: R_earth = 6371000.0
  real, parameter :: pi = 3.1415926535

contains

  subroutine calculate_advection(dt, dlon, dlat, q, w_model)
    real, intent(in)    :: dt, dlon, dlat
    real, intent(inout) :: q(:,:)
    real, intent(inout) :: w_model(:,:)

    real, allocatable   :: q_star(:,:), w_star(:,:)
    real, allocatable   :: div_q1(:,:), div_w1(:,:)
    real, allocatable   :: div_q2(:,:), div_w2(:,:)
    real, allocatable   :: r(:,:), r_one(:,:)
    integer             :: i, j

    allocate(q_star(nx, ny), w_star(nx, ny))
    allocate(div_q1(nx,ny), div_w1(nx,ny), div_q2(nx,ny), div_w2(nx,ny))  
    allocate(r(nx, ny), r_one(nx, ny))

    ! =================================================
    ! RK2 Predictor Step
    ! =================================================

    ! 1. calc the current ratio (r^n)
    do j = 1, ny
      do i = 1, nx
        r(i,j) = q(i,j) / max(w_model(i,j), 1.0e-12)
        r_one(i,j) = 1.0
      end do ! i
    end do ! j 

    ! 2. calc the current div of flux
    call calc_flux_divergence(r, u_eff, v_eff, dlon, dlat, div_q1)
    call calc_flux_divergence(r_one, u_eff, v_eff, dlon, dlat, div_w1)


    ! 3. Predict (euler forward)
    q_star = q - dt * div_q1
    w_star = w_model - dt * div_w1
    where (q_star < 0.0) q_star = 0.0
    where (w_star < 0.0) w_star = 0.0

    ! =================================================
    ! RK2 Corrector Step
    ! =================================================
   
    ! 4. recalc the ratio (r^*)
    do j = 1, ny
      do i = 1, nx
        r(i,j) = q_star(i,j) / max(w_star(i,j), 1.0e-12)
      end do ! i
    end do ! j 

    ! 5. recalc the div of flux
    call calc_flux_divergence(r, u_eff, v_eff, dlon, dlat, div_q2)
    call calc_flux_divergence(r_one, u_eff, v_eff, dlon, dlat, div_w2)

    ! 6. time integration 
    q = q - 0.5 * dt * (div_q1 + div_q2)
    w_model = w_model - 0.5 * dt * (div_w1 + div_w2)

    where (q < 0.0) q = 0.0
    where (w_model < 0.0) w_model = 0.0
    
    deallocate(q_star, w_star)
    deallocate(div_q1, div_w1, div_q2, div_w2)
    deallocate(r, r_one)
  end subroutine calculate_advection  

  subroutine calc_flux_advection(r_in, u_flux, v_flux, dlon, dlat, div_out)
    real, intent(in)  :: r_in(:,:), u_flux(:,:), v_flux(:,:)
    real, intent(in)  :: dlon, dlat
    real, intent(out) :: div_out(nx, ny)
  
    real              :: dx, dy, lat_rad
    real              :: u_face, v_face, fx_left, fx_right, fy_bot, fy_top
    integer           :: i, j, i_prev, i_next, j_prev, j_next

    dy = R_earth * (dlat * pi / 180.0)

    do j = 1, ny
      lat_rad = lat(j) * pi / 180.0
      dx      = R_earth * cos(lat_rad) * (dlon * pi / 180.0)

      j_prev  = max(1, j-1)
      j_next  = min(ny, j+1)
 
      do i = 1, nx
        i_prev = max(1, i-1)
        i_next = min(nx, i+1)
        
        ! =================================================
        ! X axis Boundary
        ! =================================================

        ! West Boundary
        u_face = 0.5 * (u_flux(i_prev, j) + u_flux(i, j))
        if (u_face > 0.0) then
          fx_left = u_face * r_in(i_prev, j)
        else
          fx_left = u_face * r_in(i, j)
        end if ! (u_face > 0.0)
        ! East Boundary
        u_face = 0.5 * (u_flux(i, j) + u_flux(i_next, j))
        if (u_face > 0.0) then
          fx_right = u_face * r_in(i, j)
        else
          fx_right = u_face * r_in(i_next, j)
        end if ! (u_face > 0.0)
  
        ! =================================================
        ! X axis Boundary
        ! =================================================
        
        ! South Boundary
        v_face = 0.5 * (v_flux(i, j_prev) + v_flux(i,j))
        if (v_face > 0.0) then
          fy_bot = v_face * r_in(i, j_prev) 
        else
          fy_bot = v_face * r_in(i, j)
        end if ! (v_face > 0.0)
        ! North Boundary
        v_face = 0.5 * (v_flux(i, j) + v_flux(i,j_next))
        if (v_face > 0.0) then
          fy_top = v_face * r_in(i, j) 
        else
          fy_top = v_face * r_in(i, j_next)
        end if ! (v_face > 0.0)

        div_out(i,j) = (fx_right - fx_left) / dx + (fy_top - fy_bot) / dy

      end do ! i
    end do  ! j
  end subroutine calc_flux_advection
end module mod_advect
