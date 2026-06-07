module mod_advect
  use mod_io, only: nx, ny, lat, tcwv, u_eff, v_eff
  implicit none

  real, parameter :: R_earth = 6371000.0
  real, parameter :: pi = 3.1415926535

contains

  subroutine calculate_advection(dt, dlon, dlat, q)
    real, intent(in)    :: dt, dlon, dlat
    real, intent(inout) :: q(:,:)

    real, allocatable   :: tcwv_star(:,:) ,tcwv_next(:,:)
    real, allocatable   :: tend1(:,:), tend2(:,:)
    real, allocatable   :: u_vel(:,:), v_vel(:,:)
    integer             :: i, j

    allocate(tcwv_star(nx, ny), tcwv_next(nx, ny))
    allocate(tend1(nx, ny), tend2(nx, ny))
    allocate(u_vel(nx, ny), v_vel(nx, ny))

    tcwv_next = tcwv

    do j = 2, ny - 1
      do i = 2, nx -1 
        if (tcwv(i,j) > 1.0e-3) then
          u_vel(i,j) = u_eff(i,j) / tcwv(i,j)
          v_vel(i,j) = v_eff(i,j) / tcwv(i,j)
        else 
          u_vel = 0.0
          v_vel = 0.0
        end if ! (tcwv(i,j) > 1.0e-3)
      end do ! i
    end do ! j

    call calc_tendency(q, u_vel, v_vel, dlon, dlat, tend1)
    tcwv_star = q + dt * tend1
    where (tcwv_star < 0.0) tcwv_star = 0.0

    call calc_tendency(tcwv_star, u_vel, v_vel, dlon, dlat, tend2)
    tcwv_next = q + 0.5 * dt * (tend1 + tend2)
    where (tcwv_next < 0.0) tcwv_next = 0.0
  
    q = tcwv_next

    deallocate(tcwv_star, tcwv_next, tend1, tend2, u_vel, v_vel)

  end subroutine calculate_advection

  subroutine calc_tendency(q_in, u_in, v_in, dlon, dlat, tend_out)
    real, intent(in)  :: q_in(:,:), u_in(:,:), v_in(:,:)
    real, intent(in)  :: dlon, dlat
    real, intent(out) :: tend_out(nx, ny)
    real              :: dx, dy, lat_rad, flux_x, flux_y
    integer           :: i, j, i_prev, i_next, j_prev, j_next

    tend_out = 0.0
    dy       = R_earth * (dlat * pi / 180.0)
    
    do j = 1, ny 
      lat_rad = lat(j) * pi / 180.0
      dx      = R_earth * (dlat * pi /180.0)
      do i = 1, nx
        i_prev = max(1,i-1); i_next = min(nx, i+1)
        j_prev = max(1,j-1); j_next = min(nx, j+1)

        if (u_in(i,j) > 0.0) then 
          flux_x = u_in(i,j) * (q_in(i,j) - q_in(i_prev, j)) / dx
        else
          flux_x = u_in(i,j) * (q_in(i_next,j) - q_in(i, j)) / dx
        end if ! (u_in(i,j) > 0.0)
        
        if (v_in(i,j) > 0.0) then 
          flux_y = v_in(i,j) * (q_in(i,j) - q_in(i, j_prev)) / dy
        else
          flux_y = v_in(i,j) * (q_in(i,j_next) - q_in(i, j)) / dy
        end if ! (v_in(i,j) > 0.0)

        tend_out(i,j) = - (flux_x + flux_y)
      end do ! i
    end do ! j 
  end subroutine calc_tendency


end module mod_advect
