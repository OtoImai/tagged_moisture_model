program main
  use mod_io 
  use mod_advect
  use mod_physics
  implicit none
  
  character(len=100) :: input_file  = "/raid/users/imai/Work/dev/tagged_moisture_model/Pre/era5/test/merged.nc"
  character(len=100) :: output_file =  "/raid/users/imai/Work/dev/tagged_moisture_model/Work/test/output_tracer.nc"

  character(len=100) :: base_time   = "hours since 2025-07-01 00:00:00"
  
  real, parameter    :: dt          = 600.0
  integer, parameter :: sub_steps   = 6
  real, parameter    :: dlon        = 1.0
  real, parameter    :: dlat        = 1.0

  real, allocatable  :: tracer(:,:)
  real, allocatable  :: w_model(:,:)
  real, allocatable  :: source_mask(:,:)

  integer            :: t, s, i, j
  real               :: total_tracer

  print *, "# ==================================================================="
  print *, "#  START the COLORED MOISTURE TRACKING MODEL"
  print *, "# ==================================================================="

  call init_grid(input_file)
  allocate(tracer(nx, ny), w_model(nx, ny), source_mask(nx, ny))
  
  call init_output(output_file, base_time)
  print *, "OK"

  tracer = 0.0
  source_mask = 0.0

  ! example: around the Sea of Japan
  do j = 1, ny
    do i = 1, nx
      if (lon(i) >= 100.0 .and. lon(i) <= 140.0 .and. lat(j) >= 15.0  .and. lat(j) <= 45.0) then 
        source_mask(i,j) = 1.0
      end if
    end do ! i
  end do ! j  

  print *, "Source mask (Sea of Japan) configured."

  do t = 1, nt
    call read_step(t)
    do s = 1, sub_steps
      call calculate_advection(dt, dlon, dlat, tracer)
      call calculate_physics(dt, tracer, source_mask)
    end do ! s
  
    call write_output_step(t, tracer)
  
    total_tracer = sum(tracer)
    if (mod(t, 10) == 0 .or. t ==1) then
      print *, " Step:", t, "/", nt, " | Total Tracer Mass index:", total_tracer
    end if
  end do ! t 

  call close_io()
  call close_output()

  deallocate(tracer, source_mask)

  print *, "# ==================================================================="
  print *, "#  CALCULATION COMPLETED"
  print *, "# ==================================================================="

end program main
