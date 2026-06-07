program main
  use mod_io 
  use mod_advect
  use mod_physics
  implicit none
  
  character(len=100) :: input_file  = "/raid/users/imai/Work/dev/tagged_moisture_model/Pre/era5/test_tohoku/merged.nc"
  character(len=100) :: output_file =  "/raid/users/imai/Work/dev/tagged_moisture_model/Work/output_tracer_tohoku.nc"

  character(len=100) :: base_time   = "hours since 2025-07-01 00:00:00"
  
  real, parameter    :: dt          = 600.0 ! Model time step [s] (10 minutes)
  integer, parameter :: sub_steps   = 6     ! Number of sub-steps per data hour (3600s / 600s = 6)
  real, parameter    :: dlon        = 1.0
  real, parameter    :: dlat        = 1.0

  real, allocatable  :: tracer(:,:)
  real, allocatable  :: w_model(:,:)
  real, allocatable  :: source_mask(:,:)

  integer            :: t, s, i, j
  real               :: ratio, total_tracer

  print *, "# ==================================================================="
  print *, "#  START the COLORED MOISTURE TRACKING MODEL"
  print *, "# ==================================================================="

  call init_grid(input_file)
  allocate(tracer(nx, ny), w_model(nx, ny), source_mask(nx, ny))
  
  call init_output(output_file, base_time)

  tracer = 0.0
  source_mask = 0.0

  ! example 
  do j = 1, ny
    do i = 1, nx
      if (lon(i) >= 118.0 .and. lon(i) <= 130.0 .and. lat(j) >= 23.0  .and. lat(j) <= 33.0) then 
        source_mask(i,j) = 1.0
      end if
    end do ! i
  end do ! j  

  print *, "Source mask configured."

  ! initialize the first step 
  call read_step(1)
  w_model = tcwv
  
  ! write out the initial state
  call write_output_step(1, tracer, w_model)

  print *, "Starting simulation loop..."  

  do t = 1, nt - 1
    do s = 1, sub_steps
      ! Advection step 
      call calculate_advection(dt, dlon, dlat, tracer, w_model)

      ! Phsics step
      call calculate_physics(dt, tracer, w_model, source_mask)
    end do ! s

    ! ============================================================
    ! Nudgingg / Mass Correction Step
    ! ============================================================
    call read_step(t + 1)

    do j = 1, ny
      do i = 1, nx
        if (w_model(i,j) > 1.0e-12) then
          ratio        = tcwv(i,j) / w_model(i,j)
          tracer(i,j)  = tracer(i,j) * ratio
          w_model(i,j) = tcwv(i,j)
        else
          tracer(i,j) = 0.0
          w_model(i,j) = 0.0
        end if ! (w_model(i,j) > 1.0e-12)
      end do ! i
    end do ! j
    
    call write_output_step(t + 1, tracer, w_model)

    total_tracer = sum(tracer)
    if (mod(t, 10) == 0 .or. t == nt - 1) then
      print '(A,I5,A,I5,A,E14.6)', " Step: ", t, " / ", nt - 1, " | Total Tagged Tracer Mass index: ", total_tracer
    end if ! (mod(t, 10) == 0 .or. t == nt - 1)
  end do ! t

  call close_io()
  call close_output()

  deallocate(tracer, w_model, source_mask)

  print *, "# ==================================================================="
  print *, "#  CALCULATION COMPLETED"
  print *, "# ==================================================================="

end program main
