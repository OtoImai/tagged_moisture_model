module mod_io
  use netcdf
  implicit none
  
  ! for input
  integer :: nx, ny, nt, ncid
  real, allocatable :: lon(:), lat(:)
  real, allocatable :: tcwv(:,:), u_eff(:,:), v_eff(:,:), evap(:,:), precip(:,:)

  ! for output
  integer :: ncid_out, varid_tracer_out, varid_w_model_out

contains

  subroutine check(status)
    integer, intent(in) :: status
      if (status /= nf90_noerr) then
        print *, "NetCDF Error: ", trim(nf90_strerror(status))
        stop "Program terminated due to NetCDF error."
    end if
  end subroutine check

  subroutine init_grid(filename)
    character(len=*), intent(in) :: filename
    integer :: dimid_lon, dimid_lat, dimid_time, varid_lon, varid_lat
    
    print *, "Initializing grid from: ", trim(filename)
    call check( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
    
    call check( nf90_inq_dimid(ncid, "longitude", dimid_lon) )
    call check( nf90_inquire_dimension(ncid, dimid_lon, len = nx) )
    call check( nf90_inq_dimid(ncid, "latitude", dimid_lat) )
    call check( nf90_inquire_dimension(ncid, dimid_lat, len = ny) )
    call check( nf90_inq_dimid(ncid, "valid_time", dimid_time) )
    call check( nf90_inquire_dimension(ncid, dimid_time, len = nt) )
    
    allocate(lon(nx), lat(ny))
    allocate(tcwv(nx, ny), u_eff(nx, ny), v_eff(nx, ny), evap(nx, ny), precip(nx, ny))
    
    call check( nf90_inq_varid(ncid, "longitude", varid_lon) )
    call check( nf90_get_var(ncid, varid_lon, lon) )
    call check( nf90_inq_varid(ncid, "latitude", varid_lat) )
    call check( nf90_get_var(ncid, varid_lat, lat) )
    
    print *, "Grid Initialization Complete."
  end subroutine init_grid

  subroutine read_step(t_index)
    integer, intent(in) :: t_index
    integer :: varid_tcwv, varid_viwve, varid_viwvn, varid_e, varid_tp
    
    call check( nf90_inq_varid(ncid, "tcwv", varid_tcwv) )
    call check( nf90_inq_varid(ncid, "viwve", varid_viwve) )
    call check( nf90_inq_varid(ncid, "viwvn", varid_viwvn) )
    call check( nf90_inq_varid(ncid, "e", varid_e) )
    call check( nf90_inq_varid(ncid, "tp", varid_tp) )
    
    call check( nf90_get_var(ncid, varid_tcwv, tcwv, start=(/1,1,t_index/), count=(/nx,ny,1/)) )
    call check( nf90_get_var(ncid, varid_viwve, u_eff, start=(/1,1,t_index/), count=(/nx,ny,1/)) )
    call check( nf90_get_var(ncid, varid_viwvn, v_eff, start=(/1,1,t_index/), count=(/nx,ny,1/)) )
    call check( nf90_get_var(ncid, varid_e, evap, start=(/1,1,t_index/), count=(/nx,ny,1/)) )
    call check( nf90_get_var(ncid, varid_tp, precip, start=(/1,1,t_index/), count=(/nx,ny,1/)) )
    
  end subroutine read_step

  subroutine close_io()
    call check( nf90_close(ncid) ) 
  end subroutine close_io
  
  subroutine init_output(filename, base_time)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: base_time 
    integer                      :: dimid_lon, dimid_lat, dimid_time
    integer                      :: varid_lon, varid_lat, varid_time
    
    print *, "Creating output file: ", trim(filename)
    call check( nf90_create(trim(filename), NF90_CLOBBER, ncid_out) )

    
    call check( nf90_def_dim(ncid_out, "lon", nx, dimid_lon) )
    call check( nf90_def_dim(ncid_out, "lat", ny, dimid_lat) )
    call check( nf90_def_dim(ncid_out, "time", nt, dimid_time) )
    
    call check( nf90_def_var(ncid_out, "lon", NF90_REAL, (/dimid_lon/), varid_lon) )
    call check( nf90_def_var(ncid_out, "lat", NF90_REAL, (/dimid_lat/), varid_lat) )
    call check( nf90_def_var(ncid_out, "time", NF90_INT, (/dimid_time/), varid_time) )

    call check( nf90_def_var(ncid_out, "tracer", NF90_REAL, (/dimid_lon, dimid_lat, dimid_time/), varid_tracer_out) )
    call check( nf90_def_var(ncid_out, "w_model", NF90_REAL, (/dimid_lon, dimid_lat, dimid_time/), varid_w_model_out) )

    call check( nf90_put_att(ncid_out, varid_lon, "units", "degrees_east") )
    call check( nf90_put_att(ncid_out, varid_lon, "long_name", "longitude") )
    
    call check( nf90_put_att(ncid_out, varid_lat, "units", "degrees_north") )
    call check( nf90_put_att(ncid_out, varid_lat, "long_name", "latitude") )    

    call check( nf90_put_att(ncid_out, varid_time, "units", trim(base_time)) )
    call check( nf90_put_att(ncid_out, varid_time, "long_name", "time") )
    call check( nf90_put_att(ncid_out, varid_time, "long_name", "time") )

    call check( nf90_put_att(ncid_out, varid_tracer_out, "units", "kg m**-2") )
    call check( nf90_put_att(ncid_out, varid_tracer_out, "long_name", "Tracer Moisture") )
    call check( nf90_put_att(ncid_out, varid_w_model_out, "units", "kg m**-2") )
    call check( nf90_put_att(ncid_out, varid_w_model_out, "long_name", "Model Total Moisture") )
    
    call check( nf90_enddef(ncid_out) )
    
    call check( nf90_put_var(ncid_out, varid_lon, lon) )
    call check( nf90_put_var(ncid_out, varid_lat, lat) )
  end subroutine init_output

  subroutine write_output_step(t_index, tracer_data, w_model_data)
    integer, intent(in) :: t_index
    real, intent(in)    :: tracer_data(nx, ny)
    real, intent(in)    :: w_model_data(nx, ny)
    integer :: varid_time
    
    call check( nf90_put_var(ncid_out, varid_tracer_out, tracer_data, start=(/1, 1, t_index/), count=(/nx, ny, 1/)) )
    call check( nf90_put_var(ncid_out, varid_w_model_out, w_model_data, start=(/1, 1, t_index/), count=(/nx, ny, 1/)) )
    
    call check( nf90_inq_varid(ncid_out, "time", varid_time) )
    call check( nf90_put_var(ncid_out, varid_time, (/t_index/), start=(/t_index/), count=(/1/)) )
  end subroutine write_output_step

  subroutine close_output()
    call check( nf90_close(ncid_out) )
    print *, "Output file closed successfully."
  end subroutine close_output

end module mod_io
