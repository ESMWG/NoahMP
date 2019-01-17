module hrldas_netcdf_test
  use netcdf
#ifdef MPP_LAND
  use mpi
  use module_mpp_land
#endif
  implicit none
contains
  subroutine hrldas_netcdf_get_var_real2d(ncid, varname, vardata, varunit, ierr)
    implicit none
    integer, intent(in) :: ncid
    real, intent(out) :: vardata(:,:)
    character(len=*), intent(in) :: varname
    character(len=*), intent(out) :: varunit
    integer, intent(out) :: ierr

#ifdef MPP_LAND
    real, allocatable :: buffer(:,:)
    integer :: data_size
    integer :: tag, i, mpi_ierr

    if (my_id == io_id) then ! io node reads data and sends them to the other nodes
       do i = 0, nproc-1
          if (i /= my_id) then ! read data into buffer and send the buffer to node i
             allocate(buffer(iolocal_allnx(i+1),iolocal_allny(i+1)))
             call netcdf_get(ncid, varname, &
                  & buffer, varunit, ierr, &
                  & iolocal_allstartx(i+1), iolocal_allstarty(i+1), &
                  & iolocal_allnx(i+1), iolocal_allny(i+1))
             tag = 101
             data_size = iolocal_allnx(i+1) * iolocal_allny(i+1)
             call mpi_send(buffer, data_size, MPI_REAL, i, &
                  & tag, MPI_COMM_WORLD, mpi_ierr)
             deallocate(buffer)
             tag = 102
             call mpi_send(varunit, len(varunit), MPI_CHARACTER, i, &
                  & tag, MPI_COMM_WORLD, mpi_ierr)
             tag = 103
             call mpi_send(ierr, 1, MPI_INTEGER, i, &
                  & tag, MPI_COMM_WORLD, mpi_ierr)
          endif
       enddo
       call netcdf_get(ncid, varname, vardata, varunit, ierr, &
            & my_startx, my_starty, my_nx, my_ny)
    else ! receive data from the io node
       tag = 101
       call mpi_recv(vardata, my_nx*my_ny, MPI_REAL, io_id, &
            & tag, MPI_COMM_WORLD, mpp_status, mpi_ierr)
       tag = 102
       call mpi_recv(varunit, len(varunit), MPI_CHARACTER, io_id, &
            & tag, MPI_COMM_WORLD, mpp_status, mpi_ierr)
       tag = 103
       call mpi_recv(ierr, 1, MPI_INTEGER, io_id, &
            & tag, MPI_COMM_WORLD, mpp_status, mpi_ierr)
    endif
#else
    call netcdf_get(ncid, varname, vardata, varunit, ierr)
#endif
  contains
    subroutine netcdf_get(ncid, varname, vardata, varunit, ierr, &
         xstart, ystart, nx, ny)
      implicit none
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: varname
      real, intent(out) :: vardata(:,:)
      character(len=*), intent(out) :: varunit
      integer, intent(out) :: ierr
      integer, optional, intent(in) :: xstart, ystart, nx, ny
      integer :: varid
      ierr = nf90_inq_varid(ncid, trim(varname), varid)
      if (ierr /= 0) then
         varunit = 'unavailable'
         return
      endif
      ierr = nf90_get_att(ncid, varid, "units", varunit)
      if (ierr /= 0) then
         varunit = 'unkown'
         return
      endif
      if (present(xstart)) then
         ierr = nf90_get_var(ncid, varid, vardata, start=(/xstart,ystart/), count=(/nx,ny/))
      else
         ierr = nf90_get_var(ncid, varid, vardata)
      endif
    end subroutine netcdf_get
  end subroutine hrldas_netcdf_get_var_real2d


  subroutine hrldas_netcdf_get_var_integer2d(ncid, varname, vardata, varunit, ierr)
    implicit none
    integer, intent(in) :: ncid
    integer, intent(out) :: vardata(:,:)
    character(len=*), intent(in) :: varname
    character(len=*), intent(out) :: varunit
    integer, intent(out) :: ierr

#ifdef MPP_LAND
    integer, allocatable :: buffer(:,:)
    integer :: data_size
    integer :: tag, i, mpi_ierr

    if (my_id == io_id) then ! io node reads data and sends them to the other nodes
       do i = 0, nproc-1
          if (i /= my_id) then ! read data into buffer and send the buffer to node i
             allocate(buffer(iolocal_allnx(i+1),iolocal_allny(i+1)))
             call netcdf_get(ncid, varname, &
                  & buffer, varunit, ierr, &
                  & iolocal_allstartx(i+1), iolocal_allstarty(i+1), &
                  & iolocal_allnx(i+1), iolocal_allny(i+1))
             tag = 101
             data_size = iolocal_allnx(i+1) * iolocal_allny(i+1)
             call mpi_send(buffer, data_size, MPI_INTEGER, i, &
                  & tag, MPI_COMM_WORLD, mpi_ierr)
             deallocate(buffer)
             tag = 102
             call mpi_send(varunit, len(varunit), MPI_CHARACTER, i, &
                  & tag, MPI_COMM_WORLD, mpi_ierr)
             tag = 103
             call mpi_send(ierr, 1, MPI_INTEGER, i, &
                  & tag, MPI_COMM_WORLD, mpi_ierr)
          endif
       enddo
       call netcdf_get(ncid, varname, vardata, varunit, ierr, &
            & my_startx, my_starty, my_nx, my_ny)
    else ! receive data from the io node
       tag = 101
       call mpi_recv(vardata, my_nx*my_ny, MPI_INTEGER, io_id, &
            & tag, MPI_COMM_WORLD, mpp_status, mpi_ierr)
       tag = 102
       call mpi_recv(varunit, len(varunit), MPI_CHARACTER, io_id, &
            & tag, MPI_COMM_WORLD, mpp_status, mpi_ierr)
       tag = 103
       call mpi_recv(ierr, 1, MPI_INTEGER, io_id, &
            & tag, MPI_COMM_WORLD, mpp_status, mpi_ierr)
    endif
#else
    call netcdf_get(ncid, varname, vardata, varunit, ierr)
#endif
  contains
    subroutine netcdf_get(ncid, varname, vardata, varunit, ierr, &
         xstart, ystart, nx, ny)
      implicit none
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: varname
      integer, intent(out) :: vardata(:,:)
      character(len=*), intent(out) :: varunit
      integer, intent(out) :: ierr
      integer, optional, intent(in) :: xstart, ystart, nx, ny
      integer :: varid
      ierr = nf90_inq_varid(ncid, trim(varname), varid)
      if (ierr /= 0) then
         varunit = 'unavailable'
         return
      endif
      ierr = nf90_get_att(ncid, varid, "units", varunit)
      if (ierr /= 0) then
         varunit = 'unkown'
         return
      endif
      if (present(xstart)) then
         ierr = nf90_get_var(ncid, varid, vardata, start=(/xstart,ystart/), count=(/nx,ny/))
      else
         ierr = nf90_get_var(ncid, varid, vardata)
      endif
    end subroutine netcdf_get
  end subroutine hrldas_netcdf_get_var_integer2d
end module hrldas_netcdf_test

program main
  use mpi
  use netcdf
  use module_mpp_land
  use hrldas_netcdf_test
  implicit none
  integer :: nx = 30
  integer :: ny = 26
  character(len=256) :: varname
  real, allocatable :: vardatar(:,:)
  integer, allocatable :: vardatai(:,:)
  character(len=256) :: varunit
  integer :: ncid, ierr
  call mpp_land_init()
  call mpp_land_partition(1, nx, ny, 1)
  allocate(vardatar(my_nx,my_ny))
  allocate(vardatai(my_nx,my_ny))

  if (my_id == io_id) then
     ierr = nf90_open('test_dongting.nc', NF90_NOWRITE, ncid)
  endif

  varname = 'XLAT'
  call hrldas_netcdf_get_var_real2d(ncid, varname, vardatar, varunit, ierr)
  print *, 'myid=', my_id, 'ix=', my_startx, 'iy=', my_starty, &
       & trim(varname), '=', vardatar(1,1:3), 'unit=',trim(varunit)

  call sleep(1)
  print *, ''
  call sleep(1)

  varname = 'IVGTYP'
  call hrldas_netcdf_get_var_integer2d(ncid, varname, vardatai, varunit, ierr)
  print *, 'myid=', my_id, 'ix=', my_startx, 'iy=', my_starty, &
       & trim(varname), '=', vardatai(2,2), 'unit=',trim(varunit)

  if (my_id == io_id) then
     ierr = nf90_close(ncid)
  end if

  deallocate(vardatar)
  deallocate(vardatai)
  call mpi_finalize(ierr)
end program main
