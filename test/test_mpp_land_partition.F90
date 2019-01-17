program main
  ! test the partitioning function of MODULE_MPP_LAND
  use mpi
  use module_mpp_land
  implicit none
  integer :: nx, ny
  integer :: ierr
  nx = 101
  ny = 101
  call mpp_land_init()
  call mpp_land_partition(1, nx, ny, 1)

  call mpp_land_sync()

  print *, 'myid=', my_id,&
       & ';startx=' ,my_startx, ';nx=', my_nx,&
       & ';starty=', my_starty, ';ny=', my_ny

  call MPI_Finalize(ierr)
end program main
