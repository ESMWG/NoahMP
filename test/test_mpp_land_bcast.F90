program main
  ! test the bcast function of MODULE_MPP_LAND
  use mpi
  use module_mpp_land
  implicit none
  integer :: nx, ny
  integer :: bufi(3)
  real :: bufr(3)
  integer :: i
  integer :: ierr
  nx = 101
  ny = 101
  call mpp_land_init()
  call mpp_land_partition(1, nx, ny, 1)

  if (my_id == io_id) then
     do i = 1, 3
        bufi(i) = 2 * i
        bufr(i) = 10.0 + i
     enddo
  endif
  print *, 'BEFORE BCAST INTEGER, myid=',my_id, 'buf=', bufi
  call mpp_land_bcast(bufi, 3, io_id)
  print *, 'AFTER BCAST INTEGER, myid=',my_id, 'buf=', bufi

  call sleep(5)

  print *, 'BEFORE BCAST REAL, myid=',my_id, 'buf=', bufr
  call mpp_land_bcast(bufr, 3, io_id)
  print *, 'AFTER BCAST real, myid=',my_id, 'buf=', bufr

  call mpi_finalize(ierr)
end program main
