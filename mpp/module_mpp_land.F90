! This is a module for parallel Land model.
module module_mpp_land
  use mpi
  use module_cpl_land
  implicit none

  integer, public :: my_id ! my process id
  integer, public :: io_id ! the IO process id. (the last process)
  integer, public :: left_id, right_id, up_id, down_id ! my neighbors
  integer, public :: nproc ! total processes, get by mpi initialization.
  integer, public :: nprocx ! total processes in the x direction
  integer, public :: nprocy ! total processes in the y direction
  integer, public :: iprocx ! the position of the current process in the 2d topography
  integer, public :: iprocy ! the position of the current process in the 2d topography
  integer, public :: global_nx, global_ny
  integer, public :: my_nx, my_ny
  integer, public :: global_rt_nx, global_rt_ny
  integer, public :: my_rtnx, my_rtny, rt_AGGFACTRT
  integer, public :: my_startx, my_starty

  integer :: mpp_status(MPI_STATUS_SIZE)

  integer :: overlap_n
  integer, allocatable, dimension(:), public :: iolocal_allnx, iolocal_allny
  integer, allocatable, dimension(:), public :: iolocal_allrtnx, iolocal_allrtny
  integer, allocatable, dimension(:), public :: iolocal_allstartx, iolocal_allstarty
  integer, allocatable, dimension(:), public :: mpp_nlinks

  interface mpp_land_bcast
   module procedure mpp_land_bcast_character
   module procedure mpp_land_bcast_integer
   module procedure mpp_land_bcast_integer_1
   module procedure mpp_land_bcast_real
   module procedure mpp_land_bcast_real_1
   module procedure mpp_land_bcast_real8
   module procedure mpp_land_bcast_real8_1
   module procedure mpp_land_bcast_logical_1
  end interface mpp_land_bcast

  interface check_land
     module procedure check_landreal1
     module procedure check_landreal1d
     module procedure check_landreal2d
     module procedure check_landreal3d
  end interface check_land

  interface write_io_land
     module procedure write_io_real3d
  end interface write_io_land

contains

  subroutine mpp_land_logicaltopo2d()
    ! setup the virtual 2D cartesian logical topography
    implicit none
    integer :: ndim, ierr
    integer, dimension(0:1) :: dims, coords
    logical :: cyclic(0:1) = (/.false., .false./) ! not cyclic
    logical :: reorder = .false.

    call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

    ! setup virtual cartesian grid topology
    call mpp_land_get_nprocsxy(nproc, nprocx, nprocy)
    ndim = 2
    dims(0) = nprocy ! rows
    dims(1) = nprocx ! columns
    call MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, &
         cyclic, reorder, cartGridComm, ierr)

#ifdef HYDRO_D
    if(my_id == IO_id) then
       write(*,*) ""
       write(*,*) "total process:", nproc
       write(*,*) "nprocx =", nprocx, "nprocy=",nprocy
    end if
#endif

    ! get the row and column of the current process in the logical topography.
    ! left --> right, 0 -->nprocx -1
    ! up --> down, 0 --> nprocy -1
    iprocx = mod(my_id, nprocx)
    iprocy = my_id / nprocx

    ! get the neighbors. -1 means no neighbor.
    if (iprocy > 0) then
       down_id = my_id - nprocx
    else
       down_id = -1
    end if
    if (iprocy < nprocy - 1) then
       up_id = my_id + nprocx
    else
       up_id = -1
    end if

    if (iprocx > 0) then
       left_id = my_id - 1
    else
       left_id = -1
    end if
    if (iprocx < nprocx - 1) then
       right_id = my_id + 1
    else
       right_id = -1
    end if

    ! the IO node is the last processor.
    IO_id = nproc - 1

    ! duplicates, np_up_down=nprocy, np_left_right=nprocx
    ! should be removed?
    call MPI_CART_GET(cartGridComm, 2, dims, cyclic, coords, ierr)
    p_up_down = coords(0)
    p_left_right = coords(1)
    np_up_down = nprocy
    np_left_right = nprocx

    call mpp_land_sync()
  end subroutine mpp_land_logicaltopo2d


  subroutine mpp_land_get_nprocsxy(nproc, nx, ny)
   ! calculate the nx and ny based on the total nproc.
   implicit none
   integer, intent(in) :: nproc
   integer, intent(out) :: nx, ny
   integer :: i, j, max
   max = nproc
   do j = 1, nproc
      if(mod(nproc,j) .eq. 0) then
         i = nproc / j
         if(abs(i-j) .lt. max) then
            max = abs(i-j)
            nx = i
            ny = j
         end if
      end if
   end do
 end subroutine mpp_land_get_nprocsxy


  subroutine mpp_land_init()
    ! initialize the parallel environment and the 2D topography
    implicit none
    integer :: ierr
    logical :: mpi_inited

    call mpi_initialized(mpi_inited, ierr)
    if (.not. mpi_inited) then
       call mpi_init(ierr)  ! stand alone land model.
       if (ierr /= MPI_SUCCESS) then
          stop 'ERROR[mpp_land_init]: unable to intialize MPI'
       endif
    endif

    ! create 2d logical mapping of the CPU.
    call mpp_land_logicaltopo2d()
  end subroutine mpp_land_init


  subroutine mpp_land_partition(overlap, in_global_nx, in_global_ny, AGGFACTRT)
    ! partition the land grid cells according to the 2d logical MPI topography
    ! assign iolocal_allnx, iolocal_allny, iolocal_allstartx, iolocal_allstarty
    ! assign iolocal_allrtnx, iolocal_allrtny, iolocal_all
    implicit none
    integer, intent(in) :: overlap ! the overlaped grid number (unused now)
    integer, intent(in) :: in_global_nx, in_global_ny
    integer, intent(in) :: AGGFACTRT
    integer :: i
    integer :: tag
    integer :: ierr
    integer :: buf(4)

    global_nx = in_global_nx
    global_ny = in_global_ny
    rt_AGGFACTRT = AGGFACTRT
    global_rt_nx = global_nx * AGGFACTRT
    global_rt_ny = global_ny * AGGFACTRT

    tag = 1
    if (my_id == io_id) then
       call mpp_land_partition_calc()
       do i = 0, nproc - 1
          if(i /= my_id) then
             !send to the other nodes
             buf(1) = iolocal_allnx(i+1)
             buf(2) = iolocal_allny(i+1)
             buf(3) = iolocal_allstartx(i+1)
             buf(4) = iolocal_allstarty(i+1)
             call mpi_send(buf, size(buf), MPI_INTEGER, i, &
                  & tag, MPI_COMM_WORLD, ierr)
          else
             my_nx = iolocal_allnx(i+1)
             my_ny = iolocal_allny(i+1)
             my_startx = iolocal_allstartx(i+1)
             my_starty = iolocal_allstarty(i+1)
          end if
       end do
    else
       ! receive from the io node
       call mpi_recv(buf, size(buf), MPI_INTEGER, io_id, &
            & tag, MPI_COMM_WORLD, mpp_status, ierr)
       my_nx = buf(1)
       my_ny = buf(2)
       my_startx = buf(3)
       my_starty = buf(4)
    end if

    my_rtnx = my_nx * AGGFACTRT + 2
    my_rtny = my_ny * AGGFACTRT + 2
    if(left_id < 0) my_rtnx = my_rtnx -1
    if(right_id < 0) my_rtnx = my_rtnx -1
    if(up_id < 0) my_rtny = my_rtny -1
    if(down_id < 0) my_rtny = my_rtny -1

#ifdef HYDRO_D
    write(*,*) "my_id=", my_id, "global_rt_nx=", global_rt_nx
    write(*,*) "my_id=", my_id, "global_rt_nx=", global_rt_ny
    write(*,*) "my_id=", my_id, "global_nx=", global_nx
    write(*,*) "my_id=", my_id, "global_nx=", global_ny
#endif
  end subroutine mpp_land_partition


  subroutine mpp_land_partition_calc()
    ! calculate the land grid cell partitioning
    implicit none
    integer :: base_nx
    integer :: base_ny
    integer :: id
    integer, allocatable :: id2d(:,:)
    integer :: iprocx, iprocy
    integer :: startx, starty

    if (my_id /= io_id) return

    allocate(iolocal_allnx(nproc))
    allocate(iolocal_allny(nproc))
    allocate(iolocal_allstartx(nproc))
    allocate(iolocal_allstarty(nproc))

    ! nx and ny
    base_nx = int(global_nx / nprocx)
    base_ny = int(global_ny / nprocy)
    do id = 0, nproc-1
       iprocx = mod(id, nprocx)
       iprocy = id / nprocx
       if (iprocx <= mod(global_nx, nprocx)-1) then
          iolocal_allnx(id+1) = base_nx + 1
       else
          iolocal_allnx(id+1) = base_nx
       end if
       if (iprocy <= mod(global_ny, nprocy)-1) then
          iolocal_allny(id+1) = base_ny + 1
       else
          iolocal_allny(id+1) = base_ny
       end if
    end do

    ! startx and starty
    allocate(id2d(nprocx,nprocy))
    do id = 0, nproc-1
       iprocx = mod(id, nprocx)
       iprocy = id / nprocx
       id2d(iprocx+1,iprocy+1) = id
    end do
    do iprocx = 0, nprocx-1
       if (iprocx == 0) then
          startx = 1
       else
          startx = startx + iolocal_allnx(id2d(iprocx,1)+1)
       end if
       do iprocy = 0, nprocy-1
          iolocal_allstartx(id2d(iprocx+1,iprocy+1)+1) = startx
       end do
    end do
    do iprocy = 0, nprocy-1
       if (iprocy == 0) then
          starty = 1
       else
          starty = starty + iolocal_allny(id2d(1,iprocy)+1)
       end if
       do iprocx = 0, nprocx-1
          iolocal_allstarty(id2d(iprocx+1,iprocy+1)+1) = starty
       end do
    end do

#ifdef HYDRO_D
    write(*,*) 'MPP_LAND_PARTITION_CALC in the IO node:'
    write(*,*) 'NODE ID:'
    do iprocy = nprocy-1, 0, -1
       do iprocx = 0, nprocx-1
          id = id2d(iprocx+1,iprocy+1)
          write(*,'(" (", I9, ") ")',advance='NO') id
       end do
       write(*,*)
    end do
    write(*,*) 'NX/NY:'
    do iprocy = nprocy-1, 0, -1
       do iprocx = 0, nprocx-1
          id = id2d(iprocx+1,iprocy+1)
          write(*,'(" (",I4, "/", I4,") ")',advance='NO') &
               & iolocal_allnx(id+1), &
               & iolocal_allny(id+1)
       end do
       write(*,*)
    end do
    write(*,*) 'STARTX/STOPX:'
    do iprocy = nprocy-1, 0, -1
       do iprocx = 0, nprocx-1
          id = id2d(iprocx+1,iprocy+1)
          write(*,'(" (",I4, "/", I4,") ")',advance='NO') &
               & iolocal_allstartx(id+1), &
               & iolocal_allstartx(id+1)+iolocal_allnx(id+1)-1
       end do
       write(*,*)
    end do
    write(*,*) 'STARTY/STOPY:'
    do iprocy = nprocy-1, 0, -1
       do iprocx = 0, nprocx-1
          id = id2d(iprocx+1,iprocy+1)
          write(*,'(" (",I4, "/", I4,") ")',advance='NO') &
               & iolocal_allstarty(id+1), &
               & iolocal_allstarty(id+1)+iolocal_allny(id+1)-1
       end do
       write(*,*)
    end do
#endif
    deallocate(id2d)
  end subroutine mpp_land_partition_calc


  subroutine mpp_land_comlr_real(in_out_data, nx, ny, flag)
    ! communicate message on left right direction.
    implicit none
    integer, intent(in) :: nx, ny
    real, intent(inout) :: in_out_data(nx,ny)
    real :: data_r(2,ny)
    integer :: data_size, tag, ierr
    integer :: flag ! 99 replace the boundary, else get the sum

    if(flag == 99) then ! replace the data
       if(right_id >= 0) then ! send to right first
          tag = 11
          data_size = ny
          call mpi_send(in_out_data(nx-1,:), data_size, MPI_REAL, &
               right_id, tag, MPI_COMM_WORLD, ierr)
       end if
       if(left_id >= 0) then ! receive from left
          tag = 11
          data_size = ny
          call mpi_recv(in_out_data(1,:), data_size, MPI_REAL, &
               left_id, tag, MPI_COMM_WORLD, mpp_status, ierr)
       endif

       if(left_id >= 0 ) then ! send to left second
          tag = 21
          data_size = ny
          call mpi_send(in_out_data(2,:), data_size, MPI_REAL, &
               left_id,tag,MPI_COMM_WORLD,ierr)
       endif
       if(right_id >= 0) then ! receive from right
          tag = 21
          data_size = ny
          call mpi_recv(in_out_data(nx,:), data_size, MPI_REAL, &
               right_id, tag, MPI_COMM_WORLD, mpp_status, ierr)
       endif
    else ! get the sum
       if(right_id >= 0) then ! send to right first
          tag = 11
          data_size = 2 * ny
          call mpi_send(in_out_data(nx-1:nx,:), data_size, MPI_REAL, &
               right_id, tag, MPI_COMM_WORLD, ierr)
       end if
       if(left_id >= 0) then ! receive from left
          tag = 11
          data_size = 2 * ny
          call mpi_recv(data_r, data_size, MPI_REAL, left_id, tag, &
               MPI_COMM_WORLD, mpp_status, ierr)
          in_out_data(1,:) = in_out_data(1,:) + data_r(1,:)
          in_out_data(2,:) = in_out_data(2,:) + data_r(2,:)
       endif

       if(left_id >= 0) then ! send to left second
          tag = 21
          data_size = 2 * ny
          call mpi_send(in_out_data(1:2,:), data_size, MPI_REAL, &
               left_id, tag, MPI_COMM_WORLD, ierr)
       endif
       if(right_id >= 0) then ! receive from right
          tag = 21
          data_size = 2 * ny
          call mpi_recv(in_out_data(nx-1:nx,:), data_size, MPI_REAL, &
               right_id, tag, MPI_COMM_WORLD, mpp_status, ierr)
       endif
    endif ! end if for flag
  end subroutine mpp_land_comlr_real


  subroutine mpp_land_comlr_real8(in_out_data, nx, ny, flag)
    ! Communicate message on left right direction.
    implicit none
    integer, intent(in) :: nx, ny
    real*8, intent(inout) :: in_out_data(nx,ny)
    real*8 :: data_r(2,ny)
    integer :: data_size, tag, ierr
    integer :: flag ! 99 replace the boundary, else get the sum

    if(flag == 99) then ! replace the data
       if(right_id >= 0) then ! send to right first
          tag = 11
          data_size = ny
          call mpi_send(in_out_data(nx-1,:), data_size, MPI_DOUBLE_PRECISION, &
               right_id, tag, MPI_COMM_WORLD, ierr)
       end if
       if(left_id >= 0) then ! receive from left
          tag = 11
          data_size = ny
          call mpi_recv(in_out_data(1,:), data_size, MPI_DOUBLE_PRECISION, &
               left_id, tag, MPI_COMM_WORLD, mpp_status, ierr)
       endif

       if(left_id >= 0) then ! send to left second
          tag = 21
          data_size = ny
          call mpi_send(in_out_data(2,:), data_size, MPI_DOUBLE_PRECISION, &
               left_id, tag, MPI_COMM_WORLD, ierr)
       endif
       if(right_id >= 0) then ! receive from right
          tag = 21
          data_size = ny
          call mpi_recv(in_out_data(nx,:), data_size, MPI_DOUBLE_PRECISION, &
               right_id, tag, MPI_COMM_WORLD, mpp_status, ierr)
       endif
    else ! get the sum
       if(right_id >= 0) then ! send to right first
          tag = 11
          data_size = 2 * ny
          call mpi_send(in_out_data(nx-1:nx,:), data_size, MPI_DOUBLE_PRECISION, &
               right_id, tag, MPI_COMM_WORLD, ierr)
       end if
       if(left_id >= 0) then ! receive from left
          tag = 11
          data_size = 2 * ny
          call mpi_recv(data_r, data_size, MPI_DOUBLE_PRECISION, left_id,tag, &
               MPI_COMM_WORLD, mpp_status, ierr)
          in_out_data(1,:) = in_out_data(1,:) + data_r(1,:)
          in_out_data(2,:) = in_out_data(2,:) + data_r(2,:)
       endif

       if(left_id >= 0) then ! send to left second
          tag = 21
          data_size = 2 * ny
          call mpi_send(in_out_data(1:2,:), data_size, MPI_DOUBLE_PRECISION, &
               left_id, tag, MPI_COMM_WORLD, ierr)
       endif
       if(right_id >= 0) then ! receive from right
          tag = 21
          data_size = 2 * ny
          call mpi_recv(in_out_data(nx-1:nx,:), data_size, MPI_DOUBLE_PRECISION, &
               right_id, tag, MPI_COMM_WORLD, mpp_status,ierr)
       endif
    endif ! end if black for flag
  end subroutine mpp_land_comlr_real8


  subroutine mpp_land_comub_real(in_out_data, nx, ny, flag)
    ! communicate message on up down direction.
    implicit none
    integer, intent(in) :: nx, ny
    real, intent(inout) :: in_out_data(nx,ny)
    real :: data_r(nx,2)
    integer :: data_size, tag, ierr
    integer :: flag ! 99 replace the boundary, else get the sum of the boundary

    if(flag == 99) then ! replace the boundary data
       if(up_id >= 0) then ! send to up first
          tag = 31
          data_size = nx
          call mpi_send(in_out_data(:,ny-1), data_size, MPI_REAL, &
               up_id, tag, MPI_COMM_WORLD, ierr)
       endif
       if(down_id >= 0) then ! receive from down
          tag = 31
          data_size = nx
          call mpi_recv(in_out_data(:,1), data_size, MPI_REAL, &
               down_id, tag, MPI_COMM_WORLD, mpp_status, ierr)
       endif

       if(down_id >= 0) then ! send down
          tag = 41
          data_size = nx
          call mpi_send(in_out_data(:,2), data_size, MPI_REAL, &
               down_id, tag, MPI_COMM_WORLD, ierr)
       endif
       if(up_id >= 0) then ! receive from upper
          tag = 41
          data_size = nx
          call mpi_recv(in_out_data(:,ny), data_size, MPI_REAL, &
               up_id, tag, MPI_COMM_WORLD, mpp_status, ierr)
       endif
    else  ! flag = 1
       if(up_id >= 0) then ! send to up first
          tag = 31
          data_size = nx * 2
          call mpi_send(in_out_data(:,ny-1:ny), data_size, MPI_REAL, &
               up_id, tag, MPI_COMM_WORLD, ierr)
       endif
       if(down_id >= 0) then ! receive from down
          tag = 31
          data_size = nx * 2
          call mpi_recv(data_r, data_size, MPI_REAL, &
               down_id, tag, MPI_COMM_WORLD, mpp_status, ierr)
          in_out_data(:,1) = in_out_data(:,1) + data_r(:,1)
          in_out_data(:,2) = in_out_data(:,2) + data_r(:,2)
       endif

       if(down_id >= 0) then ! send down
          tag = 41
          data_size = nx * 2
          call mpi_send(in_out_data(:,1:2), data_size, MPI_REAL, &
               down_id, tag, MPI_COMM_WORLD, ierr)
       endif
       if(up_id >= 0) then ! receive from upper
          tag = 41
          data_size = nx * 2
          call mpi_recv(in_out_data(:,ny-1:ny), data_size, MPI_REAL, &
               up_id, tag, MPI_COMM_WORLD, mpp_status, ierr)
       endif
    endif ! end of block flag
  end subroutine mpp_land_comub_real


  subroutine mpp_land_comub_real8(in_out_data, nx, ny, flag)
    ! communicate message on up down direction.
    implicit none
    integer, intent(in) :: nx, ny
    real*8, intent(inout) :: in_out_data(nx,ny)
    real*8 :: data_r(nx,2)
    integer :: data_size, tag, ierr
    integer :: flag ! 99 replace the boundary, else get the sum of the boundary

    if(flag == 99) then ! replace the boundary data
       if(up_id >= 0) then ! send to up first
          tag = 31
          data_size = nx
          call mpi_send(in_out_data(:,ny-1), data_size, MPI_DOUBLE_PRECISION,   &
               up_id, tag, MPI_COMM_WORLD, ierr)
       endif
       if(down_id >= 0) then ! receive from down
          tag = 31
          data_size = nx
          call mpi_recv(in_out_data(:,1), data_size, MPI_DOUBLE_PRECISION, &
               down_id, tag, MPI_COMM_WORLD, mpp_status,ierr)
       endif

       if(down_id >= 0) then ! send down
          tag = 41
          data_size = nx
          call mpi_send(in_out_data(:,2), data_size, MPI_DOUBLE_PRECISION, &
               down_id, tag, MPI_COMM_WORLD, ierr)
       endif
       if(up_id >= 0) then ! receive from upper
          tag = 41
          data_size = nx
          call mpi_recv(in_out_data(:,ny), data_size, MPI_DOUBLE_PRECISION, &
               up_id, tag, MPI_COMM_WORLD, mpp_status, ierr)
       endif
    else  ! flag = 1
       if(up_id >= 0) then ! send to up first
          tag = 31
          data_size = nx * 2
          call mpi_send(in_out_data(:,ny-1:ny), data_size, MPI_DOUBLE_PRECISION, &
               up_id, tag, MPI_COMM_WORLD, ierr)
       endif
       if(down_id >= 0) then ! receive from down
          tag = 31
          data_size = nx * 2
          call mpi_recv(data_r, data_size, MPI_DOUBLE_PRECISION, &
               down_id, tag, MPI_COMM_WORLD, mpp_status, ierr)
          in_out_data(:,1) = in_out_data(:,1) + data_r(:,1)
          in_out_data(:,2) = in_out_data(:,2) + data_r(:,2)
       endif

       if(down_id >= 0) then ! send down
          tag = 41
          data_size = nx * 2
          call mpi_send(in_out_data(:,1:2), data_size, MPI_DOUBLE_PRECISION, &
               down_id, tag, MPI_COMM_WORLD, ierr)
       endif
       if(up_id >= 0) then ! receive from upper
          tag = 41
          data_size = nx * 2
          call mpi_recv(in_out_data(:,ny-1:ny), data_size, MPI_DOUBLE_PRECISION, &
               up_id, tag, MPI_COMM_WORLD, mpp_status, ierr)
       endif
    endif ! end of block flag
  end subroutine mpp_land_comub_real8


  subroutine mpp_land_com_real(in_out_data, nx, ny, flag)
    ! Communicate message on left right and up bottom directions.
    implicit none
    integer, intent(in) :: nx, ny
    integer, intent(in) :: flag != 99  test only for land model. (replace the boundary).
    != 1   get the sum of the boundary value.
    real, intent(inout) :: in_out_data(nx,ny)

    call mpp_land_comlr_real(in_out_data, nx, ny, flag)
    call mpp_land_comub_real(in_out_data, nx, ny, flag)
  end subroutine mpp_land_com_real


  subroutine mpp_land_com_real8(in_out_data ,nx, ny, flag)
    ! Communicate message on left right and up bottom directions.
    implicit none
    integer, intent(in) :: nx, ny
    integer, intent(in) :: flag != 99  test only for land model. (replace the boundary).
    != 1   get the sum of the boundary value.
    real*8, intent(inout) :: in_out_data(nx,ny)

    call mpp_land_comlr_real8(in_out_data, nx, ny, flag)
    call mpp_land_comub_real8(in_out_data, nx, ny, flag)
  end subroutine mpp_land_com_real8


  subroutine mpp_land_com_integer(data, nx, ny, flag)
    ! Communicate message on left right and up bottom directions.
    implicit none
    integer, intent(in) :: nx, ny
    integer, intent(in) :: flag != 99  test only for land model. (replace the boundary).
    != 1   get the sum of the boundary value.
    integer, intent(inout) :: data(nx,ny)
    real :: in_out_data(nx,ny)

    in_out_data = data + 0.0
    call mpp_land_comlr_real(in_out_data, NX, NY, flag)
    call mpp_land_comub_real(in_out_data, NX, NY, flag)
    data = in_out_data + 0
  end subroutine mpp_land_com_integer


  subroutine decompose_data_real3d(in_buff, out_buff, klevel)
    implicit none
    integer, intent(in) :: klevel
    real, intent(in) :: in_buff(:,:,:)
    real, intent(out) :: out_buff(:,:,:)
    integer :: k
    do k = 1, klevel
       call decompose_data_real(in_buff(:,k,:), out_buff(:,k,:))
    end do
  end subroutine decompose_data_real3d


  subroutine decompose_data_real(in_buff, out_buff)
    ! usage: all of the cpu call this subroutine.
    ! the IO node will distribute the data to rest of the node.
    implicit none
    real, intent(in) :: in_buff(:,:)
    real, intent(out) :: out_buff(:,:)
    integer :: tag, i, ierr, data_size
    integer :: ibegin, iend, jbegin, jend

    tag = 2
    if(my_id == io_id) then
       do i = 0, nproc - 1
          ibegin = iolocal_allstartx(i+1)
          iend   = iolocal_allstartx(i+1) + iolocal_allnx(i+1) -1
          jbegin = iolocal_allstarty(i+1)
          jend   = iolocal_allstarty(i+1) + iolocal_allny(i+1) -1

          if(i == my_id) then
             out_buff = in_buff(ibegin:iend,jbegin:jend)
          else
             ! send data to the rest process.
             data_size = iolocal_allnx(i+1) * iolocal_allny(i+1)
             call mpi_send(in_buff(ibegin:iend,jbegin:jend), data_size, &
                  MPI_REAL, i, tag, MPI_COMM_WORLD, ierr)
          end if
       end do
    else
       data_size = my_nx * my_ny
       call mpi_recv(out_buff, data_size, MPI_REAL, io_id, &
            tag ,MPI_COMM_WORLD, mpp_status, ierr)
    end if
  end subroutine decompose_data_real


  subroutine decompose_data_int(in_buff, out_buff)
    ! usage: all of the cpu call this subroutine.
    ! the IO node will distribute the data to rest of the node.
    implicit none
    integer, intent(in) :: in_buff(:,:)
    integer, intent(out) :: out_buff(:,:)
    integer :: tag, i, ierr, data_size
    integer :: ibegin, iend, jbegin, jend

    tag = 2
    if(my_id == io_id) then
       do i = 0, nproc - 1
          ibegin = iolocal_allstartx(i+1)
          iend   = iolocal_allstartx(i+1) + iolocal_allnx(i+1) -1
          jbegin = iolocal_allstarty(i+1)
          jend   = iolocal_allstarty(i+1) + iolocal_allny(i+1) -1
          if(i == my_id) then
             out_buff = in_buff(ibegin:iend,jbegin:jend)
          else
             ! send data to the rest process
             data_size = iolocal_allnx(i+1) * iolocal_allny(i+1)
             call mpi_send(in_buff(ibegin:iend,jbegin:jend), data_size, &
                  MPI_INTEGER, i, tag, MPI_COMM_WORLD, ierr)
          end if
       end do
    else
       data_size = my_nx * my_ny
       call mpi_recv(out_buff, data_size, MPI_INTEGER, io_id, &
            tag, MPI_COMM_WORLD, mpp_status, ierr)
    end if
  end subroutine decompose_data_int


  subroutine write_io_int(in_buff, out_buff)
    ! the IO node will receive the data from the rest process.
    implicit none
    integer, intent(in) :: in_buff(:,:)
    integer, intent(out) :: out_buff(:,:)
    integer :: tag, i, ierr, data_size
    integer :: ibegin,iend,jbegin,jend
    tag = 2
    if(my_id /= io_id) then
       data_size = my_nx * my_ny
       call mpi_send(in_buff, data_size, MPI_INTEGER, io_id, &
            tag, MPI_COMM_WORLD, ierr)
    else
       do i = 0, nproc - 1
          ibegin = iolocal_allstartx(i+1)
          iend   = iolocal_allstartx(i+1) + iolocal_allnx(i+1) -1
          jbegin = iolocal_allstarty(i+1)
          jend   = iolocal_allstarty(i+1) + iolocal_allny(i+1) -1
          if(i == io_id) then
             out_buff(ibegin:iend,jbegin:jend) = in_buff
          else
             data_size = iolocal_allnx(i+1) * iolocal_allny(i+1)
             call mpi_recv(out_buff(ibegin:iend,jbegin:jend), data_size, &
                  MPI_INTEGER, i, tag, MPI_COMM_WORLD, mpp_status, ierr)
          end if
       end do
    end if
  end subroutine write_io_int


  subroutine write_io_real3d(in_buff, out_buff, klevel)
    ! the IO node will receive the data from the rest process.
    implicit none
    integer, intent(in) :: klevel
    real, intent(in) :: in_buff(:,:,:)
    real, intent(out) :: out_buff(:,:,:)
    integer :: k
    do k = 1, klevel
       call write_io_real(in_buff(:,k,:), out_buff(:,k,:))
    end do
  end subroutine write_io_real3d


  subroutine write_io_real(in_buff, out_buff)
    ! the IO node will receive the data from the rest process.
    implicit none
    real, intent(in) :: in_buff(:,:)
    real, intent(out) :: out_buff(:,:)
    integer :: tag, i, ierr, data_size
    integer :: ibegin, iend, jbegin, jend
    tag = 2
    if(my_id /= io_id) then
       data_size = my_nx * my_ny
       call mpi_send(in_buff, data_size, MPI_REAL, io_id, &
            tag, MPI_COMM_WORLD, ierr)
    else
       do i = 0, nproc - 1
          ibegin = iolocal_allstartx(i+1)
          iend   = iolocal_allstartx(i+1) + iolocal_allnx(i+1) -1
          jbegin = iolocal_allstarty(i+1)
          jend   = iolocal_allstarty(i+1) + iolocal_allny(i+1) -1
          if(i == io_id) then
             out_buff(ibegin:iend,jbegin:jend) = in_buff
          else
             data_size = iolocal_allnx(i+1)*iolocal_allny(i+1)
             call mpi_recv(out_buff(ibegin:iend,jbegin:jend), data_size, &
                  MPI_REAL, i, tag, MPI_COMM_WORLD, mpp_status, ierr)
          end if
       end do
    end if
  end subroutine write_io_real

  subroutine write_io_rt_real(in_buff, out_buff)
    ! the IO node will receive the data from the rest process.
    implicit none
    real, intent(in) :: in_buff(:,:)
    real, intent(out) :: out_buff(:,:)
    integer :: tag, i, ierr, data_size
    integer :: ibegin, iend, jbegin, jend
    tag = 2
    if(my_id /= io_id) then
       data_size = my_rtnx * my_rtny
       call mpi_send(in_buff, data_size, MPI_REAL, io_id, &
            tag, MPI_COMM_WORLD, ierr)
    else
       do i = 0, nproc - 1
          ibegin = iolocal_allstartx(i+1) * rt_AGGFACTRT - (rt_AGGFACTRT-1)
          if(ibegin > 1) ibegin = ibegin - 1
          iend   = ibegin + iolocal_allrtnx(i+1) -1
          jbegin = iolocal_allstarty(i+1) * rt_AGGFACTRT - (rt_AGGFACTRT-1)
          if(jbegin > 1) jbegin = jbegin - 1
          jend   = jbegin + iolocal_allrtny(i+1) -1
          if(i == io_id) then
             out_buff(ibegin:iend,jbegin:jend) = in_buff
          else
             data_size = iolocal_allrtnx(i+1) * iolocal_allrtny(i+1)
             call mpi_recv(out_buff(ibegin:iend,jbegin:jend), data_size, &
                  MPI_REAL, i, tag, MPI_COMM_WORLD, mpp_status, ierr)
          end if
       end do
    end if
  end subroutine write_io_rt_real


  subroutine write_io_rt_int(in_buff, out_buff)
    ! the IO node will receive the data from the rest process.
    implicit none
    integer, intent(in), dimension(:,:) :: in_buff
    integer, intent(out), dimension(:,:) :: out_buff
    integer :: tag, i, ierr, data_size
    integer :: ibegin, iend, jbegin, jend
    tag = 2
    if(my_id /= io_id) then
       data_size = my_rtnx * my_rtny
       call mpi_send(in_buff, data_size, MPI_INTEGER, io_id, &
            tag,MPI_COMM_WORLD, ierr)
    else
       do i = 0, nproc - 1
          ibegin = iolocal_allstartx(i+1) * rt_AGGFACTRT - (rt_AGGFACTRT-1)
          if(ibegin > 1) ibegin = ibegin - 1
          iend   = ibegin + iolocal_allrtnx(i+1) -1
          jbegin = iolocal_allstarty(i+1) * rt_AGGFACTRT - (rt_AGGFACTRT - 1)
          if(jbegin > 1) jbegin = jbegin - 1
          jend   = jbegin + iolocal_allrtny(i+1) -1
          if(i == io_id) then
             out_buff(ibegin:iend,jbegin:jend) = in_buff
          else
             data_size = iolocal_allrtnx(i+1) * iolocal_allrtny(i+1)
             call mpi_recv(out_buff(ibegin:iend,jbegin:jend), data_size, &
                  MPI_INTEGER, i, tag, MPI_COMM_WORLD, mpp_status, ierr)
          end if
       end do
    end if
  end subroutine write_io_rt_int


  subroutine mpp_land_abort()
    implicit none
    integer :: ierr
    call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
  end subroutine mpp_land_abort


  subroutine mpp_land_sync()
    implicit none
    integer :: ierr
    call MPI_barrier(MPI_COMM_WORLD ,ierr)
    if(ierr .ne. 0) call mpp_land_abort()
  end subroutine mpp_land_sync


  subroutine mpp_land_bcast_integer(buffer, count, root)
   implicit none
   integer, dimension(:) :: buffer
   integer, intent(in) :: count
   integer, intent(in) :: root
   integer :: ierr
   call mpi_bcast(buffer, count, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
   if (ierr /= MPI_SUCCESS) call mpi_abort(MPI_COMM_WORLD, 1, ierr)
 end subroutine mpp_land_bcast_integer

 subroutine mpp_land_bcast_integer_1(buffer, count, root)
   implicit none
   integer :: buffer
   integer, intent(in) :: count
   integer, intent(in) :: root
   integer :: ierr
   call mpi_bcast(buffer, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
   if (ierr /= MPI_SUCCESS) call mpi_abort(MPI_COMM_WORLD, 1, ierr)
 end subroutine mpp_land_bcast_integer_1

 subroutine mpp_land_bcast_real(buffer, count, root)
   implicit none
   real, dimension(:) :: buffer
   integer, intent(in) :: count
   integer, intent(in) :: root
   integer :: ierr
   call mpi_bcast(buffer, count, MPI_REAL, root, MPI_COMM_WORLD, ierr)
   if (ierr /= MPI_SUCCESS) call mpi_abort(MPI_COMM_WORLD, 1, ierr)
 end subroutine mpp_land_bcast_real

 subroutine mpp_land_bcast_real_1(buffer, count, root)
   implicit none
   real :: buffer
   integer, intent(in) :: count
   integer, intent(in) :: root
   integer :: ierr
   call mpi_bcast(buffer, 1, MPI_REAL, root, MPI_COMM_WORLD, ierr)
   if (ierr /= MPI_SUCCESS) call mpi_abort(MPI_COMM_WORLD, 1, ierr)
 end subroutine mpp_land_bcast_real_1

 subroutine mpp_land_bcast_real8(buffer, count, root)
   implicit none
   real*8, dimension(:) :: buffer
   integer, intent(in) :: count
   integer, intent(in) :: root
   integer :: ierr
   call mpi_bcast(buffer, count, MPI_REAL8, root, MPI_COMM_WORLD, ierr)
   if (ierr /= MPI_SUCCESS) call mpi_abort(MPI_COMM_WORLD, 1, ierr)
 end subroutine mpp_land_bcast_real8

 subroutine mpp_land_bcast_real8_1(buffer, count, root)
   implicit none
   real*8 :: buffer
   integer, intent(in) :: count
   integer, intent(in) :: root
   integer :: ierr
   call mpi_bcast(buffer, 1, MPI_REAL8, root, MPI_COMM_WORLD, ierr)
   if (ierr /= MPI_SUCCESS) call mpi_abort(MPI_COMM_WORLD, 1, ierr)
 end subroutine mpp_land_bcast_real8_1

 subroutine mpp_land_bcast_character(buffer, count, root)
   implicit none
   character(len=*) :: buffer
   integer, intent(in) :: count
   integer, intent(in) :: root
   integer :: ierr
   call mpi_bcast(buffer, count, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
   if (ierr /= MPI_SUCCESS) call mpi_abort(MPI_COMM_WORLD, 1, ierr)
 end subroutine mpp_land_bcast_character

 subroutine mpp_land_bcast_logical_1(buffer, count, root)
   implicit none
   logical :: buffer
   integer, intent(in) :: count
   integer, intent(in) :: root
   integer :: ierr
   call mpi_bcast(buffer, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
   if (ierr /= MPI_SUCCESS) call mpi_abort(MPI_COMM_WORLD, 1, ierr)
 end subroutine mpp_land_bcast_logical_1

  subroutine decompose_rt_real(in_buff, out_buff, g_nx, g_ny, nx, ny)
    ! usage: all of the cpu call this subroutine.
    ! the IO node will distribute the data to rest of the node.
    implicit none
    integer :: g_nx, g_ny, nx, ny
    real, intent(in) :: in_buff(:,:)
    real, intent(out) :: out_buff(:,:)
    integer :: tag, i, ierr, data_size
    integer :: ibegin, iend, jbegin, jend

    tag = 2
    if(my_id == io_id) then
       do i = 0, nproc - 1
          ibegin = iolocal_allstartx(i+1) * rt_AGGFACTRT - (rt_AGGFACTRT-1)
          if(ibegin .gt. 1) ibegin = ibegin - 1
          iend = ibegin + iolocal_allrtnx(i+1) -1
          jbegin = iolocal_allstarty(i+1) * rt_AGGFACTRT - (rt_AGGFACTRT-1)
          if(jbegin .gt. 1) jbegin = jbegin - 1
          jend = jbegin + iolocal_allrtny(i+1) -1

          if(i == my_id) then
             out_buff = in_buff(ibegin:iend,jbegin:jend)
          else
             ! send data to the rest process.
             data_size = (iend-ibegin+1) * (jend-jbegin+1)
             call mpi_send(in_buff(ibegin:iend,jbegin:jend), data_size, &
                  MPI_REAL, i, tag, MPI_COMM_WORLD, ierr)
          end if
       end do
    else
       data_size = nx * ny
       call mpi_recv(out_buff, data_size, MPI_REAL, io_id, &
            tag, MPI_COMM_WORLD, mpp_status, ierr)
    end if
  end subroutine decompose_rt_real


  subroutine decompose_rt_int(in_buff, out_buff, g_nx, g_ny, nx, ny)
    ! usage: all of the cpu call this subroutine.
    ! the IO node will distribute the data to rest of the node.
    implicit none
    integer :: g_nx, g_ny, nx, ny
    integer, intent(in) :: in_buff(:,:)
    integer, intent(out) :: out_buff(:,:)
    integer :: tag, i, ierr, data_size
    integer :: ibegin, iend, jbegin, jend

    tag = 2
    call mpp_land_sync()
    if(my_id == io_id) then
       do i = 0, nproc - 1
          ibegin = iolocal_allstartx(i+1) * rt_AGGFACTRT - (rt_AGGFACTRT-1)
          if(ibegin .gt. 1) ibegin = ibegin - 1
          iend = ibegin + iolocal_allrtnx(i+1) -1
          jbegin = iolocal_allstarty(i+1) * rt_AGGFACTRT - (rt_AGGFACTRT-1)
          if(jbegin .gt. 1) jbegin = jbegin - 1
          jend = jbegin + iolocal_allrtny(i+1) -1

          if(i == my_id) then
             out_buff = in_buff(ibegin:iend,jbegin:jend)
          else
             ! send data to the rest process.
             data_size = (iend-ibegin+1) * (jend-jbegin+1)
             call mpi_send(in_buff(ibegin:iend,jbegin:jend), data_size, &
                  MPI_INTEGER, i, tag, MPI_COMM_WORLD, ierr)
          end if
       end do
    else
       data_size = nx * ny
       call mpi_recv(out_buff, data_size, MPI_INTEGER, io_id, &
            tag, MPI_COMM_WORLD, mpp_status, ierr)
    end if
  end subroutine decompose_rt_int


  subroutine wrf_land_set_init(info, total_pe, AGGFACTRT)
    implicit none
    integer :: total_pe
    integer :: info(9,total_pe), AGGFACTRT
    integer :: ierr, status
    integer :: i

    call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

    if(nproc .ne. total_pe) then
       write(6,*) "Error: nproc .ne. total_pe ", nproc, total_pe
       call mpp_land_abort()
    endif

    ! get the neighbors.  -1 means no neighbor.
    left_id = info(2,my_id+1)
    right_id = info(3,my_id+1)
    up_id = info(4,my_id+1)
    down_id = info(5,my_id+1)
    IO_id = 0

    allocate(iolocal_allnx(nproc), stat=status)
    allocate(iolocal_allny(nproc), stat=status)
    allocate(iolocal_allrtnx(nproc), stat=status)
    allocate(iolocal_allrtny(nproc), stat=status)
    allocate(iolocal_allstarty(nproc), stat=status)
    allocate(iolocal_allstartx(nproc), stat=status)

    i = my_id + 1
    my_nx = info(7,i) - info(6,i) + 1
    my_ny = info(9,i) - info(8,i) + 1

    global_nx = 0
    global_ny = 0
    do i = 1, nproc
       global_nx = max(global_nx, info(7,i))
       global_ny = max(global_ny, info(9,i))
    enddo

    my_rtnx = my_nx * AGGFACTRT + 2
    my_rtny = my_ny * AGGFACTRT + 2
    if(left_id .lt. 0) my_rtnx = my_rtnx -1
    if(right_id .lt. 0) my_rtnx = my_rtnx -1
    if(up_id .lt. 0) my_rtny = my_rtny -1
    if(down_id .lt. 0) my_rtny = my_rtny -1

    global_rt_nx = global_nx * AGGFACTRT
    global_rt_ny = global_ny * AGGFACTRT
    rt_AGGFACTRT = AGGFACTRT

    do i =1,nproc
       iolocal_allnx(i) = info(7,i) - info(6,i) + 1
       iolocal_allny(i) = info(9,i) - info(8,i) + 1
       iolocal_allstartx(i)        = info(6,i)
       iolocal_allstarty(i)        = info(8,i)

       iolocal_allrtnx(i) = (info(7,i) - info(6,i) + 1) * AGGFACTRT + 2
       iolocal_allrtny(i) = (info(9,i) - info(8,i) + 1 ) * AGGFACTRT + 2
       if(info(2,i) .lt. 0) iolocal_allrtnx(i) = iolocal_allrtnx(i) -1
       if(info(3,i) .lt. 0) iolocal_allrtnx(i) = iolocal_allrtnx(i) -1
       if(info(4,i) .lt. 0) iolocal_allrtny(i) = iolocal_allrtny(i) -1
       if(info(5,i) .lt. 0) iolocal_allrtny(i) = iolocal_allrtny(i) -1
    end do
  end subroutine wrf_land_set_init

  subroutine mpp_land_get_myid()
    implicit none
    integer :: ierr
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  end subroutine mpp_land_get_myid

  subroutine mpp_channel_com_real(Link_location, ix, jy, Link_V, size, flag)
    ! communicate the data for channel routine.
    implicit none
    integer :: ix, jy, size
    integer :: Link_location(ix,jy)
    integer :: i, j, flag
    real :: Link_V(size), tmp_inout(ix,jy)

    tmp_inout = -999

    if(size .eq. 0) then
       tmp_inout = -999
    else
       !     map the Link_V data to tmp_inout(ix,jy)
       do i = 1,ix
          if(Link_location(i,1) .gt. 0) &
               tmp_inout(i,1) = Link_V(Link_location(i,1))
          if(Link_location(i,2) .gt. 0) &
               tmp_inout(i,2) = Link_V(Link_location(i,2))
          if(Link_location(i,jy-1) .gt. 0) &
               tmp_inout(i,jy-1) = Link_V(Link_location(i,jy-1))
          if(Link_location(i,jy) .gt. 0) &
               tmp_inout(i,jy) = Link_V(Link_location(i,jy))
       enddo
       do j = 1,jy
          if(Link_location(1,j) .gt. 0) &
               tmp_inout(1,j) = Link_V(Link_location(1,j))
          if(Link_location(2,j) .gt. 0) &
               tmp_inout(2,j) = Link_V(Link_location(2,j))
          if(Link_location(ix-1,j) .gt. 0) &
               tmp_inout(ix-1,j) = Link_V(Link_location(ix-1,j))
          if(Link_location(ix,j) .gt. 0) &
               tmp_inout(ix,j) = Link_V(Link_location(ix,j))
       enddo
    endif

    !   commu nicate tmp_inout
    call MPP_LAND_COM_REAL(tmp_inout, ix, jy, flag)

    !map the data back to Link_V
    if(size .eq. 0) return
    do j = 1,jy
       if( (Link_location(1,j) .gt. 0) .and. (tmp_inout(1,j) .ne. -999) ) &
            Link_V(Link_location(1,j)) = tmp_inout(1,j)
       if((Link_location(2,j) .gt. 0) .and. (tmp_inout(2,j) .ne. -999) ) &
            Link_V(Link_location(2,j)) = tmp_inout(2,j)
       if((Link_location(ix-1,j) .gt. 0) .and. (tmp_inout(ix-1,j) .ne. -999)) &
            Link_V(Link_location(ix-1,j)) = tmp_inout(ix-1,j)
       if((Link_location(ix,j) .gt. 0) .and. (tmp_inout(ix,j) .ne. -999) )&
            Link_V(Link_location(ix,j)) = tmp_inout(ix,j)
    enddo
    do i = 1,ix
       if((Link_location(i,1) .gt. 0) .and. (tmp_inout(i,1) .ne. -999) )&
            Link_V(Link_location(i,1)) = tmp_inout(i,1)
       if( (Link_location(i,2) .gt. 0) .and. (tmp_inout(i,2) .ne. -999) )&
            Link_V(Link_location(i,2)) = tmp_inout(i,2)
       if((Link_location(i,jy-1) .gt. 0) .and. (tmp_inout(i,jy-1) .ne. -999) ) &
            Link_V(Link_location(i,jy-1)) = tmp_inout(i,jy-1)
       if((Link_location(i,jy) .gt. 0) .and. (tmp_inout(i,jy) .ne. -999) ) &
            Link_V(Link_location(i,jy)) = tmp_inout(i,jy)
    enddo
  end subroutine mpp_channel_com_real


  subroutine mpp_channel_com_int(Link_location, ix, jy, Link_V, size, flag)
    ! communicate the data for channel routine.
    implicit none
    integer :: ix, jy, size
    integer :: Link_location(ix,jy)
    integer :: i, j, flag
    integer :: Link_V(size), tmp_inout(ix,jy)

    if(size .eq. 0) then
       tmp_inout = -999
    else
       !     map the Link_V data to tmp_inout(ix,jy)
       do i = 1,ix
          if(Link_location(i,1) .gt. 0) &
               tmp_inout(i,1) = Link_V(Link_location(i,1))
          if(Link_location(i,2) .gt. 0) &
               tmp_inout(i,2) = Link_V(Link_location(i,2))
          if(Link_location(i,jy-1) .gt. 0) &
               tmp_inout(i,jy-1) = Link_V(Link_location(i,jy-1))
          if(Link_location(i,jy) .gt. 0) &
               tmp_inout(i,jy) = Link_V(Link_location(i,jy))
       enddo
       do j = 1,jy
          if(Link_location(1,j) .gt. 0) &
               tmp_inout(1,j) = Link_V(Link_location(1,j))
          if(Link_location(2,j) .gt. 0) &
               tmp_inout(2,j) = Link_V(Link_location(2,j))
          if(Link_location(ix-1,j) .gt. 0) &
               tmp_inout(ix-1,j) = Link_V(Link_location(ix-1,j))
          if(Link_location(ix,j) .gt. 0) &
               tmp_inout(ix,j) = Link_V(Link_location(ix,j))
       enddo
    endif

    !   commu nicate tmp_inout
    call mpp_land_com_integer(tmp_inout, ix, jy, flag)

    !map the data back to Link_V
    if(size .eq. 0) return
    do j = 1,jy
       if( (Link_location(1,j) .gt. 0) .and. (tmp_inout(1,j) .ne. -999) ) &
            Link_V(Link_location(1,j)) = tmp_inout(1,j)
       if((Link_location(2,j) .gt. 0) .and. (tmp_inout(2,j) .ne. -999) ) &
            Link_V(Link_location(2,j)) = tmp_inout(2,j)
       if((Link_location(ix-1,j) .gt. 0) .and. (tmp_inout(ix-1,j) .ne. -999)) &
            Link_V(Link_location(ix-1,j)) = tmp_inout(ix-1,j)
       if((Link_location(ix,j) .gt. 0) .and. (tmp_inout(ix,j) .ne. -999) )&
            Link_V(Link_location(ix,j)) = tmp_inout(ix,j)
    enddo
    do i = 1,ix
       if((Link_location(i,1) .gt. 0) .and. (tmp_inout(i,1) .ne. -999) )&
            Link_V(Link_location(i,1)) = tmp_inout(i,1)
       if( (Link_location(i,2) .gt. 0) .and. (tmp_inout(i,2) .ne. -999) )&
            Link_V(Link_location(i,2)) = tmp_inout(i,2)
       if((Link_location(i,jy-1) .gt. 0) .and. (tmp_inout(i,jy-1) .ne. -999) ) &
            Link_V(Link_location(i,jy-1)) = tmp_inout(i,jy-1)
       if((Link_location(i,jy) .gt. 0) .and. (tmp_inout(i,jy) .ne. -999) ) &
            Link_V(Link_location(i,jy)) = tmp_inout(i,jy)
    enddo
  end subroutine mpp_channel_com_int


  subroutine mpp_land_max_int1(v)
    implicit none
    integer :: v, r1, max
    integer :: i, ierr, tag
    if(my_id .eq. IO_id) then
       max = v
       do i = 0, nproc - 1
          if(i .ne. my_id) then
             !block receive  from other node.
             tag = 101
             call mpi_recv(r1, 1, MPI_INTEGER,i, &
                  tag, MPI_COMM_WORLD, mpp_status, ierr)
             if(max <= r1) max = r1
          end if
       end do
    else
       tag =  101
       call mpi_send(v, 1, MPI_INTEGER, IO_id, &
            tag, MPI_COMM_WORLD, ierr)
    end if
    call mpp_land_bcast(max, 1, io_id)
    v = max
  end subroutine mpp_land_max_int1


  subroutine mpp_land_max_real1(v)
    implicit none
    real :: v, r1, max
    integer :: i, ierr, tag
    if(my_id .eq. IO_id) then
       max = v
       do i = 0, nproc - 1
          if(i .ne. my_id) then
             !block receive  from other node.
             tag = 101
             call mpi_recv(r1, 1, MPI_REAL, i, &
                  tag, MPI_COMM_WORLD, mpp_status, ierr)
             if(max <= r1) max = r1
          end if
       end do
    else
       tag =  101
       call mpi_send(v, 1, MPI_REAL, IO_id, &
            tag, MPI_COMM_WORLD, ierr)
    end if
    call mpp_land_bcast(max, 1, io_id)
    v = max
  end subroutine mpp_land_max_real1

  subroutine mpp_same_int1(v)
    implicit none
    integer :: v,r1
    integer :: i, ierr, tag
    if(my_id .eq. IO_id) then
       do i = 0, nproc - 1
          if(i .ne. my_id) then
             !block receive  from other node.
             tag = 109
             call mpi_recv(r1, 1, MPI_INTEGER, i, &
                  tag, MPI_COMM_WORLD, mpp_status, ierr)
             if(v .ne. r1) v = -99
          end if
       end do
    else
       tag =  109
       call mpi_send(v, 1, MPI_INTEGER, IO_id, &
            tag, MPI_COMM_WORLD, ierr)
    end if
    call mpp_land_bcast(v, 1, io_id)
  end subroutine mpp_same_int1


  subroutine write_chanel_real(v, map_l2g, gnlinks, nlinks, g_v)
    implicit none
    integer :: gnlinks, nlinks, map_l2g(nlinks)
    real :: recv(nlinks), v(nlinks)
    ! real g_v(gnlinks), tmp_v(gnlinks)
    integer :: i, ierr, tag, k
    integer :: length, node, message_len
    integer, allocatable, dimension(:) :: tmp_map
    real, allocatable, dimension(:) :: tmp_v
    real, dimension(:) :: g_v

    if(my_id .eq. io_id) then
       allocate(tmp_map(gnlinks))
       allocate(tmp_v(gnlinks))
       if(nlinks .le. 0) then
          tmp_map = -999
       else
          tmp_map(1:nlinks) = map_l2g(1:nlinks)
       endif
    else
       allocate(tmp_map(1))
       allocate(tmp_v(1))
    endif

    if(my_id .eq. IO_id) then
       do i = 0, nproc - 1
          message_len = mpp_nlinks(i+1)
          if(i .ne. my_id) then
             !block receive  from other node.
             tag = 109
             call mpi_recv(tmp_map(1:message_len), message_len, MPI_INTEGER, i, &
                  tag, MPI_COMM_WORLD, mpp_status, ierr)
             tag = 119

             call mpi_recv(tmp_v(1:message_len), message_len, MPI_REAL, i, &
                  tag, MPI_COMM_WORLD, mpp_status, ierr)

             do k = 1,message_len
                node = tmp_map(k)
                if(node .gt. 0) then
                   g_v(node) = tmp_v(k)
                else
                   write(6,*) "Maping infor k=",k," node=", node
                endif
             enddo
          else
             do k = 1,nlinks
                node = map_l2g(k)
                if(node .gt. 0) then
                   g_v(node) = v(k)
                else
                   write(6,*) "local Maping infor k=",k," node=",node
                endif
             enddo
          end if
       end do
    else
       tag = 109
       call mpi_send(map_l2g, nlinks, MPI_INTEGER, IO_id, &
            tag, MPI_COMM_WORLD, ierr)
       tag = 119
       call mpi_send(v, nlinks, MPI_REAL, IO_id, &
            tag, MPI_COMM_WORLD, ierr)
    end if
    deallocate(tmp_map)
    deallocate(tmp_v)
  end subroutine write_chanel_real


  subroutine write_chanel_int(v, map_l2g, gnlinks, nlinks, g_v)
    implicit none
    integer :: gnlinks, nlinks, map_l2g(nlinks)
    integer :: recv(nlinks), v(nlinks)
    integer, allocatable, dimension(:) :: tmp_map , tmp_v
    integer, dimension(:) :: g_v
    integer :: i, ierr, tag, k
    integer :: length, node, message_len

    if(my_id .eq. io_id) then
       allocate(tmp_map(gnlinks))
       allocate(tmp_v(gnlinks))
       if(nlinks .le. 0) then
          tmp_map = -999
       else
          tmp_map(1:nlinks) = map_l2g(1:nlinks)
       endif
    else
       allocate(tmp_map(1))
       allocate(tmp_v(1))
    endif

    if(my_id .eq. IO_id) then
       do i = 0, nproc - 1
          message_len = mpp_nlinks(i+1)
          if(i .ne. my_id) then
             !block receive  from other node.
             tag = 109
             call mpi_recv(tmp_map(1:message_len), message_len, MPI_INTEGER, i, &
                  tag, MPI_COMM_WORLD, mpp_status,ierr)
             tag = 119
             call mpi_recv(tmp_v(1:message_len), message_len, MPI_INTEGER, i, &
                  tag, MPI_COMM_WORLD, mpp_status, ierr)

             do k = 1,message_len
                if(tmp_map(k) .gt. 0) then
                   node = tmp_map(k)
                   g_v(node) = tmp_v(k)
                else
                   write(6,*) "Maping infor k=",k," node=",tmp_v(k)
                endif
             enddo
          else
             do k = 1,nlinks
                if(map_l2g(k) .gt. 0) then
                   node = map_l2g(k)
                   g_v(node) = v(k)
                else
                   write(6,*) "Maping infor k=",k," node=",map_l2g(k)
                endif
             enddo
          end if
       end do
    else
       tag =  109
       call mpi_send(map_l2g, nlinks, MPI_INTEGER, IO_id, &
            tag, MPI_COMM_WORLD, ierr)
       tag = 119
       call mpi_send(v, nlinks, MPI_INTEGER, IO_id, &
            tag, MPI_COMM_WORLD, ierr)
    end if
    deallocate(tmp_map)
    deallocate(tmp_v)
  end subroutine write_chanel_int


  subroutine write_lake_real(v, nodelist_in, nlakes)
    implicit none
    real :: recv(nlakes), v(nlakes)
    integer :: nodelist(nlakes), nlakes, nodelist_in(nlakes)
    integer :: i, ierr, tag, k
    integer :: node

    nodelist = nodelist_in
    if(my_id .eq. IO_id) then
       do i = 0, nproc - 1
          if(i .ne. my_id) then
             !block receive  from other node.
             tag = 129
             call mpi_recv(nodelist, nlakes, MPI_INTEGER, i, &
                  tag, MPI_COMM_WORLD, mpp_status, ierr)
             tag = 139
             call mpi_recv(recv(:), nlakes, MPI_REAL, i, &
                  tag, MPI_COMM_WORLD, mpp_status, ierr)

             do k = 1,nlakes
                if(nodelist(k) .gt. -99) then
                   node = nodelist(k)
                   v(node) = recv(node)
                endif
             enddo
          end if
       end do
    else
       tag = 129
       call mpi_send(nodelist, nlakes, MPI_INTEGER, IO_id, &
            tag, MPI_COMM_WORLD, ierr)
       tag = 139
       call mpi_send(v, nlakes, MPI_REAL, IO_id, &
            tag, MPI_COMM_WORLD, ierr)
    end if
  end subroutine write_lake_real

  subroutine read_rst_crt_r(unit, out, size)
    implicit none
    integer :: unit, size, ierr, ierr2
    real :: out(size), out1(size)
    if(my_id .eq. IO_id) then
       read(unit, IOSTAT=ierr2) out1
       if(ierr2 .eq. 0) out = out1
    endif
    call mpp_land_bcast(ierr2, 1, io_id)
    if(ierr2 .ne. 0) return
    call mpi_bcast(out, size, MPI_REAL, &
         IO_id, MPI_COMM_WORLD, ierr)
  end subroutine read_rst_crt_r

  subroutine write_rst_crt_r(unit, cd, map_l2g, gnlinks, nlinks)
    implicit none
    integer :: unit, gnlinks, nlinks, map_l2g(nlinks)
    real :: cd(nlinks)
    real :: g_cd(gnlinks)
    call write_chanel_real(cd, map_l2g, gnlinks, nlinks, g_cd)
    write(unit) g_cd
  end subroutine write_rst_crt_r

  subroutine sum_real8(vin, nsize)
    implicit none
    integer, intent(in) :: nsize
    real*8 :: vin(nsize), recv(nsize)
    real*8 :: v(nsize)
    integer :: i, tag, ierr
    tag = 319
    if(my_id .eq. IO_id) then
       do i = 0, nproc - 1
          if(i .ne. my_id) then
             call mpi_recv(recv, nsize, MPI_DOUBLE_PRECISION,i, &
                  tag, MPI_COMM_WORLD, mpp_status, ierr)
             vin(:) = vin(:) + recv(:)
          endif
       end do
       v = vin
    else
       call mpi_send(vin, nsize, MPI_DOUBLE_PRECISION, IO_id, &
            tag, MPI_COMM_WORLD, ierr)
    endif
    call mpp_land_bcast(v, nsize, io_id)
    vin = v
  end subroutine sum_real8


  subroutine gather_1d_real_tmp(vl, s_in, e_in, vg, sg)
    implicit none
    integer :: sg, s,e, size, s_in, e_in
    integer :: index_s(2)
    integer :: tag, ierr, i
    !   s: start index, e: end index
    real :: vl(e_in-s_in+1), vg(sg)
    s = s_in
    e = e_in

    if(my_id .eq. IO_id) then
       vg(s:e) = vl
    end if

    index_s(1) = s
    index_s(2) = e
    size = e - s + 1

    if(my_id .eq. IO_id) then
       do i = 0, nproc - 1
          if(i .ne. my_id) then
             !block receive  from other node.
             tag = 202
             call mpi_recv(index_s, 2, MPI_INTEGER, i, &
                  tag, MPI_COMM_WORLD, mpp_status, ierr)
             tag = 203
             e = index_s(2)
             s = index_s(1)
             size = e - s + 1
             call mpi_recv(vg(s:e), size, MPI_REAL, &
                  i, tag, MPI_COMM_WORLD, mpp_status, ierr)
          endif
       end do
    else
       tag = 202
       call mpi_send(index_s, 2, MPI_INTEGER, IO_id, &
            tag, MPI_COMM_WORLD, ierr)
       tag = 203
       call mpi_send(vl, size, MPI_REAL, IO_id, &
            tag, MPI_COMM_WORLD, ierr)
    end if
  end subroutine gather_1d_real_tmp

  subroutine sum_real1(inout)
    implicit none
    real:: inout, send
    integer :: ierr
    send = inout
    call MPI_ALLREDUCE(send, inout, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  end subroutine sum_real1

  subroutine sum_double(inout)
    implicit none
    real*8 :: inout, send
    integer :: ierr
    send = inout
    !yw CALL MPI_ALLREDUCE(send,inout,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(send, inout, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  end subroutine sum_double

  subroutine mpp_chrt_nlinks_collect(nlinks)
    ! collect the nlinks
    implicit none
    integer :: nlinks
    integer :: i, ierr, status, tag
    allocate(mpp_nlinks(nproc), stat=status)
    tag = 138
    mpp_nlinks = 0
    if(my_id .eq. IO_id) then
       do i = 0,nproc -1
          if(i .ne. my_id) then
             call mpi_recv(mpp_nlinks(i+1), 1, MPI_INTEGER, i, &
                  tag, MPI_COMM_WORLD, mpp_status, ierr)
          else
             mpp_nlinks(i+1) = 0
          end if
       end do
    else
       call mpi_send(nlinks, 1, MPI_INTEGER, IO_id, &
            tag, MPI_COMM_WORLD, ierr)
    endif
  end subroutine mpp_chrt_nlinks_collect


  subroutine mpp_land_get_localxy(startx, endx, starty, endy)
    implicit none
    integer, intent(out) :: startx, endx
    integer, intent(out) :: starty, endy
    startx = my_startx
    endx = my_startx + my_nx - 1
    starty = my_starty
    endy = my_starty + my_ny - 1
  end subroutine mpp_land_get_localxy


  subroutine check_landreal1(unit, inVar)
    implicit none
    integer :: unit
    real :: inVar
    if(my_id .eq. IO_id) then
       write(unit,*) inVar
       call flush(unit)
    endif
  end subroutine check_landreal1

  subroutine check_landreal1d(unit, inVar)
    implicit none
    integer :: unit
    real :: inVar(:)
    if(my_id .eq. IO_id) then
       write(unit,*) inVar
       call flush(unit)
    endif
  end subroutine check_landreal1d


  subroutine check_landreal2d(unit, inVar)
    implicit none
    integer :: unit
    real :: inVar(:,:)
    real :: g_var(global_nx,global_ny)
    call write_io_real(inVar,g_var)
    if(my_id .eq. IO_id) then
       write(unit,*) g_var
       call flush(unit)
    endif
  end subroutine check_landreal2d


  subroutine check_landreal3d(unit, inVar)
    implicit none
    integer :: unit, k, klevel
    real :: inVar(:,:,:)
    real :: g_var(global_nx,global_ny)
    klevel = size(inVar,2)
    do k = 1, klevel
       call write_io_real(inVar(:,k,:),g_var)
       if(my_id .eq. IO_id) then
          write(unit,*) g_var
          call flush(unit)
       endif
    end do
  end subroutine check_landreal3d

end module module_mpp_land
