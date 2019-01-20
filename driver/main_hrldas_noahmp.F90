program hrldas_noahmp
  ! this is the main program to drive the NoahMP land model.
  use module_hrldas_noahmp_namelist
  use module_hrldas_noahmp_driver, only: land_driver_init, land_driver_exe
  implicit none
  integer :: itime, ntime

  ! read configurations
  call hrldas_noahmp_namelist_read()
  call hrldas_noahmp_namelist_get_ntime(ntime)

  ! initilization
  call land_driver_init()

  ! run
  do itime = 1, ntime
     call land_driver_exe(ITIME)
  end do

end program hrldas_noahmp
