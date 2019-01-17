module module_hrldas_noahmp_namelist
  implicit none

  private
  public :: config
  public :: hrldas_noahmp_namelist_read
  public :: hrldas_noahmp_namelist_get_ntime

  integer, parameter, private :: MAX_PATH_LEN = 512
  integer, parameter, private :: MAX_SOIL_LEVELS = 100
  integer, parameter, private :: undefined_integer = -999
  real, parameter, private :: undefined_real = -1.0e20
  character(len=MAX_PATH_LEN), private :: default_namelist_filename = 'noahmp.namelist'

  type hrldas_noahmp_namelist_t
     ! temporal discretization
     integer :: start_year = 0
     integer :: start_month = 0
     integer :: start_day = 0
     integer :: start_hour = 0
     integer :: start_min = 0
     integer :: start_sec = 1

     integer :: stop_year = 0
     integer :: stop_month = 0
     integer :: stop_day =0
     integer :: stop_hour = 0
     integer :: stop_min = 0
     integer :: stop_sec = 0

     integer :: kday            !TODO
     integer :: khour           !TODO

     ! spatial discretization
     integer :: nsoil = 4 ! number of soil layers
     real :: soil_thickness(MAX_SOIL_LEVELS) = 0.1 ! [m]

     integer :: xstart = 1
     integer :: ystart = 1
     integer :: xend = 0
     integer :: yend = 0

     ! integration, input, and output
     character(len=MAX_PATH_LEN) :: const_file = 'const.nc'
     character(len=MAX_PATH_LEN) :: init_file = 'init.nc'
     character(len=MAX_PATH_LEN) :: restart_file = ' '
     logical :: from_restart = .false.
     integer :: model_timestep = 0 ! [s]

     character(len=MAX_PATH_LEN) :: indir = 'input'
     integer :: input_timestep = 3600
     real :: zlvl = 30.0
     character(len=MAX_PATH_LEN) :: outdir = 'output'
     integer :: output_timestep = 3600
     character(len=MAX_PATH_LEN) :: resdir = 'restart'
     integer :: restart_timestep = 3600

     character(len=MAX_PATH_LEN) :: mmf_runoff_file = ' '

     ! noahmp parameterization options
     integer :: opt_veg ! dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1
     integer :: opt_crs ! canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
     integer :: opt_btr ! soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
     integer :: opt_run ! runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
     integer :: opt_sfc ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
     integer :: opt_frz ! supercooled liquid water (1-> NY06; 2->Koren99)
     integer :: opt_inf ! frozen soil permeability (1-> NY06; 2->Koren99)
     integer :: opt_rad ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
     integer :: opt_alb ! snow surface albedo (1->BATS; 2->CLASS)
     integer :: opt_snf ! rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)
     integer :: opt_tbot ! lower boundary of soil temperature (1->zero-flux; 2->Noah)
     integer :: opt_stc ! snow/soil temperature time scheme
  end type hrldas_noahmp_namelist_t

  type(hrldas_noahmp_namelist_t) :: config

contains

  ! read namelist from the namelist_filename file
  subroutine hrldas_noahmp_namelist_read(filename)
    implicit none
    character(len=*), intent(in), optional :: filename
    character(len=MAX_PATH_LEN) :: namelist_filename = ' '
    ! local variables
    integer :: iostat
    ! temporal discretization
    integer :: start_year
    integer :: start_month
    integer :: start_day
    integer :: start_hour
    integer :: start_min
    integer :: start_sec
    integer :: stop_year
    integer :: stop_month
    integer :: stop_day
    integer :: stop_hour
    integer :: stop_min
    integer :: stop_sec
    ! TODO - TO BE REMOVED
    integer :: kday
    integer :: khour
    ! spatial discretization
    integer :: nsoil
    real :: soil_layer_thickness(MAX_SOIL_LEVELS) ! [m]
    integer :: xstart = 1
    integer :: ystart = 1
    integer :: xend = 0
    integer :: yend = 0
    ! integration, input, and output
    character(len=MAX_PATH_LEN) :: const_file = 'const.nc'
    character(len=MAX_PATH_LEN) :: init_file = 'init.nc'
    character(len=MAX_PATH_LEN) :: restart_file = 'restart.nc'
    logical :: from_restart = .false.
    integer :: model_timestep
    character(len=MAX_PATH_LEN) :: indir = '.'
    integer :: input_timestep
    real :: zlvl
    character(len=MAX_PATH_LEN) :: outdir = '.'
    integer :: output_timestep
    character(len=MAX_PATH_LEN) :: resdir = '.'
    integer :: restart_timestep
    character(len=MAX_PATH_LEN) :: mmf_runoff_file = ' '

    ! noahmp parameterization options
    integer :: dynamic_veg_option = undefined_integer
    integer :: canopy_stomatal_resistance_option = undefined_integer
    integer :: btr_option = undefined_integer
    integer :: runoff_option = undefined_integer
    integer :: surface_drag_option = undefined_integer
    integer :: supercooled_water_option = undefined_integer
    integer :: frozen_soil_option = undefined_integer
    integer :: radiative_transfer_option = undefined_integer
    integer :: snow_albedo_option = undefined_integer
    integer :: precipitation_partition_option = undefined_integer
    integer :: tbot_option = undefined_integer
    integer :: temp_time_scheme_option

    namelist / HRLDAS_NOAHMP / &
         & start_year, start_month, start_day, &
         & start_hour, start_min, start_sec, &
         & stop_year, stop_month, stop_day, &
         & stop_hour, stop_min, stop_sec, &
         & kday, khour, &
         & const_file, init_file, restart_file, &
         & from_restart, model_timestep, &
         & indir, input_timestep, zlvl, &
         & outdir, output_timestep, resdir, restart_timestep, &

         & nsoil, soil_layer_thickness, &
         & xstart, ystart, xend, yend, &
         & mmf_runoff_file, &

         & dynamic_veg_option, &
         & canopy_stomatal_resistance_option, &
         & btr_option, &
         & runoff_option, &
         & surface_drag_option, &
         & supercooled_water_option, &
         & frozen_soil_option, &
         & radiative_transfer_option, &
         & snow_albedo_option, &
         & precipitation_partition_option, &
         & tbot_option, &
         & temp_time_scheme_option

    ! dummy values
    nsoil = -999
    soil_layer_thickness = -999.0
    start_year = -999
    start_month = -999
    start_day = -999
    start_hour = -999
    start_min = -999
    start_sec = -999
    stop_year = -999
    stop_month = -999
    stop_day = -999
    stop_hour = -999
    stop_min = -999
    stop_sec = -999
    zlvl = -999.0

    kday = -999
    khour = -999
    input_timestep = 0
    output_timestep = 0
    model_timestep = 0

    restart_file = ' '

    if (present(filename)) then
       namelist_filename = filename
    else
       namelist_filename = default_namelist_filename
    endif
    open(30, file=namelist_filename, form="FORMATTED", status="OLD", action="READ")
    read(30, HRLDAS_NOAHMP, iostat=iostat)
    if (iostat /= 0) then
       write(*,'(/, "ERROR: Problem reading namelist HRLDAS_NOAHMP from file ", A/)') namelist_filename
       rewind(30)
       read(30, HRLDAS_NOAHMP)
       stop "ERROR: Problem reading namelist HRLDAS_NOAHMP"
    endif
    close(30)

    ! NAMELIST check begin
    if (nsoil < 0) then
       write(*, '("ERROR[namelist]: NSOIL must be set.")')
       stop
    endif

    if ((khour < 0) .and. (kday < 0)) then
       write(*, '("ERROR[namelist]: Either KHOUR or KDAY must be defined.")')
       stop
    else if (( khour < 0 ) .and. (kday > 0)) then
       khour = kday * 24
    else if ((khour > 0) .and. (kday > 0)) then
       write(*, '("WARNING[namelist]: KHOUR and KDAY are both defined.")')
    else
       ! all is well.  KHOUR defined
    endif

    if (input_timestep <= 0) then
       write(*, '("ERROR[namelist]: INPUT_TIMESTEP must be greater than zero.")')
       stop
    endif

    if (model_timestep <= 0) then
       write(*, '("ERROR[namelist]: MODEL_TIMESTEP must be greater than zero.")')
       write(*, '("                                900 seconds is recommended.")')
       stop
    endif

    ! Check that OUTPUT_TIMESTEP fits into MODEL_TIMESTEP:
    if (output_timestep > 0) then
       if (mod(output_timestep, model_timestep) > 0) then
          write(*, '("ERROR[namelist]: OUTPUT_TIMESTEP shoud be set to an integer multiple of MODEL_TIMESTEP.")')
          write(*, '("INFO[namelist]:  OUTPUT_TIMESTEP = ", I12, " seconds")') output_timestep
          write(*, '("INFO[namelist]:  MODEL_TIMESTEP  = ", I12, " seconds")') model_timestep
          stop
       endif
    endif

    ! Check that RESTART_TIMESTEP fits into MODEL_TIMESTEP:
    if (restart_timestep /= 0) then
       if (mod(restart_timestep, model_timestep) > 0) then
          write(*, '("ERROR[namelist]: RESTART_TIMESTEP shoud be set to an integer multiple of MODEL_TIMESTEP.")')
          write(*, '("INFO[namelist]:  RESTART_TIMESTEP = ", I12, " seconds")') restart_timestep
          write(*, '("INFO[namelist]:  MODEL_TIMESTEP   = ", I12, " seconds")') model_timestep

          stop
       endif
    endif

    if (dynamic_veg_option == 2) then
       if (canopy_stomatal_resistance_option /= 1) then
          write(*, '("ERROR[namelist]: CANOPY_STOMATAL_RESISTANCE_OPTION must be 1 when DYNAMIC_VEG_OPTION == 2")')
          stop
       endif
    endif
    ! NAMELIST check end

    config%start_year = start_year
    config%start_month = start_month
    config%start_day = start_day
    config%start_hour = start_hour
    config%start_min = start_min
    config%start_sec = start_sec
    config%stop_year = stop_year
    config%stop_month = stop_month
    config%stop_day = stop_day
    config%stop_hour = stop_hour
    config%stop_min = stop_min
    config%stop_sec = stop_sec
    config%kday = kday
    config%khour = khour

    config%nsoil = nsoil
    config%soil_thickness(1:config%nsoil) = soil_layer_thickness(1:config%nsoil)
    config%xstart = 1
    config%ystart = 1
    config%xend = 0
    config%yend = 0

    config%const_file = const_file
    config%init_file = init_file
    config%restart_file = restart_file
    config%from_restart = from_restart
    config%model_timestep = model_timestep

    config%indir = indir
    config%outdir = outdir
    config%resdir = resdir
    config%input_timestep = input_timestep
    config%output_timestep = output_timestep
    config%restart_timestep = restart_timestep

    config%opt_veg = dynamic_veg_option
    config%opt_crs = canopy_stomatal_resistance_option
    config%opt_btr = btr_option
    config%opt_run = runoff_option
    config%opt_sfc = surface_drag_option
    config%opt_frz = supercooled_water_option
    config%opt_inf = frozen_soil_option
    config%opt_rad = radiative_transfer_option
    config%opt_alb = snow_albedo_option
    config%opt_snf = precipitation_partition_option
    config%opt_tbot = tbot_option
    config%opt_stc = temp_time_scheme_option
  end subroutine hrldas_noahmp_namelist_read


  subroutine hrldas_noahmp_namelist_get_ntime(ntime)
    implicit none
    integer, intent(out) :: ntime
    ntime = config%khour * 3600. / config%model_timestep
  end subroutine hrldas_noahmp_namelist_get_ntime

end module module_hrldas_noahmp_namelist
