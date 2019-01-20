module module_hrldas_noahmp_driver
  use module_hrldas_noahmp_namelist
  use module_hrldas_noahmp_vars
  use module_hrldas_netcdf_io
  use module_sf_noahmp_groundwater
  use module_sf_noahmpdrv, only: noahmp_init, noahmplsm, soil_veg_gen_parm
  use module_sf_noahmplsm, only: read_mp_veg_parameters
  use module_date_utilities
#ifdef MPP_LAND
  use mpi
  use module_mpp_land, only: mpp_land_init, mpp_land_partition
  use module_mpp_land, only: my_id, io_id, nproc, mpp_land_get_localxy, mpp_land_bcast
  use module_cpl_land, only: cpl_land_init
#endif
  implicit none

  integer, parameter :: MAX_PATH_LEN = 512
  integer, parameter :: DATETIME_STR_LEN = 19
  integer, parameter :: MAX_STR_LEN = 256

  ! model timer
  integer            :: ITIMESTEP ! timestep number
  integer            :: YR        ! 4-digit year
  real               :: JULIAN_IN ! Julian day
  ! parameter table
  character(LEN=256) :: MMINSL = 'STAS' ! soil classification
  character(LEN=256) :: LLANDUSE ! (=USGS, using USGS landuse classification)
  real               :: XICE_THRESHOLD ! fraction of grid determining seaice
  integer            :: ISICE     ! land cover category for ice
  integer            :: ISURBAN   ! land cover category for urban
  integer            :: ISWATER   ! land cover category for water
  integer            :: ISOILWATER! soil category for water
  integer            :: IZ0TLND   ! option of Chen adjustment of Czil (not used)

  ! Needed for initializing NoahMP
  integer, allocatable, dimension(:,:) :: landmask ! 1 for land
  logical                           :: FNDSOILW    ! soil water present in input
  logical                           :: FNDSNOWH    ! snow depth present in input
  real, allocatable, dimension(:,:) :: CHSTARXY    ! for consistency with MP_init; delete later
  real, allocatable, dimension(:,:) :: SEAICE      ! seaice fraction

  character(len=9), parameter :: version = "v20140513"
  integer :: LDASIN_VERSION

  ! execution timer:
  integer :: clock_count_1 = 0
  integer :: clock_count_2 = 0
  integer :: clock_rate    = 0
  real    :: timing_sum    = 0.0

  integer :: sflx_count_sum
  integer :: count_before_sflx
  integer :: count_after_sflx

  !  DECLARE/Initialize constants
  integer                             :: I
  integer                             :: J, K
  integer                             :: SLOPETYP
  integer                             :: YEARLEN

  !  File naming, parallel
  character(len=19)  :: olddate, newdate, startdate
  character(len=19)  :: tmpdate
  logical            :: lexist
  integer            :: imode
  integer            :: ixfull
  integer            :: jxfull
  integer            :: ixpar
  integer            :: jxpar
  integer            :: xstartpar
  integer            :: ystartpar
  integer            :: rank = 0
  character(len=256) :: inflnm, outflnm, inflnm_template
  character(len=256) :: restart_flnm
  integer            :: ierr

  ! Attributes from CONST_FILE
  integer :: IX
  integer :: JX
  real :: dx
  real :: DY
  real :: TRUELAT1
  real :: TRUELAT2
  real :: CEN_LON
  integer :: MAPPROJ

contains

  subroutine land_driver_init(wrfits, wrfite, wrfjts, wrfjte)
    implicit none
    integer, optional, intent(in) :: wrfits, wrfite, wrfjts, wrfjte

    ! initilization for stand alone parallel code.
#ifdef MPP_LAND
    call mpp_land_init()
#endif

    ! Initialize gridded domain
    call hrldas_diminfo(config%const_file, ix, jx)
#ifdef MPP_LAND
    call mpp_land_partition(1, ix, jx, 1)
    call mpp_land_get_localxy(xstart, xend, ystart, yend)
    write(*,*) 'INFO[DRIVER]: Parallel run using ', nproc, ' nodes'
#else
    write(*,*) 'INFO[DRIVER]: Sequential run'
    xstart = 1
    xend = ix
    ystart = 1
    yend = jx
#endif

    ids = xstart
    ide = xend
    jds = ystart
    jde = yend
    kds = 1
    kde = 2
    its = xstart
    ite = xend
    jts = ystart
    jte = yend
    kts = 1
    kte = 2
    ims = xstart
    ime = xend
    jms = ystart
    jme = yend
    kms = 1
    kme = 2

    !  Allocate multi-dimension fields for subwindow calculation
    ixfull = xend-xstart+1
    jxfull = yend-ystart+1

    ixpar = ixfull
    jxpar = jxfull
    xstartpar = 1
    ystartpar = 1

    call hrldas_hdrinfo(config%const_file, &
         xstart, xend, ystart, yend, &
         llanduse, iswater, isurban, isice, isoilwater, &
         & dx, dy, truelat1, truelat2, cen_lon, mapproj)

    write(startdate,'(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2)') &
         config%start_year, config%start_month, config%start_day, config%start_hour, config%start_min, config%start_sec

    olddate = startdate

    ! initialize timer, domain, parameter, and initial states
    ! TODO
    ! initialize state
    call hrldas_noahmp_vars_init()

    ITIMESTEP = 1

    allocate(seaice(xstart:xend,ystart:yend)) ! seaice fraction
    seaice = undefined_real
    XLAND          = 1.0   ! water = 2.0, land = 1.0
    XICE           = 0.0   ! fraction of grid that is seaice
    XICE_THRESHOLD = 0.5   ! fraction of grid determining seaice

    ! Read Landuse Type and Soil Texture and Other Information
    call hrldas_const_read(config%const_file, xstart, xend, ystart, yend, &
         ISWATER, lat2d, lon2d, XLAND, IVGTYP, ISLTYP, TERRAIN, TMN, SEAICE, MSFTX, MSFTY)

    where(SEAICE > 1.0e-6) XICE = 1.0

    ! For OPT_RUN = 5 (MMF groundwater), read in necessary extra fields
    if (config%opt_run == 5) then
       !!yw not tested for this if block
       call READ_MMF_RUNOFF(config%mmf_runoff_file, xstart, xend, ystart, yend, &
            ZWTXY, EQZWT, RIVERBEDXY, RIVERCONDXY, PEXPXY, FDEPTHXY)
    end if

    ! Initialize Model State
    SLOPETYP = 2
    DZS = config%soil_thickness(1:config%nsoil)

    if (config%from_restart) then
       call hrldas_restart_read(trim(config%restart_file), &
            xstart, xend, ystart, yend, config%nsoil, olddate)

       ITIMESTEP = 2
       call hrldas_restart_get(xstart, xend, ystart, yend, "SOIL_T"  , TSLB     )
       call hrldas_restart_get(xstart, xend, ystart, yend, "SNOW_T"  , TSNOXY   )
       call hrldas_restart_get(xstart, xend, ystart, yend, "SMC"     , SMOIS    )
       call hrldas_restart_get(xstart, xend, ystart, yend, "SH2O"    , SH2O     )
       call hrldas_restart_get(xstart, xend, ystart, yend, "ZSNSO"   , ZSNSOXY  )
       call hrldas_restart_get(xstart, xend, ystart, yend, "SNICE"   , SNICEXY  )
       call hrldas_restart_get(xstart, xend, ystart, yend, "SNLIQ"   , SNLIQXY  )
       call hrldas_restart_get(xstart, xend, ystart, yend, "QSNOW"   , QSNOWXY  )
       call hrldas_restart_get(xstart, xend, ystart, yend, "FWET"    , FWETXY   )
       call hrldas_restart_get(xstart, xend, ystart, yend, "SNEQVO"  , SNEQVOXY )
       call hrldas_restart_get(xstart, xend, ystart, yend, "EAH"     , EAHXY    )
       call hrldas_restart_get(xstart, xend, ystart, yend, "TAH"     , TAHXY    )
       call hrldas_restart_get(xstart, xend, ystart, yend, "ALBOLD"  , ALBOLDXY )
       call hrldas_restart_get(xstart, xend, ystart, yend, "CM"      , CMXY     )
       call hrldas_restart_get(xstart, xend, ystart, yend, "CH"      , CHXY     )
       call hrldas_restart_get(xstart, xend, ystart, yend, "ISNOW"   , ISNOWXY  )
       call hrldas_restart_get(xstart, xend, ystart, yend, "CANLIQ"  , CANLIQXY )
       call hrldas_restart_get(xstart, xend, ystart, yend, "CANICE"  , CANICEXY )
       call hrldas_restart_get(xstart, xend, ystart, yend, "SNEQV"   , SNOW     )
       call hrldas_restart_get(xstart, xend, ystart, yend, "SNOWH"   , SNOWH    )
       call hrldas_restart_get(xstart, xend, ystart, yend, "TV"      , TVXY     )
       call hrldas_restart_get(xstart, xend, ystart, yend, "TG"      , TGXY     )
       call hrldas_restart_get(xstart, xend, ystart, yend, "ZWT"     , ZWTXY    )
       call hrldas_restart_get(xstart, xend, ystart, yend, "WA"      , WAXY     )
       call hrldas_restart_get(xstart, xend, ystart, yend, "WT"      , WTXY     )
       call hrldas_restart_get(xstart, xend, ystart, yend, "WSLAKE"  , WSLAKEXY )
       call hrldas_restart_get(xstart, xend, ystart, yend, "LFMASS"  , LFMASSXY )
       call hrldas_restart_get(xstart, xend, ystart, yend, "RTMASS"  , RTMASSXY )
       call hrldas_restart_get(xstart, xend, ystart, yend, "STMASS"  , STMASSXY )
       call hrldas_restart_get(xstart, xend, ystart, yend, "WOOD"    , WOODXY   )
       call hrldas_restart_get(xstart, xend, ystart, yend, "STBLCP"  , STBLCPXY )
       call hrldas_restart_get(xstart, xend, ystart, yend, "FASTCP"  , FASTCPXY )
       call hrldas_restart_get(xstart, xend, ystart, yend, "LAI"     , LAI      )
       call hrldas_restart_get(xstart, xend, ystart, yend, "SAI"     , XSAIXY   )
       call hrldas_restart_get(xstart, xend, ystart, yend, "FPAR"    , VEGFRA   )
       call hrldas_restart_get(xstart, xend, ystart, yend, "GVFMIN"  , GVFMIN   )
       call hrldas_restart_get(xstart, xend, ystart, yend, "SHDMAX"  , SHDMAX   )
       call hrldas_restart_get(xstart, xend, ystart, yend, "GVFMAX"  , GVFMAX   )
       call hrldas_restart_get(xstart, xend, ystart, yend, "ACMELT"  , ACSNOM   )
       call hrldas_restart_get(xstart, xend, ystart, yend, "ACSNOW"  , ACSNOW   )
       call hrldas_restart_get(xstart, xend, ystart, yend, "TAUSS"   , TAUSSXY  )
       call hrldas_restart_get(xstart, xend, ystart, yend, "QSFC"    , QSFC     )
       call hrldas_restart_get(xstart, xend, ystart, yend, "SFCRUNOFF",SFCRUNOFF   )
       call hrldas_restart_get(xstart, xend, ystart, yend, "UDRUNOFF" ,UDRUNOFF    )

       ! below for opt_run = 5
       call hrldas_restart_get(xstart, xend, ystart, yend, "SMOISEQ"   , SMOISEQ    )
       call hrldas_restart_get(xstart, xend, ystart, yend, "AREAXY"    , AREAXY     )
       call hrldas_restart_get(xstart, xend, ystart, yend, "SMCWTDXY"  , SMCWTDXY   )
       call hrldas_restart_get(xstart, xend, ystart, yend, "QRFXY"     , QRFXY      )
       call hrldas_restart_get(xstart, xend, ystart, yend, "DEEPRECHXY", DEEPRECHXY )
       call hrldas_restart_get(xstart, xend, ystart, yend, "QSPRINGXY" , QSPRINGXY  )
       call hrldas_restart_get(xstart, xend, ystart, yend, "QSLATXY"   , QSLATXY    )
       call hrldas_restart_get(xstart, xend, ystart, yend, "QRFSXY"    , QRFSXY     )
       call hrldas_restart_get(xstart, xend, ystart, yend, "QSPRINGSXY", QSPRINGSXY )
       call hrldas_restart_get(xstart, xend, ystart, yend, "RECHXY"    , RECHXY     )
       call hrldas_restart_get(xstart, xend, ystart, yend, "FDEPTHXY"   ,FDEPTHXY   )
       call hrldas_restart_get(xstart, xend, ystart, yend, "RIVERCONDXY",RIVERCONDXY)
       call hrldas_restart_get(xstart, xend, ystart, yend, "RIVERBEDXY" ,RIVERBEDXY )
       call hrldas_restart_get(xstart, xend, ystart, yend, "EQZWT"      ,EQZWT      )
       call hrldas_restart_get(xstart, xend, ystart, yend, "PEXPXY"     ,PEXPXY     )

       STEPWTD = nint(WTDDT*60./config%model_timestep)
       STEPWTD = max(STEPWTD,1)

       call read_mp_veg_parameters(trim(LLANDUSE))
       call SOIL_VEG_GEN_PARM(LLANDUSE, MMINSL)

    else
       SMOIS     =  undefined_real
       TSLB      =  undefined_real
       SH2O      =  undefined_real
       CANLIQXY  =  undefined_real
       TSK       =  undefined_real
       RAINBL_tmp = undefined_real
       SNOW      =  undefined_real
       SNOWH     =  undefined_real

       call hrldas_init_read(config%init_file, xstart, xend, ystart, yend, &
            config%nsoil, SMOIS, TSLB, &
            CANWAT, TSK, SNOW, SNOWH, FNDSNOWH)
       SNOW = SNOW * 1000. ! Convert snow water equivalent to mm.

       VEGFRA = undefined_real
       LAI    = undefined_real
       GVFMIN = undefined_real
       SHDMAX = undefined_real
       !  call hrldas_parm_veg_read(inflnm, xstart, xend, ystart, yend, &
       !       OLDDATE, IVGTYP, VEGFRA, LAI, GVFMIN, SHDMAX)
       !  VEGFRA = VEGFRA * 100.
       !  GVFMIN = GVFMIN * 100.
       !  SHDMAX = SHDMAX * 100.

       allocate(CHSTARXY(xstart:xend,ystart:yend)) ! for consistency with MP_init; delete later
       CHSTARXY = undefined_real
       FNDSOILW = .false.
       call NOAHMP_INIT(LLANDUSE, SNOW, SNOWH, CANWAT, ISLTYP, IVGTYP, ISURBAN, & ! call from WRF phys_init
            TSLB, SMOIS, SH2O, DZS, FNDSOILW, FNDSNOWH, ISICE, iswater, &
            TSK, ISNOWXY, TVXY, TGXY, CANICEXY, TMN, XICE, &
            CANLIQXY, EAHXY, TAHXY, CMXY, CHXY, &
            FWETXY, SNEQVOXY, ALBOLDXY, QSNOWXY, WSLAKEXY, ZWTXY, WAXY, &
            WTXY, TSNOXY, ZSNSOXY, SNICEXY, SNLIQXY, LFMASSXY, RTMASSXY, &
            STMASSXY, WOODXY, STBLCPXY, FASTCPXY, XSAIXY, &
            T2MVXY, T2MBXY, CHSTARXY, &
            config%nsoil, .false., &
            .true., config%opt_run, &
            ids,ide+1, jds,jde+1, kds,kde, & ! domain
            ims,ime, jms,jme, kms,kme, & ! memory
            its,ite, jts,jte, kts,kte, & ! tile
            smoiseq, smcwtdxy, rechxy, deeprechxy, areaxy ,dx, dy, msftx, msfty, &
            wtddt, stepwtd, real(config%model_timestep), qrfsxy, qspringsxy, qslatxy, &
            fdepthxy, terrain, riverbedxy, eqzwt, rivercondxy, pexpxy &
            )

       TAUSSXY = 0.0   ! Need to be added to _INIT later
    endif

    ! Begin Time Loop

    call system_clock(count=clock_count_1)   ! Start a timer

  end subroutine land_driver_init


  subroutine land_driver_exe(ITIME)
    implicit none
    integer :: ITIME ! timestep loop
    ! TIMELOOP : DO ITIME=1,NTIME
    ! Output for restart BEFORE the processing for this particular time has begun,
    ! so that upon restart, we're ready to go on this time step.
    !---------------------------------------------------------------------------------
    ! Read the forcing data.
    !---------------------------------------------------------------------------------
    ! For HRLDAS, we're assuming (for now) that each time period is in a
    ! separate file.  So we can open a new one right now.

    inflnm = trim(config%indir) // '/input.' &
         // startdate(1:4) // startdate(6:7) // startdate(9:10) // 'T' &
         // startdate(12:13) // startdate(15:16) // startdate(18:19)

    ! Build a filename template
    inflnm_template = trim(config%indir)//"/input.<date>"

#ifdef MPP_LAND
    call mpp_land_bcast(olddate(1:19), 19, io_id)
#endif

    call hrldas_input_read(INFLNM_TEMPLATE, config%input_timestep, OLDDATE, &
         xstart, xend, ystart, yend, &
         T_PHY(:,1,:), QV_CURR(:,1,:),U_PHY(:,1,:),V_PHY(:,1,:), &
         P8W(:,1,:), GLW, SWDOWN, RAINBL_tmp, LAI, VEGFRA)
    VEGFRA = VEGFRA * 100.0

    P8W(:,2,:)     = P8W(:,1,:)   ! WRF uses lowest two layers
    T_PHY(:,2,:)   = T_PHY(:,1,:) ! Only pressure is needed in two layer but fill the rest
    U_PHY(:,2,:)   = U_PHY(:,1,:) !
    V_PHY(:,2,:)   = V_PHY(:,1,:) !
    QV_CURR(:,2,:) = QV_CURR(:,1,:) !
    RAINBL = RAINBL_tmp * real(config%model_timestep) ! RAINBL in WRF is [mm]
    DZ8W = 2.0 * config%zlvl           ! 2* to be consistent with WRF model level

    ! Noah-MP updates we can do before spatial loop.

    ! create a few fields that are IN in WRF - coszen, julian_in,yr

    do J = ystart,yend
       do I = xstart,xend
          call CALC_DECLIN(OLDDATE(1:19), lat2d(I,J), lon2d(I,J), COSZEN(I,J), JULIAN_IN)
       end do
    end do

    read(OLDDATE(1:4),*) YR
    YEARLEN = 365 ! find length of year for phenology (also S Hemisphere)
    if (mod(YR,4) == 0) then
       YEARLEN = 366
       if (mod(YR,100) == 0) then
          YEARLEN = 365
          if (mod(YR,400) == 0) then
             YEARLEN = 366
          endif
       endif
    endif

    ! Call to Noah-MP driver same as surface_driver
    sflx_count_sum = 0 ! Timing

    ! Timing information for SFLX:
    call system_clock(count=count_before_sflx, count_rate=clock_rate)

    if (itime == 1 .and. .not. config%from_restart) then
       ! Initial guess only.  EAH gets updated in SFLX for later time steps.
       EAHXY = (P8W(:,1,:)*QV_CURR(:,1,:))/(0.622+QV_CURR(:,1,:))

       ! TAH: Canopy air temperature (K)
       ! Initial guess only. TAH gets updated in SFLX for later time steps.
       TAHXY = T_PHY(:,1,:)

       CHXY = 0.1
       CMXY = 0.1
    endif

    call noahmplsm(ITIMESTEP, YR, JULIAN_IN, COSZEN, lat2d, &
         DZ8W, real(config%model_timestep), DZS, config%nsoil, DX, &
         IVGTYP, ISLTYP, VEGFRA, SHDMAX, TMN, &
         XLAND, XICE, XICE_THRESHOLD, ISICE, ISURBAN, &
         config%opt_veg, config%opt_crs, config%opt_btr, config%opt_run, config%opt_sfc, config%opt_frz, &
         config%opt_inf, config%opt_rad, config%opt_alb, config%opt_snf, config%opt_tbot, config%opt_stc, &
         IZ0TLND, &
         T_PHY, QV_CURR, U_PHY, V_PHY, SWDOWN, GLW, &
         P8W, RAINBL, &
         TSK, HFX, QFX, LH, GRDFLX, SMSTAV, &
         SMSTOT,SFCRUNOFF, UDRUNOFF, ALBEDO, SNOWC, SMOIS, &
         SH2O, TSLB, SNOW, SNOWH, CANWAT, ACSNOM, &
         ACSNOW, EMISS, QSFC, &
         ISNOWXY, TVXY, TGXY, CANICEXY, CANLIQXY, EAHXY, &
         TAHXY, CMXY, CHXY, FWETXY, SNEQVOXY, ALBOLDXY, &
         QSNOWXY, WSLAKEXY, ZWTXY, WAXY, WTXY, TSNOXY, &
         ZSNSOXY,  SNICEXY, SNLIQXY, LFMASSXY, RTMASSXY, STMASSXY, &
         WOODXY, STBLCPXY, FASTCPXY, LAI, XSAIXY, TAUSSXY, &
         SMOISEQ, SMCWTDXY, DEEPRECHXY, RECHXY, &
         T2MVXY, T2MBXY, Q2MVXY, Q2MBXY, &
         TRADXY, NEEXY, GPPXY, NPPXY, FVEGXY, RUNSFXY, &
         RUNSBXY, ECANXY, EDIRXY, ETRANXY, FSAXY, FIRAXY, &
         APARXY, PSNXY, SAVXY, SAGXY,   RSSUNXY, RSSHAXY, &
         BGAPXY, WGAPXY, TGVXY, TGBXY, CHVXY, CHBXY, &
         SHGXY, SHCXY, SHBXY,  EVGXY, EVBXY, GHVXY, &
         GHBXY, IRGXY, IRCXY, IRBXY, TRXY, EVCXY, &
         CHLEAFXY, CHUCXY, CHV2XY, CHB2XY, &
         ids,ide, jds,jde, kds,kde, &
         ims,ime, jms,jme, kms,kme, &
         its,ite, jts,jte, kts,kte)

    call system_clock(count=count_after_sflx, count_rate=clock_rate)
    sflx_count_sum = sflx_count_sum + ( count_after_sflx - count_before_sflx )

    if(config%opt_run .eq. 5 .and. mod(ITIME,STEPWTD) .eq. 0)then
       call wrf_message('calling WTABLE' )

       !gmm update wtable from lateral flow and shed water to rivers
       call WTABLE_MMF_NOAHMP(&
            config%nsoil,  XLAND, XICE,       XICE_THRESHOLD, ISICE,    &
            ISLTYP,      SMOISEQ,    DZS,        WTDDT,                &
            FDEPTHXY,    AREAXY,     TERRAIN,    ISURBAN,    IVGTYP,   &
            RIVERCONDXY, RIVERBEDXY, EQZWT,      PEXPXY,               &
            SMOIS,       SH2O,       SMCWTDXY,   ZWTXY,                &
            QRFXY,       DEEPRECHXY, QSPRINGXY,                        &
            QSLATXY,     QRFSXY,     QSPRINGSXY, RECHXY,               &
            IDS,IDE, JDS,JDE, KDS,KDE,                                 &
            IMS,IME, JMS,JME, KMS,KME,                                 &
            ITS,ITE, JTS,JTE, KTS,KTE )

    endif
    ! END of surface_driver consistent code

    ! Output for history
    OUTPUT_FOR_HISTORY: if (config%output_timestep > 0) then
       if (mod((ITIME-1) * config%model_timestep, config%output_timestep) == 0) then

          call hrldas_output_prepare(trim(config%outdir), version, &
               config%output_timestep, llanduse, &
               ixfull, jxfull, ixpar, jxpar, xstartpar, ystartpar, &
               iswater, mapproj, dx, dy, truelat1, truelat2, cen_lon, &
               config%nsoil, nsnow, dzs, startdate, olddate, IVGTYP, ISLTYP)

          DEFINE_MODE_LOOP : do imode = 1, 2
             call hrldas_output_definemode(imode)
             ! For 3D arrays, we need to know whether the Z dimension is snow layers, or soil layers.
             ! Properties - Assigned or predicted
             call hrldas_output_add(IVGTYP     , "IVGTYP"  , "Dominant vegetation category"         , "category"              )
             call hrldas_output_add(ISLTYP     , "ISLTYP"  , "Dominant soil category"               , "category"              )
             call hrldas_output_add(FVEGXY     , "FVEG"    , "Green Vegetation Fraction"              , "-"                   )
             call hrldas_output_add(LAI        , "LAI"     , "Leaf area index"                      , "-"                     )
             call hrldas_output_add(XSAIXY     , "SAI"     , "Stem area index"                      , "-"                     )
             ! Forcing
             call hrldas_output_add(SWDOWN     , "SWFORC"  , "Shortwave forcing"                    , "W m{-2}"               )
             call hrldas_output_add(COSZEN     , "COSZ"    , "Cosine of zenith angle"                    , "W m{-2}"               )
             call hrldas_output_add(GLW        , "LWFORC"  , "Longwave forcing"                    , "W m{-2}"               )
             call hrldas_output_add(RAINBL_tmp     , "RAINRATE", "Precipitation rate"                   , "kg m{-2} s{-1}"        )
             ! Grid energy budget terms
             call hrldas_output_add(EMISS      , "EMISS"   , "Grid emissivity"                    , ""               )
             call hrldas_output_add(FSAXY      , "FSA"     , "Total absorbed SW radiation"          , "W m{-2}"               )
             call hrldas_output_add(FIRAXY     , "FIRA"    , "Total net LW radiation to atmosphere" , "W m{-2}"               )
             call hrldas_output_add(GRDFLX     , "GRDFLX"  , "Heat flux into the soil"              , "W m{-2}"               )
             call hrldas_output_add(HFX        , "HFX"     , "Total sensible heat to atmosphere"    , "W m{-2}"               )
             call hrldas_output_add(LH         , "LH"      , "Total latent heat to atmosphere"    , "W m{-2}"               )
             call hrldas_output_add(ECANXY     , "ECAN"    , "Canopy water evaporation rate"        , "kg m{-2} s{-1}"        )
             call hrldas_output_add(ETRANXY    , "ETRAN"   , "Transpiration rate"                   , "kg m{-2} s{-1}"        )
             call hrldas_output_add(EDIRXY     , "EDIR"    , "Direct from soil evaporation rate"    , "kg m{-2} s{-1}"        )
             call hrldas_output_add(ALBEDO     , "ALBEDO"  , "Surface albedo"                         , "-"                   )
             ! Grid water budget terms - in addition to above
             call hrldas_output_add(UDRUNOFF   , "UGDRNOFF", "Accumulated underground runoff"       , "mm"                    )
             call hrldas_output_add(SFCRUNOFF  , "SFCRNOFF", "Accumulatetd surface runoff"          , "mm"                    )
             call hrldas_output_add(CANLIQXY   , "CANLIQ"  , "Canopy liquid water content"          , "mm"                    )
             call hrldas_output_add(CANICEXY   , "CANICE"  , "Canopy ice water content"             , "mm"                    )
             call hrldas_output_add(ZWTXY      , "ZWT"     , "Depth to water table"                 , "m"                     )
             call hrldas_output_add(WAXY       , "WA"      , "Water in aquifer"                     , "kg m{-2}"              )
             call hrldas_output_add(WTXY       , "WT"      , "Water in aquifer and saturated soil"  , "kg m{-2}"              )
             ! Additional needed to close the canopy energy budget
             call hrldas_output_add(SAVXY      , "SAV"     , "Solar radiative heat flux absorbed by vegetation", "W m{-2}"    )
             call hrldas_output_add(TRXY       , "TR"      , "Transpiration heat"                     , "W m{-2}"             )
             call hrldas_output_add(EVCXY      , "EVC"     , "Canopy evap heat"                       , "W m{-2}"             )
             call hrldas_output_add(IRCXY      , "IRC"     , "Canopy net LW rad"                      , "W m{-2}"             )
             call hrldas_output_add(SHCXY      , "SHC"     , "Canopy sensible heat"                   , "W m{-2}"             )
             ! Additional needed to close the under canopy ground energy budget
             call hrldas_output_add(IRGXY      , "IRG"     , "Ground net LW rad"                      , "W m{-2}"             )
             call hrldas_output_add(SHGXY      , "SHG"     , "Ground sensible heat"                   , "W m{-2}"             )
             call hrldas_output_add(EVGXY      , "EVG"     , "Ground evap heat"                       , "W m{-2}"             )
             call hrldas_output_add(GHVXY      , "GHV"     , "Ground heat flux + to soil vegetated"   , "W m{-2}"             )
             ! Needed to close the bare ground energy budget
             call hrldas_output_add(SAGXY      , "SAG"     , "Solar radiative heat flux absorbed by ground", "W m{-2}"        )
             call hrldas_output_add(IRBXY      , "IRB"     , "Net LW rad to atm bare"                 , "W m{-2}"             )
             call hrldas_output_add(SHBXY      , "SHB"     , "Sensible heat to atm bare"              , "W m{-2}"             )
             call hrldas_output_add(EVBXY      , "EVB"     , "Evaporation heat to atm bare"           , "W m{-2}"             )
             call hrldas_output_add(GHBXY      , "GHB"     , "Ground heat flux + to soil bare"        , "W m{-2}"             )
             ! Above-soil temperatures
             call hrldas_output_add(TRADXY     , "TRAD"    , "Surface radiative temperature"        , "K"                     )
             call hrldas_output_add(TGXY       , "TG"      , "Ground temperature"                   , "K"                     )
             call hrldas_output_add(TVXY       , "TV"      , "Vegetation temperature"               , "K"                     )
             call hrldas_output_add(TAHXY      , "TAH"     , "Canopy air temperature"               , "K"                     )
             call hrldas_output_add(TGVXY      , "TGV"     , "Ground surface Temp vegetated"          , "K"                   )
             call hrldas_output_add(TGBXY      , "TGB"     , "Ground surface Temp bare"               , "K"                   )
             call hrldas_output_add(T2MVXY     , "T2MV"    , "2m Air Temp vegetated"                  , "K"                   )
             call hrldas_output_add(T2MBXY     , "T2MB"    , "2m Air Temp bare"                       , "K"                   )
             ! Above-soil moisture
             call hrldas_output_add(Q2MVXY     , "Q2MV"    , "2m mixing ratio vegetated"              , "kg/kg"               )
             call hrldas_output_add(Q2MBXY     , "Q2MB"    , "2m mixing ratio bare"                   , "kg/kg"               )
             call hrldas_output_add(EAHXY      , "EAH"     , "Canopy air vapor pressure"            , "Pa"                    )
             call hrldas_output_add(FWETXY     , "FWET"    , "Wetted or snowed fraction of canopy"  , "fraction"              )
             ! Snow and soil - 3D terms
             call hrldas_output_add(ZSNSOXY(:,-nsnow+1:0,:),  "ZSNSO_SN" , "Snow layer depths from snow surface", "m", "SNOW")
             call hrldas_output_add(SNICEXY    , "SNICE"   , "Snow layer ice"                       , "mm"             , "SNOW")
             call hrldas_output_add(SNLIQXY    , "SNLIQ"   , "Snow layer liquid water"              , "mm"             , "SNOW")
             call hrldas_output_add(TSLB       , "SOIL_T"  , "soil temperature"                     , "K"              , "SOIL")
             call hrldas_output_add(SMOIS      , "SOIL_M"  , "volumetric soil moisture"             , "m{3} m{-3}"     , "SOIL")
             call hrldas_output_add(SH2O       , "SOIL_W"  , "liquid volumetric soil moisture"      , "m3 m-3"         , "SOIL")
             call hrldas_output_add(TSNOXY     , "SNOW_T"  , "snow temperature"                     , "K"              , "SNOW")
             ! Snow - 2D terms
             call hrldas_output_add(SNOWH      , "SNOWH"   , "Snow depth"                           , "m"                     )
             call hrldas_output_add(SNOW       , "SNEQV"   , "Snow water equivalent"                , "kg m{-2}"              )
             call hrldas_output_add(QSNOWXY    , "QSNOW"   , "Snowfall rate"                        , "mm s{-1}"              )
             call hrldas_output_add(ISNOWXY    , "ISNOW"   , "Number of snow layers"                , "count"                 )
             call hrldas_output_add(SNOWC      , "FSNO"    , "Snow-cover fraction on the ground"      , ""                    )
             call hrldas_output_add(ACSNOW     , "ACSNOW"  , "accumulated snow fall"                  , "mm"                  )
             call hrldas_output_add(ACSNOM     , "ACSNOM"  , "accumulated melting water out of snow bottom" , "mm"            )
             ! Exchange coefficients
             call hrldas_output_add(CMXY       , "CM"      , "Momentum drag coefficient"            , ""                      )
             call hrldas_output_add(CHXY       , "CH"      , "Sensible heat exchange coefficient"   , ""                      )
             call hrldas_output_add(CHVXY      , "CHV"     , "Exchange coefficient vegetated"         , "m s{-1}"             )
             call hrldas_output_add(CHBXY      , "CHB"     , "Exchange coefficient bare"              , "m s{-1}"             )
             call hrldas_output_add(CHLEAFXY   , "CHLEAF"  , "Exchange coefficient leaf"              , "m s{-1}"             )
             call hrldas_output_add(CHUCXY     , "CHUC"    , "Exchange coefficient bare"              , "m s{-1}"             )
             call hrldas_output_add(CHV2XY     , "CHV2"    , "Exchange coefficient 2-meter vegetated" , "m s{-1}"             )
             call hrldas_output_add(CHB2XY     , "CHB2"    , "Exchange coefficient 2-meter bare"      , "m s{-1}"             )
             ! Carbon allocation model
             call hrldas_output_add(LFMASSXY   , "LFMASS"  , "Leaf mass"                            , "g m{-2}"               )
             call hrldas_output_add(RTMASSXY   , "RTMASS"  , "Mass of fine roots"                   , "g m{-2}"               )
             call hrldas_output_add(STMASSXY   , "STMASS"  , "Stem mass"                            , "g m{-2}"               )
             call hrldas_output_add(WOODXY     , "WOOD"    , "Mass of wood and woody roots"         , "g m{-2}"               )
             call hrldas_output_add(STBLCPXY   , "STBLCP"  , "Stable carbon in deep soil"           , "g m{-2}"               )
             call hrldas_output_add(FASTCPXY   , "FASTCP"  , "Short-lived carbon in shallow soil"   , "g m{-2}"               )
             call hrldas_output_add(NEEXY      , "NEE"     , "Net ecosystem exchange"                 , "g m{-2} s{-1} CO2"   )
             call hrldas_output_add(GPPXY      , "GPP"     , "Net instantaneous assimilation"         , "g m{-2} s{-1} C"     )
             call hrldas_output_add(NPPXY      , "NPP"     , "Net primary productivity"               , "g m{-2} s{-1} C"     )
             call hrldas_output_add(PSNXY      , "PSN"     , "Total photosynthesis"                   , "umol CO@ m{-2} s{-1}")
             call hrldas_output_add(APARXY     , "APAR"    , "Photosynthesis active energy by canopy" , "W m{-2}"             )

             ! Carbon allocation model
             if(config%opt_run == 5) then
                call hrldas_output_add(SMCWTDXY   , "SMCWTD"   , "Leaf mass"                            , "g m{-2}"               )
                call hrldas_output_add(RECHXY     , "RECH"     , "Mass of fine roots"                   , "g m{-2}"               )
                call hrldas_output_add(QRFSXY     , "QRFS"     , "Stem mass"                            , "g m{-2}"               )
                call hrldas_output_add(QSPRINGSXY , "QSPRINGS" , "Mass of wood and woody roots"         , "g m{-2}"               )
                call hrldas_output_add(QSLATXY    , "QSLAT"    , "Stable carbon in deep soil"           , "g m{-2}"               )
             endif

          enddo DEFINE_MODE_LOOP

          call hrldas_output_finalize(1)

       endif
    endif OUTPUT_FOR_HISTORY

    if (IVGTYP(xstart,ystart)==ISWATER) then
       write(*,'(" ***DATE=", A19)', advance="NO") olddate
    else
       write(*,'(" ***DATE=", A19, 6F10.5)', advance="NO") olddate, &
            TSLB(xstart,1,ystart), LAI(xstart,ystart)
    endif

    ! Update the time
    call geth_newdate(newdate, olddate, config%model_timestep)
    olddate = newdate


    ! update the timer
    call system_clock(count=clock_count_2, count_rate=clock_rate)
    timing_sum = timing_sum + float(clock_count_2-clock_count_1)/float(clock_rate)
    write(*,'("    Timing: ",f6.2," Cumulative:  ", f10.2, "  SFLX: ", f6.2 )') &
         float(clock_count_2-clock_count_1)/float(clock_rate), &
         timing_sum, real(sflx_count_sum) / real(clock_rate)
    clock_count_1 = clock_count_2

    ! jlm - this is the propsed, new code chunk essentially copied from above
    if ((config%restart_timestep .gt. 0) .and. (mod(ITIME, int(config%restart_timestep/config%model_timestep)) == 0)) then

       call hrldas_noahmp_vars_write_restart()

    endif

    ITIMESTEP = ITIMESTEP + 1
  end subroutine land_driver_exe

  subroutine hrldas_noahmp_vars_write_restart()
    implicit none
    ! character(len=*),intent(in)  :: filename
    call hrldas_restart_prepare(trim(config%resdir), version, llanduse, olddate, startdate, &
         ixfull, jxfull, ixpar, jxpar, xstartpar, ystartpar, &
         config%nsoil, nsnow, dx, dy, truelat1, truelat2, mapproj, &
         cen_lon, iswater, ivgtyp)

    ! Since the restart files are not really for user consumption, the units and description are
    ! a little superfluous.  I've made them optional arguments.  If not present, they both default to "-".
    ! jlm - for data assimilation, the restarts are key (at least until we do the
    ! jlm - assimilation in memory). The restarts are also the main part of the DA
    ! jlm - diagnostics. Having units would be nice and having missing values for each
    ! jlm - variable is also desirable.

    call hrldas_restart_add(TSLB      , "SOIL_T", layers="SOIL")
    call hrldas_restart_add(TSNOXY    , "SNOW_T", layers="SNOW")
    call hrldas_restart_add(SMOIS     , "SMC"   , layers="SOIL")
    call hrldas_restart_add(SH2O      , "SH2O"  , layers="SOIL")
    call hrldas_restart_add(ZSNSOXY   , "ZSNSO" , layers="SOSN")
    call hrldas_restart_add(SNICEXY   , "SNICE" , layers="SNOW")
    call hrldas_restart_add(SNLIQXY   , "SNLIQ" , layers="SNOW")
    call hrldas_restart_add(QSNOWXY   , "QSNOW" )
    call hrldas_restart_add(FWETXY    , "FWET"  )
    call hrldas_restart_add(SNEQVOXY  , "SNEQVO")
    call hrldas_restart_add(EAHXY     , "EAH"   )
    call hrldas_restart_add(TAHXY     , "TAH"   )
    call hrldas_restart_add(ALBOLDXY  , "ALBOLD")
    call hrldas_restart_add(CMXY      , "CM"    )
    call hrldas_restart_add(CHXY      , "CH"    )
    call hrldas_restart_add(ISNOWXY   , "ISNOW" )
    call hrldas_restart_add(CANLIQXY  , "CANLIQ")
    call hrldas_restart_add(CANICEXY  , "CANICE")
    call hrldas_restart_add(SNOW      , "SNEQV" )
    call hrldas_restart_add(SNOWH     , "SNOWH" )
    call hrldas_restart_add(TVXY      , "TV"    )
    call hrldas_restart_add(TGXY      , "TG"    )
    call hrldas_restart_add(ZWTXY     , "ZWT"   )
    call hrldas_restart_add(WAXY      , "WA"    )
    call hrldas_restart_add(WTXY      , "WT"    )
    call hrldas_restart_add(WSLAKEXY  , "WSLAKE")
    call hrldas_restart_add(LFMASSXY  , "LFMASS")
    call hrldas_restart_add(RTMASSXY  , "RTMASS")
    call hrldas_restart_add(STMASSXY  , "STMASS")
    call hrldas_restart_add(WOODXY    , "WOOD"  )
    call hrldas_restart_add(STBLCPXY  , "STBLCP")
    call hrldas_restart_add(FASTCPXY  , "FASTCP")
    call hrldas_restart_add(LAI       , "LAI"   )
    call hrldas_restart_add(XSAIXY    , "SAI"   )
    call hrldas_restart_add(VEGFRA    , "FPAR"  )
    call hrldas_restart_add(GVFMIN    , "GVFMIN")
    call hrldas_restart_add(GVFMAX    , "GVFMAX")
    call hrldas_restart_add(SHDMAX    , "SHDMAX")
    call hrldas_restart_add(ACSNOM    , "ACMELT")
    call hrldas_restart_add(ACSNOW    , "ACSNOW")
    call hrldas_restart_add(TAUSSXY   , "TAUSS" )
    call hrldas_restart_add(QSFC      , "QSFC"  )
    call hrldas_restart_add(SFCRUNOFF , "SFCRUNOFF")
    call hrldas_restart_add(UDRUNOFF  , "UDRUNOFF" )

    ! below for opt_run = 5
    call hrldas_restart_add(SMOISEQ   , "SMOISEQ"  , layers="SOIL"  )
    call hrldas_restart_add(AREAXY    , "AREAXY"     )
    call hrldas_restart_add(SMCWTDXY  , "SMCWTDXY"   )
    call hrldas_restart_add(DEEPRECHXY, "DEEPRECHXY" )
    call hrldas_restart_add(QSLATXY   , "QSLATXY"    )
    call hrldas_restart_add(QRFSXY    , "QRFSXY"     )
    call hrldas_restart_add(QSPRINGSXY, "QSPRINGSXY" )
    call hrldas_restart_add(RECHXY    , "RECHXY"     )
    call hrldas_restart_add(QRFXY     , "QRFXY"      )
    call hrldas_restart_add(QSPRINGXY , "QSPRINGXY"  )
    call hrldas_restart_add(FDEPTHXY , "FDEPTHXY"  )
    call hrldas_restart_add(RIVERCONDXY , "RIVERCONDXY"  )
    call hrldas_restart_add(RIVERBEDXY , "RIVERBEDXY"  )
    call hrldas_restart_add(EQZWT , "EQZWT"  )
    call hrldas_restart_add(PEXPXY , "PEXPXY"  )
    call hrldas_restart_finalize()

  end subroutine hrldas_noahmp_vars_write_restart

  subroutine hrldas_noahmp_vars_write_output(filename)
    character(len=*), intent(in) :: filename

    call hrldas_output_prepare(trim(config%outdir), version, &
         config%output_timestep, llanduse, &
         ixfull, jxfull, ixpar, jxpar, xstartpar, ystartpar, &
         iswater, mapproj, dx, dy, truelat1, truelat2, cen_lon, &
         config%nsoil, nsnow, dzs, startdate, olddate, IVGTYP, ISLTYP)

    do imode = 1, 2
       call hrldas_output_definemode(imode)
       ! For 3D arrays, we need to know whether the Z dimension is snow layers, or soil layers.
       ! Properties - Assigned or predicted
       call hrldas_output_add(IVGTYP     , "IVGTYP"  , "Dominant vegetation category"         , "category"              )
       call hrldas_output_add(ISLTYP     , "ISLTYP"  , "Dominant soil category"               , "category"              )
       call hrldas_output_add(FVEGXY     , "FVEG"    , "Green Vegetation Fraction"              , "-"                   )
       call hrldas_output_add(LAI        , "LAI"     , "Leaf area index"                      , "-"                     )
       call hrldas_output_add(XSAIXY     , "SAI"     , "Stem area index"                      , "-"                     )
       ! Forcing
       call hrldas_output_add(SWDOWN     , "SWFORC"  , "Shortwave forcing"                    , "W m{-2}"               )
       call hrldas_output_add(COSZEN     , "COSZ"    , "Cosine of zenith angle"                    , "W m{-2}"               )
       call hrldas_output_add(GLW        , "LWFORC"  , "Longwave forcing"                    , "W m{-2}"               )
       call hrldas_output_add(RAINBL_tmp     , "RAINRATE", "Precipitation rate"                   , "kg m{-2} s{-1}"        )
       ! Grid energy budget terms
       call hrldas_output_add(EMISS      , "EMISS"   , "Grid emissivity"                    , ""               )
       call hrldas_output_add(FSAXY      , "FSA"     , "Total absorbed SW radiation"          , "W m{-2}"               )
       call hrldas_output_add(FIRAXY     , "FIRA"    , "Total net LW radiation to atmosphere" , "W m{-2}"               )
       call hrldas_output_add(GRDFLX     , "GRDFLX"  , "Heat flux into the soil"              , "W m{-2}"               )
       call hrldas_output_add(HFX        , "HFX"     , "Total sensible heat to atmosphere"    , "W m{-2}"               )
       call hrldas_output_add(LH         , "LH"      , "Total latent heat to atmosphere"    , "W m{-2}"               )
       call hrldas_output_add(ECANXY     , "ECAN"    , "Canopy water evaporation rate"        , "kg m{-2} s{-1}"        )
       call hrldas_output_add(ETRANXY    , "ETRAN"   , "Transpiration rate"                   , "kg m{-2} s{-1}"        )
       call hrldas_output_add(EDIRXY     , "EDIR"    , "Direct from soil evaporation rate"    , "kg m{-2} s{-1}"        )
       call hrldas_output_add(ALBEDO     , "ALBEDO"  , "Surface albedo"                         , "-"                   )
       ! Grid water budget terms - in addition to above
       call hrldas_output_add(UDRUNOFF   , "UGDRNOFF", "Accumulated underground runoff"       , "mm"                    )
       call hrldas_output_add(SFCRUNOFF  , "SFCRNOFF", "Accumulatetd surface runoff"          , "mm"                    )
       call hrldas_output_add(CANLIQXY   , "CANLIQ"  , "Canopy liquid water content"          , "mm"                    )
       call hrldas_output_add(CANICEXY   , "CANICE"  , "Canopy ice water content"             , "mm"                    )
       call hrldas_output_add(ZWTXY      , "ZWT"     , "Depth to water table"                 , "m"                     )
       call hrldas_output_add(WAXY       , "WA"      , "Water in aquifer"                     , "kg m{-2}"              )
       call hrldas_output_add(WTXY       , "WT"      , "Water in aquifer and saturated soil"  , "kg m{-2}"              )
       ! Additional needed to close the canopy energy budget
       call hrldas_output_add(SAVXY      , "SAV"     , "Solar radiative heat flux absorbed by vegetation", "W m{-2}"    )
       call hrldas_output_add(TRXY       , "TR"      , "Transpiration heat"                     , "W m{-2}"             )
       call hrldas_output_add(EVCXY      , "EVC"     , "Canopy evap heat"                       , "W m{-2}"             )
       call hrldas_output_add(IRCXY      , "IRC"     , "Canopy net LW rad"                      , "W m{-2}"             )
       call hrldas_output_add(SHCXY      , "SHC"     , "Canopy sensible heat"                   , "W m{-2}"             )
       ! Additional needed to close the under canopy ground energy budget
       call hrldas_output_add(IRGXY      , "IRG"     , "Ground net LW rad"                      , "W m{-2}"             )
       call hrldas_output_add(SHGXY      , "SHG"     , "Ground sensible heat"                   , "W m{-2}"             )
       call hrldas_output_add(EVGXY      , "EVG"     , "Ground evap heat"                       , "W m{-2}"             )
       call hrldas_output_add(GHVXY      , "GHV"     , "Ground heat flux + to soil vegetated"   , "W m{-2}"             )
       ! Needed to close the bare ground energy budget
       call hrldas_output_add(SAGXY      , "SAG"     , "Solar radiative heat flux absorbed by ground", "W m{-2}"        )
       call hrldas_output_add(IRBXY      , "IRB"     , "Net LW rad to atm bare"                 , "W m{-2}"             )
       call hrldas_output_add(SHBXY      , "SHB"     , "Sensible heat to atm bare"              , "W m{-2}"             )
       call hrldas_output_add(EVBXY      , "EVB"     , "Evaporation heat to atm bare"           , "W m{-2}"             )
       call hrldas_output_add(GHBXY      , "GHB"     , "Ground heat flux + to soil bare"        , "W m{-2}"             )
       ! Above-soil temperatures
       call hrldas_output_add(TRADXY     , "TRAD"    , "Surface radiative temperature"        , "K"                     )
       call hrldas_output_add(TGXY       , "TG"      , "Ground temperature"                   , "K"                     )
       call hrldas_output_add(TVXY       , "TV"      , "Vegetation temperature"               , "K"                     )
       call hrldas_output_add(TAHXY      , "TAH"     , "Canopy air temperature"               , "K"                     )
       call hrldas_output_add(TGVXY      , "TGV"     , "Ground surface Temp vegetated"          , "K"                   )
       call hrldas_output_add(TGBXY      , "TGB"     , "Ground surface Temp bare"               , "K"                   )
       call hrldas_output_add(T2MVXY     , "T2MV"    , "2m Air Temp vegetated"                  , "K"                   )
       call hrldas_output_add(T2MBXY     , "T2MB"    , "2m Air Temp bare"                       , "K"                   )
       ! Above-soil moisture
       call hrldas_output_add(Q2MVXY     , "Q2MV"    , "2m mixing ratio vegetated"              , "kg/kg"               )
       call hrldas_output_add(Q2MBXY     , "Q2MB"    , "2m mixing ratio bare"                   , "kg/kg"               )
       call hrldas_output_add(EAHXY      , "EAH"     , "Canopy air vapor pressure"            , "Pa"                    )
       call hrldas_output_add(FWETXY     , "FWET"    , "Wetted or snowed fraction of canopy"  , "fraction"              )
       ! Snow and soil - 3D terms
       call hrldas_output_add(ZSNSOXY(:,-nsnow+1:0,:),  "ZSNSO_SN" , "Snow layer depths from snow surface", "m", "SNOW")
       call hrldas_output_add(SNICEXY    , "SNICE"   , "Snow layer ice"                       , "mm"             , "SNOW")
       call hrldas_output_add(SNLIQXY    , "SNLIQ"   , "Snow layer liquid water"              , "mm"             , "SNOW")
       call hrldas_output_add(TSLB       , "SOIL_T"  , "soil temperature"                     , "K"              , "SOIL")
       call hrldas_output_add(SMOIS      , "SOIL_M"  , "volumetric soil moisture"             , "m{3} m{-3}"     , "SOIL")
       call hrldas_output_add(SH2O       , "SOIL_W"  , "liquid volumetric soil moisture"      , "m3 m-3"         , "SOIL")
       call hrldas_output_add(TSNOXY     , "SNOW_T"  , "snow temperature"                     , "K"              , "SNOW")
       ! Snow - 2D terms
       call hrldas_output_add(SNOWH      , "SNOWH"   , "Snow depth"                           , "m"                     )
       call hrldas_output_add(SNOW       , "SNEQV"   , "Snow water equivalent"                , "kg m{-2}"              )
       call hrldas_output_add(QSNOWXY    , "QSNOW"   , "Snowfall rate"                        , "mm s{-1}"              )
       call hrldas_output_add(ISNOWXY    , "ISNOW"   , "Number of snow layers"                , "count"                 )
       call hrldas_output_add(SNOWC      , "FSNO"    , "Snow-cover fraction on the ground"      , ""                    )
       call hrldas_output_add(ACSNOW     , "ACSNOW"  , "accumulated snow fall"                  , "mm"                  )
       call hrldas_output_add(ACSNOM     , "ACSNOM"  , "accumulated melting water out of snow bottom" , "mm"            )
       ! Exchange coefficients
       call hrldas_output_add(CMXY       , "CM"      , "Momentum drag coefficient"            , ""                      )
       call hrldas_output_add(CHXY       , "CH"      , "Sensible heat exchange coefficient"   , ""                      )
       call hrldas_output_add(CHVXY      , "CHV"     , "Exchange coefficient vegetated"         , "m s{-1}"             )
       call hrldas_output_add(CHBXY      , "CHB"     , "Exchange coefficient bare"              , "m s{-1}"             )
       call hrldas_output_add(CHLEAFXY   , "CHLEAF"  , "Exchange coefficient leaf"              , "m s{-1}"             )
       call hrldas_output_add(CHUCXY     , "CHUC"    , "Exchange coefficient bare"              , "m s{-1}"             )
       call hrldas_output_add(CHV2XY     , "CHV2"    , "Exchange coefficient 2-meter vegetated" , "m s{-1}"             )
       call hrldas_output_add(CHB2XY     , "CHB2"    , "Exchange coefficient 2-meter bare"      , "m s{-1}"             )
       ! Carbon allocation model
       call hrldas_output_add(LFMASSXY   , "LFMASS"  , "Leaf mass"                            , "g m{-2}"               )
       call hrldas_output_add(RTMASSXY   , "RTMASS"  , "Mass of fine roots"                   , "g m{-2}"               )
       call hrldas_output_add(STMASSXY   , "STMASS"  , "Stem mass"                            , "g m{-2}"               )
       call hrldas_output_add(WOODXY     , "WOOD"    , "Mass of wood and woody roots"         , "g m{-2}"               )
       call hrldas_output_add(STBLCPXY   , "STBLCP"  , "Stable carbon in deep soil"           , "g m{-2}"               )
       call hrldas_output_add(FASTCPXY   , "FASTCP"  , "Short-lived carbon in shallow soil"   , "g m{-2}"               )
       call hrldas_output_add(NEEXY      , "NEE"     , "Net ecosystem exchange"                 , "g m{-2} s{-1} CO2"   )
       call hrldas_output_add(GPPXY      , "GPP"     , "Net instantaneous assimilation"         , "g m{-2} s{-1} C"     )
       call hrldas_output_add(NPPXY      , "NPP"     , "Net primary productivity"               , "g m{-2} s{-1} C"     )
       call hrldas_output_add(PSNXY      , "PSN"     , "Total photosynthesis"                   , "umol CO@ m{-2} s{-1}")
       call hrldas_output_add(APARXY     , "APAR"    , "Photosynthesis active energy by canopy" , "W m{-2}"             )

       ! Carbon allocation model
       if(config%opt_run == 5) then
          call hrldas_output_add(SMCWTDXY   , "SMCWTD"   , "Leaf mass"                            , "g m{-2}"               )
          call hrldas_output_add(RECHXY     , "RECH"     , "Mass of fine roots"                   , "g m{-2}"               )
          call hrldas_output_add(QRFSXY     , "QRFS"     , "Stem mass"                            , "g m{-2}"               )
          call hrldas_output_add(QSPRINGSXY , "QSPRINGS" , "Mass of wood and woody roots"         , "g m{-2}"               )
          call hrldas_output_add(QSLATXY    , "QSLAT"    , "Stable carbon in deep soil"           , "g m{-2}"               )
       endif

    enddo

    call hrldas_output_finalize(1)

  end subroutine hrldas_noahmp_vars_write_output

end module module_hrldas_noahmp_driver


function wrf_dm_on_monitor()
  implicit none
  logical :: wrf_dm_on_monitor
  wrf_dm_on_monitor = .true.
end function wrf_dm_on_monitor


subroutine CALC_DECLIN(NOWDATE, LATITUDE, LONGITUDE, COSZ, JULIAN)
  use MODULE_DATE_UTILITIES
  implicit none
  real, parameter :: DEGRAD = 3.14159265/180.
  real, parameter :: DPD    = 360./365.
  ! !ARGUMENTS:
  character(LEN=19), intent(IN)  :: NOWDATE    ! YYYY-MM-DD_HH:MM:SS
  real,              intent(IN)  :: LATITUDE
  real,              intent(IN)  :: LONGITUDE
  real,              intent(OUT) :: COSZ
  real,              intent(OUT) :: JULIAN
  real                           :: HRANG
  real                           :: DECLIN
  real                           :: OBECL
  real                           :: SINOB
  real                           :: SXLONG
  real                           :: ARG
  real                           :: TLOCTIM
  integer                        :: IDAY
  integer                        :: IHOUR
  integer                        :: IMINUTE
  integer                        :: ISECOND

  call GETH_IDTS(NOWDATE(1:10), NOWDATE(1:4)//"-01-01", IDAY)
  read(NOWDATE(12:13), *) IHOUR
  read(NOWDATE(15:16), *) IMINUTE
  read(NOWDATE(18:19), *) ISECOND
  JULIAN = real(IDAY) + real(IHOUR)/24.

  ! FOR SHORT WAVE RADIATION

  DECLIN=0.

  !-----OBECL : OBLIQUITY = 23.5 DEGREE.

  OBECL=23.5*DEGRAD
  SINOB=sin(OBECL)

  !-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:

  if(JULIAN.ge.80.)SXLONG=DPD*(JULIAN-80.)*DEGRAD
  if(JULIAN.lt.80.)SXLONG=DPD*(JULIAN+285.)*DEGRAD
  ARG=SINOB*sin(SXLONG)
  DECLIN=asin(ARG)

  TLOCTIM = real(IHOUR) + real(IMINUTE)/60.0 + real(ISECOND)/3600.0 + LONGITUDE/15.0 ! LOCAL TIME IN HOURS
  TLOCTIM = AMOD(TLOCTIM+24.0, 24.0)
  HRANG=15.*(TLOCTIM-12.)*DEGRAD
  COSZ=sin(LATITUDE*DEGRAD)*sin(DECLIN)+cos(LATITUDE*DEGRAD)*cos(DECLIN)*cos(HRANG)

end subroutine CALC_DECLIN
