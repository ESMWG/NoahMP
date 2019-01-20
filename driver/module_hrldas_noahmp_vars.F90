module module_hrldas_noahmp_vars
  use module_hrldas_noahmp_namelist, only: config
  use module_hrldas_netcdf_io
  implicit none

  real, parameter :: undefined_real = -1.0e20
  integer, parameter :: undefined_integer = -999

  integer :: xstart, xend
  integer :: ystart, yend
  integer :: ids,ide, jds,jde, kds,kde ! d -> domain
  integer :: ims,ime, jms,jme, kms,kme ! m -> memory
  integer :: its,ite, jts,jte, kts,kte ! t -> tile

  ! domain
  integer, parameter :: NSNOW = 3    ! number of snow layers fixed to 3
  real, allocatable, dimension(:,:) :: lat2d ! latitude [degree]
  real, allocatable, dimension(:,:) :: lon2d ! longitude [degree]
  real, allocatable, dimension(:,:,:)  :: DZ8W      ! thickness of atmo layers [m]
  real, allocatable, dimension(:)      :: DZS       ! thickness of soil layers [m]

  ! IN only (as defined in WRF)
  real,    allocatable, dimension(:,:)    :: COSZEN    ! cosine zenith angle
  integer, allocatable, dimension(:,:)    :: IVGTYP    ! vegetation type
  integer, allocatable, dimension(:,:)    :: ISLTYP    ! soil type
  real,    allocatable, dimension(:,:)    :: VEGFRA    ! vegetation fraction []
  real,    allocatable, dimension(:,:)    :: SHDMAX    ! annual max vegetation fraction []
  real,    allocatable, dimension(:,:)    :: TMN       ! deep soil temperature [K]
  real,    allocatable, dimension(:,:)    :: XLAND     ! =2 ocean; =1 land/seaice
  real,    allocatable, dimension(:,:)    :: XICE      ! fraction of grid that is seaice
  real,    allocatable, dimension(:,:,:)  :: T_PHY     ! 3D atmospheric temperature valid at mid-levels [K]
  real,    allocatable, dimension(:,:,:)  :: QV_CURR   ! 3D water vapor mixing ratio [kg/kg_dry]
  real,    allocatable, dimension(:,:,:)  :: U_PHY     ! 3D U wind component [m/s]
  real,    allocatable, dimension(:,:,:)  :: V_PHY     ! 3D V wind component [m/s]
  real,    allocatable, dimension(:,:)    :: SWDOWN    ! solar down at surface [W m-2]
  real,    allocatable, dimension(:,:)    :: GLW       ! longwave down at surface [W m-2]
  real,    allocatable, dimension(:,:,:)  :: P8W       ! 3D pressure, valid at interface [Pa]
  real,    allocatable, dimension(:,:)    :: RAINBL, RAINBL_tmp    ! precipitation entering land model [mm]

  ! INOUT (with generic LSM equivalent) (as defined in WRF)
  real,    allocatable, dimension(:,:)    :: TSK       ! surface radiative temperature [K]
  real,    allocatable, dimension(:,:)    :: HFX       ! sensible heat flux [W m-2]
  real,    allocatable, dimension(:,:)    :: QFX       ! latent heat flux [kg s-1 m-2]
  real,    allocatable, dimension(:,:)    :: LH        ! latent heat flux [W m-2]
  real,    allocatable, dimension(:,:)    :: GRDFLX    ! ground/snow heat flux [W m-2]
  real,    allocatable, dimension(:,:)    :: SMSTAV    ! soil moisture avail. [not used]
  real,    allocatable, dimension(:,:)    :: SMSTOT    ! total soil water [mm][not used]
  real,    allocatable, dimension(:,:)    :: SFCRUNOFF ! accumulated surface runoff [m]
  real,    allocatable, dimension(:,:)    :: UDRUNOFF  ! accumulated sub-surface runoff [m]
  real,    allocatable, dimension(:,:)    :: ALBEDO    ! total grid albedo []
  real,    allocatable, dimension(:,:)    :: SNOWC     ! snow cover fraction []
  real,    allocatable, dimension(:,:,:)  :: SMOISEQ   ! volumetric soil moisture [m3/m3]
  real,    allocatable, dimension(:,:,:)  :: SMOIS     ! volumetric soil moisture [m3/m3]
  real,    allocatable, dimension(:,:,:)  :: SH2O      ! volumetric liquid soil moisture [m3/m3]
  real,    allocatable, dimension(:,:,:)  :: TSLB      ! soil temperature [K]
  real,    allocatable, dimension(:,:)    :: SNOW      ! snow water equivalent [mm]
  real,    allocatable, dimension(:,:)    :: SNOWH     ! physical snow depth [m]
  real,    allocatable, dimension(:,:)    :: CANWAT    ! total canopy water + ice [mm]
  real,    allocatable, dimension(:,:)    :: ACSNOM    ! accumulated snow melt leaving pack
  real,    allocatable, dimension(:,:)    :: ACSNOW    ! accumulated snow on grid
  real,    allocatable, dimension(:,:)    :: EMISS     ! surface bulk emissivity
  real,    allocatable, dimension(:,:)    :: QSFC      ! bulk surface specific humidity

  ! INOUT (with no Noah LSM equivalent) (as defined in WRF)
  integer, allocatable, dimension(:,:)    :: ISNOWXY   ! actual no. of snow layers
  real,    allocatable, dimension(:,:)    :: TVXY      ! vegetation leaf temperature
  real,    allocatable, dimension(:,:)    :: TGXY      ! bulk ground surface temperature
  real,    allocatable, dimension(:,:)    :: CANICEXY  ! canopy-intercepted ice (mm)
  real,    allocatable, dimension(:,:)    :: CANLIQXY  ! canopy-intercepted liquid water (mm)
  real,    allocatable, dimension(:,:)    :: EAHXY     ! canopy air vapor pressure (pa)
  real,    allocatable, dimension(:,:)    :: TAHXY     ! canopy air temperature (k)
  real,    allocatable, dimension(:,:)    :: CMXY      ! bulk momentum drag coefficient
  real,    allocatable, dimension(:,:)    :: CHXY      ! bulk sensible heat exchange coefficient
  real,    allocatable, dimension(:,:)    :: FWETXY    ! wetted or snowed fraction of the canopy (-)
  real,    allocatable, dimension(:,:)    :: SNEQVOXY  ! snow mass at last time step(mm h2o)
  real,    allocatable, dimension(:,:)    :: ALBOLDXY  ! snow albedo at last time step (-)
  real,    allocatable, dimension(:,:)    :: QSNOWXY   ! snowfall on the ground [mm/s]
  real,    allocatable, dimension(:,:)    :: WSLAKEXY  ! lake water storage [mm]
  real,    allocatable, dimension(:,:)    :: ZWTXY     ! water table depth [m]
  real,    allocatable, dimension(:,:)    :: WAXY      ! water in the "aquifer" [mm]
  real,    allocatable, dimension(:,:)    :: WTXY      ! groundwater storage [mm]
  real,    allocatable, dimension(:,:)    :: SMCWTDXY  ! groundwater storage [mm]
  real,    allocatable, dimension(:,:)    :: DEEPRECHXY! groundwater storage [mm]
  real,    allocatable, dimension(:,:)    :: RECHXY    ! groundwater storage [mm]
  real,    allocatable, dimension(:,:,:)  :: TSNOXY    ! snow temperature [K]
  real,    allocatable, dimension(:,:,:)  :: ZSNSOXY   ! snow layer depth [m]
  real,    allocatable, dimension(:,:,:)  :: SNICEXY   ! snow layer ice [mm]
  real,    allocatable, dimension(:,:,:)  :: SNLIQXY   ! snow layer liquid water [mm]
  real,    allocatable, dimension(:,:)    :: LFMASSXY  ! leaf mass [g/m2]
  real,    allocatable, dimension(:,:)    :: RTMASSXY  ! mass of fine roots [g/m2]
  real,    allocatable, dimension(:,:)    :: STMASSXY  ! stem mass [g/m2]
  real,    allocatable, dimension(:,:)    :: WOODXY    ! mass of wood (incl. woody roots) [g/m2]
  real,    allocatable, dimension(:,:)    :: STBLCPXY  ! stable carbon in deep soil [g/m2]
  real,    allocatable, dimension(:,:)    :: FASTCPXY  ! short-lived carbon, shallow soil [g/m2]
  real,    allocatable, dimension(:,:)    :: LAI       ! leaf area index
  real,    allocatable, dimension(:,:)    :: LAI_tmp   ! leaf area index
  real,    allocatable, dimension(:,:)    :: XSAIXY    ! stem area index
  real,    allocatable, dimension(:,:)    :: TAUSSXY   ! snow age factor

  ! OUT (with no Noah LSM equivalent) (as defined in WRF)
  real,    allocatable, dimension(:,:)    :: T2MVXY    ! 2m temperature of vegetation part
  real,    allocatable, dimension(:,:)    :: T2MBXY    ! 2m temperature of bare ground part
  real,    allocatable, dimension(:,:)    :: Q2MVXY    ! 2m mixing ratio of vegetation part
  real,    allocatable, dimension(:,:)    :: Q2MBXY    ! 2m mixing ratio of bare ground part
  real,    allocatable, dimension(:,:)    :: TRADXY    ! surface radiative temperature (k)
  real,    allocatable, dimension(:,:)    :: NEEXY     ! net ecosys exchange (g/m2/s CO2)
  real,    allocatable, dimension(:,:)    :: GPPXY     ! gross primary assimilation [g/m2/s C]
  real,    allocatable, dimension(:,:)    :: NPPXY     ! net primary productivity [g/m2/s C]
  real,    allocatable, dimension(:,:)    :: FVEGXY    ! Noah-MP vegetation fraction [-]
  real,    allocatable, dimension(:,:)    :: RUNSFXY   ! surface runoff [mm/s]
  real,    allocatable, dimension(:,:)    :: RUNSBXY   ! subsurface runoff [mm/s]
  real,    allocatable, dimension(:,:)    :: ECANXY    ! evaporation of intercepted water (mm/s)
  real,    allocatable, dimension(:,:)    :: EDIRXY    ! soil surface evaporation rate (mm/s]
  real,    allocatable, dimension(:,:)    :: ETRANXY   ! transpiration rate (mm/s)
  real,    allocatable, dimension(:,:)    :: FSAXY     ! total absorbed solar radiation (w/m2)
  real,    allocatable, dimension(:,:)    :: FIRAXY    ! total net longwave rad (w/m2) [+ to atm]
  real,    allocatable, dimension(:,:)    :: APARXY    ! photosyn active energy by canopy (w/m2)
  real,    allocatable, dimension(:,:)    :: PSNXY     ! total photosynthesis (umol co2/m2/s) [+]
  real,    allocatable, dimension(:,:)    :: SAVXY     ! solar rad absorbed by veg. (w/m2)
  real,    allocatable, dimension(:,:)    :: SAGXY     ! solar rad absorbed by ground (w/m2)
  real,    allocatable, dimension(:,:)    :: RSSUNXY   ! sunlit leaf stomatal resistance (s/m)
  real,    allocatable, dimension(:,:)    :: RSSHAXY   ! shaded leaf stomatal resistance (s/m)
  real,    allocatable, dimension(:,:)    :: BGAPXY    ! between gap fraction
  real,    allocatable, dimension(:,:)    :: WGAPXY    ! within gap fraction
  real,    allocatable, dimension(:,:)    :: TGVXY     ! under canopy ground temperature[K]
  real,    allocatable, dimension(:,:)    :: TGBXY     ! bare ground temperature [K]
  real,    allocatable, dimension(:,:)    :: CHVXY     ! sensible heat exchange coefficient vegetated
  real,    allocatable, dimension(:,:)    :: CHBXY     ! sensible heat exchange coefficient bare-ground
  real,    allocatable, dimension(:,:)    :: SHGXY     ! veg ground sen. heat [w/m2]   [+ to atm]
  real,    allocatable, dimension(:,:)    :: SHCXY     ! canopy sen. heat [w/m2]   [+ to atm]
  real,    allocatable, dimension(:,:)    :: SHBXY     ! bare sensible heat [w/m2]  [+ to atm]
  real,    allocatable, dimension(:,:)    :: EVGXY     ! veg ground evap. heat [w/m2]  [+ to atm]
  real,    allocatable, dimension(:,:)    :: EVBXY     ! bare soil evaporation [w/m2]  [+ to atm]
  real,    allocatable, dimension(:,:)    :: GHVXY     ! veg ground heat flux [w/m2]  [+ to soil]
  real,    allocatable, dimension(:,:)    :: GHBXY     ! bare ground heat flux [w/m2] [+ to soil]
  real,    allocatable, dimension(:,:)    :: IRGXY     ! veg ground net LW rad. [w/m2] [+ to atm]
  real,    allocatable, dimension(:,:)    :: IRCXY     ! canopy net LW rad. [w/m2] [+ to atm]
  real,    allocatable, dimension(:,:)    :: IRBXY     ! bare net longwave rad. [w/m2] [+ to atm]
  real,    allocatable, dimension(:,:)    :: TRXY      ! transpiration [w/m2]  [+ to atm]
  real,    allocatable, dimension(:,:)    :: EVCXY     ! canopy evaporation heat [w/m2]  [+ to atm]
  real,    allocatable, dimension(:,:)    :: CHLEAFXY  ! leaf exchange coefficient
  real,    allocatable, dimension(:,:)    :: CHUCXY    ! under canopy exchange coefficient
  real,    allocatable, dimension(:,:)    :: CHV2XY    ! veg 2m exchange coefficient
  real,    allocatable, dimension(:,:)    :: CHB2XY    ! bare 2m exchange coefficient

  ! 2D variables not used in WRF - should be removed?
  real,    allocatable, dimension(:,:)    :: TERRAIN     ! terrain height
  real,    allocatable, dimension(:,:)    :: GVFMIN      ! annual minimum in vegetation fraction
  real,    allocatable, dimension(:,:)    :: GVFMAX      ! annual maximum in vegetation fraction

  ! Needed for MMF_RUNOFF (IOPT_RUN = 5); not part of MP driver in WRF
  real,    allocatable, dimension(:,:)    :: MSFTX
  real,    allocatable, dimension(:,:)    :: MSFTY
  real,    allocatable, dimension(:,:)    :: EQZWT
  real,    allocatable, dimension(:,:)    :: RIVERBEDXY
  real,    allocatable, dimension(:,:)    :: RIVERCONDXY
  real,    allocatable, dimension(:,:)    :: PEXPXY
  real,    allocatable, dimension(:,:)    :: FDEPTHXY
  real,    allocatable, dimension(:,:)    :: AREAXY
  real,    allocatable, dimension(:,:)    :: QRFSXY
  real,    allocatable, dimension(:,:)    :: QSPRINGSXY
  real,    allocatable, dimension(:,:)    :: QRFXY
  real,    allocatable, dimension(:,:)    :: QSPRINGXY
  real,    allocatable, dimension(:,:)    :: QSLATXY
  real                                    :: WTDDT  = 30.0    ! frequency of groundwater call [minutes]
  integer                                 :: STEPWTD          ! step of groundwater call

contains

  subroutine hrldas_noahmp_vars_init()
    implicit none
    ! domain
    allocate(lat2d(xstart:xend,ystart:yend))
    allocate(lon2d(xstart:xend,ystart:yend))
    allocate(DZS(1:config%nsoil))
    allocate(DZ8W(xstart:xend,kds:kde,ystart:yend))
    lat2d = undefined_real
    lon2d = undefined_real
    DZ8W = undefined_real
    DZS = undefined_real

    allocate(COSZEN    (xstart:xend,ystart:yend))    ! cosine zenith angle
    allocate(IVGTYP    (xstart:xend,ystart:yend))    ! vegetation type
    allocate(ISLTYP    (xstart:xend,ystart:yend))    ! soil type
    allocate(VEGFRA    (xstart:xend,ystart:yend))    ! vegetation fraction []
    allocate(SHDMAX    (xstart:xend,ystart:yend))    ! annual max vegetation fraction []
    allocate(TMN       (xstart:xend,ystart:yend))    ! deep soil temperature [K]
    allocate(XLAND     (xstart:xend,ystart:yend))    ! =2 ocean; =1 land/seaice
    allocate(XICE      (xstart:xend,ystart:yend))    ! fraction of grid that is seaice
    allocate(T_PHY     (xstart:xend,kds:kde,ystart:yend))  ! 3D atmospheric temperature valid at mid-levels [K]
    allocate(QV_CURR   (xstart:xend,kds:kde,ystart:yend))  ! 3D water vapor mixing ratio [kg/kg_dry]
    allocate(U_PHY     (xstart:xend,kds:kde,ystart:yend))  ! 3D U wind component [m/s]
    allocate(V_PHY     (xstart:xend,kds:kde,ystart:yend))  ! 3D V wind component [m/s]
    allocate(SWDOWN    (xstart:xend,ystart:yend))    ! solar down at surface [W m-2]
    allocate(GLW       (xstart:xend,ystart:yend))    ! longwave down at surface [W m-2]
    allocate(P8W       (xstart:xend,kds:kde,ystart:yend))  ! 3D pressure, valid at interface [Pa]
    allocate(RAINBL    (xstart:xend,ystart:yend))    ! precipitation entering land model [mm]
    allocate(RAINBL_tmp(xstart:xend,ystart:yend))    ! precipitation entering land model [mm]

    ! INOUT (with generic LSM equivalent) (as defined in WRF)
    allocate(TSK       (xstart:xend,ystart:yend))  ! surface radiative temperature [K]
    allocate(HFX       (xstart:xend,ystart:yend))  ! sensible heat flux [W m-2]
    allocate(QFX       (xstart:xend,ystart:yend))  ! latent heat flux [kg s-1 m-2]
    allocate(LH        (xstart:xend,ystart:yend))  ! latent heat flux [W m-2]
    allocate(GRDFLX    (xstart:xend,ystart:yend))  ! ground/snow heat flux [W m-2]
    allocate(SMSTAV    (xstart:xend,ystart:yend))  ! soil moisture avail. [not used]
    allocate(SMSTOT    (xstart:xend,ystart:yend))  ! total soil water [mm][not used]
    allocate(SFCRUNOFF (xstart:xend,ystart:yend))  ! accumulated surface runoff [m]
    allocate(UDRUNOFF  (xstart:xend,ystart:yend))  ! accumulated sub-surface runoff [m]
    allocate(ALBEDO    (xstart:xend,ystart:yend))  ! total grid albedo []
    allocate(SNOWC     (xstart:xend,ystart:yend))  ! snow cover fraction []
    allocate(SMOISEQ   (xstart:xend,1:config%nsoil,ystart:yend))     ! eq volumetric soil moisture [m3/m3]
    allocate(SMOIS     (xstart:xend,1:config%nsoil,ystart:yend))     ! volumetric soil moisture [m3/m3]
    allocate(SH2O      (xstart:xend,1:config%nsoil,ystart:yend))     ! volumetric liquid soil moisture [m3/m3]
    allocate(TSLB      (xstart:xend,1:config%nsoil,ystart:yend))     ! soil temperature [K]
    allocate(SNOW      (xstart:xend,ystart:yend))  ! snow water equivalent [mm]
    allocate(SNOWH     (xstart:xend,ystart:yend))  ! physical snow depth [m]
    allocate(CANWAT    (xstart:xend,ystart:yend))  ! total canopy water + ice [mm]
    allocate(ACSNOM    (xstart:xend,ystart:yend))  ! accumulated snow melt leaving pack
    allocate(ACSNOW    (xstart:xend,ystart:yend))  ! accumulated snow on grid
    allocate(EMISS     (xstart:xend,ystart:yend))  ! surface bulk emissivity
    allocate(QSFC      (xstart:xend,ystart:yend))  ! bulk surface specific humidity

    ! INOUT (with no Noah LSM equivalent) (as defined in WRF)
    allocate(ISNOWXY   (xstart:xend,ystart:yend)) ! actual no. of snow layers
    allocate(TVXY      (xstart:xend,ystart:yend)) ! vegetation leaf temperature
    allocate(TGXY      (xstart:xend,ystart:yend)) ! bulk ground surface temperature
    allocate(CANICEXY  (xstart:xend,ystart:yend)) ! canopy-intercepted ice (mm)
    allocate(CANLIQXY  (xstart:xend,ystart:yend)) ! canopy-intercepted liquid water (mm)
    allocate(EAHXY     (xstart:xend,ystart:yend)) ! canopy air vapor pressure (pa)
    allocate(TAHXY     (xstart:xend,ystart:yend)) ! canopy air temperature (k)
    allocate(CMXY      (xstart:xend,ystart:yend)) ! bulk momentum drag coefficient
    allocate(CHXY      (xstart:xend,ystart:yend)) ! bulk sensible heat exchange coefficient
    allocate(FWETXY    (xstart:xend,ystart:yend)) ! wetted or snowed fraction of the canopy (-)
    allocate(SNEQVOXY  (xstart:xend,ystart:yend)) ! snow mass at last time step(mm h2o)
    allocate(ALBOLDXY  (xstart:xend,ystart:yend)) ! snow albedo at last time step (-)
    allocate(QSNOWXY   (xstart:xend,ystart:yend)) ! snowfall on the ground [mm/s]
    allocate(WSLAKEXY  (xstart:xend,ystart:yend)) ! lake water storage [mm]
    allocate(ZWTXY     (xstart:xend,ystart:yend)) ! water table depth [m]
    allocate(WAXY      (xstart:xend,ystart:yend)) ! water in the "aquifer" [mm]
    allocate(WTXY      (xstart:xend,ystart:yend)) ! groundwater storage [mm]
    allocate(SMCWTDXY  (xstart:xend,ystart:yend)) ! soil moisture below the bottom of the column (m3m-3)
    allocate(DEEPRECHXY(xstart:xend,ystart:yend)) ! recharge to the water table when deep (m)
    allocate(RECHXY    (xstart:xend,ystart:yend)) ! recharge to the water table (diagnostic) (m)
    allocate(TSNOXY    (xstart:xend,-NSNOW+1:0,    ystart:yend)) ! snow temperature [K]
    allocate(ZSNSOXY   (xstart:xend,-NSNOW+1:config%nsoil,ystart:yend)) ! snow layer depth [m]
    allocate(SNICEXY   (xstart:xend,-NSNOW+1:0,    ystart:yend)) ! snow layer ice [mm]
    allocate(SNLIQXY   (xstart:xend,-NSNOW+1:0,    ystart:yend)) ! snow layer liquid water [mm]
    allocate(LFMASSXY  (xstart:xend,ystart:yend)) ! leaf mass [g/m2]
    allocate(RTMASSXY  (xstart:xend,ystart:yend)) ! mass of fine roots [g/m2]
    allocate(STMASSXY  (xstart:xend,ystart:yend)) ! stem mass [g/m2]
    allocate(WOODXY    (xstart:xend,ystart:yend)) ! mass of wood (incl. woody roots) [g/m2]
    allocate(STBLCPXY  (xstart:xend,ystart:yend)) ! stable carbon in deep soil [g/m2]
    allocate(FASTCPXY  (xstart:xend,ystart:yend)) ! short-lived carbon, shallow soil [g/m2]
    allocate(LAI       (xstart:xend,ystart:yend)) ! leaf area index
    allocate(LAI_tmp   (xstart:xend,ystart:yend)) ! leaf area index
    allocate(XSAIXY    (xstart:xend,ystart:yend)) ! stem area index
    allocate(TAUSSXY   (xstart:xend,ystart:yend)) ! snow age factor

    ! OUT (with no Noah LSM equivalent) (as defined in WRF)
    allocate(T2MVXY    (xstart:xend,ystart:yend)) ! 2m temperature of vegetation part
    allocate(T2MBXY    (xstart:xend,ystart:yend)) ! 2m temperature of bare ground part
    allocate(Q2MVXY    (xstart:xend,ystart:yend)) ! 2m mixing ratio of vegetation part
    allocate(Q2MBXY    (xstart:xend,ystart:yend)) ! 2m mixing ratio of bare ground part
    allocate(TRADXY    (xstart:xend,ystart:yend)) ! surface radiative temperature (k)
    allocate(NEEXY     (xstart:xend,ystart:yend)) ! net ecosys exchange (g/m2/s CO2)
    allocate(GPPXY     (xstart:xend,ystart:yend)) ! gross primary assimilation [g/m2/s C]
    allocate(NPPXY     (xstart:xend,ystart:yend)) ! net primary productivity [g/m2/s C]
    allocate(FVEGXY    (xstart:xend,ystart:yend)) ! Noah-MP vegetation fraction [-]
    allocate(RUNSFXY   (xstart:xend,ystart:yend)) ! surface runoff [mm/s]
    allocate(RUNSBXY   (xstart:xend,ystart:yend)) ! subsurface runoff [mm/s]
    allocate(ECANXY    (xstart:xend,ystart:yend)) ! evaporation of intercepted water (mm/s)
    allocate(EDIRXY    (xstart:xend,ystart:yend)) ! soil surface evaporation rate (mm/s]
    allocate(ETRANXY   (xstart:xend,ystart:yend)) ! transpiration rate (mm/s)
    allocate(FSAXY     (xstart:xend,ystart:yend)) ! total absorbed solar radiation (w/m2)
    allocate(FIRAXY    (xstart:xend,ystart:yend)) ! total net longwave rad (w/m2) [+ to atm]
    allocate(APARXY    (xstart:xend,ystart:yend)) ! photosyn active energy by canopy (w/m2)
    allocate(PSNXY     (xstart:xend,ystart:yend)) ! total photosynthesis (umol co2/m2/s) [+]
    allocate(SAVXY     (xstart:xend,ystart:yend)) ! solar rad absorbed by veg. (w/m2)
    allocate(SAGXY     (xstart:xend,ystart:yend)) ! solar rad absorbed by ground (w/m2)
    allocate(RSSUNXY   (xstart:xend,ystart:yend)) ! sunlit leaf stomatal resistance (s/m)
    allocate(RSSHAXY   (xstart:xend,ystart:yend)) ! shaded leaf stomatal resistance (s/m)
    allocate(BGAPXY    (xstart:xend,ystart:yend)) ! between gap fraction
    allocate(WGAPXY    (xstart:xend,ystart:yend)) ! within gap fraction
    allocate(TGVXY     (xstart:xend,ystart:yend)) ! under canopy ground temperature[K]
    allocate(TGBXY     (xstart:xend,ystart:yend)) ! bare ground temperature [K]
    allocate(CHVXY     (xstart:xend,ystart:yend)) ! sensible heat exchange coefficient vegetated
    allocate(CHBXY     (xstart:xend,ystart:yend)) ! sensible heat exchange coefficient bare-ground
    allocate(SHGXY     (xstart:xend,ystart:yend)) ! veg ground sen. heat [w/m2]   [+ to atm]
    allocate(SHCXY     (xstart:xend,ystart:yend)) ! canopy sen. heat [w/m2]   [+ to atm]
    allocate(SHBXY     (xstart:xend,ystart:yend)) ! bare sensible heat [w/m2]  [+ to atm]
    allocate(EVGXY     (xstart:xend,ystart:yend)) ! veg ground evap. heat [w/m2]  [+ to atm]
    allocate(EVBXY     (xstart:xend,ystart:yend)) ! bare soil evaporation [w/m2]  [+ to atm]
    allocate(GHVXY     (xstart:xend,ystart:yend)) ! veg ground heat flux [w/m2]  [+ to soil]
    allocate(GHBXY     (xstart:xend,ystart:yend)) ! bare ground heat flux [w/m2] [+ to soil]
    allocate(IRGXY     (xstart:xend,ystart:yend)) ! veg ground net LW rad. [w/m2] [+ to atm]
    allocate(IRCXY     (xstart:xend,ystart:yend)) ! canopy net LW rad. [w/m2] [+ to atm]
    allocate(IRBXY     (xstart:xend,ystart:yend)) ! bare net longwave rad. [w/m2] [+ to atm]
    allocate(TRXY      (xstart:xend,ystart:yend)) ! transpiration [w/m2]  [+ to atm]
    allocate(EVCXY     (xstart:xend,ystart:yend)) ! canopy evaporation heat [w/m2]  [+ to atm]
    allocate(CHLEAFXY  (xstart:xend,ystart:yend)) ! leaf exchange coefficient
    allocate(CHUCXY    (xstart:xend,ystart:yend)) ! under canopy exchange coefficient
    allocate(CHV2XY    (xstart:xend,ystart:yend)) ! veg 2m exchange coefficient
    allocate(CHB2XY    (xstart:xend,ystart:yend)) ! bare 2m exchange coefficient

    ! 2D variables not used in WRF - should be removed?
    allocate(TERRAIN   (xstart:xend,ystart:yend)) ! terrain height
    allocate(GVFMIN    (xstart:xend,ystart:yend)) ! annual minimum in vegetation fraction
    allocate(GVFMAX    (xstart:xend,ystart:yend)) ! annual maximum in vegetation fraction

    ! Needed for MMF_RUNOFF (IOPT_RUN = 5); not part of MP driver in WRF
    allocate(MSFTX      (xstart:xend,ystart:yend)) !
    allocate(MSFTY      (xstart:xend,ystart:yend)) !
    allocate(EQZWT      (xstart:xend,ystart:yend)) !
    allocate(RIVERBEDXY (xstart:xend,ystart:yend)) !
    allocate(RIVERCONDXY(xstart:xend,ystart:yend)) !
    allocate(PEXPXY     (xstart:xend,ystart:yend)) !
    allocate(FDEPTHXY   (xstart:xend,ystart:yend)) !
    allocate(AREAXY     (xstart:xend,ystart:yend)) !
    allocate(QRFSXY     (xstart:xend,ystart:yend)) !
    allocate(QSPRINGSXY (xstart:xend,ystart:yend)) !
    allocate(QRFXY      (xstart:xend,ystart:yend)) !
    allocate(QSPRINGXY  (xstart:xend,ystart:yend)) !
    allocate(QSLATXY    (xstart:xend,ystart:yend)) !

    !initialize the value
    COSZEN = undefined_real
    IVGTYP = undefined_integer
    ISLTYP = undefined_integer
    VEGFRA = undefined_real
    SHDMAX = undefined_real
    TMN = undefined_real
    XLAND = undefined_real
    XICE = undefined_real
    T_PHY = undefined_real
    QV_CURR = undefined_real
    U_PHY = undefined_real
    V_PHY = undefined_real
    SWDOWN = undefined_real
    GLW = undefined_real
    P8W = undefined_real
    RAINBL = undefined_real
    RAINBL_tmp = undefined_real
    TSK = undefined_real
    QFX = undefined_real
    SMSTAV = undefined_real
    SMSTOT = undefined_real
    SMOIS = undefined_real
    SH2O = undefined_real
    TSLB = undefined_real
    SNOW = undefined_real
    SNOWH = undefined_real
    CANWAT = undefined_real
    ACSNOM = undefined_real
    ACSNOW = undefined_real
    QSFC = undefined_real
    SFCRUNOFF = undefined_real
    UDRUNOFF = undefined_real
    SMOISEQ = undefined_real
    ALBEDO = undefined_real
    ISNOWXY = undefined_integer
    TVXY = undefined_real
    TGXY = undefined_real
    CANICEXY = undefined_real
    CANLIQXY = undefined_real
    EAHXY = undefined_real
    TAHXY = undefined_real
    CMXY = undefined_real
    CHXY = undefined_real
    FWETXY = undefined_real
    SNEQVOXY = undefined_real
    ALBOLDXY = undefined_real
    QSNOWXY = undefined_real
    WSLAKEXY = undefined_real
    ZWTXY = undefined_real
    WAXY = undefined_real
    WTXY = undefined_real
    TSNOXY = undefined_real
    SNICEXY = undefined_real
    SNLIQXY = undefined_real
    LFMASSXY = undefined_real
    RTMASSXY = undefined_real
    STMASSXY = undefined_real
    WOODXY =undefined_real
    STBLCPXY = undefined_real
    FASTCPXY = undefined_real
    LAI =undefined_real
    LAI_tmp =undefined_real
    XSAIXY = undefined_real
    TAUSSXY = undefined_real

    SMCWTDXY = undefined_real
    DEEPRECHXY = undefined_real
    RECHXY = undefined_real
    ZSNSOXY = undefined_real
    GRDFLX = undefined_real
    HFX = undefined_real
    LH = undefined_real
    EMISS = undefined_real
    SNOWC = undefined_real
    T2MVXY = undefined_real
    T2MBXY =undefined_real
    Q2MVXY = undefined_real
    Q2MBXY = undefined_real
    TRADXY = undefined_real
    NEEXY = undefined_real
    GPPXY = undefined_real
    NPPXY = undefined_real
    FVEGXY = undefined_real
    RUNSFXY = undefined_real
    RUNSBXY = undefined_real
    ECANXY = undefined_real
    EDIRXY = undefined_real
    ETRANXY = undefined_real
    FSAXY = undefined_real
    FIRAXY = undefined_real
    APARXY = undefined_real
    PSNXY = undefined_real
    SAVXY = undefined_real
    FIRAXY = undefined_real
    SAGXY = undefined_real
    RSSUNXY = undefined_real
    RSSHAXY = undefined_real
    BGAPXY = undefined_real
    WGAPXY = undefined_real
    TGVXY = undefined_real
    TGBXY = undefined_real
    CHVXY = undefined_real
    CHBXY = undefined_real
    SHGXY = undefined_real
    SHCXY = undefined_real
    SHBXY = undefined_real
    EVGXY = undefined_real
    EVBXY = undefined_real
    GHVXY = undefined_real
    GHBXY = undefined_real
    IRGXY = undefined_real
    IRCXY = undefined_real
    IRBXY = undefined_real
    TRXY = undefined_real
    EVCXY = undefined_real
    CHLEAFXY = undefined_real
    CHUCXY = undefined_real
    CHV2XY = undefined_real
    CHB2XY = undefined_real
    TERRAIN = undefined_real
    GVFMIN = undefined_real
    GVFMAX = undefined_real
    MSFTX = undefined_real
    MSFTY = undefined_real
    EQZWT = undefined_real
    RIVERBEDXY = undefined_real
    RIVERCONDXY = undefined_real
    PEXPXY = undefined_real
    FDEPTHXY = undefined_real
    AREAXY = undefined_real
    QRFSXY = undefined_real
    QSPRINGSXY = undefined_real
    QRFXY = undefined_real
    QSPRINGXY = undefined_real
    QSLATXY = undefined_real

  end subroutine hrldas_noahmp_vars_init

  subroutine hrldas_noahmp_vars_probe()
    implicit none
    integer :: i
    integer :: j

    i = xstart
    j = ystart

    print *, ''
    print *, 'DEBUG - BEGIN'
    print *, 'LAT2D = ', lat2d(i,j)
    print *, 'LON2D = ', lon2d(i,j)
    print *, 'DZ8W = ', DZ8W(i,1,j)
    print *, 'DZS = ', DZS(1)
    print *, 'COSZEN = ', COSZEN(i,j)
    print *, 'IVGTYP = ', IVGTYP(i,j)
    print *, 'ISLTYP = ', ISLTYP(i,j)
    print *, 'VEGFRA = ', VEGFRA(i,j)
    print *, 'SHDMAX = ', SHDMAX(i,j)
    print *, 'TMN = ', TMN(i,j)
    print *, 'XLAND = ', XLAND(i,j)
    print *, 'XICE = ', XICE(i,j)
    print *, 'T_PHY = ', T_PHY(i,1,j)
    print *, 'QV_CURR = ', QV_CURR(i,1,j)
    print *, 'U_PHY = ', U_PHY(i,1,j)
    print *, 'V_PHY = ', V_PHY(i,1,j)
    print *, 'SWDOWN = ', SWDOWN(i,j)
    print *, 'GLW = ', GLW(i,j)
    print *, 'P8W = ', P8W(i,1,j)
    print *, 'RAINBL = ', RAINBL(i,j)
    print *, 'RAINBL_tmp = ', RAINBL_tmp(i,j)
    print *, 'TSK = ', TSK(i,j)
    print *, 'HFX = ', HFX(i,j)
    print *, 'QFX = ', QFX(i,j)
    print *, 'LH = ', LH(i,j)
    print *, 'GRDFLX = ', GRDFLX(i,j)
    print *, 'SMSTAV = ', SMSTAV(i,j)
    print *, 'SMSTOT = ', SMSTOT(i,j)
    print *, 'SFCRUNOFF = ', SFCRUNOFF(i,j)
    print *, 'UDRUNOFF = ', UDRUNOFF(i,j)
    print *, 'ALBEDO = ', ALBEDO(i,j)
    print *, 'SNOWC = ', SNOWC(i,j)
    print *, 'SMOISEQ = ', SMOISEQ(i,1,j)
    print *, 'SMOIS = ', SMOIS(i,1,j)
    print *, 'SH2O = ', SH2O(i,1,j)
    print *, 'TSLB = ', TSLB(i,1,j)
    print *, 'SNOW = ', SNOW(i,j)
    print *, 'SNOWH = ', SNOWH(i,j)
    print *, 'CANWAT = ', CANWAT(i,j)
    print *, 'ACSNOM = ', ACSNOM(i,j)
    print *, 'ACSNOW = ', ACSNOW(i,j)
    print *, 'EMISS = ', EMISS(i,j)
    print *, 'QSFC = ', QSFC(i,j)
    print *, 'ISNOWXY = ', ISNOWXY(i,j)
    print *, 'TVXY = ', TVXY(i,j)
    print *, 'TGXY = ', TGXY(i,j)
    print *, 'CANICEXY = ', CANICEXY(i,j)
    print *, 'CANLIQXY = ', CANLIQXY(i,j)
    print *, 'EAHXY = ', EAHXY(i,j)
    print *, 'TAHXY = ', TAHXY(i,j)
    print *, 'CMXY = ', CMXY(i,j)
    print *, 'CHXY = ', CHXY(i,j)
    print *, 'FWETXY = ', FWETXY(i,j)
    print *, 'SNEQVOXY = ', SNEQVOXY(i,j)
    print *, 'ALBOLDXY = ', ALBOLDXY(i,j)
    print *, 'QSNOWXY = ', QSNOWXY(i,j)
    print *, 'WSLAKEXY = ', WSLAKEXY(i,j)
    print *, 'ZWTXY = ', ZWTXY(i,j)
    print *, 'WAXY = ', WAXY(i,j)
    print *, 'WTXY = ', WTXY(i,j)
    print *, 'SMCWTDXY = ', SMCWTDXY(i,j)
    print *, 'DEEPRECHXY = ', DEEPRECHXY(i,j)
    print *, 'RECHXY = ', RECHXY(i,j)
    print *, 'TSNOXY = ', TSNOXY(i,1,j)
    print *, 'ZSNSOXY = ', ZSNSOXY(i,1,j)
    print *, 'SNICEXY = ', SNICEXY(i,1,j)
    print *, 'SNLIQXY = ', SNLIQXY(i,1,j)
    print *, 'SNLIQXY = ', LFMASSXY(i,j)
    print *, 'RTMASSXY = ', RTMASSXY(i,j)
    print *, 'STMASSXY = ', STMASSXY(i,j)
    print *, 'WOODXY = ', WOODXY(i,j)
    print *, 'STBLCPXY = ', STBLCPXY(i,j)
    print *, 'FASTCPXY = ', FASTCPXY(i,j)
    print *, 'LAI = ', LAI(i,j)
    print *, 'LAI_tmp = ', LAI_tmp(i,j)
    print *, 'XSAIXY = ', XSAIXY(i,j)
    print *, 'TAUSSXY = ', TAUSSXY(i,j)
    print *, 'T2MVXY = ', T2MVXY(i,j)
    print *, 'T2MBXY = ', T2MBXY(i,j)
    print *, 'Q2MVXY = ', Q2MVXY(i,j)
    print *, 'Q2MBXY = ', Q2MBXY(i,j)
    print *, 'TRADXY = ', TRADXY(i,j)
    print *, 'NEEXY = ', NEEXY(i,j)
    print *, 'GPPXY = ', GPPXY(i,j)
    print *, 'NPPXY = ', NPPXY(i,j)
    print *, 'FVEGXY = ', FVEGXY(i,j)
    print *, 'RUNSFXY = ', RUNSFXY(i,j)
    print *, 'RUNSBXY = ', RUNSBXY(i,j)
    print *, 'ECANXY = ', ECANXY(i,j)
    print *, 'EDIRXY = ', EDIRXY(i,j)
    print *, 'ETRANXY = ', ETRANXY(i,j)
    print *, 'FSAXY = ', FSAXY(i,j)
    print *, 'FIRAXY = ', FIRAXY(i,j)
    print *, 'APARXY = ', APARXY(i,j)
    print *, 'PSNXY = ', PSNXY(i,j)
    print *, 'SAVXY = ', SAVXY(i,j)
    print *, 'SAGXY = ', SAGXY(i,j)
    print *, 'RSSUNXY = ', RSSUNXY(i,j)
    print *, 'RSSHAXY = ', RSSHAXY(i,j)
    print *, 'BGAPXY = ', BGAPXY(i,j)
    print *, 'WGAPXY = ', WGAPXY(i,j)
    print *, 'TGVXY = ', TGVXY(i,j)
    print *, 'TGBXY = ', TGBXY(i,j)
    print *, 'CHVXY = ', CHVXY(i,j)
    print *, 'CHBXY = ', CHBXY(i,j)
    print *, 'SHGXY = ', SHGXY(i,j)
    print *, 'SHCXY = ', SHCXY(i,j)
    print *, 'SHBXY = ', SHBXY(i,j)
    print *, 'EVGXY = ', EVGXY(i,j)
    print *, 'EVBXY = ', EVBXY(i,j)
    print *, 'GHVXY = ', GHVXY(i,j)
    print *, 'GHBXY = ', GHBXY(i,j)
    print *, 'IRGXY = ', IRGXY(i,j)
    print *, 'IRCXY = ', IRCXY(i,j)
    print *, 'IRBXY = ', IRBXY(i,j)
    print *, 'TRXY = ', TRXY(i,j)
    print *, 'EVCXY = ', EVCXY(i,j)
    print *, 'CHLEAFXY = ', CHLEAFXY(i,j)
    print *, 'CHUCXY = ', CHUCXY(i,j)
    print *, 'CHV2XY = ', CHV2XY(i,j)
    print *, 'CHB2XY = ', CHB2XY(i,j)
    print *, 'TERRAIN = ', TERRAIN(i,j)
    print *, 'GVFMIN = ', GVFMIN(i,j)
    print *, 'GVFMAX = ', GVFMAX(i,j)

    print *, 'DEBUG - END'
  end subroutine hrldas_noahmp_vars_probe

end module module_hrldas_noahmp_vars
