.SUFFIXES:
.SUFFIXES: .o .F90

ifeq ($(WRF_HYDRO_RAPID),1)
include ${TAO_DIR}/bmake/tao_common
endif

include ../makefile.in

ifeq ($(WRF_HYDRO_RAPID),1)
PHDF5_INC=-I ${TACC_HDF5_INC}
RAPID_MACRO = ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_LIB} ${PHDF5_INC} \
        -Wl,-rpath,${TACC_HDF5_LIB} -L${TACC_HDF5_LIB} -lhdf5 -lz 
RAPID_LIB =  -lrapid
else
RAPID_MACRO = 
RAPID_LIB =  
endif

EXE = noahmp.exe

OBJS = \
       ../driver/main_hrldas_noahmp.o \
       ../driver/module_hrldas_noahmp_driver.o \
       ../driver/module_hrldas_noahmp_vars.o \
       ../driver/module_hrldas_noahmp_namelist.o \
       ../driver/module_hrldas_netcdf_io.o \
       ../phys/module_sf_noahmpdrv.o \
       ../phys/module_sf_noahmplsm.o \
       ../phys/module_sf_noahmp_glacier.o \
       ../phys/module_sf_noahmp_groundwater.o \
       ../phys/module_sf_myjsfc.o \
       ../phys/module_sf_sfclay.o \
       ../util/module_wrf_utilities.o \
       ../util/module_model_constants.o \
       ../util/module_date_utilities.o \
       ../util/kwm_string_utilities.o
ifeq ($(MPPFLAGS), -DMPP_LAND)
    OBJS += ../mpp/module_cpl_land.o\
	    ../mpp/module_mpp_land.o
endif

all: $(EXE)

$(EXE): $(OBJS)
	@echo ""
	echo "${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_LIB} ${PHDF5_INC} -Wl,-rpath,${TACC_HDF5_LIB} -L${TACC_HDF5_LIB} -lhdf5 -lz"
	$(F90) -o $@ $(OBJS) $(RAPID_LIB) $(NETCDFLIB) $(RAPID_MACRO)

clean:
	$(RMD) *~ $(EXE)

