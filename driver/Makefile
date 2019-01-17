# Makefile 
.SUFFIXES:
.SUFFIXES: .o .F90

include ../makefile.in

OBJS = \
       main_hrldas_noahmp.o \
       module_hrldas_noahmp_driver.o \
       module_hrldas_noahmp_vars.o \
       module_hrldas_noahmp_namelist.o \
       module_hrldas_netcdf_io.o

MODFLAGS += -D_HRLDAS_OFFLINE_ -I../util -I../phys
ifeq ($(MPPFLAGS), -DMPP_LAND)
    MODFLAGS += -I../mpp $(MPPFLAGS)
endif

.F90.o:
	@echo ""
	$(F90) -o $@ -c -cpp $(CPPFLAGS) $(F90FLAGS) $(NETCDFINC) $(MODFLAGS) $<

# Dependency
all: $(OBJS)

main_hrldas_noahmp.o: module_hrldas_noahmp_driver.o

module_hrldas_noahmp_driver.o: module_hrldas_noahmp_vars.o \
                               module_hrldas_noahmp_namelist.o \
                               module_hrldas_netcdf_io.o

module_hrldas_noahmp_vars.o: module_hrldas_noahmp_namelist.o \
                             module_hrldas_netcdf_io.o

clean:
	$(RMD) *.o *.mod *.stb *~
