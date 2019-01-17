# Makefile 
.SUFFIXES:
.SUFFIXES: .o .F90

include ../makefile.in

OBJS = \
       module_sf_noahmpdrv.o \
       module_sf_noahmplsm.o \
       module_sf_noahmp_glacier.o \
       module_sf_noahmp_groundwater.o \
       module_sf_myjsfc.o \
       module_sf_sfclay.o

MODFLAGS += -D_HRLDAS_OFFLINE_ -I../util

.F90.o:
	@echo ""
	$(F90) -o $@ -c -cpp $(CPPFLAGS) $(F90FLAGS) $(MODFLAGS) $<

# Dependencies:
all: $(OBJS)

module_sf_noahmplsm.o: module_sf_myjsfc.o module_sf_sfclay.o

module_sf_noahmpdrv.o: module_sf_noahmplsm.o module_sf_noahmp_glacier.o module_sf_noahmp_groundwater.o

clean:
	$(RM) *.o *.mod *.stb *~

