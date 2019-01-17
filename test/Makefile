.SUFFIXES:
.SUFFIXES: .o .F90 .exe

include ../makefile.in

MODFLAGS += -I../mpp -DMPP_LAND

EXES = test_hrldas_netcdf_io.exe \
       test_mpp_land_partition.exe \
       test_mpp_land_bcast.exe

.F90.o:
	@echo ""
	$(F90) -o $@ -Wall -c -cpp $(CPPFLAGS) $(F90FLAGS) $(NETCDFINC) $(MODFLAGS) $<

.o.exe:
	$(F90) -o $@ $^ $(NETCDFLIB)

# Dependencies:
all: $(EXES)

test_hrldas_netcdf_io.exe: test_hrldas_netcdf_io.o \
                           ../mpp/module_mpp_land.o \
                           ../mpp/module_cpl_land.o

test_mpp_land_partition.exe: test_mpp_land_partition.o \
                             ../mpp/module_mpp_land.o \
                             ../mpp/module_cpl_land.o

test_mpp_land_bcast.exe: test_mpp_land_bcast.o \
                         ../mpp/module_mpp_land.o \
                         ../mpp/module_cpl_land.o

clean:
	$(RMD) $(EXES) *.o *.mod *.stb *~ a.out
