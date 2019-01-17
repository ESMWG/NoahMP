# Makefile 
.SUFFIXES:
.SUFFIXES: .o .F90

include ../makefile.in

OBJS = module_cpl_land.o module_mpp_land.o 

.F90.o:
	@echo ""
	$(F90) -o $@ -c -cpp $(CPPFLAGS) $(F90FLAGS) $(NETCDFINC) $(MODFLAGS) $<

# Dependencies:
all: $(OBJS)

module_mpp_land.o: module_cpl_land.o

clean:
	$(RMD) *.o *.mod *.stb *~
