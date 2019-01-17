# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .F

include ../makefile.in

OBJS = \
	module_date_utilities.o \
	module_model_constants.o \
	module_wrf_utilities.o \
	kwm_string_utilities.o
	
MODFLAGS += -D_HRLDAS_OFFLINE_

.F.o:
	@echo ""
	$(F90) -o $@ -c -cpp $(CPPFLAGS) $(F90FLAGS) $(NETCDFINC) $(MODFLAGS) $<

# Dependencies:
all:	$(OBJS)

clean:
	$(RM) *.o *.mod *.stb *~

