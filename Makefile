include makefile.in

ifeq ($(MPPFLAGS), -DMPP_LAND)
    SUBDIRS = mpp
endif
SUBDIRS += \
	   util \
	   phys \
	   driver \
	   run

TOPTARGETS := all clean

$(TOPTARGETS): $(SUBDIRS)

$(SUBDIRS):
	make -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS)
