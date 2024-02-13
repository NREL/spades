DEBUG = FALSE
USE_MPI  = FALSE
USE_OMP  = FALSE

COMP = llvm
BL_NO_FORT=TRUE

DIM = 2

AMREX_HOME ?= ../../../amrex

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
EBASE := vac_hop
CEXE_sources += vac_hop.cpp

INCLUDE_LOCATIONS += ./
VPATH_LOCATIONS += ./

Pdirs := Base Boundary

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
