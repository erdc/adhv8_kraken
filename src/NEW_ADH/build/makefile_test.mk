
# C Compiler
#CC = /usr/local/opt/open-mpi/bin/mpicc
petsc.pc := $(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig/petsc.pc
# Additional libraries that support pkg-config can be added to the list of PACKAGES below.
PACKAGES := $(petsc.pc)
#CC=/opt/homebrew/bin/mpicc
CC := $(shell pkg-config --variable=ccompiler $(PACKAGES))
CXX := $(shell pkg-config --variable=cxxcompiler $(PACKAGES))
FC := $(shell pkg-config --variable=fcompiler $(PACKAGES))
CFLAGS_OTHER := $(shell pkg-config --cflags-only-other $(PACKAGES))
CFLAGS := $(shell pkg-config --variable=cflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
CXXFLAGS := $(shell pkg-config --variable=cxxflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
FFLAGS := $(shell pkg-config --variable=fflags_extra $(PACKAGES))
CPPFLAGS := $(shell pkg-config --cflags-only-I $(PACKAGES))
LDFLAGS := $(shell pkg-config --libs-only-L --libs-only-other $(PACKAGES))
LDFLAGS += $(patsubst -L%, $(shell pkg-config --variable=ldflag_rpath $(PACKAGES))%, $(shell pkg-config --libs-only-L $(PACKAGES)))
LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm
CUDAC := $(shell pkg-config --variable=cudacompiler $(PACKAGES))
CUDAC_FLAGS := $(shell pkg-config --variable=cudaflags_extra $(PACKAGES))
CUDA_LIB := $(shell pkg-config --variable=cudalib $(PACKAGES))
CUDA_INCLUDE := $(shell pkg-config --variable=cudainclude $(PACKAGES))


OBJECTS_DIR = ./obj


#------------- Directories
SOURCE_DIR = \
$(wildcard $(MASTER_SRC_DIR)../debug) \
$(wildcard $(MASTER_SRC_DIR)../tools) \
$(wildcard $(MASTER_SRC_DIR)../messg) \
$(wildcard $(MASTER_SRC_DIR)../friction) \
$(wildcard $(MASTER_SRC_DIR)../xdmf) \
$(wildcard $(MASTER_SRC_DIR)../dofmaps) \
$(wildcard $(MASTER_SRC_DIR)../tokens) \
$(wildcard $(MASTER_SRC_DIR)../structs/svect) \
$(wildcard $(MASTER_SRC_DIR)../structs/snode) \
$(wildcard $(MASTER_SRC_DIR)../structs/stensor) \
$(wildcard $(MASTER_SRC_DIR)../structs/selem) \
$(wildcard $(MASTER_SRC_DIR)../structs/squad) \
$(wildcard $(MASTER_SRC_DIR)../structs/smodel) \
$(wildcard $(MASTER_SRC_DIR)../structs/slist_items) \
$(wildcard $(MASTER_SRC_DIR)../structs/smpi) \
$(wildcard $(MASTER_SRC_DIR)../structs/smeteor) \
$(wildcard $(MASTER_SRC_DIR)../structs/sflags) \
$(wildcard $(MASTER_SRC_DIR)../structs/sstr_value) \
$(wildcard $(MASTER_SRC_DIR)../structs/sgrid) \
$(wildcard $(MASTER_SRC_DIR)../structs/sarray) \
$(wildcard $(MASTER_SRC_DIR)../structs/smat) \
$(wildcard $(MASTER_SRC_DIR)../structs/sdvar) \
$(wildcard $(MASTER_SRC_DIR)../structs/ssw) \
$(wildcard $(MASTER_SRC_DIR)../structs/slin_sys) \
$(wildcard $(MASTER_SRC_DIR)../structs/sio) \
$(wildcard $(MASTER_SRC_DIR)../structs/sseries) \
$(wildcard $(MASTER_SRC_DIR)../structs/smodel_super) \
$(wildcard $(MASTER_SRC_DIR)../structs/smodel_design) \
$(wildcard $(MASTER_SRC_DIR)../residual) \
$(wildcard $(MASTER_SRC_DIR)../jacobian)\
$(wildcard $(MASTER_SRC_DIR)../la) \
$(wildcard $(MASTER_SRC_DIR)../newton) \
$(wildcard $(MASTER_SRC_DIR)../fe) \
$(wildcard $(MASTER_SRC_DIR)../timeloop) \
$(wildcard $(MASTER_SRC_DIR)../models) \
$(wildcard $(MASTER_SRC_DIR)../models/sw2) \
$(wildcard $(MASTER_SRC_DIR)../models/poisson) \
$(wildcard $(MASTER_SRC_DIR)../models/heat) \
$(wildcard $(MASTER_SRC_DIR)../test/la) \
$(wildcard $(MASTER_SRC_DIR)../test/residual) \
$(wildcard $(MASTER_SRC_DIR)../test/jacobian) \
$(wildcard $(MASTER_SRC_DIR)../test/newton) \
$(wildcard $(MASTER_SRC_DIR)../test/nonlinear_newton) \
$(wildcard $(MASTER_SRC_DIR)../test/timeloop) \
$(wildcard $(MASTER_SRC_DIR)../test/sw2_wd) \
$(wildcard $(MASTER_SRC_DIR)../main)







INCLUDE_DIR = $(MASTER_SRC_DIR)../include \
$(MASTER_SRC_DIR)../debug/include \
$(MASTER_SRC_DIR)../tools \
$(MASTER_SRC_DIR)../messg \
$(MASTER_SRC_DIR)../friction \
$(MASTER_SRC_DIR)../xdmf \
$(MASTER_SRC_DIR)../dofmaps \
$(MASTER_SRC_DIR)../tokens \
$(MASTER_SRC_DIR)../structs/svect \
$(MASTER_SRC_DIR)../structs/snode \
$(MASTER_SRC_DIR)../structs/stensor \
$(MASTER_SRC_DIR)../structs/selem \
$(MASTER_SRC_DIR)../structs/squad \
$(MASTER_SRC_DIR)../structs/smodel \
$(MASTER_SRC_DIR)../structs/slist_items \
$(MASTER_SRC_DIR)../structs/smpi \
$(MASTER_SRC_DIR)../structs/smeteor \
$(MASTER_SRC_DIR)../structs/sflags \
$(MASTER_SRC_DIR)../structs/sstr_value \
$(MASTER_SRC_DIR)../structs/sgrid \
$(MASTER_SRC_DIR)../structs/sarray \
$(MASTER_SRC_DIR)../structs/smat \
$(MASTER_SRC_DIR)../structs/sdvar \
$(MASTER_SRC_DIR)../structs/ssw \
$(MASTER_SRC_DIR)../structs/slin_sys\
$(MASTER_SRC_DIR)../structs/sio \
$(MASTER_SRC_DIR)../structs/sseries \
$(MASTER_SRC_DIR)../structs/smodel_super \
$(MASTER_SRC_DIR)../structs/smodel_design \
$(MASTER_SRC_DIR)../residual \
$(MASTER_SRC_DIR)../jacobian \
$(MASTER_SRC_DIR)../la \
$(MASTER_SRC_DIR)../newton \
$(MASTER_SRC_DIR)../fe \
$(MASTER_SRC_DIR)../timeloop \
$(MASTER_SRC_DIR)../models \
$(MASTER_SRC_DIR)../models/sw2 \
$(MASTER_SRC_DIR)../models/poisson \
$(MASTER_SRC_DIR)../models/heat \
$(MASTER_SRC_DIR)../test/la \
$(MASTER_SRC_DIR)../test/residual \
$(MASTER_SRC_DIR)../test/jacobian \
$(MASTER_SRC_DIR)../test/newton \
$(MASTER_SRC_DIR)../test/nonlinear_newton \
$(MASTER_SRC_DIR)../test/timeloop \
$(MASTER_SRC_DIR)../test/sw2_wd




#/opt/homebrew/include/suitesparse
#/opt/homebrew/include


OBJ_MK = $(addprefix $(OBJECTS_DIR)/, $(SOURCE_DIR))
#----------------------------

VPATH               = $(SOURCE_DIR)

#------------- Files
SOURCE              = $(foreach dir,    $(SOURCE_DIR),  $(wildcard  $(dir)/*.c))
FOBJECTS            = $(addprefix $(OBJECTS_DIR)/, $(SOURCE:.c=.o))
OBJECTS             = $(addprefix $(OBJECTS_DIR)/, $(notdir $(FOBJECTS)))
DEPS                = $(foreach dir,    $(INCLUDE_DIR), $(wildcard  $(dir)/*.h))
#----------------------------


#------------- Flags
OPT                 =
IFLAGS              += $(foreach dir,    $(INCLUDE_DIR), -I$(dir))
LFLAGS              += 
CFLAGS              += -O3 -D_ADH_HDF5 -Wno-incompatible-pointer-types-discards-qualifiers #-D_PETSC #-D_DEBUG #-D_MPI #-D_MESSG #-Wall
CFLAGS              += -L/opt/homebrew/lib -lumfpack -pedantic -std=c99 -I/opt/homebrew/include/suitesparse -I/opt/homebrew/include -lhdf5 -lscotch -lscotcherr -lm  #-L/opt/homebrew/lib
FLAGS               = $(OPT) $(IFLAGS) $(LFLAGS) $(CFLAGS) $(LDLIBS) 
#----------------------------


#$(info $(wildcard $(MASTER_SRC_DIR)) = "$(wildcard $(MASTER_SRC_DIR))")
#$(info SOURCE = "$(SOURCE)")
#$(info OBJECTS = "$(OBJECTS)")
#$(info DEPS = "$(DEPS)")
#$(info IFLAGS = "$(IFLAGS)")
#$(info OBJ_MK = "$(OBJ_MK)")
#$(info FOBJECTS = "$(FOBJECTS)")

BINARY = adh

latest : $(BINARY)

$(BINARY) : $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LFLAGS) $(CFLAGS) $(LDLIBS)
	#$(CC) -o $@ $^  $(OBJECTS) $(LDLIBS) $(CFLAGS)

$(OBJECTS_DIR)/%.o  : %.c $(DEPS)
	#$(LINK.C) -o $@ $^ $(LDLIBS)
	#$(CC) $(FLAGS)  -c -o $@ $< 
	$(CC) $(FLAGS)  -c -o $@ $<
	
	
.PHONY: clean
local_Cclean:
	rm -f $(OBJECTS_DIR)/*
	rm -f ./adh
	clear

#include ${PETSC_DIR}/lib/petsc/conf/variables
#include variables
