
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
$(wildcard $(MASTER_SRC_DIR)../friction) \
$(wildcard $(MASTER_SRC_DIR)../xdmf) \
$(wildcard $(MASTER_SRC_DIR)../dofmaps) \
$(wildcard $(MASTER_SRC_DIR)../tokens) \
$(wildcard $(MASTER_SRC_DIR)../structs/svect) \
$(wildcard $(MASTER_SRC_DIR)../structs/snode) \
$(wildcard $(MASTER_SRC_DIR)../structs/stensor) \
$(wildcard $(MASTER_SRC_DIR)../structs/selem) \
$(wildcard $(MASTER_SRC_DIR)../structs/squad) \
$(wildcard $(MASTER_SRC_DIR)../structs/selem_physics) \
$(wildcard $(MASTER_SRC_DIR)../structs/slist_items) \
$(wildcard $(MASTER_SRC_DIR)../structs/smpi) \
$(wildcard $(MASTER_SRC_DIR)../structs/smeteor) \
$(wildcard $(MASTER_SRC_DIR)../structs/sflags) \
$(wildcard $(MASTER_SRC_DIR)../structs/sstr_value) \
$(wildcard $(MASTER_SRC_DIR)../structs/sgrid) \
$(wildcard $(MASTER_SRC_DIR)../structs/sarray) \
$(wildcard $(MASTER_SRC_DIR)../structs/smat) \
$(wildcard $(MASTER_SRC_DIR)../structs/sio) \
$(wildcard $(MASTER_SRC_DIR)../structs/sseries) \
$(wildcard $(MASTER_SRC_DIR)../structs/smodel_super) \
$(wildcard $(MASTER_SRC_DIR).../structs/smodel_design) \
$(wildcard $(MASTER_SRC_DIR)../la) \
$(wildcard $(MASTER_SRC_DIR)../test/la) \
$(wildcard $(MASTER_SRC_DIR)../test/residual) \
$(wildcard $(MASTER_SRC_DIR)../main)

INCLUDE_DIR = $(MASTER_SRC_DIR)../include \
$(MASTER_SRC_DIR)../debug/include \
$(MASTER_SRC_DIR)../tools \
$(MASTER_SRC_DIR)../friction \
$(MASTER_SRC_DIR)../xdmf \
$(MASTER_SRC_DIR)../dofmaps \
$(MASTER_SRC_DIR)../tokens \
$(MASTER_SRC_DIR)../structs/svect \
$(MASTER_SRC_DIR)../structs/snode \
$(MASTER_SRC_DIR)../structs/stensor \
$(MASTER_SRC_DIR)../structs/selem \
$(MASTER_SRC_DIR)../structs/squad \
$(MASTER_SRC_DIR)../structs/selem_physics \
$(MASTER_SRC_DIR)../structs/slist_items \
$(MASTER_SRC_DIR)../structs/smpi \
$(MASTER_SRC_DIR)../structs/smeteor \
$(MASTER_SRC_DIR)../structs/sflags \
$(MASTER_SRC_DIR)../structs/sstr_value \
$(MASTER_SRC_DIR)../structs/sgrid \
$(MASTER_SRC_DIR)../structs/sarray \
$(MASTER_SRC_DIR)../structs/smat \
$(MASTER_SRC_DIR)../structs/sio \
$(MASTER_SRC_DIR)../structs/sseries \
$(MASTER_SRC_DIR)../structs/smodel_super \
$(MASTER_SRC_DIR)../structs/smodel_design \
$(MASTER_SRC_DIR)../la \
$(MASTER_SRC_DIR)../test/la \
$(MASTER_SRC_DIR)../test/residual
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
CFLAGS              += -D_PETSC -O3 -D_MPI #-Wall
CFLAGS              += -L/opt/homebrew/lib -lumfpack -pedantic -std=c99 -I/opt/homebrew/include/suitesparse -I/opt/homebrew/opt/petsc/include #-L/opt/homebrew/lib
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