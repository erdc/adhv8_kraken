# C Compiler
CC                  = /usr/local/opt/open-mpi/bin/mpicc

#------------- Directories
#SOURCE_DIR          = src lib
SOURCE_DIR = \
$(wildcard $(MASTER_SRC_DIR)../../../debug) \
$(wildcard $(MASTER_SRC_DIR)../../../tools) \
$(wildcard $(MASTER_SRC_DIR)../../../friction) \
$(wildcard $(MASTER_SRC_DIR)../../../xdmf) \
$(wildcard $(MASTER_SRC_DIR)../../svect) \
$(wildcard $(MASTER_SRC_DIR)../../snode) \
$(wildcard $(MASTER_SRC_DIR)../../stensor) \
$(wildcard $(MASTER_SRC_DIR)../../selem) \
$(wildcard $(MASTER_SRC_DIR)../../squad) \
$(wildcard $(MASTER_SRC_DIR)../../selem_physics) \
$(wildcard $(MASTER_SRC_DIR)../../slist_items) \
$(wildcard $(MASTER_SRC_DIR)../../smpi) \
$(wildcard $(MASTER_SRC_DIR)../../smeteor) \
$(wildcard $(MASTER_SRC_DIR)../../sflags) \
$(wildcard $(MASTER_SRC_DIR)../../sstr_value) \
$(wildcard $(MASTER_SRC_DIR)../../sgrid) \
$(wildcard $(MASTER_SRC_DIR)../../sarray) \
$(wildcard $(MASTER_SRC_DIR)../../sio) \
$(wildcard $(MASTER_SRC_DIR)../../sseries) \
$(wildcard $(MASTER_SRC_DIR)../../smodel_super) \
$(wildcard $(MASTER_SRC_DIR)../../smodel_design) \
$(wildcard $(MASTER_SRC_DIR)./)


OBJECTS_DIR         = ./obj

#INCLUDE_DIR         = . ./include
INCLUDE_DIR = $(MASTER_SRC_DIR)../../../include \
$(MASTER_SRC_DIR)../../../debug/include \
$(MASTER_SRC_DIR)../../../tools \
$(MASTER_SRC_DIR)../../../friction \
$(MASTER_SRC_DIR)../../../xdmf \
$(MASTER_SRC_DIR)../../svect \
$(MASTER_SRC_DIR)../../snode \
$(MASTER_SRC_DIR)../../stensor \
$(MASTER_SRC_DIR)../../selem \
$(MASTER_SRC_DIR)../../squad \
$(MASTER_SRC_DIR)../../selem_physics \
$(MASTER_SRC_DIR)../../slist_items \
$(MASTER_SRC_DIR)../../smpi \
$(MASTER_SRC_DIR)../../smeteor \
$(MASTER_SRC_DIR)../../sflags \
$(MASTER_SRC_DIR)../../sstr_value \
$(MASTER_SRC_DIR)../../sgrid \
$(MASTER_SRC_DIR)../../sarray \
$(MASTER_SRC_DIR)../../sio \
$(MASTER_SRC_DIR)../../sseries \
$(MASTER_SRC_DIR)../../smodel_super \
$(MASTER_SRC_DIR)../../smodel_design


OBJ_MK              = $(addprefix $(OBJECTS_DIR)/, $(SOURCE_DIR))
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
IFLAGS              = $(foreach dir,    $(INCLUDE_DIR), -I$(dir))
LFLAGS              =
CFLAGS              = -Wall #-O3
FLAGS               = $(OPT) $(IFLAGS) $(LFLAGS) $(CFLAGS)
#----------------------------


#$(info $(wildcard $(MASTER_SRC_DIR)) = "$(wildcard $(MASTER_SRC_DIR))")
#$(info SOURCE = "$(SOURCE)")
#$(info OBJECTS = "$(OBJECTS)")
#$(info DEPS = "$(DEPS)")
#$(info IFLAGS = "$(IFLAGS)")
#$(info OBJ_MK = "$(OBJ_MK)")
#$(info FOBJECTS = "$(FOBJECTS)")

BINARY = test

all : $(BINARY)

$(BINARY) : $(OBJECTS)
	$(CC) -o $@ $(OBJECTS)

$(OBJECTS_DIR)/%.o  : %.c $(DEPS)
	$(CC) $(FLAGS)  -c -o $@ $<
	
	
.PHONY: clean
clean:
	rm -f $(OBJECTS_DIR)/*
	rm -f ./test
	clear
