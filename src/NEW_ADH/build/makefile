# C Compiler
CC = /usr/local/opt/open-mpi/bin/mpicc
OBJECTS_DIR = ./obj

#------------- Directories
SOURCE_DIR = \
$(wildcard $(MASTER_SRC_DIR)../debug) \
$(wildcard $(MASTER_SRC_DIR)../tools) \
$(wildcard $(MASTER_SRC_DIR)../friction) \
$(wildcard $(MASTER_SRC_DIR)../xdmf) \
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
$(wildcard $(MASTER_SRC_DIR)../structs/sio) \
$(wildcard $(MASTER_SRC_DIR)../structs/sseries) \
$(wildcard $(MASTER_SRC_DIR)../structs/smodel_super) \
$(wildcard $(MASTER_SRC_DIR).../structs/smodel_design) \
$(wildcard $(MASTER_SRC_DIR)../main)

INCLUDE_DIR = $(MASTER_SRC_DIR)../include \
$(MASTER_SRC_DIR)../debug/include \
$(MASTER_SRC_DIR)../tools \
$(MASTER_SRC_DIR)../friction \
$(MASTER_SRC_DIR)../xdmf \
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
$(MASTER_SRC_DIR)../structs/sio \
$(MASTER_SRC_DIR)../structs/sseries \
$(MASTER_SRC_DIR)../structs/smodel_super \
$(MASTER_SRC_DIR)../structs/smodel_design

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

BINARY = adh

all : $(BINARY)

$(BINARY) : $(OBJECTS)
	$(CC) -o $@ $(OBJECTS)

$(OBJECTS_DIR)/%.o  : %.c $(DEPS)
	$(CC) $(FLAGS)  -c -o $@ $<
	
	
.PHONY: clean
clean:
	rm -f $(OBJECTS_DIR)/*
	rm -f ./adh
	clear
