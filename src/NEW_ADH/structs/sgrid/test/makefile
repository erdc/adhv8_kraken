

#works for MPI with H5 enabled
#I also have to export LIBRARY_PATH=/usr/local/opt/hdf5_parallel/lib/
CC = /usr/local/opt/open-mpi/bin/mpicc
CFLAGS = -g -D_DEBUG -D_MESSG -D_MPI -D_ADH_HDF5 -Wno-implicit-function-declaration -lhdf5 -L/usr/local/opt/hdf5_parallel/lib/
INC_MPI = /usr/local/opt/open-mpi/include/
#Parallel hdf5
INC_HDF5= /usr/local/opt/hdf5_parallel/include/ 


#comment these lines out for parallel compilation
####
#serial and hdf5
#CC = /usr/local/bin/gcc-13
#CFLAGS = -g -D_DEBUG -D_ADH_HDF5 -Wno-implicit-function-declaration -lhdf5
#serial hdf5
#INC_HDF5=/usr/local/include/
#something garbage
#INC_MPI = /usr/local/icnlude/
####


OBJ := ./obj

MASTER_SRC_DIR = ./
MASTER_SRC_FILES = $(wildcard $(MASTER_SRC_DIR)/*.c)
DEBUG_SRC_DIR = ../debug
DEBUG_SRC_FILES = $(wildcard $(DEBUG_SRC_DIR)/*.c)
XDMF_SRC_DIR = ../xdmf
XDMF_SRC_FILES = $(wildcard $(XDMF_SRC_DIR)/*.c)
MPI_SRC_DIR = ../messg
MPI_SRC_FILES = $(wildcard $(MPI_SRC_DIR)/*.c)

INC_DIR = $(MASTER_SRC_DIR)../debug/include/ 
INC_DIR2 = $(MASTER_SRC_DIR)../include/

SOURCES = $(subst $(MASTER_SRC_DIR),$(SOURCE_DIR),$(MASTER_SRC_FILES)) \
$(subst $(DEBUG_SRC_DIR),$(SOURCE_DIR),$(DEBUG_SRC_FILES)) \
$(subst $(MPI_SRC_DIR),$(SOURCE_DIR),$(MPI_SRC_FILES)) \
$(subst $(XDMF_SRC_DIR),$(SOURCE_DIR),$(XDMF_SRC_FILES))
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SOURCES))

#$(info SOURCES="$(SOURCES)")
#$(info OBJECTS="$(OBJECTS)")

test: $(OBJECTS)
	$(CC) $^ -o $@ $(CFLAGS) -I$(INC_DIR) -I$(INC_DIR2) -I$(MASTER_SRC_DIR) -I$(INC_MPI) -I$(INC_HDF5)

$(OBJ)/%.o: $(MASTER_SRC_DIR)/%.c
	$(CC) -c $< -o $@ $(CFLAGS) -I$(INC_DIR) -I$(INC_DIR2) -I$(MASTER_SRC_DIR) -I$(INC_MPI) -I$(INC_HDF5)

$(OBJ)/%.o: $(DEBUG_SRC_DIR)/%.c
	$(CC) -c $< -o $@ $(CFLAGS) -I$(INC_DIR) 

$(OBJ)/%.o: $(MPI_SRC_DIR)/%.c
	$(CC) -c $< -o $@ $(CFLAGS) -I$(INC_DIR) -I$(INC_DIR2) -I$(MASTER_SRC_DIR) -I$(INC_MPI) -I$(INC_HDF5)
	
$(OBJ)/%.o: $(XDMF_SRC_DIR)/%.c
	$(CC) -c $< -o $@ $(CFLAGS) -I$(INC_DIR) -I$(INC_DIR2) -I$(MASTER_SRC_DIR) -I$(INC_MPI) -I$(INC_HDF5)

.PHONY: clean
clean:
	rm -f $(OBJ)/*
	rm -f ./test
	clear
