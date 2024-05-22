
CC = /usr/local/Cellar/open-mpi/5.0.2/bin/mpicc
CFLAGS = -g -D_DEBUG -D_MESSG -D_MPI #-lmpi -L/usr/local/Cellar/open-mpi/5.0.2/lib/libmpi.dylib

OBJ := ./obj

MASTER_SRC_DIR = ./
MASTER_SRC_FILES = $(wildcard $(MASTER_SRC_DIR)/*.c)
DEBUG_SRC_DIR = ../debug
DEBUG_SRC_FILES = $(wildcard $(DEBUG_SRC_DIR)/*.c)
MPI_SRC_DIR = ../messg
MPI_SRC_FILES = $(wildcard $(MPI_SRC_DIR)/*.c)

INC_DIR = $(MASTER_SRC_DIR)../debug/include/ 
INC_DIR2 = $(MASTER_SRC_DIR)../include/
INC_MPI = /usr/local/Cellar/open-mpi/5.0.2/include/

SOURCES = $(subst $(MASTER_SRC_DIR),$(SOURCE_DIR),$(MASTER_SRC_FILES)) \
$(subst $(DEBUG_SRC_DIR),$(SOURCE_DIR),$(DEBUG_SRC_FILES)) \
$(subst $(MPI_SRC_DIR),$(SOURCE_DIR),$(MPI_SRC_FILES))
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SOURCES))

#$(info SOURCES="$(SOURCES)")
#$(info OBJECTS="$(OBJECTS)")

test: $(OBJECTS)
	$(CC) $^ -o $@ $(CFLAGS) -I$(INC_DIR) -I$(INC_DIR2) -I$(MASTER_SRC_DIR) -I$(INC_MPI)

$(OBJ)/%.o: $(MASTER_SRC_DIR)/%.c
	$(CC) -c $< -o $@ $(CFLAGS) -I$(INC_DIR) -I$(INC_DIR2) -I$(MASTER_SRC_DIR) -I$(INC_MPI)

$(OBJ)/%.o: $(DEBUG_SRC_DIR)/%.c
	$(CC) -c $< -o $@ $(CFLAGS) -I$(INC_DIR) 

$(OBJ)/%.o: $(MPI_SRC_DIR)/%.c
	$(CC) -c $< -o $@ $(CFLAGS) -I$(INC_DIR) -I$(INC_DIR2) -I$(MASTER_SRC_DIR) -I$(INC_MPI)

.PHONY: clean
clean:
	rm -f $(OBJ)/*
	rm -f ./test
	clear
