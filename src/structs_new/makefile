
CC = gcc
CFLAGS = -g 

OBJ := ./obj

MASTER_SRC_DIR = ./
MASTER_SRC_FILES = $(wildcard $(MASTER_SRC_DIR)/*.c)
DEBUG_SRC_DIR = ../debug
DEBUG_SRC_FILES = $(wildcard $(DEBUG_SRC_DIR)/*.c)

INC_DIR = $(MASTER_SRC_DIR)../debug/include

SOURCES = $(subst $(MASTER_SRC_DIR),$(SOURCE_DIR),$(MASTER_SRC_FILES)) \
$(subst $(DEBUG_SRC_DIR),$(SOURCE_DIR),$(DEBUG_SRC_FILES))
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SOURCES))

#$(info SOURCES="$(SOURCES)")
#$(info OBJECTS="$(OBJECTS)")

test: $(OBJECTS)
	$(CC) $^ -o $@ $(CFLAGS) -I$(INC_DIR)/

$(OBJ)/%.o: $(MASTER_SRC_DIR)/%.c
	$(CC) -c $< -o $@ $(CFLAGS) -I$(INC_DIR)/

$(OBJ)/%.o: $(DEBUG_SRC_DIR)/%.c
	$(CC) -c $< -o $@ $(CFLAGS) -I$(INC_DIR)/

.PHONY: clean
clean:
	rm -f $(OBJ)/*
	rm ./test
