
CC = gcc
CFLAGS = -I../debug/include -g

OBJ := ./obj

MASTER_SRC_DIR = ./
MASTER_SRC_FILES = $(wildcard $(MASTER_SRC_DIR)/*.c)
DEBUG_SRC_DIR = ../debug
DEBUG_SRC_FILES = $(wildcard $(DEBUG_SRC_DIR)/*.c)

SOURCES = $(subst $(MASTER_SRC_DIR),$(SOURCE_DIR),$(MASTER_SRC_FILES)) \
$(subst $(DEBUG_SRC_DIR),$(SOURCE_DIR),$(DEBUG_SRC_FILES))
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SOURCES))

#$(info SOURCES="$(SOURCES)")
#$(info OBJECTS="$(OBJECTS)")

test: $(OBJECTS)
	$(CC) $^ -o $@ $(CFLAGS)

$(OBJ)/%.o: $(MASTER_SRC_DIR)/%.c
	$(CC) -c $< -o $@ $(CFLAGS)

$(OBJ)/%.o: $(DEBUG_SRC_DIR)/%.c
	$(CC) -c $< -o $@ $(CFLAGS)

.PHONY: clean
clean:
	rm ./test
	rm -f $(OBJ)/*
