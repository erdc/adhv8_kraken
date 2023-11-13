# created by cjt
SHELL=/bin/sh

.SUFFIXES: .o .c .f .F

# compiler options
CC = gcc     # <------ change this for your compiler
CFLAGS = -g -Wall -Wextra -Wno-unused -Wno-unused-parameter -Wuninitialized -O3  #GNU options

# directory locations
CURRENT_DIR = .
INCD= ../../../include
STRUCTSD= ../../../structs
INITIOD= ../../../initio
DEBUGD= ../../../debug
FED= ../../../fe
TOOLSD= ../../../tools
ELEMD= ../../../elem
COLD= ../../../columns
NODED=../../../node
METEORD=../../../meteor
LIB= -lm  #GNU option

INCLUDE= $(INCD)/global_header.h 

SRC= $(CURRENT_DIR)/scale_and_rotate_3dm.c \
	$(TOOLSD)/math_simple.c \
	$(TOOLSD)/root_print.c \
	$(TOOLSD)/createChildMesh.c \
	$(TOOLSD)/math_blas_wrappers.c \
	$(TOOLSD)/tl_grid_mass.c \
	$(TOOLSD)/tc_timeunits.c \
	$(TOOLSD)/tc_eval_series.c \
	$(TOOLSD)/tc_trap_area.c \
	$(STRUCTSD)/svect.c \
	$(STRUCTSD)/svect2d.c \
	$(STRUCTSD)/selem_1d.c \
	$(STRUCTSD)/selem_2d.c \
	$(STRUCTSD)/selem_3d.c \
	$(STRUCTSD)/stensor.c \
	$(STRUCTSD)/snode.c \
	$(STRUCTSD)/sgrid.c \
	$(STRUCTSD)/sio.c \
	$(STRUCTSD)/sarray.c \
	$(STRUCTSD)/sdof.c \
	$(STRUCTSD)/smpi.c \
	$(STRUCTSD)/sseries.c \
	$(STRUCTSD)/smeteor.c \
	$(STRUCTSD)/squad.c \
	$(INITIOD)/read_geo.c \
	$(INITIOD)/get_node_id.c \
	$(INITIOD)/get_material_id.c \
	$(INITIOD)/read_int_field.c \
	$(INITIOD)/io_main.c \
	$(INITIOD)/print_header.c \
	$(INITIOD)/read_dbl_field.c \
	$(INITIOD)/read_text_field.c \
	$(INITIOD)/open_run_files.c \
	$(INITIOD)/parse_card.c \
	$(DEBUGD)/tl_debug.c \
	$(DEBUGD)/tl_debug_alloc.c \
	$(DEBUGD)/tl_debug_error.c \
	$(DEBUGD)/tl_debug_free.c \
	$(DEBUGD)/tl_debug_init.c \
	$(DEBUGD)/tl_debug_realloc.c \
	$(ELEMD)/element_djac_nrml_gradPhi.c \
	$(COLD)/tl_columns.c \
	$(FED)/fe_elem_utilities.c \
	$(FED)/fe_elem1d_integrations.c \
	$(FED)/fe_elem2d_integrations.c \
	$(FED)/fe_elem3d_integrations.c \
	$(NODED)/node_renumber_double.c \
	$(NODED)/node_renumber_vect2d.c \
	$(METEORD)/station_node_contrib.c \
	$(METEORD)/winds_proc_teeter.c \
	$(METEORD)/winds_proc_wu.c \
	$(DEBUGD)/tag.c 


LOBJ=   $(SRC:%.c=%.o)

all:	$(LOBJ)
	@echo "Building ..."
	$(CC) $(CFLAGS) -I$(INCD) -I$(STRUCTSD) -I$(DEBUGD/include) -I$(FED) -o scaleRotate3dm *.o
	rm -f *.o

 
%.o: %.c
	${CC} ${CFLAGS} -c $*.c -I$(INCD) -I$(STRUCTSD) -I$(DEBUGD)/include -I$(FED) 

clean:
	rm -f *.o ./scaleRotate3dm


