################################################################################
# User input, if needed.

SHELL := /bin/bash

EXE = pt
BUILDDIR = build

LIBINCLUDES = 
LIBS = -lm 
CC = /usr/local/Cellar/open-mpi/4.0.3/bin/mpicc
CFLAGS = -D_XDMF -D_MPI
CDEBUGFLAGS = -g -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter


################################################################################
# Do not edit below this.

TARGET = $(BUILDDIR)/$(EXE)
HEADERS = $(wildcard *.h)

OBJDIR = $(BUILDDIR)/obj
OBJS = $(patsubst %.c, %.o, $(wildcard *.c))
OBJECTS = $(addprefix $(OBJDIR)/, $(OBJS))

TAGFILE = $(OBJDIR)/.build

#ifeq ($(XDMF),ON)
#	CFLAGS += -D_ADH_HDF5  -D_ADH_DEFINE_COPY
#	# corey-tag : Need to set paths to HDF5, SZIP, ZLIB, and XML2 libraries that AdH generally uses.
#	LIBINCLUDES += -I/usr/local/Cellar/hdf5/1.12.0/include
#   	LIBINCLUDES += -I/usr/local/include
#   	LIBINCLUDES += -I/opt/local/include
#   	LIBINCLUDES += -I/opt/local/include/libxml2/libxml
#	LIBINCLUDES += -I/usr/local/Cellar/open-mpi/4.0.3/include
#   	LIBS += -L/usr/local/Cellar/hdf5/1.12.0/lib -lhdf5
#   	LIBS += -L/usr/local/lib -lsz
#   	LIBS += -L/opt/local/lib -lz
#   	LIBS += -L/usr/lib -lxml2
#	LIBS += -ldl
#endif

CFLAGS += -D_ADH_HDF5  -D_ADH_DEFINE_COPY
LIBINCLUDES += -I/usr/local/Cellar/hdf5/1.12.0/include
LIBINCLUDES += -I/usr/local/include
LIBINCLUDES += -I/opt/local/include
LIBINCLUDES += -I/opt/local/include/libxml2/
#LIBINCLUDES += -I/opt/local/include/libxml2/libxml/
LIBINCLUDES += -I/usr/local/Cellar/open-mpi/4.0.3/include
LIBS += -L/usr/local/Cellar/hdf5/1.12.0/lib -lhdf5
LIBS += -L/usr/local/lib -lsz
LIBS += -L/opt/local/lib -lz
LIBS += -L/usr/lib -lxml2
#LIBS += -L/opt/local/lib/libxml2.dylib
LIBS += -ldl

################################################################################
# Makefile for Oden Institute computers (only tested on hostname "opal")

ifeq (${LMHOST},lic.ices.utexas.edu)
	CC = /opt/apps/ossw/libraries/mpich2/mpich2-3.2/ubt18/gcc-8.2/bin/mpicc
	LIBINCLUDES += -I/usr/include/hdf5/serial
	LIBINCLUDES += -I/usr/include/libxml2
	LIBINCLUDES += -I/opt/apps/ossw/libraries/mpich2/mpich2-3.2/ubt18/gcc-8.2/include
	LIBS += -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5
	LIBS += -lxml2 -lsz -lz -ldl
endif

################################################################################
.PHONY: all clean debug

all: BUILDTXT += release
all: prep $(TARGET)

debug: CFLAGS += $(CDEBUGFLAGS)
debug: BUILDTXT += debug
debug: prep $(TARGET)

prep:
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(OBJDIR)
	$(shell [[ `cat $(TAGFILE) 2>/dev/null` == $(BUILDTXT) ]] || echo $(BUILDTXT) > $(TAGFILE))

$(OBJDIR)/%.o: %.c $(HEADERS) $(TAGFILE)
	$(CC) $(CFLAGS) $(LIBINCLUDES) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) $(CFLAGS) $(LIBS) -o $@

clean:
	-rm -f $(OBJECTS)
	-rm -f $(TARGET)
	-rm -f $(TAGFILE)
