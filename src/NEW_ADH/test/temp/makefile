petsc.pc := $(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig/PETSc.pc

# Additional libraries that support pkg-config can be added to the list of PACKAGES below.
PACKAGES := $(petsc.pc)

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
#CUDAC := $(shell pkg-config --variable=cudacompiler $(PACKAGES))
#CUDAC_FLAGS := $(shell pkg-config --variable=cudaflags_extra $(PACKAGES))
#CUDA_LIB := $(shell pkg-config --variable=cudalib $(PACKAGES))
#CUDA_INCLUDE := $(shell pkg-config --variable=cudainclude $(PACKAGES))

print:
	@echo CC=$(CC)
	@echo CXX=$(CXX)
	@echo FC=$(FC)
	@echo CFLAGS=$(CFLAGS)
	@echo CXXFLAGS=$(CXXFLAGS)
	@echo FFLAGS=$(FFLAGS)
	@echo CPPFLAGS=$(CPPFLAGS)
	@echo LDFLAGS=$(LDFLAGS)
	@echo LDLIBS=$(LDLIBS)
	@echo CUDAC=$(CUDAC)
	@echo CUDAC_FLAGS=$(CUDAC_FLAGS)
	@echo CUDA_LIB=$(CUDA_LIB)
	@echo CUDA_INCLUDE=$(CUDA_INCLUDE)
	@echo LINK.C=$(LINK.C)
	@echo INCLUDEDIR=$(INCLUDE_DIR)



petsctest : petsctest.o
	$(CC) -o $@ $^ $(LDLIBS) $(CFLAGS)

CFLAGS += -lumfpack -pedantic -std=c99 -I/opt/homebrew/include/suitesparse -L/opt/homebrew/lib
inhouse: inhouse.o 
	$(CC) -o inhouse inhouse.o $(CFLAGS)



clean:
	rm -f *.o
	rm -f ./petsctest
	rm -f ./inhouse
	clear
