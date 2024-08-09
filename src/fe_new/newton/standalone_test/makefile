include variables
include rules
CFLAGS += -pedantic -std=c99

test: petsctest.o
	 -${CLINKER} -o petsctest petsctest.o ${PETSC_LIB}
	${RM} petsctest.o

advec: advec.o
	-${CLINKER} -o advec advec.o ${PETSC_LIB}
	${RM} advec.o
