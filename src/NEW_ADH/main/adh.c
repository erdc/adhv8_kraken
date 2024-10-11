#include "adh.h"
int main(int argc, char **argv) {
    SMODEL_SUPER sm;
#ifdef _MPI
    MPI_Init(NULL, NULL);
    
    int i,j;
    i=inhouse_test(argc, argv);
    j=petsc_test(argc, argv);
    MPI_Finalize();
#endif

    //not returning 0 will result in error
    return 0;
}
