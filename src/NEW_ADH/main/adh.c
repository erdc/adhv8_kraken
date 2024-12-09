#include "adh.h"
int main(int argc, char **argv) {

    SMODEL_SUPER sm;

#ifdef _MPI
    printf("HERE\n");
    MPI_Init(NULL, NULL);
    printf("MPI INITIALIZED");
    int i,j;
    i=inhouse_test(argc, argv);
    j=petsc_test(argc, argv);
    MPI_Finalize();
#endif
#ifdef _DEBUG
        debug_initialize();
#endif


    int ierr;
    //residual_test(argc,argv);
    //try an assembly
    //jacobian_test(argc,argv);
    //try a newton solve
    ierr = newton_test(argc,argv);
    assert(ierr==0);

    
    //not returning 0 will result in MPI error
    return 0;
}
