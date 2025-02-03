#include "adh.h"
int main(int argc, char **argv) {

    //SMODEL_SUPER sm;

    //set resid routines and inc routines
    set_models(fe_resid, fe_init);

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
#ifdef _PETSC
    printf("Calling PETSC Initialize\n");
    PetscCall(PetscInitialize(&argc,&argv,NULL,
    "Compute e in parallel with PETSc.\n\n"));
     printf("Called PETSC Initialize\n");
#endif

    int ierr = 0;
    residual_test(argc,argv);
    printf(">>>>>>>>>>>>>>>Residual test complete<<<<<<<<<<<<<<<<<<<<<\n");
    //try an assembly
    jacobian_test(argc,argv);
    printf(">>>>>>>>>>>>>>>Jacobian test complete<<<<<<<<<<<<<<<<<<<<<\n");
    //try a newton solve
    ierr = newton_test(argc,argv);
    assert(ierr==0);
    printf(">>>>>>>>>>>>>>>Newton test complete<<<<<<<<<<<<<<<<<<<<<<<<\n");
    ////try a nonlinear newton solve
    ierr = nonlinear_newton_test(argc,argv);
    assert(ierr==0);
    printf(">>>>>>>>>>>>>>>Nonlinear Newton test complete<<<<<<<<<<<<<<<<<<<<<<<<\n");
    //try a time step
    ierr = timeloop_test(argc,argv);
    assert(ierr==0);
   printf(">>>>>>>>>>>>>>>TIMELOOP test complete<<<<<<<<<<<<<<<<<<<<<<<<\n");
    //try a time step
    ierr = sw2_wd_test(argc, argv);
    assert(ierr==0);
    printf(">>>>>>>>>>>>>>>SW2 WD test complete<<<<<<<<<<<<<<<<<<<<<<<<\n");
    

    free_bcgstab();
    //not returning 0 will result in MPI error
    return 0;
}
