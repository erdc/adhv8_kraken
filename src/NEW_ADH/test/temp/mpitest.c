#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv) {

    //MPI_Init(NULL, NULL);
//    SMODEL_SUPER sm;
//    int i,j;
//    i=inhouse_test(argc, argv);
//    j=petsc_test(argc, argv);
//    MPI_Finalize();


    MPI_Init(NULL, NULL);

    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("The number of spawned processes are %d and this is the process %d\n",size,rank);

    MPI_Finalize();

    return 0;
}
