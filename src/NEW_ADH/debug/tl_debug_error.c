/*!
 \file tl_error.c
 \brief Write Error Message and Exit
 */

#include "header_tl_alloc.h"

#ifndef _DEBUG
void tl_error(char *error_mssg) {
    int ierr = 0;
#ifdef _MESSG
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    char location[150];
    sprintf(location, "MPI_COMM_WOLRD ID: %d ", myid);
    ierr += fputs(location, stderr);
#endif
    
    ierr += fputs(error_mssg, stderr);
    fprintf(stderr, "\n");
    fflush(stdout);
    exit(7);
}

#else
void tl_error_debug(char *error_mssg, int line, char *file) {
    int ierr = 0;
    char location[150];
#ifdef _MESSG
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    char id[150];
    sprintf(id, "MPI_COMM_WORLD ID: %d ", myid);
    ierr += fputs(id, stderr);
#endif
    sprintf(location, "ERROR: (file:line) %s:%d :: ", file, line);
    ierr += fputs(location, stderr);
    ierr += fputs("MESSAGE: ", stderr);
    ierr += fputs(error_mssg, stderr);
    fprintf(stderr, "\n");
    fflush(stdout);
    exit(7);
}
#endif
