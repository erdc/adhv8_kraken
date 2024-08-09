#include "header_tl_alloc.h"
#ifdef _MESSG
#include "mpi.h"
#endif
static int count = 0;

#ifndef _DEBUG
void tag() {
 count ++;
    printf("checkpoint: %d \n",count);
}
#else
void tag_debug(int linenumber,char *filename
#ifdef _MESSG
    ,MPI_Comm comm
#endif
) {
    int ierr = 0, myid = 0;
    count ++;
#ifdef _MESSG
    ierr = MPI_Comm_rank(comm, &myid);
    fflush(stdout); // fflush the screen output
    for (ierr=0;ierr<100000;ierr++); // a little pause to write everything to screen
    //ierr = MPI_Barrier(comm); // wait for other PEs (this may not work all the time)
#endif
    printf("MYID %d at checkpoint: %d (file:line) %s:%d\n",myid,count,filename, linenumber);
}
#endif
