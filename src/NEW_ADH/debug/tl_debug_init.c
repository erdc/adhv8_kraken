#define DEBUG_INIT
/*!
 \file tl_debug_init.c
 \brief Initialize Debug Information
 */
#include "header_tl_alloc.h"
#include "debug.h"

#ifndef _DEBUG
int allocated_memory = 0;
int  max_allocated_memory = 0;
#endif

#ifdef _MESSG
// cjt :: since COMM_WORLD cannot be used for CSTORM :: only good for 1 grid in CSTORM
extern MPI_Comm debug_comm;

//*****************************************//
//*****************************************//
void debug_initialize(MPI_Comm tcomm)
#else
void debug_initialize()
#endif
{
    
#ifdef _MESSG
    MPI_Comm_dup(tcomm, &debug_comm);
#endif
    
#ifdef _DEBUG
    int i = 0;
    picket_value = PICKET;
    picket_pointer = (char *)&picket_value;
    
    nan_value = MY_NAN;
    nan_pointer = (char *)&nan_value;
    
    /* allocate storage for linked lists for memory tables */
    mem_hash = (DEBUG_MEM_ENTRY **) malloc(DEBUG_HASHSIZE * sizeof(DEBUG_MEM_ENTRY *));
    for(i = 0; i < DEBUG_HASHSIZE; i++){
        mem_hash[i] = NULL;
    }
#endif
    allocated_memory = 0;
    max_allocated_memory = 0;
    return;
}

//*****************************************//
//*****************************************//

void debug_finalize()
{
#ifdef _MESSG
    int myid, ierr_code;
    ierr_code = MPI_Comm_rank(debug_comm, &myid); // only works for 1 grid in CSTORM :: MPI_COMM_WORLD,&myid);
#else
    int myid=0;
#endif
#ifdef _DEBUG
    int ierr = 0;
    
    /* Check for Unfreed and then free mem_hash */
    tl_check_unfreed();
    free(mem_hash);
#endif
    if(myid >= 0) {
        printf("\n");
        printf("**********************************************************************\n");
        printf("Memory/Debug Info:\n");
        printf("MYID = %4d, Exiting, (unfreed) Allocated Memory of %zu Bytes\n", myid, allocated_memory);
        printf("MYID = %4d, Maximum Allocated Memory %zu B = %zu kB = %zu MB = %zu GB\n",
               myid, max_allocated_memory, max_allocated_memory / 1024,
               max_allocated_memory / 1048576, max_allocated_memory / 1073741824);
        printf("**********************************************************************\n");
    }
    return;
}
