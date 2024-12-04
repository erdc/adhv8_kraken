/*!
 \file tl_alloc.c
 \brief Allocate Memory
 */

/*!
 \brief Allocate memory
 
 The purpose of this routine is to allocate memory.  An attempt
 was made to do some elementary memory checking.  Normally, a pointer
 X would point to these bytes:
 X = VVVVVVVVVV
 If compiled with "-D_DEBUG" then X points to these bytes:
 X = LL_DDDD_VVVVVVVVVV_DDDD
 Then on either side of useful memory 'V' are "pickets", called here
 'DDDD'.  In the pickets are placed NaN values, so that if they are
 touched, we will know.   Also, we will place 'LL' at the very
 beginning.  In this memory, we will store the exact size originally
 allocated, that is, the size of the 'V' storage.  Also, we will
 initialize all useful memory to NaN values, so that if it is used
 before it is initialized by the calling program, we will know.
 
 Possible issues, which may/may not be important:
 * It is possible that the picket NaN values are used but meaningless
 because of the particular error. The main point of this exercise
 was to put something in there with a distinct bit pattern and if
 we get lucky that it was grabbed and used in a certain way, all
 the better.
 * The same thing goes for initializing all useful memory to NaN.
 
 Note: mypntr/pntr below are "char".  The thinking is that chars
 are 1-byte.  So all pointer arithmetic will be valid and
 we do not need to worry about adding/subtracting a specific
 number of bytes from it.
 
 \param obj_size sizeof(double), sizeof(int), etc.
 \param number Number of Units needed
 */

// macros
#define MAX(a,b) (a > b ? a : b)


#ifdef _MESSG
//#include "mpi.h"
#else
static int myid = 0;
#endif

/************************************************************/
/************************************************************/
#ifndef _DEBUG

#include "header_tl_alloc.h"
#include "debug.h"
void *tl_alloc(size_t obj_size, /* the size of the object to be allocated */
               int number       /* the number of objects to be allocated */
)
{
    size_t isize = 0;             /* the size of the space to be allocated */
    char *pntr = NULL;            /* pointer for memory allocation purposes */
    char *mypntr = NULL;          /* pointer for memory allocation purposes */
#ifdef _MESSG
    int myid, ierr_code;
    ierr_code = MPI_Comm_rank(debug_comm, &myid); // only works for 1 grid in CSTORM :: MPI_COMM_WORLD, &myid);
#else
    int myid = 0;
#endif
    
    /* trap bad memory request */
    if ((obj_size == 0) || (number <= 0)) {
        //   printf("MYID %4d, Bad Memory Request, Run with _DEBUG\n", myid);
        exit(EXIT_FAILURE);
    }
    
    /* compute the size of the desired space */
    isize = obj_size * ((size_t) (number));
    
    /* allocate the array, point to beginning of useful memory */
    mypntr = (char *) malloc(isize);
    if (mypntr == NULL) {
        printf("MYID %4d, Memory Request Failing (tl_alloc), Current Allocated_Memory = %d MB\n", myid, allocated_memory / 1048576);
        printf("MYID %4d, Asking for = %lu MB\n", myid, isize / 1048576);
        printf("MYID %4d, Run with USE_PACKAGE_DEBUG turned ON, to find FILE:LINE of request.\n", myid);
        exit(EXIT_FAILURE);
    }
    pntr = mypntr;
    
    allocated_memory += isize;
    max_allocated_memory = MAX(max_allocated_memory, allocated_memory);
    return (pntr);
}

/************************************************************/
/************************************************************/
#else
#include "header_tl_alloc.h"
#include "debug.h"

/*!
 \brief Allocate memory (Debug Version)
 
 \param obj_size sizeof(double), sizeof(int), etc.
 \param number Number of Units needed
 \param linenumber Should be macro __LINE__
 \param *filename Should be macro __FILE__
 */
void *tl_alloc_debug(size_t obj_size,   /* the size of the object to be allocated */
                     int number,    /* the number of objects to be allocated */
                     int linenumber,    /* line number */
                     char *filename /* Calling routine filename */
)
{   
    size_t isize = 0;             /* the size of the space to be allocated */
    size_t picket_size = 0;       /* the size of the picket for memory debugging */
    size_t debug_size = 0;        /* the size of the array for memory debugging */
    char *pntr = NULL;            /* pointer for memory allocation purposes */
    char *mypntr = NULL;          /* pointer for memory allocation purposes */
#ifdef _MESSG
    int myid, ierr_code;
    
    ierr_code = MPI_Comm_rank(debug_comm, &myid); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD, &myid);
#endif
    /* trap bad memory request */
    if ((obj_size == 0) || (number <= 0)) {
        printf("MYID %4d, Bad Memory Request, (file:line) %s:%d\n", myid, filename, linenumber);
        exit(EXIT_FAILURE);
    }
    
    /* compute the size of the desired space */
    isize = obj_size * ((size_t) (number));
    
    /* set the picket size */
    debug_size = sizeof(unsigned long);
    picket_size = sizeof(int) * ((size_t) (BWIDTH));
    
    /* allocate the array, point to beginning of useful memory */
    //printf("picket_size: %lu\n",picket_size);
    //printf("debug_size: %lu\n", debug_size);
    //printf("isize: %lu\n", isize);
    //printf("size: %lu\n",isize + 2 * picket_size + debug_size);
    mypntr = (char *) malloc(isize + 2 * picket_size + debug_size);
    if (mypntr == NULL) {
        printf("MYID %4d, Error (file:line) %s:%d\n", myid, filename, linenumber);
        exit(EXIT_FAILURE);
    }
    pntr = mypntr + picket_size + debug_size;
    //tag();
    //printf("In tl_alloc_debug\n");
    /* Set Pickets */
    //printf("mypntr: %s\n",mypntr);
    
    tl_set_picket(mypntr + debug_size, isize);
    //tag();
    /* Store Size at very beginning of memory */
    *((unsigned long *) (mypntr)) = (unsigned long) (isize);
    
    /* Initialize useful memory to NaN */
    tl_set_nan(pntr, isize);
    
    /* add entry to hash table */
    tl_add_hash_entry(pntr, isize, linenumber, filename);
    
    allocated_memory += isize;
    max_allocated_memory = MAX(max_allocated_memory, allocated_memory);
    return (pntr);
}
#endif
/************************************************************/
/************************************************************/
