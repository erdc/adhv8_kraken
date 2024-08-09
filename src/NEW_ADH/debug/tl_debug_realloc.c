/*!
 \file tl_realloc.c
 \brief Re-Allocate Memory
 */

// macros
#define MAX(a,b) (a > b ? a : b)

#include "header_tl_alloc.h"
#include "debug.h"

/*!
 \brief Re-Allocate memory
 
 Re-allocate memory.  In case of memory debugging, check pickets
 and size first.  See tl_alloc.c for description.
 
 Note: mypntr below is a "char".  The thinking is that chars
 are 1-byte.  So all pointer arithmetic will be valid and
 we do not need to worry about adding/subtracting a specific
 number of bytes from it.
 
 \param obj_size sizeof(double), sizeof(int), etc.
 \param number Number of Units needed
 \param old_number Previous number of Units needed
 \param pntr Pointer to memory we will realloc
 */

/************************************************************/
/************************************************************/
#ifndef _DEBUG
void *tl_realloc(
                 size_t obj_size,		/* the size of the object to be allocated */
                 int number,			/* the number of objects to be allocated */
                 int old_number,		/* the number of objects previously allocated */
                 void *pntr			/* the pointer */
)
{
    size_t isize = 0;		/* the size of the space to be allocated */
    size_t old_size = 0;		/* the previous amount of space allocated */
#ifdef _MESSG
    int myid, ierr_code;
    ierr_code = MPI_Comm_rank(debug_comm, &myid); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD, &myid);
#else
    int myid = 0;
#endif
    
    /* compute the size of the desired space */
    isize = obj_size * ((size_t) (number));
    old_size = obj_size * ((size_t) (old_number));
    
    /* trap bad memory request */
    if((obj_size == 0) || (number <= 0) || (old_number < 0)) {
        printf("MYID %4d, Bad Memory Request, Run with _DEBUG\n", myid);
        exit(EXIT_FAILURE);
    }
    
    /* allocate the arrays */
    if(old_size == 0) {
        pntr = malloc(isize);
    }
    else {
        pntr = realloc(pntr, isize);
    }
    
    /* Check for failure */
    if(pntr == NULL) {
        //printf
        //("MYID %4d, Memory Request Failing (tl_realloc), Current Allocated_Memory = %lu GB \n",
        // myid, allocated_memory / 1048576);
        printf("MYID %4d, Asking for = %lu MB\n", myid, isize / 1048576);
        printf
        ("MYID %4d, Run with USE_PACKAGE_DEBUG turned ON, to find FILE:LINE of request.\n",
         myid);
        exit(EXIT_FAILURE);
    }
    /* Update allocated_memory count */
    allocated_memory -= old_size;
    allocated_memory += isize;
    max_allocated_memory = MAX(max_allocated_memory, allocated_memory);
    return (pntr);
}

/************************************************************/
/************************************************************/
#else
/*!
 \brief Re-Allocate memory (Debug Version)
 
 \param obj_size sizeof(double), sizeof(int), etc.
 \param number Number of Units needed
 \param old_number Previous number of Units needed
 \param linenumber Should be macro __LINE__
 \param filename Should be macro __FILE__
 \param pntr Pointer to memory we will realloc
 */
void *tl_realloc_debug(
                       size_t obj_size,		/* the size of the object to be allocated */
                       int number,			/* the number of objects to be allocated */
                       int old_number,		/* the number of objects previously allocated */
                       int linenumber,		/* line number */
                       char *filename,		/* calling routine filename */
                       void *pntr			/* the pointer */
)
{
    size_t isize = 0;		/* the size of the space to be allocated */
    size_t old_size = 0;		/* the previous amount of space allocated */
    size_t picket_size = 0;	/* picket size for memory debugging */
    size_t debug_size = 0;	/* the size of the array for memory debugging */
    char *mypntr = NULL;		/* pointer for memory allocation purposes */
    
#ifdef _MESSG
    int myid, ierr_code;
    ierr_code = MPI_Comm_rank(debug_comm, &myid); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD, &myid);
#else
    int myid = 0;
#endif
    
    /* compute the size of the desired space */
    isize = obj_size * ((size_t) (number));
    old_size = obj_size * ((size_t) (old_number));
    
    /* trap bad memory request */
    if((obj_size == 0) || (number <= 0) || (old_number < 0)) {
        printf("MYID %4d, Bad Memory Request, (file:line) %s:%d\n", myid, filename,
               linenumber);
        exit(EXIT_FAILURE);
    }
    
    picket_size = sizeof(int) * ((size_t) (BWIDTH));
    debug_size = sizeof(unsigned long);
    
    if(old_size == 0) {
        /* allocate enough memory (useful + debug info) */
        mypntr = (char *)malloc(isize + 2 * picket_size + debug_size);
        if(mypntr == NULL)
        {
            printf("MYID %4d, Error (file:line) %s:%d\n", myid, filename, linenumber);
            printf("Malloc Failed\n");
            exit(EXIT_FAILURE);
        }
        
        /* Set Picket on either side of useful memory */
        tl_set_picket(mypntr + debug_size, isize);
        
        /* Store Size at very beginning of memory */
        *((unsigned long *)(mypntr)) = (unsigned long)(isize);
        
        /* Point to beginning of useful memory */
        pntr = mypntr + picket_size + debug_size;
        
        /* Initialize useful memory to NaN */
        mypntr = (char *)(pntr);
        tl_set_nan(mypntr, isize);
        
        /* add entry to hash table */
        tl_add_hash_entry(pntr, isize, linenumber, filename);
        
    } else {
        
        /* Check allocated vs. free size */
        mypntr = (char *)(pntr);
        mypntr = mypntr - picket_size - debug_size;
        if((unsigned long)(old_size) != *((unsigned long *)(mypntr))) {
            printf("MYID %4d, Error (file:line) %s:%d\n", myid, filename, linenumber);
            printf("Freeing-> %lu , %lu <-Originally Bytes\n", old_size,
                   *((unsigned long *)(mypntr)));
            exit(EXIT_FAILURE);
        }
        
        /* Check pickets */
        mypntr = (char *)(pntr);
        tl_check_picket(mypntr, old_size, linenumber, filename);
        
        /* remove old entry from hash table */
        tl_remove_hash_entry(pntr);
        
        /* reallocate, point to beginning of useful memory */
        mypntr = (char *)(pntr);
        mypntr = mypntr - picket_size - debug_size;
        mypntr = (char *)realloc(mypntr, isize + 2 * picket_size + debug_size);
        if(mypntr == NULL) {
            printf("MYID %4d, Error (file:line) %s:%d\n", myid, filename, linenumber);
            printf("Realloc Failed\n");
            exit(EXIT_FAILURE);
        }
        pntr = mypntr + picket_size + debug_size;
        
        /* add new entry to hash table */
        tl_add_hash_entry(pntr, isize, linenumber, filename);
        
        /* Set Picket on either side of useful memory */
        tl_set_picket(mypntr + debug_size, isize);
        
        /* Store Size at very beginning of memory */
        *((unsigned long *)(mypntr)) = (unsigned long)(isize);
    }
    
    /* Update allocated_memory count */
    allocated_memory -= old_size;
    allocated_memory += isize;
    max_allocated_memory = MAX(max_allocated_memory, allocated_memory);
    return (pntr);
}
#endif
/************************************************************/
/************************************************************/
