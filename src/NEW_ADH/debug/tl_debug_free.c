/*!
 \file tl_free.c
 \brief Free Allocated Memory
 */

#include "header_tl_alloc.h"
#include "debug.h"

/*!
 \brief Free allocated memory
 
 This routine frees allocated memory.  When compiled with
 "-D_DEBUG" it checks the original size allocated and
 examines the pickets on either side.  See description in
 tl_alloc.c
 
 On return, give a NULL value so that the pointer to the
 memory just freed can be set to NULL.
 
 Note: mypntr below is a "char".  The thinking is that chars
 are 1-byte.  So all pointer arithmetic will be valid and
 we do not need to worry about adding/subtracting a specific
 number of bytes from it.
 
 \param obj_size sizeof(double), sizeof(int), etc.
 \param number Number of Units needed
 \param pntr Pointer to Memory being freed
 */

/************************************************************/
/************************************************************/
#ifndef _DEBUG
void *tl_free(
              size_t obj_size,		/* the size of the object to be allocated */
              int number,			/* the number of objects to be allocated */
              void *pntr			/* the pointer */
)
{
    size_t isize = 0;		/* the size of the space to be allocated */
    if((obj_size != 0) && (number > 0)) {
        isize = obj_size * ((size_t) (number));
        free(pntr);
        allocated_memory -= isize;
    }
    return (NULL);
}

/************************************************************/
/************************************************************/
#else
/*!
 \brief Free allocated memory (Debug Version)
 
 \param obj_size sizeof(double), sizeof(int), etc.
 \param number Number of Units needed
 \param linenumber Should be macro __LINE__
 \param *filename Should be macro __FILE__
 \param pntr Pointer to Memory being freed
 */
void *tl_free_debug(
                    size_t obj_size,    /* the size of the object to be allocated */
                    int number,			/* the number of objects to be allocated */
                    int linenumber,		/* line number */
                    char *filename,		/* calling routine filename */
                    void *pntr			/* the pointer */
)
{
    size_t isize = 0;		/* the size of the space to be allocated */
    size_t picket_size = 0;	/* the size of the picket for memory debugging */
    size_t debug_size = 0;	/* the size of the array for memory debugging */
    char *mypntr = NULL;		/* pointer for memory purposes */
#ifdef _MESSG
    int myid, ierr_code;
    ierr_code = MPI_Comm_rank(debug_comm, &myid); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD, &myid);
#else
    int myid = 0;
#endif
    
    if((obj_size == 0) || (number < 0)) {
        printf("MYID %4d, Error (file:line) %s:%d\n", myid, filename, linenumber);
        printf("Error in isize, this is an error\n");
        exit(8);
    } else {
        isize = obj_size * ((size_t) (number));
        picket_size = sizeof(int) * ((size_t) (BWIDTH));
        debug_size = sizeof(unsigned long);
        
        /* Check allocated vs. free size */
        mypntr = (char *)(pntr);
        mypntr = mypntr - picket_size - debug_size;
        if((unsigned long)(isize) != *((unsigned long *)(mypntr))) {
            printf("MYID %4d,Error (file:line) %s:%d\n", myid, filename, linenumber);
            printf("Freeing-> %lu , %lu <-Originally Bytes\n", isize, *((unsigned long *)(mypntr)));
            exit(9);
        }
        
        /* Check pickets */
        mypntr = (char *)(pntr);
        tl_check_picket(mypntr, isize, linenumber, filename);
        
        /* Finally free it */
        mypntr = (char *)(pntr);
        free(mypntr - picket_size - debug_size);
        allocated_memory -= isize;
        
        /* Remove entry in hash table */
        tl_remove_hash_entry(pntr);
    }
    return (NULL);
}
#endif

/************************************************************/
/************************************************************/
/************************************************************/

#ifndef _DEBUG
/* double array free */
void *tl_free2(
               size_t obj1_size,		/* the size of the first object to be allocated */
               int number1,			/* the number of first objects to be allocated */
               size_t obj2_size,      /* the size of the second object to be allocated */
               int number2,			/* the number of second objects to be allocated */
               void **pntr			/* the pointer */
)
{
    int i=0;
    for (i=0; i<number1; i++) {
        if (pntr[i] != NULL) {
            pntr[i] = (void *) tl_free(obj2_size, number2, pntr[i]);
        }
    }
    pntr = (void **) tl_free(obj1_size, number1, pntr);
    return (NULL);
}

#else
void *tl_free2_debug(
                     size_t obj1_size,      /* the size of the first object to be allocated */
                     int number1,           /* the number of first objects to be allocated */
                     size_t obj2_size,      /* the size of the second object to be allocated */
                     int number2,           /* the number of second objects to be allocated */
                     int line,
                     char * file,
                     void **pntr            /* the pointer */
)
{
    int i=0;
    for (i=0; i<number1; i++) {
        if (pntr[i] != NULL) {
            pntr[i] = (void *) tl_free_debug(obj2_size, number2, line, file, pntr[i]);
        }
    }
    pntr = (void **) tl_free_debug(obj1_size, number1, line, file, pntr);
    return (NULL);
}
#endif


/************************************************************/
/************************************************************/
/************************************************************/

#ifndef _DEBUG
/* triple array free */
void *tl_free3(
               size_t obj1_size,		/* the size of the first object to be allocated */
               int number1,			/* the number of first objects to be allocated */
               size_t obj2_size,      /* the size of the second object to be allocated */
               int number2,			/* the number of second objects to be allocated */
               size_t obj3_size,      /* the size of the second object to be allocated */
               int number3,			/* the number of second objects to be allocated */
               void ***pntr			/* the pointer */
)
{
    int i=0;
    for (i=0; i<number1; i++) {
        if (pntr[i] != NULL) {
            pntr[i] = (void **) tl_free2(obj2_size, number2, obj3_size, number3, pntr[i]);
        }
    }
    pntr = (void ***) tl_free(obj1_size, number1, pntr);
    return (NULL);
}

#else
/* triple array free */
void *tl_free3_debug(
                     size_t obj1_size,      /* the size of the first object to be allocated */
                     int number1,           /* the number of first objects to be allocated */
                     size_t obj2_size,      /* the size of the second object to be allocated */
                     int number2,           /* the number of second objects to be allocated */
                     size_t obj3_size,      /* the size of the second object to be allocated */
                     int number3,           /* the number of second objects to be allocated */
                     int line,
                     char * file,
                     void ***pntr           /* the pointer */
)
{
    int i=0;
    for (i=0; i<number1; i++) {
        if (pntr[i] != NULL) {
            pntr[i] = (void **) tl_free2_debug(obj2_size, number2, obj3_size, number3, line, file, pntr[i]);
        }
    }
    pntr = (void  ***) tl_free_debug(obj1_size, number1, line, file, pntr);
    return (NULL);
}
#endif
