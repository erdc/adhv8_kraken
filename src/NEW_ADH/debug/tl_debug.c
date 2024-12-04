/*!
 \file tl_debug.c
 \brief Checking for allocations when compiled with "-D_DEBUG"
 
 tl_debug and associated files:
 ??.200? - Written by Jeff Henley
 01.2006 - Attempted to make general for all machines, Jeff Hensley
 02.2006 - Added Hash Tables, Jeff Hensley
 
 What we've tried to do: Put a "picket" on either side of allocated
 memory.  Inside this picket, place NaN, so that if something grabs
 that memory, and then uses it, a.out will die.  At time of releasing
 allocated memory, this also checks the pickets to make sure they
 weren't altered in any way.  Also, we can place NaN's inside of
 the useful memory, so that if it is used before initialization, you
 should get an execption.
 
 Note: Here, we've tried to take "picket_value", stored as an int.
 Then use a pointer to char to step through each byte of
 picket_value.  The thinking is that chars are 1-byte.  So all
 pointer arithmetic will be valid and we do not need to worry about
 adding/subtracting a specific number of bytes from it.
 */
#include "header_tl_alloc.h"
#include "debug.h"
#ifdef _DEBUG
/*!
 \brief Set Pickets Before/After Allocate Memory
 */
void tl_set_picket(char *pntr, size_t isize) {

#ifdef _MESSG
    int myid, ierr_code;
    ierr_code = MPI_Comm_rank(debug_comm, &myid); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD, &myid);
#else
    int myid = 0;
#endif

    size_t ii = 0, jj = 0;      /* loop counters */
    char *mypntr = NULL;        /* Pointer to 1 Byte */
    
    /* Fill picket fence in before array */
    mypntr = pntr;
    for (ii = 0; ii < BWIDTH; ii++) {
        for (jj = 0; jj < sizeof(int); jj++) {
            *(mypntr++) = picket_pointer[jj];
        }
    }

    
    /* Fill picket fence in after array */
    mypntr = pntr + isize + sizeof(int) * ((size_t) (BWIDTH));
    for (ii = 0; ii < BWIDTH; ii++) {
        for (jj = 0; jj < sizeof(int); jj++) {
            *(mypntr++) = picket_pointer[jj];
        }
    }
    return;
}

/*!
 \brief Initialize memory to NaN
 */
void tl_set_nan(char *pntr, size_t isize)
{
    size_t ii = 0, jj = 0;      /* loop counters */
    size_t num = 0;             /* num. of int-size words allocated */
    char *mypntr = NULL;        /* Pointer to 1 Byte */
    
    mypntr = pntr;
    num = isize / sizeof(int);
    for (ii = 0; ii < num; ii++) {
        for (jj = 0; jj < sizeof(int); jj++) {
            *(mypntr++) = nan_pointer[jj];
        }
    }
    return;
}

/*!
 \brief Check Picket, Each Side of Allocated Memory
 */
void tl_check_picket(char *pntr, size_t isize, int linenumber, char *filename)
{
    size_t ii = 0, jj = 0;      /* loop counters */
    char *mypntr = NULL;        /* Pointer to 1 Byte */
#ifdef _MESSG
    int myid, ierr_code;
    ierr_code = MPI_Comm_rank(debug_comm, &myid); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD, &myid);
#else
    int myid = 0;
#endif
    /* Check Under Index */
    mypntr = pntr - sizeof(int) * ((size_t) (BWIDTH));
    for (ii = 0; ii < BWIDTH; ii++) {
        for (jj = 0; jj < sizeof(int); jj++) {
            if (*(mypntr++) != picket_pointer[jj]) {
                printf("Corrupted Picket - Lower (file:line) %s:%d MYID %d \n", filename, linenumber, myid);
                exit(EXIT_FAILURE);
            }
        }
    }
    
    /* Check Over Index */
    mypntr = pntr + isize;
    for (ii = 0; ii < BWIDTH; ii++) {
        for (jj = 0; jj < sizeof(int); jj++) {
            if (*(mypntr++) != picket_pointer[jj]) {
                printf("Corrupted Picket - Upper (file:line) %s:%d MYID %d \n", filename, linenumber, myid);
                exit(EXIT_FAILURE);
            }
        }
    }
    return;
}

/*!
 \brief Check All the Pickets
 */
void tl_check_all_pickets(char *filename, int linenumber)
{
    size_t ii = 0;              /* loop counter */
    
#ifdef _MESSG
    int myid, ierr_code;
    ierr_code = MPI_Comm_rank(debug_comm, &myid); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD, &myid);
#else
    int myid = 0;
#endif
    DEBUG_MEM_ENTRY *entry = NULL;
    
    printf("     myid %4d, checking all pickets at %s:%d, ", myid, filename, linenumber);
    for (ii = 0; ii < DEBUG_HASHSIZE; ii++) {
        entry = mem_hash[ii];
        
        while (entry != NULL) {
            tl_check_picket((char *) entry->address, entry->isize, entry->linenumber, entry->filename);
            entry = entry->next;
        }
    }
    printf("Finished\n");
    return;
}

/*!
 \brief Get Allocated memory size
 */
unsigned long tl_get_allocated_memory_size(void *pntr)
{
    size_t picket_size = 0;     /* the size of the picket for memory debugging */
    size_t debug_size = 0;      /* the size of the array for memory debugging */
    char *mypntr = NULL;        /* pointer for memory purposes */
    
    /* Set up Pieces */
    picket_size = sizeof(int) * ((size_t) (BWIDTH));
    debug_size = sizeof(unsigned long);
    
    /* Find True Beginning of Memory, where Info is stored */
    mypntr = (char *) pntr;
    mypntr = mypntr - picket_size - debug_size;
    
    return (*((unsigned long *) (mypntr)));
}

/*!
 \brief Check for Unfreed Memory
 */
void tl_check_unfreed(void)
{
    int ii = 0;
    DEBUG_MEM_ENTRY *entry = NULL;
    
#ifdef _MESSG
    int myid, ierr_code;
    ierr_code = MPI_Comm_rank(debug_comm, &myid); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD, &myid);
#else
    int myid = 0;
#endif
    
    for (ii = 0; ii < DEBUG_HASHSIZE; ii++) {
        entry = mem_hash[ii];
        
        while (entry != NULL) {
            printf("WARNING: MYID: %d :: Unfreed Memory from (file:line): %s:%d\n", myid, entry->filename, entry->linenumber);
            entry = entry->next;
        }
    }
    return;
}

/*!
 \brief Add Hash Entry
 */
void tl_add_hash_entry(void *pntr, size_t isize, int linenumber, char *filename)
{
    int index = 0;
    unsigned long address = 0;
    
    address = ((unsigned long) (pntr));
    index = hash_function(address);
    Push_Hash_Entry(&mem_hash[index], address, isize, linenumber, filename);
    return;
}

/*!
 \brief Remove Hash Entry
 */
void tl_remove_hash_entry(void *pntr)
{
    int index = 0;
    unsigned long address = 0;
    
    address = ((unsigned long) (pntr));
    index = hash_function(address);
    Delete_Hash_Entry(&mem_hash[index], address);
    return;
}

/*!
 \brief Our Hash Function
 */
int hash_function(unsigned long address)
{
    return (address % DEBUG_HASHSIZE);
}

/*!
 \brief Push an Entry onto Hash Table
 */
void Push_Hash_Entry(DEBUG_MEM_ENTRY ** headRef, unsigned long address, size_t isize, int linenumber, char *filename)
{
    size_t len = 0;
    
    DEBUG_MEM_ENTRY *newNode = (DEBUG_MEM_ENTRY *) malloc(sizeof(DEBUG_MEM_ENTRY));
    
    newNode->address = address;
    newNode->isize = isize;
    newNode->linenumber = linenumber;
    
    len = strlen(filename) + 1;
    newNode->filename = (char *) malloc(len);
    memcpy(newNode->filename, filename, len);
    
    newNode->next = *headRef;
    *headRef = newNode;
    
    return;
}

/*!
 \brief Delete an Entry in Hash Table
 */
void Delete_Hash_Entry(DEBUG_MEM_ENTRY ** head, unsigned long address)
{
    DEBUG_MEM_ENTRY *entry = NULL, *next = NULL, *prev = NULL;
    entry = *head;
    
    while (entry != NULL) {
        next = entry->next;
        if (entry->address == address) {    /* found it! */
            if (prev == NULL) {
                *head = next;
                free(entry->filename);
                free(entry);
                return;
            }
            else {
                prev->next = next;
                free(entry->filename);
                free(entry);
                return;
            }
        }
        else {
            prev = entry;
            entry = next;
        }
    }
    
    printf("ERROR: Did not find entry in hash table\n");
    printf("Possibly already freed this memory\n");
    exit(EXIT_FAILURE);
}

/*!
 \brief Return the length of linked list
 */
int LengthList(DEBUG_MEM_ENTRY * head)
{
    int count = 0;
    DEBUG_MEM_ENTRY *current = head;
    
    while (current != NULL) {
        count++;
        current = current->next;
    }
    
    return count;
}

/* iwidth = width of bins for histogram
 istart = starting size for bin
 inum   = number of bins; **there will be two more bins, one for
 < istart, another for > (istart+ inum*iwidth)
 
 For each range (lenght of lists) is printed the number of lists that
 have lenght within that range.
 Lists of length zero are not included in the count.
 
 Example: CreateHistogram(1,2,5) would produce a count of lists with
 1-2 elements
 3-4 elements
 5-6 elements
 7-8 elements
 9-10 elements
 >10 elements
 */
void CreateHistogram(int istart, int iwidth, int inum)
{
    int length = 0;             /* Length of List */
    int ii = 0;                 /* Loop Counter */
    int *bins = NULL;           /* Hash Table Bins */
    
    if (inum == -1) {
        for (ii = 0; ii < DEBUG_HASHSIZE; ii++)
            length += LengthList(mem_hash[ii]);
        printf("Total Number of Debug Hash Entries: %d\n", length);
        return;
    }
    
    bins = (int *) malloc(sizeof(int) * (inum + 2));
    for (ii = 0; ii < inum + 2; ii++)
        bins[ii] = 0;
    
    for (ii = 0; ii < DEBUG_HASHSIZE; ii++) {
        length = LengthList(mem_hash[ii]);
        if (length == 0)
            continue;
        else if (length < istart)
            (bins[0])++;
        else if (length > (istart + iwidth * inum)) {
            (bins[inum + 1])++;
        }
        else {
            (bins[(length - istart) / iwidth + 1])++;
        }
    }
    
    /* Print out the histogram */
    printf("HISTOGRAM for DEBUG HASH TABLE\n\n");
    printf("Range    No. of Entries\n");
    if (iwidth == 1) {
        printf(" <%d      %d\n", istart, bins[0]);
        for (ii = 1; ii < inum + 1; ii++) {
            printf(" %d       %d\n", istart + ii - 1, bins[ii]);
        }
        printf(" >%d      %d\n", istart + iwidth * inum - 1, bins[inum + 1]);
    }
    else {
        printf(" <%d      %d\n", istart, bins[0]);
        for (ii = 1; ii < inum + 1; ii++) {
            printf(" %d-%d       %d\n", istart + (ii - 1) * iwidth - 1, istart + ii * iwidth - 1, bins[ii]);
        }
        printf(" >%d      %d\n", istart + iwidth * inum - 1, bins[inum + 1]);
    }
    return;
}

#endif
