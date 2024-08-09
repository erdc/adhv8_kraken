/*!
   \file debug.h
   \brief Include file for -D_DEBUG work

   These variables are used to set pickets around allocated memory.  
   The pickets are filled with nans that may be trapped when
   accessed.  
   ??.200? Written by Hensley.  Added by Howington. 
   02.2006 Updated by Eslinger
 */

//****************************************************************
//****************************************************************
//****************************************************************

#ifndef H_DEBUG_
#define H_DEBUG_


#ifdef _DEBUG
#define BWIDTH 8		/* number of int-size entities for buffer on each side of data */
#define PICKET 0xfffb5b5a	/* guardpost value */
#define MY_NAN 0xfffa5a5a	/* NaN to use to initialize memory buffers */

/* HASHSIZE used for hash function */
#define DEBUG_HASHSIZE 65537

/* structure for linked lists of allocated memory locations */
typedef struct debug_mem_entry {
  struct debug_mem_entry *next;
  unsigned long address;
  size_t isize;
  int linenumber;
  char *filename;
} DEBUG_MEM_ENTRY;

void tl_check_picket(char *, size_t, int, char *);
void tl_set_nan(char *, size_t);
void tl_set_picket(char *, size_t);
void tl_add_hash_entry(void *pntr, size_t isize, int linenumber, char *filename);
void tl_remove_hash_entry(void *pntr);
void tl_check_all_pickets(char *,int);
void tl_check_unfreed(void);
void Push_Hash_Entry(DEBUG_MEM_ENTRY **,unsigned long,size_t,int,char *);
void Delete_Hash_Entry(DEBUG_MEM_ENTRY **,unsigned long);
int hash_function(unsigned long);
int LengthList(DEBUG_MEM_ENTRY *);
void CreateHistogram(int, int, int);
unsigned long tl_get_allocated_memory_size(void *);
#endif

#endif /// ifndef H_DEBUG_

//****************************************************************
//****************************************************************
//****************************************************************


#ifdef _DEBUG

#ifdef DEBUG_INIT
#define EXTERN_DEBUG
#else
#define EXTERN_DEBUG extern
#endif

EXTERN_DEBUG size_t allocated_memory;
EXTERN_DEBUG size_t max_allocated_memory;
EXTERN_DEBUG DEBUG_MEM_ENTRY **mem_hash;
EXTERN_DEBUG int picket_value;
EXTERN_DEBUG char *picket_pointer;
EXTERN_DEBUG int nan_value;
EXTERN_DEBUG char *nan_pointer;
EXTERN_DEBUG FILE *outalloc;
EXTERN_DEBUG FILE *outfree;

#else

extern int allocated_memory;
extern int max_allocated_memory;

#endif
