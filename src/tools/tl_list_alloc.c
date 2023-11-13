/* this routine allocates the storage for linked lists */
/* cjt :: reconsider this file when multiple models are being ran concurrently */

#include "global_header.h"

/* local definitions */
#define NMEM_ROW_INC 100    /* the number of rows */
#define ROW_SIZE 1000   /* the size of the rows */

/* variables shared by tl_list_alloc and tl_list_free_all */

static int imem[NLIST];         /* the current memory row */
static int jmem[NLIST];         /* the current memory column -
                                   note that jmem is initialized to ROW_SIZE so first row will be
                                   allocated on the first pass */
static int current_size[NLIST]; /* the number of row pointers allocated */
static char **mem[NLIST];       /* the rows and columns of memory for the linked list */

static int object_size[NLIST] = {
  sizeof(EDGE_LIST_ITEM),
  sizeof(FACE_LIST_ITEM),
  sizeof(ELEM1D_LIST_ITEM),
  sizeof(ELEM2D_LIST_ITEM),
  sizeof(ELEM3D_LIST_ITEM),
#ifdef _MESSG
  sizeof(NODE_LIST_ITEM),
  sizeof(ELEM_REF_LIST_ITEM)
#endif
}; /* the size of the different list items */

void tl_list_setup(void)
{
  int i;                        /* loop counter */

  for (i = 0; i < NLIST; i++) {
    imem[i] = -1;
    jmem[i] = ROW_SIZE;
    current_size[i] = 0;
    mem[i] = NULL;
  }
}

void *tl_list_alloc(int list_flag   /* indicates which list is being set */
  )
{
  int new_current_size;         /* the new current size */
  int jmem_return;              /* the return column number for the memory */

  /* allocate the row pointers on first pass */
  if (imem[list_flag] == -1) {
    current_size[list_flag] = NMEM_ROW_INC;
    mem[list_flag] = (char **) tl_alloc(sizeof(char *), current_size[list_flag]);
  }

  /* allocate the next row if previous row is full */
  if (jmem[list_flag] == ROW_SIZE) {
    imem[list_flag]++;
    jmem[list_flag] = 0;
    if (imem[list_flag] == current_size[list_flag]) {
      new_current_size = current_size[list_flag] + NMEM_ROW_INC;
      mem[list_flag] = (char **) tl_realloc(sizeof(char *), new_current_size, current_size[list_flag], (char **) mem[list_flag]);
      current_size[list_flag] = new_current_size;
    }
    mem[list_flag][imem[list_flag]] = (char *) tl_alloc(object_size[list_flag], ROW_SIZE);
  }

  /* increment the column counter and return the pointer to the list item */
  jmem_return = jmem[list_flag];
  jmem[list_flag]++;
  return ((void *) &mem[list_flag][imem[list_flag]][jmem_return * object_size[list_flag]]);
}

void tl_list_free_all(int list_flag /* indicates which list is being set */
  )
{
  int i;                        /* loop counter */
  for (i = 0; i <= imem[list_flag]; i++) {
    mem[list_flag][i] = (char *) tl_free(object_size[list_flag], ROW_SIZE, mem[list_flag][i]);
  }
  if (current_size[list_flag] != 0) {
    mem[list_flag] = (char **) tl_free(sizeof(char *), current_size[list_flag], mem[list_flag]);
  }
  imem[list_flag] = -1;
  jmem[list_flag] = ROW_SIZE;
  current_size[list_flag] = 0;
}
