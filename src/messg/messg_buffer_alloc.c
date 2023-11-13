#include "global_header.h"

/*! 
   \brief Allocate more space for a communication buffer

   See type.h for a better description
 */
void messg_buffer_alloc(int nitem,  /* the number of items to be stored in the buffer */
                        size_t sizeof_item, /* the size of the items */
                        MESSG_BUFFER * buffer   /* the buffer to be packed */
  )
{
#ifdef _MPI
  int new_allocated_size = 0;   /* the new size of the buffer */

  /* Do some checking */
  assert(nitem >= 0);
  assert(sizeof_item != 0);
  assert(buffer != NULL);

  /* Allocate the buffer and set the number of items */
  buffer->nitem = nitem;
  new_allocated_size = nitem * ((int) (sizeof_item));
  if (new_allocated_size > buffer->size)
    {
      buffer->buffer = tl_realloc(ONE, new_allocated_size, buffer->size, buffer->buffer);
      buffer->size = new_allocated_size;
    }
#endif
  return;
}
