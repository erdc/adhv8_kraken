/* ADH Version 2.0.0 6-04 */
/* allocate more space for a communication buffer */

#include "global_header.h"

#ifdef _MESSG
void messg_alloc(int nitem,     /* the number of items to be stored in the buffer */
                 int sizeof_item,   /* the size of the items */
                 MESSG_BUFFER * buffer  /* the buffer to be packed */
  )
{
  int new_size;                 /* the new size of the buffer */

  /* allocate the buffer and set the number of items */
  buffer->nitem = nitem;
  new_size = nitem * sizeof_item;
  if (new_size > buffer->size)

    {
      buffer->buffer = tl_realloc(ONE, new_size, buffer->size, buffer->buffer);
      buffer->size = new_size;
    }
}

#else /*  */
void messg_alloc(void)
{
}
#endif /*  */
