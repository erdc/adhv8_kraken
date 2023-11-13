#include "global_header.h"

/*! 
   \brief Free the given communication buffer 

   See type.h for a better description
 */
void messg_buffer_free(MESSG_BUFFER * buffer    /* the buffer to be packed */
  )
{
  /* free the buffer and reset the size */
  if ((buffer->buffer != NULL) && (buffer->size > 0))
    {
      buffer->buffer = tl_free(ONE, buffer->size, buffer->buffer);
    }

  buffer->size = 0;
  buffer->nitem = 0;
  buffer->type = UNSET_INT;
  buffer->sd = UNSET_INT;
  buffer->pos = 0;
  return;
}
