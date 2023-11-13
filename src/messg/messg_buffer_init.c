#include "global_header.h"

/*!
   \brief Initializes a buffer 

   See type.h for a better description
 */
void messg_buffer_init(MESSG_BUFFER * buffer,   /* the message buffer */
                       int i_processor  /* the target processor */
  )
{
  buffer->size = 0;
  buffer->nitem = 0;
  buffer->type = UNSET_INT;
  buffer->sd = i_processor;
  buffer->pos = 0;
  buffer->buffer = NULL;
  return;
}
