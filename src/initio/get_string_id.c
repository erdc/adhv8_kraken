/* gets and checks the string number */

#include "global_header.h"

int get_string_id(SIO io, char **data, int nstring)
{

  /* gets the string and checks it */
  int istring = read_int_field_custom(io, data, NULL, "boundary string ID", UNSET_INT, TRUE);
  if (istring > nstring || istring < 1) {
    io_read_error(io, "The specified boundary string ID does not exist.", TRUE);
  }
  istring--;
  return (istring);
}
