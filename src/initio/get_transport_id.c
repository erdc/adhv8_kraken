/* gets and checks the string number */

#include "global_header.h"

int get_transport_id(SIO io, char **data, int ntransport)
{

  /* gets the string and checks it */
  int itrns = read_int_field_custom(io, data, NULL, "transport constituent ID", UNSET_INT, TRUE);
  if (itrns > ntransport || itrns < 1) {
    io_read_error(io, "The transport consistituent ID is either too big or too small.", TRUE);
  }
  itrns--;
  return (itrns);
}
