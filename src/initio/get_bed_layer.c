/* gets and checks the layer number */

#include "global_header.h"

int get_bed_layer(SIO io, char **data, int nlayers) {
  /* gets the string and checks it */
  int ilayer = read_int_field_custom(io, data, NULL, "bed layer ID", UNSET_INT, TRUE);
  if (ilayer > nlayers || ilayer < 1) {
    io_read_error(io, "The specified bed layer does not exist.", TRUE);
  }
  ilayer--;
  return (ilayer);
}
