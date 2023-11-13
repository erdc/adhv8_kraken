/* gets and checks the material number */

#include "global_header.h"

int get_material_id(SIO io, char **data, int nmat)
{

  /* gets and checks the material number */
  int imat = read_int_field_custom(io, data, NULL, "material ID", UNSET_INT, TRUE);
  if (imat > nmat || imat < 1) {
    io_read_error(io, "The specified material does not exist.", TRUE);
  }
  imat--;
  return imat;
}
