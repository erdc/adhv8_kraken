#include "global_header.h"

void read_bc_OFF(SMODEL * mod, char *data)
{

  char line[MAXLINE];           /* the input line */
  char *subdata = NULL;         /* the data after the second card   is read */
  char *subsubdata = NULL;

  int ibc = 0;

  ibc = get_string_id(*(mod->io), &data, mod->nstring);
  mod->str_values[ibc].phys_flag = OFF;
}
