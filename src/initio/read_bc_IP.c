
#include "global_header.h"

void read_bc_IP(SMODEL *mod, char *data) {

  char line[MAXLINE];           /* the input line */
  char *subdata = NULL;         /* the data after the second card   is read */
  char *subsubdata = NULL;
  SIO io = *(mod->io);   // alias

  switch (parse_card(data, &subdata)) {

  case CARD_MIT:
    mod->solver_info.max_lin_it = read_int_field(io,&subdata);
    mod->solver_info.force_lin_it = NO;
    break;

  case CARD_NIT:
    mod->max_nonlin_it = read_int_field(io,&subdata);
    mod->solver_info.force_nonlin_it = NO;
    break;

  case CARD_ITL:
    mod->inc_nonlin = read_dbl_field(io,&subdata);
    break;

  case CARD_NTL:
    mod->tol_nonlin = read_dbl_field(io,&subdata);
    break;

  default:
    break;

  }
}
