#include "global_header.h"

void read_bc_TC(SMODEL *mod, char *data) {

  char line[MAXLINE];           /* the input line */
  char *subdata = NULL;         /* the data after the second card   is read */
  char *subsubdata = NULL;
  int query, units;

  SIO info = *(mod->io);        // alias

  switch (parse_card(data, &subdata)) {

  case CARD_T0:
    mod->t_prev = mod->t_init = read_dbl_field(info, &subdata);
    /* get optional units flag */
    units = read_query_int_field(info, &subdata, &query);
    switch (query) {
    case READ_SUCCESS:
      mod->t_init *= tc_conversion_factor(units, TO);
      break;
    case READ_FAIL:
      printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
      io_read_error(info, "Unable to read optional units card. Encountered " "an invalid integer input.", TRUE);
      break;
    case READ_NONE:
    default:
      break;
    }
    break;

  case CARD_TF:
    mod->t_final = read_dbl_field(info, &subdata);
    /* get optional units flag */
    units = read_query_int_field(info, &subdata, &query);
    switch (query) {
    case READ_SUCCESS:
      mod->t_final *= tc_conversion_factor(units, TO);
      break;
    case READ_FAIL:
      io_read_error(info, "Unable to read optional units card. Encountered " "an invalid integer input.", TRUE);
      break;
    case READ_NONE:
    default:
      break;
    }
    break;

  case CARD_NDP:
    mod->t_adpt_flag = OFF;
    break;

  default:
    break;

  }
}
