#include "global_header.h"

void read_bc_FR(SMODEL *mod, char *data)
{

  char line[MAXLINE];           /* the input line */
  char *subdata = NULL;         /* the data after the second card   is read */
  char *subsubdata = NULL;
  int ibc = 0;

  SIO info = *(mod->io);   // alias
  int nstring = mod->nstring;           // alias

  switch (parse_card(data, &subdata)) {

  case CARD_MNG:
    ibc = get_string_id(info, &subdata, nstring); 
    mod->str_values[ibc].roughness = read_dbl_field(info, &subdata);
    mod->str_values[ibc].fterms.manningsn = mod->str_values[ibc].roughness;
    mod->str_values[ibc].fterms.mng_flag = YES;
    break;

  case CARD_ERH:
    ibc = get_string_id(info, &subdata, nstring);
    if (mod->str_values[ibc].fterms.eqrheight == UNSET_FLT) {
      mod->str_values[ibc].fterms.eqrheight = read_dbl_field(info, &subdata);
    }
    mod->str_values[ibc].fterms.erh_flag = YES;
    break;

  default:
    break;

  }
}
