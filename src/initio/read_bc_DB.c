#include "global_header.h"

void read_bc_DB(SMODEL *mod, char *data)
{

  char line[MAXLINE];           /* the input line */
  char *subdata = NULL;         /* the data after the second card   is read */
  char *subsubdata = NULL;
  int ibc = UNSET_INT, iseries = UNSET_INT;
  int itrn;

  SIO info = *(mod->io);        // alias
  int nstring = mod->nstring;   // alias    
  int nseries = mod->nseries;   // alias
  switch (parse_card(data, &subdata)) {

  case CARD_VEL:
    ibc = get_string_id(info, &subdata, nstring);
    mod->str_values[ibc].flow.bc_flag = BCT_VEL_DIR;
    iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].flow.ivx = iseries;
    iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].flow.ivy = iseries;
    iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].flow.ivz = iseries;
    //printf("DB VEL READ :: ibc : %d nseries: %d nstring: %d iseries: %d\n",ibc,nseries,nstring,iseries);
    break;
  case CARD_PRS:
    ibc = get_string_id(info, &subdata, nstring);
    mod->str_values[ibc].pressure.bc_flag = BCT_PRS_DIR;
    iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].pressure.iu_0 = iseries;
    break;
  case CARD_FRS:
    ibc = get_string_id(info, &subdata, nstring);
    mod->str_values[ibc].displacement.bc_flag = BCT_FREE_DIR;
    iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].displacement.iu_0 = iseries;

    break;
  case CARD_OVL:
    ibc = get_string_id(info, &subdata, nstring);
    mod->str_values[ibc].ol_flow.bc_flag = BCT_VEL_DIR;
    iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].ol_flow.ivx = iseries;
    iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].ol_flow.ivy = iseries;
    break;
  case CARD_OVH:
    ibc = get_string_id(info, &subdata, nstring);
    mod->str_values[ibc].ol_flow.bc_flag = BCT_VEL_PRS_DIR;
    iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].ol_flow.ivx = iseries;
    iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].ol_flow.ivy = iseries;
	iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].ol_flow.iu_0 = iseries;
    break;
  case CARD_TRN:
    ibc = get_string_id(info, &subdata, nstring);
    itrn = get_transport_id(info, &subdata, mod->ntransport);
    mod->str_values[ibc].trans[itrn].bc_flag = BCT_DIR;
    iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].trans[itrn].iu_0 = iseries;
    break;
#ifdef _ADH_GROUNDWATER
  case CARD_THD:
    ibc = get_string_id(info, &subdata, nstring);
    mod->str_values[ibc].flow.bc_flag = BCT_DIR;
    iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].flow.iu_0 = iseries;
    break;
  case CARD_FLW:
    ibc = get_string_id(info, &subdata, nstring);
    mod->str_values[ibc].flow.bc_flag = BCT_DIR;
    iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].flow.iu_0 = iseries;
    break;
  case CARD_PSI:
    ibc = get_string_id(info, &subdata, nstring);
    mod->str_values[ibc].flow.bc_flag = BCT_PRS_DIR;
    iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
    mod->str_values[ibc].flow.iu_0 = iseries;
    break;
#endif
  default:

    break;

  }
}
