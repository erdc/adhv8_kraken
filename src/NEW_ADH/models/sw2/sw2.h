#ifndef H_SW2_
#define H_SW2_

/***********************************************************/
/***********************************************************/
/***********************************************************/

int fe_sw2_body_resid(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
double fe_sw2_get_roughness(SMODEL_SUPER *mod, double depth, double velocity, int id, int which);
void fe_sw2_wdflag_legacy(SMODEL_SUPER *sm);
int fe_sw2_init(SMODEL_SUPER *sm);
/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
