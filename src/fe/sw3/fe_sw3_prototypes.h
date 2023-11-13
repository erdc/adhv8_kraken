#ifndef H_FE_SW3_PROTOTYPES_
#define H_FE_SW3_PROTOTYPES_

/***********************************************************/
/***********************************************************/
/***********************************************************/

void fe_wvel_body_resid(SMODEL *mod, double *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_hvel_body_resid(SMODEL *mod,  DOF_3 *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int *node_in_column_flag, int DEBUG);
void fe_wvel_boundary_resid(SMODEL *mod, double *elem_rhs, double *elem_rhs_noBed, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_hvel_boundary_resid(SMODEL *mod,  DOF_3 *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int *node_in_column_flag, int DEBUG);
void fe_wvel_inc(SSUPER_MODEL *sm, int imod);
void fe_hvel_inc(SSUPER_MODEL *sm, int imod);
void fe_wvel_init(SSUPER_MODEL *sm, int imod);
void fe_hvel_init(SSUPER_MODEL *sm, int imod);
void fe_wvel_load(SSUPER_MODEL *sm, int imod);
void fe_hvel_load(SSUPER_MODEL *sm, int imod);
void fe_wvel_resid(SSUPER_MODEL *sm, int imod);
void fe_hvel_resid(SSUPER_MODEL *sm, int imod);
void fe_wvel_update(SSUPER_MODEL *sm, int imod);
void fe_hvel_update(SSUPER_MODEL *sm, int imod);
int fe_sw3(SMODEL *mod);
int fe_sw3_solve(SSUPER_MODEL *sm, int imod);
double fe_sw3_get_roughness(SMODEL *mod, double depth, int id);

void fe_weighted_normals(SMODEL *mod); /* Gajanan gkc proposing new addition */

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
