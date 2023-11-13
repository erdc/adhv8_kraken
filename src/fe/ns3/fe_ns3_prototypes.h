#ifndef H_FE_NS3_PROTOTYPES_
#define H_FE_NS3_PROTOTYPES_

/***********************************************************/
/***********************************************************/
/***********************************************************/

void fe_ns3_body_resid(SMODEL *mod,  DOF_4 *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_ns3_boundary_resid(SMODEL *mod,  DOF_4 *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_ns3_inc(SSUPER_MODEL *sm, int imod);
void fe_ns3_init(SSUPER_MODEL *sm, int imod);
void fe_ns3_load(SSUPER_MODEL *sm, int imod);
void fe_ns3_resid(SSUPER_MODEL *sm, int imod);
void fe_ns3_update(SSUPER_MODEL *sm, int imod);
int fe_ns3_solve(SSUPER_MODEL *sm, int imod);
int fe_ns3(SMODEL *mod);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
