#ifndef H_SW2_
#define H_SW2_

/***********************************************************/
/***********************************************************/
/***********************************************************/
//routines
int fe_sw2_body_resid(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
double fe_sw2_get_roughness(SMODEL_SUPER *mod, double depth, double velocity, int id, int which);
void fe_sw2_wdflag_legacy(SMODEL_SUPER *sm);
int fe_sw2_init(SMODEL_SUPER *sm);
//all boundary integrals (all weakly enforced boudaries)
//go into generic assembler so unified interface required
int fe_sw2_bc_discharge(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_ele(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_flapd(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_flapu(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_flux(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_h(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_hybrid(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_outflow(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_sluiced(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_sluiceu(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_tailwater(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_vel(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_vel_ele(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_weird(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
int fe_sw2_bc_weiru(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
/***********************************************************/
/***********************************************************/

#endif
