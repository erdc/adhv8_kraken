#ifndef H_FE_TRANSPORT_PROTOTYPES_
#define H_FE_TRANSPORT_PROTOTYPES_

/***********************************************************/
/***********************************************************/
/***********************************************************/

int  fe_transport(SSUPER_MODEL *sm, int imod);
void fe_transport_load(SSUPER_MODEL *sm, int imod);
void fe_transport_resid(SSUPER_MODEL *sm, int imod);
void fe_transport_update(SSUPER_MODEL *sm, int imod);
void fe_transport_inc(SSUPER_MODEL *sm, int imod);
void fe_transport_init(SSUPER_MODEL *sm, int imod);

int fe_transport_hybrid(SSUPER_MODEL *sm);
void fe_transport_hybrid_inc(SSUPER_MODEL *sm, int imod);
void fe_transport_hybrid_init(SSUPER_MODEL *sm, int imod);
void fe_transport_hybrid_load(SSUPER_MODEL *sm, int imod);
void fe_transport_hybrid_resid(SSUPER_MODEL *sm, int imod);
void fe_transport_hybrid_update(SSUPER_MODEL *sm, int imod);

void fe_2d_transport_body_resid(SMODEL *mod, double *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_2d_transport_boundary_resid(SMODEL *mod, double *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_2d_transport_wd_convection_triangle(SVECT2D *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *elem_rhs);
void fe_2d_transport_wd_temporal_triangle(SVECT2D *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *elem_rhs);
void fe_3d_transport_body_resid(SMODEL *mod, double *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_3d_transport_boundary_resid(SMODEL *mod, double *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
