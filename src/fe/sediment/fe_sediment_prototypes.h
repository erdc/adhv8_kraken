#ifndef H_FE_SEDIMENT_TRANSPORT_PROTOTYPES_
#define H_FE_SEDIMENT_TRANSPORT_PROTOTYPES_

/***********************************************************/
/***********************************************************/
/***********************************************************/

int fe_sediment_transport(SMODEL *mod, int ndim);

// bed load transport functions
void fe_bedload_boundary_resid(SMODEL *mod, double *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_bedload_body_resid(SMODEL *mod, double *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_bedload_inc(SMODEL *mod);
void fe_bedload_init(SMODEL *mod);
void fe_bedload_load(SMODEL *mod);
void fe_bedload_resid(SMODEL *mod);
void fe_bedload_update(SMODEL *mod);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
