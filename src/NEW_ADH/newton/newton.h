#ifndef _H_NEWTON_
#define _H_NEWTON_


//from newton_tools.c
void initialize_system(SMODEL_SUPER *sm);
void initialize_dirichlet_bc(SMODEL_SUPER *sm);
void update_function(SMODEL_SUPER *sm);
void increment_function(SMODEL_SUPER *sm);
void get_residual_norms(double *resid_max_norm, double *resid_l2_norm, double *inc_max_norm,
         int *imax_dof, int *iinc_dof, int *include_dof,
         int my_ndofs, int ndofs, int macro_ndofs, double *residual, double *dsol, int *bc_mask);
int fe_newton(SMODEL_SUPER* sm);

#endif
