/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  heat_residual.c This file collections functions responsible for
 *          the 2D Heat equation with constant RHS, used for testing
 *          https://jsdokken.com/dolfinx-tutorial/chapter2/heat_code.html             */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static double alpha = 3.0;
static double beta = 1.2;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 2D Heat residual with constant RHS on an element, used for testing.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] mod (SMODEL_SUPER*) - pointer to SMODEL_SUPER struct
 * @param[in,out] elem_rhs (double*) - array of doubles that will store elemental residual
 * @param[in] ie (int) - the elemental id
 * @param[in] pertubation (double) - the Newton pertubation
 * @param[in] perturb_node (int) - the node to be pertubed
 * @param[in] perturb_var (int) - the variable code to be perturbed
 * @param[in] perturb_sign (int) - the direction of Newton perturbation (-1 or +1)
 * @param[in] DEBUG (int) - a debug flag
 * \returns integer code
 *  \details Solves the body integals of the following weak, discrete body terms of the 2D Poisson equation: \n
 *  \f$
 *  \int_{\Omega} \frac{\partial u_i}{\partial t} + \nabla u_i \cdot \nabla v_i dx - \int_{\Omega} f dx = 0 \\
 *  f= \beta - 2 -2\alpha \\
 *  u_{exact} = 1 + x^2 + \alpha y^2 + \beta t\\
 *  \f$
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int heat_residual(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    SELEM_2D *elem2d = &(mod->grid->elem2d[ie]); // be careful not to shallow copy here
    //call a simple integral
    int nnodes = mod->grid->elem2d[ie].nnodes;
    //get the solution on this element, for now just use U
    double elem_u[nnodes], elem_u_old[nnodes];
    double dt = *(mod->dt);
    //this map only works for CG, want to generalize to DG in future
    //global_to_local_dbl_cg_2(elem_u, mod->sol, elem2d->nodes, nnodes, PERTURB_U, mod->node_physics_mat, mod->node_physics_mat_id);
    //global_to_local_dbl_cg(elem_u, mod->sol, elem2d->nodes, nnodes, PERTURB_U, mod->dof_map_local, mod->node_physics_mat, mod->node_physics_mat_id);
    //global_to_local_dbl_cg(elem_u_old, mod->sol_old, elem2d->nodes, nnodes, PERTURB_U, mod->dof_map_local, mod->node_physics_mat, mod->node_physics_mat_id);
    global_to_local_dbl_cg(elem_u, mod->sol, elem2d->nodes, nnodes, PERTURB_U, mod->dof_map_local, mod->node_physics_mat);
    global_to_local_dbl_cg(elem_u_old, mod->sol_old, elem2d->nodes, nnodes, PERTURB_U, mod->dof_map_local, mod->node_physics_mat);
    //for now let's let f be a constant so we have analytic solution to compare to
    double f[nnodes], dhdt[nnodes];
    int i;
    //perturb solution variable only, let's say we are using U
    if (perturb_var == PERTURB_U){
        elem_u[perturb_node] += perturb_sign * perturbation;
    }
    for (i =0;i<nnodes;i++){
        f[i] = beta - 2.0 - 2.0*alpha;
        dhdt[i] = (elem_u[i] - elem_u_old[i])/dt;
    }
    //this would be mass matrix
    //integrate_triangle_phi_f(mod->grid->elem2d[ie].djac, 1, elem_u, elem_rhs);
    //compute the Laplacian
    SVECT2D grad_u;
    svect2d_init(&grad_u);
    SVECT2D *grad_shp = elem2d->grad_shp;
    //compute diffusion
    for (i=0; i<nnodes; i++) {
        grad_u.x += elem_u[i] * grad_shp[i].x;
        grad_u.y += elem_u[i] * grad_shp[i].y;
    }
    integrate_triangle_gradPhi_dot_vbar(grad_shp, mod->grid->elem2d[ie].djac, 1.0,  grad_u, elem_rhs);
    //time derivative, implicit Euler
    integrate_triangle_phi_f(mod->grid->elem2d[ie].djac, 1.0, dhdt, elem_rhs);  
    //rhs of linear system (will be 0 when Jacobian is computed)
    integrate_triangle_phi_f(mod->grid->elem2d[ie].djac, -1.0, f, elem_rhs);

    return 0;
}
