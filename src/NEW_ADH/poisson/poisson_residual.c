/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_sw2_2d_body_resid.c This file collections functions responsible for
 *          the 2D shallow water body contributions to the elemental residual.              */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//#include "global_header.h"
#include "adh.h"
// prototypes
//void fe_sw2_temporal(int ie, SELEM_2D *elem2d, int nnodes, SVECT *elem_nds, double djac, double drying_lower_limit, double *elem_head, SVECT2D *elem_vel, double wd_factor, double dt_factor, DOF_3 *elem_rhs, char *string, int DEBUG, int DEBUG_LOCAL, int wd_flag);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 2D shallow water elemental residual.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] elem_rhs      the 2D elemental residual array
 * @param[in]  mod           a pointer to the model struct
 * @param[in]  ie            the elemental id
 * @param[in]  pertubation   the Newton pertubation
 * @param[in]  perturb_node  the node to be perturbed
 * @param[in]  perturb_var   the variable to be perturbed
 * @param[in]  perturb_sign  the direction of Newton perturbation
 * @param[in]  DEBUG         a debug flag
 * @param[in]  PRESSURE_FLAG  a flag turning pressure contributions on/off
 *
 *  \details Solves the body integals of the following weak, discrete body terms of the 2D shallow water equation: \n
 *  \f{eqnarray*}{
 *  \weakSWDAcont{e}{i}{h} \\
 *  \weakSWMxDD{e}{i}{h} \\
 *  \weakSWMyDD{e}{i}{h}
 *  \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int poisson_residual(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    SELEM_2D *elem2d = &(mod->grid->elem2d[ie]); // be careful not to shallow copy here
    //call a simple integral
    int nnodes = mod->grid->elem2d[ie].nnodes;
    //get the solution on this element, for now just use U
    double elem_u[nnodes];
    //this map only works for CG, want to generalize to DG in future
    global_to_local_dbl_cg(mod->sol, elem_u, elem2d->nodes, nnodes, PERTURB_U, mod->node_physics_mat, mod->node_physics_mat_id);

    //for now let's let f be a constant so we have analytic solution to compare to
    double f[nnodes];
    int i;
    for (i =0;i<nnodes;i++){
        f[i] = -2;
    }
    //perturb solution variable only, let's say we are using U
    if (perturb_var == PERTURB_U){
        elem_u[perturb_node] += perturb_sign * perturbation;
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
    
    integrate_triangle_gradPhi_dot_vbar(grad_shp, mod->grid->elem2d[ie].djac, -1,  grad_u, elem_rhs);

    //rhs of linear system (will be 0 when Jacobian is computed)
    integrate_triangle_phi_f(mod->grid->elem2d[ie].djac, -1, f, elem_rhs);

    return 0;
}