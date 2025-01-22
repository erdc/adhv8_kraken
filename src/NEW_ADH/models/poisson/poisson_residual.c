/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  poisson_residual.c This file collections functions responsible for
 *          the 2D Poisson equation with constant RHS, used for testing 
 *          and 2D nonlinear Poisson equation for testing            */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 2D Poisson residual (linear or nonlinear), used for testing.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
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
 *  \details If mod->LINEAR then solves the body integals of the following weak, 
 *  discrete body terms of the 2D Poisson equation: \n
 *  \f$
 *  \int_{\Omega} \nabla u_i \cdot \nabla v_i dx + \int_{\Omega} f dx = 0 \\
 *  f=6 
 *  u_{exact} = 1 + x^2 + 2y^2\\
 *  \f$
 *  If it is nonlinear then solves:
 *  \f$
 *                -\nabla \cdot (q(u) \nabla u) + f = 0\\
 *                q(u) = 1 + u^2\\
 *                f = 10 + 10x + 20y\\
 *                Omega = (0,1) x (0,1)\\
 *                u = u_{exact} on \partial \Omega\\
 *                u_{exact} = 1 + x + 2y \\
 *  \f$
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
    //global_to_local_dbl_cg_2(elem_u, mod->sol, elem2d->nodes, nnodes, PERTURB_U, mod->node_physics_mat, mod->node_physics_mat_id);
    //global_to_local_dbl_cg(elem_u, mod->sol, elem2d->nodes, nnodes, PERTURB_U, mod->dof_map_local, mod->node_physics_mat, mod->node_physics_mat_id);
    global_to_local_dbl_cg(elem_u, mod->sol, elem2d->nodes, nnodes, PERTURB_U, mod->dof_map_local, mod->node_physics_mat);
    //for now let's let f be a constant so we have analytic solution to compare to
    double f[nnodes];
    double q[nnodes];
    int i;
    SNODE nodes[nnodes];
    for (i=0; i<nnodes; i++) {
        snode_copy(&(nodes[i]), mod->grid->node[elem2d->nodes[i]]);
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
    SVECT2D grad_u_nonlinear[nnodes];



    //linear case
    if (mod->LINEAR_PROBLEM == YES){
        //compute diffusion
        for (i=0; i<nnodes; i++) {
            grad_u.x += elem_u[i] * grad_shp[i].x;
            grad_u.y += elem_u[i] * grad_shp[i].y;
        }
        for (i =0;i<nnodes;i++){
            f[i] = 6;
        }
    
        integrate_triangle_gradPhi_dot_vbar(grad_shp, mod->grid->elem2d[ie].djac, 1,  grad_u, elem_rhs);
        //rhs of linear system (will be 0 when Jacobian is computed)
        integrate_triangle_phi_f(mod->grid->elem2d[ie].djac, 1, f, elem_rhs);
    }else if(mod->LINEAR_PROBLEM == NO){
        svect2d_init_array(grad_u_nonlinear, nnodes);
        //compute diffusion
        for (i=0; i<nnodes; i++) {
            grad_u_nonlinear[0].x += elem_u[i] * grad_shp[i].x;
            grad_u_nonlinear[0].y += elem_u[i] * grad_shp[i].y;
            grad_u_nonlinear[1].x += elem_u[i] * grad_shp[i].x;
            grad_u_nonlinear[1].y += elem_u[i] * grad_shp[i].y;
            grad_u_nonlinear[2].x += elem_u[i] * grad_shp[i].x;
            grad_u_nonlinear[2].y += elem_u[i] * grad_shp[i].y;
            q[i] = 1 + elem_u[i]*elem_u[i];
            f[i] = 10.0 + 10.0*nodes[i].x + 20.0*nodes[i].y;
        }


        integrate_triangle_gradPhi_dot_f_v(grad_shp, mod->grid->elem2d[ie].djac, 1, q, grad_u_nonlinear, elem_rhs);
        //rhs of linear system (will be 0 when Jacobian is computed)
        integrate_triangle_phi_f(mod->grid->elem2d[ie].djac, 1, f, elem_rhs);
    }
    return 0;
}
