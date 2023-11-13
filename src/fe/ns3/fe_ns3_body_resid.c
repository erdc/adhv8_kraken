/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_ns3_body_resid.c This file collections functions responsible for
 *         the 3D NS p,u,v,w equation contributions to the elemental residual.              */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

// debug print precision
static int printFieldWidth = 15;
static int printPrecision  =  9;


static int DEBUG_LOCAL = OFF;
static int DEBUG_PICKETS = OFF;
static int DEBUG_NODE_ID = UNSET_INT;
static int DEBUG_EXIT_ON_NAN = ON;

static SGRID *grid = NULL;
static SQUAD *quad = NULL;
static int nnodes = 0, nnodes_quad = 0, isElementTetrahedron = TRUE, imat = UNSET_INT;
static SELEM_3D *elem3d = NULL;
static SNODE elem_nodes[NDONPRISM];
static double g = 0., dt = 0., alpha = 0.;
static SQUAD *quad;

static int INCLUDE_PRESSURE = ON;
static double elem_prs[NDONPRISMQUAD];
static double elem_dpl[NDONPRISM],elem_dpl_old[NDONPRISM], elem_dpl_older[NDONPRISM];
static double grad_phi_x[NDONPRISM], grad_phi_y[NDONPRISM], grad_phi_z[NDONPRISM];
static double elem_u[NDONPRISM], elem_v[NDONPRISM], elem_w[NDONPRISM];
static double elem_u_old[NDONPRISM], elem_v_old[NDONPRISM], elem_w_old[NDONPRISM];
static double elem_u_older[NDONPRISM], elem_v_older[NDONPRISM], elem_w_older[NDONPRISM];
static double elem_ur[NDONPRISM], elem_vr[NDONPRISM], elem_wr[NDONPRISM];
static double elem_rho[NDONPRISM];

static SVECT elem_vel[NDONPRISM], elem_vel_old[NDONPRISM], elem_vel_older[NDONPRISM];
static SVECT elem_rel_vel[NDONPRISM], elem_grid_vel[NDONPRISM];
static SVECT elem_nds[NDONPRISMQUAD], elem_nds_old[NDONPRISMQUAD], elem_nds_older[NDONPRISMQUAD];
static SVECT grad_phi[NDONPRISM], grad_phi_old[NDONPRISM], grad_phi_older[NDONPRISM];

static double elem_avg_strong_continuity = 0.;
static double elem_volume = 0., elem_volume_old = 0., elem_volume_older = 0., elem_djac = 0., elem_djac_old = 0., elem_djac_older = 0.;
static double elem_avg_density = 0., elem_avg_prs = 0.;
static double elem_avg_ur = 0., elem_avg_vr = 0., elem_avg_wr = 0.;
static double elem_avg_u = 0., elem_avg_v = 0., elem_avg_w = 0.;
static double elem_avg_u_old = 0., elem_avg_v_old = 0., elem_avg_w_old = 0.;
static double elem_avg_u_older = 0., elem_avg_v_older = 0., elem_avg_w_older = 0.;
static SVECT elem_avg_vel, elem_avg_vel_old, elem_avg_vel_older, elem_avg_vel_rel, elem_grad_u, elem_grad_v, elem_grad_w, elem_grad_p;
static SVECT elem_avg_grad_u, elem_avg_grad_v, elem_avg_grad_w, elem_avg_grad_p, elem_avg_strong_mom;
static SVECT diffusive_flux_x, diffusive_flux_y, diffusive_flux_z;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// file prototypes
int fe_ns3_body_test_for_NAN_or_INF();
void fe_ns3_body_debug(double perturbation, int perturb_node, int perturb_var, int perturb_sign);
void fe_ns3_temporal(int ie, int nnodes, SELEM_3D *elem3d, SVECT *elem_nds, double djac, double *u, double *v, double *w,
                     double dt_factor, DOF_4 *elem_rhs, char *string, double perturbation, int perturb_node, int perturb_var,
                     int perturb_sign, int DEBUG);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 3D Navier Stokes h,u,v,w equation contributions to the elemental residual.
 *  \author    Gary Brown
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] elem_rhs      the 2D elemental residual array
 * @param[in]  mod           a pointer to the model struct
 * @param[in]  ie            the elemental id
 * @param[in]  perturbation   the Newton perturbation
 * @param[in]  perturb_node  the node to be perturbed
 * @param[in]  perturb_var   the variable to be perturbed
 * @param[in]  perturb_sign  the direction of Newton perturbation
 * @param[in]  DEBUG         a debug flag
 *
 *
 * \details Solves the following weak, discrete body terms of the 3D SW equation: \n
 *  \f{eqnarray*}{ \weakSWDaContReducedKinematicBody{i} \\
 *                 \weakSWMxDDDBody{e}{i}{h}            \\
 *                 \weakSWMyDDDBody{e}{i}{h}            \f}
 * \n where \f$ \diffTensor{mx,y}{}{\diffFluxSW}{} \f$ are diffusive fluxes.
 *
 * \note cjt \:: all residuals are multiplied by dt / density
 * \note cjt \::
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_ns3_body_resid(SMODEL *mod, DOF_4 *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    
    int i;
    int print_debug_info = OFF; // OFF/ON
    INCLUDE_PRESSURE = ON;
    
    // DO NOT TAMPER WITH +++++++++++++++++++++++++++
    // for time-dependent debugging
    if (mod->t_prev > DEBUG_TIME) DEBUG_LOCAL = ON;
    
    // for node-dependent debugging
    for (i=0; i<mod->grid->elem3d[ie].nnodes; i++) {
        if (mod->grid->elem3d[ie].nodes[i] == DEBUG_NODE_ID) {
            DEBUG_LOCAL = ON;
        }
    }
    //+++++++++++++++++++++++++++++++++++++++++++++++
    
#ifdef _DEBUG
    if (DEBUG == ON) DEBUG_LOCAL = ON;
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    time_t time1;  time(&time1);
#endif
    
    // aliases
    SNS_3D *ns = mod->ns->d3;
    g = mod->gravity;
    dt = mod->dt;
    alpha = mod->tau_temporal;
    
    double tau_pg = mod->tau_pg;
    STR_VALUE *str_values = mod->str_values;
    double dtTerm1 = (1. + alpha / 2.);
    double dtTerm2 = (-1.) * (1. + alpha);
    double dtTerm3 = (alpha / 2.) ;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GRID VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    grid = mod->grid;
    quad = grid->quad_prism;
    elem3d = &(grid->elem3d[ie]);
    nnodes = elem3d->nnodes;
    nnodes_quad = elem3d->nnodes_quad;
    imat = elem3d->mat;
    SMAT_NS mat = *mod->mat[imat].ns;
    
    isElementTetrahedron = TRUE;
    if (nnodes != NDONTET) isElementTetrahedron = FALSE;
    for(i=0; i<nnodes;i++) {
        snode_copy(&(elem_nodes[i]), grid->node[elem3d->nodes[i]]);
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // pressure perturbation
    global_to_local_dbl(ns->prs, elem_prs, elem3d->nodes, nnodes);
    if (perturb_var == PERTURB_P) {elem_prs[perturb_node] += perturb_sign * perturbation;}
    
    // displacement (for moving grid nodes)
    global_to_local_dbl(ns->displacement, elem_dpl, elem3d->nodes, nnodes);
    if (perturb_var == PERTURB_D) {elem_dpl[perturb_node] += perturb_sign * perturbation;}
    
    // velocity perturbation
    global_to_local_svect(ns->vel, elem_vel, elem3d->nodes, nnodes);
    if (perturb_var == PERTURB_U) {
        elem_vel[perturb_node].x += perturb_sign * perturbation;
        //INCLUDE_PRESSURE = OFF;
    } else if (perturb_var == PERTURB_V) {
        elem_vel[perturb_node].y += perturb_sign * perturbation;
        //INCLUDE_PRESSURE = OFF;
    } else if (perturb_var == PERTURB_W) {
        elem_vel[perturb_node].z += perturb_sign * perturbation;
        //INCLUDE_PRESSURE = OFF;
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // velocities
    global_to_local_svect(ns->old_vel, elem_vel_old, elem3d->nodes, nnodes);
    global_to_local_svect(ns->older_vel, elem_vel_older, elem3d->nodes, nnodes);
    
    // moving grid variables
    sarray_init_dbl(elem_dpl_old,nnodes);
    sarray_init_dbl(elem_dpl_older,nnodes);
    svect_init_array(elem_grid_vel, nnodes);
    svect_copy_array(elem_rel_vel,elem_vel, nnodes);
    if (mod->flag.MG == ON) {
        int gnodeID = UNSET_INT;
        for (i=0; i<nnodes; i++) {
            gnodeID = elem3d->nodes[i];
            if (mod->grid->node[gnodeID].bflag == 0) {
                //
                elem_dpl_old[i]   = ns->old_displacement[gnodeID];
                elem_dpl_older[i] = ns->older_displacement[gnodeID];
                elem_grid_vel[i].z = elem_vel[i].z;
                elem_rel_vel[i].z = 0.0;
                
                // CJT :: add w perturbation here
                if (perturb_var == PERTURB_W) elem_dpl[i] += (perturb_sign * perturbation * dt);
            }
            
        }
    } else {
        global_to_local_dbl(ns->old_displacement, elem_dpl_old, elem3d->nodes, nnodes);
        global_to_local_dbl(ns->older_displacement, elem_dpl_older, elem3d->nodes, nnodes);
    }
    
    // moving grid displacements only on surface and if turned on (add w perturbation)
    if (mod->flag.MG == ON) {
        if (INCLUDE_PRESSURE == ON) {
            int gnodeID = UNSET_INT;
            
            for (i=0; i<nnodes; i++) {
                gnodeID = elem3d->nodes[i];
                if (grid->node[perturb_node].bflag == 0 && grid->node[gnodeID].bflag != 0) {
                    if (perturb_var == PERTURB_D) {
                        if (fabs(elem_nodes[i].x - elem_nodes[perturb_node].x) < SMALL && fabs(elem_nodes[i].y - elem_nodes[perturb_node].y) < SMALL) {
                            elem_prs[i] += mod->density * mod->gravity * (perturb_sign * perturbation);
                        }
                    }
                    if (perturb_var == PERTURB_W){
                        if (fabs(elem_nodes[i].x - elem_nodes[perturb_node].x) < SMALL && fabs(elem_nodes[i].y - elem_nodes[perturb_node].y) < SMALL) {
                            elem_prs[i] += mod->density * mod->gravity * (perturb_sign * perturbation * dt);
                        }
                    }
                }
            }
        }
    }
    
    
    // nodal vectors positions w/ displacements if MG
    elem_get_midpt_locations2(elem_nodes, elem_dpl,       elem_nds,       nnodes, nnodes_quad, elem3d->edges);
    elem_get_midpt_locations2(elem_nodes, elem_dpl_old,   elem_nds_old,   nnodes, nnodes_quad, elem3d->edges);
    elem_get_midpt_locations2(elem_nodes, elem_dpl_older, elem_nds_older, nnodes, nnodes_quad, elem3d->edges);
    
    // density
    global_to_local_dbl(ns->density, elem_rho, elem3d->nodes, nnodes);
    
    // calculate tetrehedral jacobians if needed
    elem_djac = 0., elem_djac_old = 0., elem_djac_older = 0.;
    elem_volume = 0., elem_volume_old = 0., elem_volume_older = 0.;
    svect_init_array(grad_phi,NDONPRISM); svect_init_array(grad_phi_old,NDONPRISM); svect_init_array(grad_phi_older,NDONPRISM);
    svect_init(&elem_grad_u); svect_init(&elem_grad_v); svect_init(&elem_grad_w); svect_init(&elem_grad_p);
    if (isElementTetrahedron == TRUE) {
        elem_djac = get_tet_linear_djac_gradPhi2(NULL, elem_nds, grad_phi);
        elem_djac_old = get_tet_linear_djac_gradPhi2(NULL, elem_nds_old, grad_phi_old);
        elem_djac_older = get_tet_linear_djac_gradPhi2(NULL, elem_nds_older, grad_phi_older);
        
        grad_phi_dot_v(grad_phi, elem_vel, &elem_grad_u, &elem_grad_v, &elem_grad_w, nnodes);
        grad_phi_f(grad_phi, elem_prs, &elem_grad_p, nnodes); // should be a quadratic gradient
        
        elem_volume = elem_djac;
        elem_volume_old = elem_djac_old;
        elem_volume_older = elem_djac_older;
    } else {
        elem_volume = get_triprism_volume(elem_nds);
        elem_volume_old = get_triprism_volume(elem_nds_old);
        elem_volume_older = get_triprism_volume(elem_nds_older);
    }
    
    // dump velocities into doubles for function calls
    dumpVector(elem_vel,       nnodes, elem_u,       elem_v,       elem_w);
    dumpVector(elem_vel_old,   nnodes, elem_u_old,   elem_v_old,   elem_w_old);
    dumpVector(elem_vel_older, nnodes, elem_u_older, elem_v_older, elem_w_older);
    dumpVector(elem_rel_vel,   nnodes, elem_ur, elem_vr, elem_wr);
    if (isElementTetrahedron == TRUE) {
        dumpVector(grad_phi, nnodes, grad_phi_x, grad_phi_y, grad_phi_z);
    }
    
    // calculate elemental averages
    if (isElementTetrahedron == TRUE) {
        
        // elementally averaged gradients :: note :: since gradients on tets are constant
        svect_copy(&elem_avg_grad_u,elem_grad_u);
        svect_copy(&elem_avg_grad_v,elem_grad_v);
        svect_copy(&elem_avg_grad_w,elem_grad_w);
        svect_copy(&elem_avg_grad_p,elem_grad_p);
        
        // elementally averaged quadratic pressure and densities
        //double sum1 = elem_prs[0] + elem_prs[1] + elem_prs[2] + elem_prs[3];
        //double sum2 = elem_prs[4] + elem_prs[5] + elem_prs[6] + elem_prs[7] + elem_prs[8] + elem_prs[9];
        //elem_avg_prs = one_20 * (-sum1 + 4. * sum2);
        
        elem_avg_prs =     one_4 * (elem_prs[0] + elem_prs[1] + elem_prs[2] + elem_prs[3]); // linear pressure average
        elem_avg_density = one_4 * (elem_rho[0] + elem_rho[1] + elem_rho[2] + elem_rho[3]);
        
        // elementally averaged velocities
        elem_avg_vel       = svect_average_array(elem_vel, nnodes);
        elem_avg_vel_old   = svect_average_array(elem_vel_old, nnodes);
        elem_avg_vel_older = svect_average_array(elem_vel_older, nnodes);
        elem_avg_vel_rel   = svect_average_array(elem_rel_vel, nnodes);
        
        // elementally averaged momentum 2nd order temporal terms
        double elem_avg_dudt = elem_avg_density*(dtTerm1 * elem_avg_vel.x + dtTerm2 * elem_avg_vel_old.x + dtTerm3 * elem_avg_vel_older.x)/dt;
        double elem_avg_dvdt = elem_avg_density*(dtTerm1 * elem_avg_vel.y + dtTerm2 * elem_avg_vel_old.y + dtTerm3 * elem_avg_vel_older.y)/dt;
        double elem_avg_dwdt = elem_avg_density*(dtTerm1 * elem_avg_vel.z + dtTerm2 * elem_avg_vel_old.z + dtTerm3 * elem_avg_vel_older.z)/dt;
        
        // elementally averaged strong 3D continuity and momentum equations (for SUPG)
        elem_avg_strong_continuity = elem_grad_u.x + elem_grad_v.y + elem_grad_w.z;
        elem_avg_strong_mom.x = elem_avg_dudt + elem_avg_density * svect_dotp(elem_avg_vel_rel, elem_grad_u);
        elem_avg_strong_mom.y = elem_avg_dvdt + elem_avg_density * svect_dotp(elem_avg_vel_rel, elem_grad_v);
        elem_avg_strong_mom.z = elem_avg_dwdt + elem_avg_density * svect_dotp(elem_avg_vel_rel, elem_grad_w);
        if (INCLUDE_PRESSURE == ON) {
            elem_avg_strong_mom.x += elem_avg_grad_p.x;
            elem_avg_strong_mom.y += elem_avg_grad_p.y;
            elem_avg_strong_mom.z += elem_avg_grad_p.z;
        }
        //if (perturb_var == PERTURB_NONE) {
            elem_avg_strong_mom.z += elem_avg_density * g; // assume density isn't perturbed or a function of a perturbed variable
        //}
        
        
        // scale this shit
        elem_avg_dudt = (dtTerm1 * elem_avg_vel.x + dtTerm2 * elem_avg_vel_old.x + dtTerm3 * elem_avg_vel_older.x)/dt;
        elem_avg_dvdt = (dtTerm1 * elem_avg_vel.y + dtTerm2 * elem_avg_vel_old.y + dtTerm3 * elem_avg_vel_older.y)/dt;
        elem_avg_dwdt = (dtTerm1 * elem_avg_vel.z + dtTerm2 * elem_avg_vel_old.z + dtTerm3 * elem_avg_vel_older.z)/dt;
        
        elem_avg_strong_continuity = (elem_grad_u.x + elem_grad_v.y + elem_grad_w.z) * (1./mod->density);
        elem_avg_strong_mom.x = elem_avg_dudt + svect_dotp(elem_avg_vel_rel, elem_grad_u);
        elem_avg_strong_mom.y = elem_avg_dvdt + svect_dotp(elem_avg_vel_rel, elem_grad_v);
        elem_avg_strong_mom.z = elem_avg_dwdt + svect_dotp(elem_avg_vel_rel, elem_grad_w);
        if (INCLUDE_PRESSURE == ON) {
            elem_avg_strong_mom.x += elem_avg_grad_p.x * (1./mod->density);
            elem_avg_strong_mom.y += elem_avg_grad_p.y * (1./mod->density);
            elem_avg_strong_mom.z += elem_avg_grad_p.z * (1./mod->density);
        }
        //if (perturb_var == PERTURB_NONE) {
            elem_avg_strong_mom.z += g; // assume density isn't perturbed or a function of a perturbed variable
        //}
        
    } else {
        
        elem_avg_density = integrate_triPrism_f(elem_nds, 1./elem_volume, elem_rho);
        elem_avg_prs =     integrate_triPrism_f(elem_nds, 1./elem_volume, elem_prs); // linear pressure average
        elem_avg_u =       integrate_triPrism_f(elem_nds, 1./elem_volume, elem_u);
        elem_avg_v =       integrate_triPrism_f(elem_nds, 1./elem_volume, elem_v);
        elem_avg_w =       integrate_triPrism_f(elem_nds, 1./elem_volume, elem_w);
        elem_avg_ur =      elem_avg_u; // no lateral grid movement
        elem_avg_vr =      elem_avg_v; // no lateral grid movement
        elem_avg_wr =      integrate_triPrism_f(elem_nds, 1./elem_volume, elem_wr);
        elem_avg_u_old =   integrate_triPrism_f(elem_nds_old, 1./elem_volume_old, elem_u_old);
        elem_avg_v_old =   integrate_triPrism_f(elem_nds_old, 1./elem_volume_old, elem_v_old);
        elem_avg_w_old =   integrate_triPrism_f(elem_nds_old, 1./elem_volume_old, elem_w_old);
        elem_avg_u_older = integrate_triPrism_f(elem_nds_older, 1./elem_volume_older, elem_u_older);
        elem_avg_v_older = integrate_triPrism_f(elem_nds_older, 1./elem_volume_older, elem_v_older);
        elem_avg_w_older = integrate_triPrism_f(elem_nds_older, 1./elem_volume_older, elem_w_older);
        elem_avg_vel.x = elem_avg_u;
        elem_avg_vel.y = elem_avg_v;
        elem_avg_vel.z = elem_avg_w;
        elem_avg_vel_rel.x = elem_avg_ur;
        elem_avg_vel_rel.y = elem_avg_vr;
        elem_avg_vel_rel.z = elem_avg_wr;
        elem_avg_grad_u = integrate_triPrism_df_full(elem_nds, 1./elem_volume, elem_u);
        elem_avg_grad_v = integrate_triPrism_df_full(elem_nds, 1./elem_volume, elem_v);
        elem_avg_grad_w = integrate_triPrism_df_full(elem_nds, 1./elem_volume, elem_w);
        
        // elementally averaged momentum 2nd order temporal terms
        double elem_avg_dudt = (dtTerm1 * elem_avg_u + dtTerm2 * elem_avg_u_old + dtTerm3 * elem_avg_u_older)/dt;
        double elem_avg_dvdt = (dtTerm1 * elem_avg_v + dtTerm2 * elem_avg_v_old + dtTerm3 * elem_avg_v_older)/dt;
        double elem_avg_dwdt = (dtTerm1 * elem_avg_w + dtTerm2 * elem_avg_w_old + dtTerm3 * elem_avg_w_older)/dt;
        
        // elementally averaged strong 3D continuity and momentum equations (for SUPG)
        elem_avg_strong_continuity = elem_avg_grad_u.x + elem_avg_grad_v.y + elem_avg_grad_w.z;
        elem_avg_strong_mom.x = elem_avg_dudt + svect_dotp(elem_avg_vel_rel, elem_avg_grad_u);
        elem_avg_strong_mom.y = elem_avg_dvdt + svect_dotp(elem_avg_vel_rel, elem_avg_grad_v);
        elem_avg_strong_mom.z = elem_avg_dwdt + svect_dotp(elem_avg_vel_rel, elem_avg_grad_w) + elem_avg_density * g;
        //if (PRESSURE_FLAG == ON) {
        elem_avg_strong_mom.x += integrate_triPrism_df(elem_nds, 1./elem_volume, elem_prs, 1);
        elem_avg_strong_mom.y += integrate_triPrism_df(elem_nds, 1./elem_volume, elem_prs, 2);
        elem_avg_strong_mom.z += integrate_triPrism_df(elem_nds, 1./elem_volume, elem_prs, 3);
        //}
    }
    
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    assert(ie == elem3d->id);
    assert(alpha > -1E-6 && alpha < 1.0);
    assert(imat >= 0);
    assert(perturb_var == PERTURB_NONE || perturb_var == PERTURB_P || perturb_var == PERTURB_U ||
           perturb_var == PERTURB_V || perturb_var == PERTURB_W || perturb_var == PERTURB_D);
    
    if (isElementTetrahedron == TRUE) {
        if (elem_djac < SMALL) {print_debug_info = ON; DEBUG_EXIT_ON_NAN = ON;}
        if (elem_djac_old < SMALL) {print_debug_info = ON; DEBUG_EXIT_ON_NAN = ON;}
        if (elem_djac_older < SMALL) {print_debug_info = ON; DEBUG_EXIT_ON_NAN = ON;}
    } else {
        if (elem_volume < SMALL) {print_debug_info = ON; DEBUG_EXIT_ON_NAN = ON;}
        if (elem_volume_old < SMALL) {print_debug_info = ON; DEBUG_EXIT_ON_NAN = ON;}
        if (elem_volume_older < SMALL) {print_debug_info = ON; DEBUG_EXIT_ON_NAN = ON;}
    }
    
    // if there is a NAN or INF, print arrays and possibly exit
    print_debug_info = fe_ns3_body_test_for_NAN_or_INF();
    if (print_debug_info != OFF) {
        tl_check_all_pickets(__FILE__,__LINE__);
        fe_ns3_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
    
    if (mod->flag.MG == OFF) {
        for (i=0; i<nnodes; i++) {
            
            if (fabs(elem_dpl[i]) > SMALL || fabs(elem_dpl_old[i]) > SMALL || fabs(elem_dpl_old[i]) > SMALL) {
                printf("ERROR :: elem_dpl[i] :: %20.10f \t elem_dpl_old[i] :: %20.10f\n",elem_dpl[i],elem_dpl_old[i]);
                tl_error("elem_dpl should be 0 for a fixed grid!\n");
            }
            if (fabs(elem_rel_vel[i].x - elem_vel[i].x) > SMALL ||
                fabs(elem_rel_vel[i].y - elem_vel[i].y) > SMALL ||
                fabs(elem_rel_vel[i].z - elem_vel[i].z) > SMALL) {
                printf("ERROR :: elem_vel[%d] :: {%20.10f, %20.10f, %20.10f} \n",i,elem_vel[i].x,elem_vel[i].y,elem_vel[i].z);
                printf("ERROR :: elem_rel_vel[%d] :: {%20.10f, %20.10f, %20.10f} \n",i,elem_rel_vel[i].x,elem_rel_vel[i].y,elem_rel_vel[i].z);
                tl_error("elem_rel_vel should equal elem_vel for a fixed grid!\n");
            }
            if (fabs(elem_grid_vel[i].x) > SMALL ||
                fabs(elem_grid_vel[i].y) > SMALL ||
                fabs(elem_grid_vel[i].z) > SMALL) {
                printf("ERROR :: elem_grid_vel[%d] :: {%20.10f, %20.10f, %20.10f} \n",i,elem_grid_vel[i].x,elem_grid_vel[i].y,elem_grid_vel[i].z);
                tl_error("elem_grid_vel should be 0 for a fixed grid!\n");
            }
            if (fabs(elem_djac - elem_djac_old) > SMALL || fabs(elem_djac - elem_djac_older) > SMALL) {
                printf("elem_djac: %20.10f \t elem_djac_old: %20.10f \t elem_djac_older: %20.10f\n",
                       elem_djac, elem_djac_old, elem_djac_older);
                tl_error("elem_djac should be = elem_djac_old = elem_djac_older for a fixed grid!\n");
            }
        }
    }
    
    if (DEBUG_LOCAL == ON) fe_ns3_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    
#endif
    
//    printf("tet: %d :: %10.5f %10.5f %10.5f \t %10.5f %10.5f %10.5f \t %10.5f %10.5f %10.5f \t %10.5f %10.5f %10.5f\n",ie,
//           grid->node[grid->elem3d[ie].nodes[0]].x,grid->node[grid->elem3d[ie].nodes[0]].y, grid->node[grid->elem3d[ie].nodes[0]].z,
//           grid->node[grid->elem3d[ie].nodes[1]].x,grid->node[grid->elem3d[ie].nodes[1]].y, grid->node[grid->elem3d[ie].nodes[1]].z,
//           grid->node[grid->elem3d[ie].nodes[2]].x,grid->node[grid->elem3d[ie].nodes[2]].y, grid->node[grid->elem3d[ie].nodes[2]].z,
//           grid->node[grid->elem3d[ie].nodes[3]].x,grid->node[grid->elem3d[ie].nodes[3]].y, grid->node[grid->elem3d[ie].nodes[3]].z);
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                  FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    
    // initialize right hand side (to be returned)
    dof4_init_array(elem_rhs, nnodes);
    
    // temporary rhs
    double integral[nnodes], integral_X[nnodes], integral_Y[nnodes], integral_Z[nnodes];
    DOF_4 rhs[nnodes]; dof4_init_array(rhs, nnodes);
    
    // zero SUPG da-cont terms for wvel addition
    //for (i=0; i<nnodes; i++) mod->sw->d3->elem_rhs_supg_dacont[i][ie] = 0.;
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                     CONTINUITY CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the continuity elemental residual. \n
     * \note
     *********************************************************************************************/
    
    sarray_init_dbl(integral,nnodes);
    if (isElementTetrahedron == TRUE) {
        // use linear velocities
        integrate_tetrahedron_gradPhi_dot_v(grad_phi, elem_djac, 1., elem_vel, integral);
        // use constant velocities
        //integrate_tetrahedron_gradPhi_dot_vcon(grad_phi, elem_djac, 1., elem_avg_vel, integral);
    } else {
        // user linear velocities
        integrate_triPrism_gradPhi_dot_v(elem_nds, 1., elem_vel, integral);
    }
    for (i=0; i<nnodes; i++) {
        elem_rhs[i].c_eq -= integral[i]; // * (1./mod->density);
    }
    
#ifdef _DEBUG
    dof4_init_array(rhs,nnodes);
    for (i=0; i<nnodes; i++) {
        rhs[i].c_eq -= integral[i];
    }
    
    if (DEBUG_LOCAL == ON) {
        rhs_4dof("3D NS || CONTNUITY: ",nnodes, ie, elem3d->nodes, rhs, printFieldWidth, printPrecision);
    }
    if (dof4_debug_continue(elem_rhs,nnodes) == 1) {
        tl_check_all_pickets(__FILE__,__LINE__);
        rhs_4dof("3D NS || CONTINUITY: ", nnodes, elem3d->id, elem3d->nodes, rhs, printFieldWidth, printPrecision);
        printf("\n");
        fe_ns3_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   MOMENTUM TEMPORAL CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the temporal addition to the 3D NS momentum elemental residual. \n
     * \note
     *********************************************************************************************/
    
    //    // ++ t(n+1) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //    fe_ns3_temporal(ie, nnodes, elem3d, elem_nds, elem_djac, elem_u, elem_v, elem_w,
    //                    elem_avg_density * dtTerm1 / dt, elem_rhs, "3D NS || TEMPORAL T(N+1) ", perturbation, perturb_node, perturb_var, perturb_sign, DEBUG);
    //
    //    // ++ t(n) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //    fe_ns3_temporal(ie, nnodes, elem3d, elem_nds_old, elem_djac_old, elem_u_old, elem_v_old, elem_w_old,
    //                    elem_avg_density * dtTerm2 / dt, elem_rhs, "3D NS || TEMPORAL T(N) ", perturbation, perturb_node, perturb_var, perturb_sign, DEBUG);
    //
    //    // ++ t(n-1) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //    fe_ns3_temporal(ie, nnodes, elem3d, elem_nds_older, elem_djac_older, elem_u_older, elem_v_older, elem_w_older,
    //                    elem_avg_density * dtTerm3 / dt, elem_rhs, "3D NS || TEMPORAL T(N-1) ", perturbation, perturb_node, perturb_var, perturb_sign, DEBUG);
    
    
    // ++ t(n+1) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_ns3_temporal(ie, nnodes, elem3d, elem_nds, elem_djac, elem_u, elem_v, elem_w,
                    dtTerm1 / dt, elem_rhs, "3D NS || TEMPORAL T(N+1) ", perturbation, perturb_node, perturb_var, perturb_sign, DEBUG);
    
    // ++ t(n) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_ns3_temporal(ie, nnodes, elem3d, elem_nds_old, elem_djac_old, elem_u_old, elem_v_old, elem_w_old,
                    dtTerm2 / dt, elem_rhs, "3D NS || TEMPORAL T(N) ", perturbation, perturb_node, perturb_var, perturb_sign, DEBUG);
    
    // ++ t(n-1) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_ns3_temporal(ie, nnodes, elem3d, elem_nds_older, elem_djac_older, elem_u_older, elem_v_older, elem_w_older,
                    dtTerm3 / dt, elem_rhs, "3D NS || TEMPORAL T(N-1) ", perturbation, perturb_node, perturb_var, perturb_sign, DEBUG);
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    ADVECTION CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the advection addition to the 3D NS momentum elemental residual. \n
     * \note
     *********************************************************************************************/
    
#ifdef _DEBUG
    if (debug.no_advection == OFF) {
#endif
        sarray_init_dbl(integral_X,nnodes);
        sarray_init_dbl(integral_Y,nnodes);
        sarray_init_dbl(integral_Z,nnodes);
        if (isElementTetrahedron) {
            integrate_tetrahedron_gradPhi_dot_f_v(grad_phi, elem_djac, 1., elem_u, elem_rel_vel, integral_X);
            integrate_tetrahedron_gradPhi_dot_f_v(grad_phi, elem_djac, 1., elem_v, elem_rel_vel, integral_Y);
            integrate_tetrahedron_gradPhi_dot_f_v(grad_phi, elem_djac, 1., elem_w, elem_rel_vel, integral_Z);
        } else {
            integrate_triPrism_gradPhi_dot_f_v(elem_nds, 1., elem_u, elem_rel_vel, integral_X);
            integrate_triPrism_gradPhi_dot_f_v(elem_nds, 1., elem_v, elem_rel_vel, integral_Y);
            integrate_triPrism_gradPhi_dot_f_v(elem_nds, 1., elem_w, elem_rel_vel, integral_Z);
        }
        for (i=0; i<nnodes; i++) {
            //elem_rhs[i].x_eq -= elem_avg_density * integral_X[i];
            //elem_rhs[i].y_eq -= elem_avg_density * integral_Y[i];
            //elem_rhs[i].z_eq -= elem_avg_density * integral_Z[i];
            
            
            elem_rhs[i].x_eq -= integral_X[i];
            elem_rhs[i].y_eq -= integral_Y[i];
            elem_rhs[i].z_eq -= integral_Z[i];
        }
        
        
         //grid divergence (questionable)
//                sarray_init_dbl(integral_X,nnodes);
//                sarray_init_dbl(integral_Y,nnodes);
//                sarray_init_dbl(integral_Z,nnodes);
//                double grid_divergence = 0.;
//                for (i=0; i<nnodes; i++) grid_divergence += elem_rel_vel[i].x * grad_phi_x[i] + elem_rel_vel[i].y * grad_phi_y[i] + elem_rel_vel[i].z * grad_phi_z[i];
//                integrate_tetrahedron_phi_f(elem_djac, grid_divergence, elem_u, integral_X);
//                integrate_tetrahedron_phi_f(elem_djac, grid_divergence, elem_v, integral_Y);
//                integrate_tetrahedron_phi_f(elem_djac, grid_divergence, elem_w, integral_Z);
//                for (i=0; i<nnodes; i++) {
//                    elem_rhs[i].x_eq -= elem_avg_density * integral_X[i];
//                    elem_rhs[i].y_eq -= elem_avg_density * integral_Y[i];
//                    elem_rhs[i].z_eq -= elem_avg_density * integral_Z[i];
//                }
        
        
#ifdef _DEBUG
    }
    
    dof4_init_array(rhs,nnodes);
    for (i=0; i<nnodes; i++) {
        rhs[i].x_eq -= elem_avg_density * integral_X[i];
        rhs[i].y_eq -= elem_avg_density * integral_Y[i];
        rhs[i].z_eq -= elem_avg_density * integral_Z[i];
    }
    
    if (DEBUG_LOCAL == ON) {
        rhs_4dof("3D NS || ADVECTION: ",nnodes, ie, elem3d->nodes, rhs, printFieldWidth, printPrecision);
    }
    if (dof4_debug_continue(elem_rhs,nnodes) == 1) {
        tl_check_all_pickets(__FILE__,__LINE__);
        rhs_4dof("3D NS || ADVECTIONS: ", nnodes, elem3d->id, elem3d->nodes, rhs, printFieldWidth, printPrecision);
        printf("\n");
        fe_ns3_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    PRESSURE CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the pressure addition to the 3D NS momentum elemental residual. \n
     * \note
     *********************************************************************************************/
    
    if (INCLUDE_PRESSURE == ON) {
        
        sarray_init_dbl(integral_X,nnodes);
        sarray_init_dbl(integral_Y,nnodes);
        sarray_init_dbl(integral_Z,nnodes);
        if (isElementTetrahedron) {
            // linear pressure
            //integrate_tetrahedron_fi(elem_djac, elem_avg_prs, grad_phi_x, integral_X);
            //integrate_tetrahedron_fi(elem_djac, elem_avg_prs, grad_phi_y, integral_Y);
            //integrate_tetrahedron_fi(elem_djac, elem_avg_prs, grad_phi_z, integral_Z);
            
            // linear pressure2
            for (i=0; i<nnodes; i++) {
                //integral_X[i] = integrate_tetrahedron_f_g(elem_djac, 1., elem_prs, grad_phi_x);
                //integral_Y[i] = integrate_tetrahedron_f_g(elem_djac, 1., elem_prs, grad_phi_y);
                //integral_Z[i] = integrate_tetrahedron_f_g(elem_djac, 1., elem_prs, grad_phi_z);
                
                integral_X[i] = integrate_tetrahedron_f(elem_djac, grad_phi_x[i], elem_prs);
                integral_Y[i] = integrate_tetrahedron_f(elem_djac, grad_phi_y[i], elem_prs);
                integral_Z[i] = integrate_tetrahedron_f(elem_djac, grad_phi_z[i], elem_prs);
            }
        } else {
            // constant pressure || note ::  is not stable for some reason ...
            //integrate_triPrism_dphi(elem_nds, elem_avg_prs, integral_X, 1);
            //integrate_triPrism_dphi(elem_nds, elem_avg_prs, integral_Y, 2);
            
            // linear pressure || note :: works, but seems to cap out early on NTL convergence
            integrate_triPrism_dphi_f(elem_nds, 1., elem_prs, integral_X, 1);
            integrate_triPrism_dphi_f(elem_nds, 1., elem_prs, integral_Y, 2);
            integrate_triPrism_dphi_f(elem_nds, 1., elem_prs, integral_Z, 3);
            
            // linear pressure using quadrature ||note ::  slow, but the convergence rates seem better than analytic, probably due to precision
            //integrate_triPrism_dphi_f_quadrature(elem_nds, grid->quad_prism, 1., elem_prs, integral_X, integral_Y, 3);
            //integrate_triPrism_dphi_f_quadrature(elem_nds, grid->quad_prism, 1., elem_prs, integral_X, integral_Y, 3);
        }
        for (i=0; i<nnodes; i++) {
            elem_rhs[i].x_eq -= integral_X[i] * (1./mod->density);
            elem_rhs[i].y_eq -= integral_Y[i] * (1./mod->density);
            elem_rhs[i].z_eq -= integral_Z[i] * (1./mod->density);
            
            //if (elem_nodes[i].id==5) printf("    BODY PRESSURE :: integral_X: %20.10f integral_Y: %20.10f integral_Z: %20.10f\n",
            //                                -integral_X[i] * (dt/elem_avg_density),
            //                                -integral_Y[i] * (dt/elem_avg_density),
            //                                -integral_Z[i] * (dt/elem_avg_density));
        }
#ifdef _DEBUG
        dof4_init_array(rhs,nnodes);
        for (i=0; i<nnodes; i++) {
            rhs[i].x_eq -= integral_X[i];
            rhs[i].y_eq -= integral_Y[i];
            rhs[i].z_eq -= integral_Z[i];
        }
        
        if (DEBUG_LOCAL == ON) {
            rhs_4dof("3D NS || PRESSURE: ",nnodes, ie, elem3d->nodes, rhs, printFieldWidth, printPrecision);
        }
        if (dof4_debug_continue(elem_rhs,nnodes) == 1) {
            tl_check_all_pickets(__FILE__,__LINE__);
            rhs_4dof("3D NS || PRESSURE: ", nnodes, elem3d->id, elem3d->nodes, rhs, printFieldWidth, printPrecision);
            printf("\n");
            fe_ns3_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
            if (DEBUG_EXIT_ON_NAN) exit(-1);
        }
#endif
        
        //if (elem_nodes[0].id==1 || elem_nodes[1].id==1 || elem_nodes[2].id==1 || elem_nodes[3].id==1)
        //printf("DPDX :: gradP: %20.10f %20.10f %20.10f\n",elem_grad_p.x,elem_grad_p.y,elem_grad_p.z);
    }
    
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                      GRAVITY CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the gravity addition to the 3D NS momentum elemental residual. \n
     * \note
     *********************************************************************************************/
    
    //if (perturb_var == PERTURB_NONE) {
        
        //sarray_init_dbl(integral_Z,nnodes);
        if (isElementTetrahedron) {
            //integrate_tetrahedron_phi_f(elem_djac, g, elem_rho, integral_Z);
            for (i=0; i<nnodes; i++) {
                elem_rhs[i].z_eq += (0.25 * elem_djac * g);   //integral_Z[i]* (dt);
                //elem_rhs[i].z_eq += integral_Z[i]* (dt);
            }
        } else {
            
        }
#ifdef _DEBUG
        dof4_init_array(rhs,nnodes);
        for (i=0; i<nnodes; i++) {
            rhs[i].z_eq += integral_Z[i];
            //if (elem_nodes[i].id==5) printf(" GRAVITY PRESSURE :: integral_X: %20.10f integral_Y: %20.10f integral_Z: %20.10f\n"
            //                                ,0.,0.,(0.25 * elem_djac * g * dt));
        }
        
        if (DEBUG_LOCAL == ON) {
            rhs_4dof("3D NS || GRAVITY: ",nnodes, ie, elem3d->nodes, rhs, printFieldWidth, printPrecision);
        }
        if (dof4_debug_continue(elem_rhs,nnodes) == 1) {
            tl_check_all_pickets(__FILE__,__LINE__);
            rhs_4dof("3D NS || GRAVITY: ", nnodes, elem3d->id, elem3d->nodes, rhs, printFieldWidth, printPrecision);
            printf("\n");
            fe_ns3_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
            if (DEBUG_EXIT_ON_NAN) exit(-1);
        }
#endif
    //}
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    DIFFUSION CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the diffusion addition to the 3D NS momentum elemental residual. \n
     * \note CJT :: OLD ADH multiplies this term by elem_avg_density
     *********************************************************************************************/
    
    STENSOR3D ev;
    double multiply = 1.;
    ev.xx = multiply * (mat.ev.xx + mod->viscosity);
    ev.xy = multiply * (mat.ev.xy + mod->viscosity);
    ev.xz = multiply * (mat.ev.xz + mod->viscosity);
    ev.yx = multiply * (mat.ev.xy + mod->viscosity);
    ev.yy = multiply * (mat.ev.yy + mod->viscosity);
    ev.yz = multiply * (mat.ev.yz + mod->viscosity);
    ev.zx = multiply * (mat.ev.xz + mod->viscosity);
    ev.zy = multiply * (mat.ev.yz + mod->viscosity);
    ev.zz = multiply * (mat.ev.zz + mod->viscosity);
    
    sarray_init_dbl(integral_X,nnodes);
    sarray_init_dbl(integral_Y,nnodes);
    sarray_init_dbl(integral_Z,nnodes);
    if (isElementTetrahedron) {
        diffusive_flux_x.x = ev.xx * elem_grad_u.x;  diffusive_flux_x.y = ev.xy * elem_grad_u.y;  diffusive_flux_x.z = ev.xz * elem_grad_u.z;
        diffusive_flux_y.x = ev.yx * elem_grad_v.x;  diffusive_flux_y.y = ev.yy * elem_grad_v.y;  diffusive_flux_y.z = ev.yz * elem_grad_v.z;
        diffusive_flux_z.x = ev.zx * elem_grad_w.x;  diffusive_flux_z.y = ev.zy * elem_grad_w.y;  diffusive_flux_z.z = ev.zz * elem_grad_w.z;
        
        integrate_tetrahedron_gradPhi_dot_vcon(grad_phi, elem_djac, 1., diffusive_flux_x, integral_X);
        integrate_tetrahedron_gradPhi_dot_vcon(grad_phi, elem_djac, 1., diffusive_flux_y, integral_Y);
        integrate_tetrahedron_gradPhi_dot_vcon(grad_phi, elem_djac, 1., diffusive_flux_z, integral_Z);
    }
    else {
        
    }
    for (i=0; i<nnodes; i++) {
        elem_rhs[i].x_eq += integral_X[i]* (1./mod->density);
        elem_rhs[i].y_eq += integral_Y[i]* (1./mod->density);
        elem_rhs[i].z_eq += integral_Z[i]* (1./mod->density);
    }
#ifdef _DEBUG
    dof4_init_array(rhs,nnodes);
    for (i=0; i<nnodes; i++) {
        rhs[i].x_eq += integral_X[i];
        rhs[i].y_eq += integral_Y[i];
        rhs[i].z_eq += integral_Z[i];
    }
    
    if (DEBUG_LOCAL == ON) {
        rhs_4dof("3D NS || DIFFUSION: ",nnodes, ie, elem3d->nodes, rhs, printFieldWidth, printPrecision);
    }
    if (dof4_debug_continue(elem_rhs,nnodes) == 1) {
        tl_check_all_pickets(__FILE__,__LINE__);
        rhs_4dof("3D NS || DIFFUSION: ", nnodes, elem3d->id, elem3d->nodes, rhs, printFieldWidth, printPrecision);
        printf("\n");
        fe_ns3_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
#endif
    
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                      SUPG CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates the SUPG  addition to the elemental residual. \n
     *  \f{eqnarray*}{
     *  \weakSwDACsupgDDD{2d}{e}{i}{c}{\velh} \\
     *  \weakSwMXYsupgDDD{3d}{e}{i}{mx}{\velrh}{g\,\depth{h}}{x} \\
     *  \weakSwMXYsupgDDD{3d}{e}{i}{my}{\velrh}{g\,\depth{h}}{y}
     *  \f}
     *
     * \note CJT \:: use elementally averaged strong residuals to simplify the integrals
     * \note CJT \:: use elementally averaged velocities to also simplify the integrals
     * \note CJT \:: found out here that, for prisms, the strong continuity equation cannot be elem averaged for best stability
     *********************************************************************************************/
    
    DOF_4 rhs_mom_supg[nnodes]; dof4_init_array(rhs_mom_supg, nnodes);
    DOF_4 rhs_con_supg[nnodes]; dof4_init_array(rhs_con_supg, nnodes);
    double tau_SUPG = 0., tau_PSPG = 0.;
    if (isElementTetrahedron == TRUE) {
        
        double upwind = 0., constant = 0.;
        double elem_vmag_rel = svect_mag(elem_avg_vel_rel); // relative velocity magnitude
        double ev_max = stensor3d_max(ev);                  // maximum diffusion coefficient
        
        // get element length
        double cell_length = get_element_length(nnodes, NULL, elem_avg_vel.x, elem_avg_vel.y, elem_avg_vel.z,
                                                grad_phi_x, grad_phi_y, grad_phi_z, elem_djac, 3, UNSET_INT);
        
        cell_length = pow(5. * elem_djac, 0.33333333333333333333);;
        
        // SUPG Tau taken from Tezduyar 1992
        tau_SUPG  = 4. * pow(elem_vmag_rel,2) / pow(cell_length,2);
        tau_SUPG += 4. / pow(dt,2);
        tau_SUPG += 16. * pow(ev_max,2) / pow(cell_length, 4.);
        tau_SUPG  = 1. / sqrt(tau_SUPG + 2. * NOT_QUITE_SMALL);
        tau_SUPG *= mod->tau_pg;
        tau_PSPG = tau_SUPG;
        
        //printf("cell length: %20.10f Charlie's cell length: %20.10f \t ev_max: %20.10f\n",cell_length,ev_max,pow(5. * elem_djac, 0.33333333333333333333));
        
        // use CHARLIE's dimen con
//        double reynolds_factor = elem_vmag_rel * cell_length / (ev_max + NOT_QUITE_SMALL);
//        if(reynolds_factor < 3) {
//            reynolds_factor /= 3.;
//        } else{
//            reynolds_factor = 1.;
//        }
//        tau_PSPG = reynolds_factor * elem_vmag_rel * cell_length;
        
        
        // PSPG contribution (since velocity and pressure basis are the same)
        constant = elem_djac * tau_PSPG / elem_avg_density;
        for (i=0; i<nnodes; i++) {
            rhs_con_supg[i].c_eq = constant *  (elem_avg_strong_mom.x * grad_phi_x[i] +
                                                elem_avg_strong_mom.y * grad_phi_y[i] +
                                                elem_avg_strong_mom.z * grad_phi_z[i]);
        }
        
        // SUPG contribution (for advection wiggles)
        constant = elem_djac * tau_SUPG;
        for (i=0; i<nnodes; i++) {
            upwind = svect_dotp(elem_avg_vel_rel, grad_phi[i]);
            rhs_mom_supg[i].x_eq = constant * (elem_avg_strong_mom.x * upwind);
            rhs_mom_supg[i].y_eq = constant * (elem_avg_strong_mom.y * upwind);
            rhs_mom_supg[i].z_eq = constant * (elem_avg_strong_mom.z * upwind);
        }
        
    }  else {
        
    }
    
    for (i=0; i<nnodes; i++) {
        //mod->sw->d3->elem_rhs_supg_dacont[i][ie] = rhs_dacont_supg[i].c_eq + rhs_mom_supg[i].c_eq;
        elem_rhs[i].c_eq += rhs_con_supg[i].c_eq;
        elem_rhs[i].x_eq += rhs_mom_supg[i].x_eq;
        elem_rhs[i].y_eq += rhs_mom_supg[i].y_eq;
        elem_rhs[i].z_eq += rhs_mom_supg[i].z_eq;
    }
    
#ifdef _DEBUG
    dof4_init_array(rhs,nnodes);
    for (i=0; i<nnodes; i++) {
        rhs[i].c_eq += rhs_con_supg[i].c_eq;
        rhs[i].x_eq += rhs_mom_supg[i].x_eq;
        rhs[i].y_eq += rhs_mom_supg[i].y_eq;
        rhs[i].z_eq += rhs_mom_supg[i].z_eq;
    }
    if (DEBUG_LOCAL == ON) {
        rhs_4dof("3D NS || SUPG/PSPG: ",nnodes, ie, elem3d->nodes, rhs, printFieldWidth, printPrecision);
    }
    if (dof4_debug_continue(elem_rhs,nnodes) == 1) {
        tl_check_all_pickets(__FILE__,__LINE__);
        rhs_4dof("3D NS || SUPG/PSPG: ", nnodes, elem3d->id, elem3d->nodes, rhs, printFieldWidth, printPrecision);
        printf("\n");
        fe_ns3_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
#endif
    
    
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *   SETS DIRICHLET BOUNDARY CONDITIONS AND EXCHANGES MOMENTUM RESIDUAL IF NB VEL IS GIVEN
     *-------------------------------------------------------------------------------------------
     * Assigned DB conditions. \n
     *********************************************************************************************/
    
    //    for (i=0; i<nnodes; i++) {
    //        j = elem3d->nodes[i] * mod->nsys; // the global DOF
    //        if (mod->bc_mask[j] == YES) {
    //            elem_rhs[i].x_eq = 0.;
    //        } else if (mod->bc_mask[j] == 2) {
    //            elem_rhs[i].y_eq = elem_rhs[i].x_eq * sw->tanvec[elem3d->nodes[i]].x + elem_rhs[i].y_eq * sw->tanvec[elem3d->nodes[i]].y;
    //            elem_rhs[i].x_eq = 0.;
    //        }
    //        if (mod->bc_mask[j + 1] == YES) {
    //            elem_rhs[i].y_eq = 0.;
    //        }
    //        if (elem_surface_nodeID[i] == elem3d->nodes[i]) {
    //            if (mod->bc_mask[j + 2] == YES) {
    //                elem_rhs[i].c_eq = 0.;
    //            }
    //        }
    //    }
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        rhs_4dof("3D NS BODY || TOTAL RHS: ",nnodes, ie, elem3d->nodes, elem_rhs, printFieldWidth, printPrecision);
        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    }
    dof4_debug(elem_rhs, nnodes, __FILE__, __LINE__);
    
    time_t time2;  time(&time2);
    TIME_IN_NS3_BODY_RESID += difftime(time2,time1);
    //printf("\n");
#endif
    
    //exit(-1);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 3D NS body temporal contributions to the shallow water equations.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_ns3_temporal(int ie, int nnodes, SELEM_3D *elem3d, SVECT *elem_nds, double djac, double *u, double *v, double *w,
                     double dt_factor, DOF_4 *elem_rhs, char *string, double perturbation, int perturb_node, int perturb_var,
                     int perturb_sign, int DEBUG) {
    int i;
    double integral_X[nnodes];  sarray_init_dbl(integral_X,nnodes);
    double integral_Y[nnodes];  sarray_init_dbl(integral_Y,nnodes);
    double integral_Z[nnodes];  sarray_init_dbl(integral_Z,nnodes);
    
    if (nnodes == NDONTET) {
        integrate_tetrahedron_phi_f(djac, dt_factor, u, integral_X);
        integrate_tetrahedron_phi_f(djac, dt_factor, v, integral_Y);
        integrate_tetrahedron_phi_f(djac, dt_factor, w, integral_Z);
    }
    else {
        integrate_triPrism_phi_f(elem_nds, dt_factor, u, integral_X);
        integrate_triPrism_phi_f(elem_nds, dt_factor, v, integral_Y);
        integrate_triPrism_phi_f(elem_nds, dt_factor, w, integral_Z);
    }
    for (i=0; i<nnodes; i++) {
        elem_rhs[i].x_eq += integral_X[i];
        elem_rhs[i].y_eq += integral_Y[i];
        elem_rhs[i].z_eq += integral_Z[i];
    }
#ifdef _DEBUG
    DOF_4 rhs[nnodes]; dof4_init_array(rhs,nnodes);
    for (i=0; i<nnodes; i++) {
        rhs[i].x_eq += integral_X[i];
        rhs[i].y_eq += integral_Y[i];
        rhs[i].z_eq += integral_Z[i];
    }
    
    if (DEBUG_LOCAL == ON) {
        rhs_4dof("3D NS || TEMPORAL: ",nnodes, ie, elem3d->nodes, rhs, printFieldWidth, printPrecision);
    }
    if (dof4_debug_continue(elem_rhs,nnodes) == 1) {
        tl_check_all_pickets(__FILE__,__LINE__);
        rhs_4dof("3D NS || TEMPORAL: ", nnodes, elem3d->id, elem3d->nodes, rhs, printFieldWidth, printPrecision);
        fe_ns3_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Checks 3D transport body variables for NANS or INFS.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int fe_ns3_body_test_for_NAN_or_INF() {
    
    int print_debug_info = OFF;
    
    
    //    int i, STOP_CODE = 0;
    //    for (i=0; i<nnodes; i++) {
    //        if (elem_prs[i] < -0.01) STOP_CODE = 1;
    //    }
    //    if (STOP_CODE == 1) {
    //        printf("\n");
    //        printScreen_debug2_dbl("elem_prs", elem_prs, nnodes, elem3d->nodes);
    //        tl_error(">> Pressures should NOT be negative.");
    //    }
    
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_prs, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_prs has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_dpl, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_dpl has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_dpl_old, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_dpl_old has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_dpl_older, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_dpl_older has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_vel, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_vel has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_vel_old, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_vel_old has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_vel_older, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_vel_older has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_grid_vel, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_grid_vel has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_rel_vel, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_rel_vel has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_nds, nnodes_quad ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_nds has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_nds_old, nnodes_quad ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_nds_old has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_nds_older, nnodes_quad ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_nds_older has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_prs,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_prs has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_density,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_density has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_u,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_u has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_v,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_v has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_w,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_w has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_u_old,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_u_old has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_v_old,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_v_old has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_w_old,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_w_old has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_u_older,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_u_older has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_v_older,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_v_older has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_w_older,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_w_older has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_Vector_Inf_or_NaN_noExit(elem_avg_grad_p,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_grad_p has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_Vector_Inf_or_NaN_noExit(elem_avg_grad_u,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_grad_u has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Vector_Inf_or_NaN_noExit(elem_avg_grad_v,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_grad_v has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Vector_Inf_or_NaN_noExit(elem_avg_grad_w,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_grad_w has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_strong_continuity,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_strong_continuity has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Vector_Inf_or_NaN_noExit(elem_avg_strong_mom,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_strong_mom has either a NAN or INF, printing debug info \n\n");
    
    //print_debug_info = Is_Double_Inf_or_NaN_noExit(molecular_diffusion,__FILE__ ,__LINE__);
    //if (print_debug_info != NO) printf("molecular_diffusion has either a NAN or INF, printing debug info \n\n");
    //print_debug_info = Is_Double_Inf_or_NaN_noExit(turbulent_diffusion,__FILE__ ,__LINE__);
    //if (print_debug_info != NO) printf("turbulent_diffusion has either a NAN or INF, printing debug info \n\n");
    //print_debug_info = Is_Double_Inf_or_NaN_noExit(vertical_diffusion,__FILE__ ,__LINE__);
    //if (print_debug_info != NO) printf("vertical_diffusion has either a NAN or INF, printing debug info \n\n");
    //print_debug_info = Is_Double_Inf_or_NaN_noExit(total_diffusion,__FILE__ ,__LINE__);
    //if (print_debug_info != NO) printf("total_diffusion has either a NAN or INF, printing debug info \n\n");
    
    if (print_debug_info == NO) print_debug_info = OFF;
    
    return print_debug_info;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Prints to screen DEBUG info for the 3D transport body residual.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_ns3_body_debug(double perturbation, int perturb_node, int perturb_var, int perturb_sign) {
    
    printf("3D NS BODY ELEM RESID :: ie: %d \t dt: %-*.*f \t volume: %-*.*f \t old volume: %-*.*f \t older volume: %-*.*f",
           elem3d->id,
           printFieldWidth,printPrecision,dt,
           printFieldWidth,printPrecision,elem_volume,
           printFieldWidth,printPrecision,elem_volume_old,
           printFieldWidth,printPrecision,elem_volume_older);
    
    if (perturb_var == PERTURB_P) {
        printf("\t perturbing P || node: %d || perturbation: %-*.*e\n",elem3d->nodes[perturb_node],printFieldWidth,printPrecision,perturb_sign*perturbation);
    } else if (perturb_var == PERTURB_U) {
        printf("\t perturbing U || node: %d || perturbation: %-*.*e\n",elem3d->nodes[perturb_node],printFieldWidth,printPrecision,perturb_sign*perturbation);
    } else if (perturb_var == PERTURB_V) {
        printf("\t perturbing V || node: %d || perturbation: %-*.*e\n",elem3d->nodes[perturb_node],printFieldWidth,printPrecision,perturb_sign*perturbation);
    } else if (perturb_var == PERTURB_W) {
        printf("\t perturbing W || node: %d || perturbation: %-*.*e\n",elem3d->nodes[perturb_node],printFieldWidth,printPrecision,perturb_sign*perturbation);
    } else {
        printf("\n");
    }
    
    //selem3d_printScreen(elem3d);
    printScreen_debug_vec("node locations: ",elem_nds, nnodes);
    printScreen_debug_vec("node locations old: ",elem_nds_old, nnodes);
    printScreen_debug_vec("node locations older: ",elem_nds_older, nnodes);
    printScreen_debug_vec("grad_phi: ",grad_phi, nnodes);
    printScreen_debug_vec("grad_phi_old: ",grad_phi_old, nnodes);
    printScreen_debug_vec("grad_phi_older: ",grad_phi_older, nnodes);
    printScreen_debug2_dbl("elem_prs", elem_prs, nnodes, elem3d->nodes);
    printScreen_debug2_dbl("elem_displacement", elem_dpl, nnodes, elem3d->nodes);
    printScreen_debug2_dbl("elem_displacement_old", elem_dpl_old, nnodes, elem3d->nodes);
    printScreen_debug2_dbl("elem_displacement_older", elem_dpl_older, nnodes, elem3d->nodes);
    printScreen_debug_svect("elem_vel", elem_vel, nnodes, elem3d->nodes);
    printScreen_debug_svect("elem_vel_old", elem_vel_old, nnodes, elem3d->nodes);
    printScreen_debug_svect("elem_vel_older", elem_vel_older, nnodes, elem3d->nodes);
    printScreen_debug_svect("elem_grid_vel", elem_grid_vel, nnodes, elem3d->nodes);
    printScreen_debug_svect("elem_rel_vel", elem_rel_vel, nnodes, elem3d->nodes);
    printScreen_debug2_dbl("elem_rho", elem_rho, nnodes, elem3d->nodes);
    printf("elem_avg_density: %-*.*f\n",printFieldWidth,printPrecision,elem_avg_density);
    printf("elem_avg_prs: %-*.*f\n",printFieldWidth,printPrecision,elem_avg_prs);
    printf("elem_avg_ur: %-*.*f \t elem_avg_vr: %-*.*f \t elem_avg_wr: %-*.*f\n",
           printFieldWidth,printPrecision,elem_avg_ur,
           printFieldWidth,printPrecision,elem_avg_vr,
           printFieldWidth,printPrecision,elem_avg_wr);
    printf("elem_avg_u: %-*.*f \t elem_avg_v: %-*.*f \t elem_avg_w: %-*.*f\n",
           printFieldWidth,printPrecision,elem_avg_u,
           printFieldWidth,printPrecision,elem_avg_v,
           printFieldWidth,printPrecision,elem_avg_w);
    printf("elem_avg_u_old: %-*.*f   \t elem_avg_v_old:   %-*.*f \t elem_avg_w_old:   %-*.*f\n",
           printFieldWidth,printPrecision,elem_avg_u_old,
           printFieldWidth,printPrecision,elem_avg_v_old,
           printFieldWidth,printPrecision,elem_avg_w_old);
    printf("elem_avg_u_older: %-*.*f \t elem_avg_v_older: %-*.*f \t elem_avg_w_older: %-*.*f\n",
           printFieldWidth,printPrecision,elem_avg_u_older,
           printFieldWidth,printPrecision,elem_avg_v_older,
           printFieldWidth,printPrecision,elem_avg_w_older);
    printf("elem_avg_grad_u: ddx = %-*.*f \t ddy = %-*.*f \t ddz = %-*.*f\n",
           printFieldWidth,printPrecision,elem_avg_grad_u.x,
           printFieldWidth,printPrecision,elem_avg_grad_u.y,
           printFieldWidth,printPrecision,elem_avg_grad_u.z);
    printf("elem_avg_grad_v: ddx = %-*.*f \t ddy = %-*.*f \t ddz = %-*.*f\n",
           printFieldWidth,printPrecision,elem_avg_grad_v.x,
           printFieldWidth,printPrecision,elem_avg_grad_v.y,
           printFieldWidth,printPrecision,elem_avg_grad_v.z);
    printf("elem_avg_grad_w: ddx = %-*.*f \t ddy = %-*.*f \t ddz = %-*.*f\n",
           printFieldWidth,printPrecision,elem_avg_grad_w.x,
           printFieldWidth,printPrecision,elem_avg_grad_w.y,
           printFieldWidth,printPrecision,elem_avg_grad_w.z);
    printf("elem_avg_grad_p: ddx = %-*.*f \t ddy = %-*.*f \t ddz = %-*.*f\n",
           printFieldWidth,printPrecision,elem_avg_grad_p.x,
           printFieldWidth,printPrecision,elem_avg_grad_p.y,
           printFieldWidth,printPrecision,elem_avg_grad_p.z);
    printf("SUPG :: elem_avg_strong_continuity: %-*.*f\n",printFieldWidth,printPrecision,elem_avg_strong_continuity);
    printf("SUPG :: elem_avg_strong_momentum: x :: %-*.*f\n",printFieldWidth,printPrecision,elem_avg_strong_mom.x);
    printf("SUPG :: elem_avg_strong_momentum: y :: %-*.*f\n",printFieldWidth,printPrecision,elem_avg_strong_mom.y);
    printf("SUPG :: elem_avg_strong_momentum: z :: %-*.*f\n",printFieldWidth,printPrecision,elem_avg_strong_mom.z);
    
    printf("diffusive_flux_x = %*.*e \t %*.*e \t %*.*e\n",
           printFieldWidth,printPrecision,diffusive_flux_x.x,
           printFieldWidth,printPrecision,diffusive_flux_x.y,
           printFieldWidth,printPrecision,diffusive_flux_x.z);
    printf("diffusive_flux_y = %*.*e \t %*.*e \t %*.*e\n",
           printFieldWidth,printPrecision,diffusive_flux_y.x,
           printFieldWidth,printPrecision,diffusive_flux_y.y,
           printFieldWidth,printPrecision,diffusive_flux_y.z);
    printf("diffusive_flux_z = %*.*e \t %*.*e \t %*.*e\n",
           printFieldWidth,printPrecision,diffusive_flux_z.x,
           printFieldWidth,printPrecision,diffusive_flux_z.y,
           printFieldWidth,printPrecision,diffusive_flux_z.z);
    
    printf("\n");
    
}



















