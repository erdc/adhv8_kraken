/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes the body residual contributions for the SW-3D continuity model.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] elem_rhs      the elemental residual array
 * @param[in]  mod           a pointer to the model struct
 * @param[in]  ie            the elemental id
 * @param[in]  perturbation   the Newton perturbation
 * @param[in]  perturb_node  the node to be perturbed
 * @param[in]  perturb_var   the variable to be perturbed
 * @param[in]  perturb_sign  the direction of Newton perturbation
 * @param[in]  DEBUG         a debug flag
 *
 *  \details Integrates the weak body terms: \n
 *  \f{eqnarray*}{\resid{i}{}{c}=dt \,*\, \sum\limits_{e} \bigg[-\bodyConv{\,3d}{r}{\phiddd{i}}{\vel^{h}} + \bodySUPG{\supg{i}{c}{e}} \bigg] \f}
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_wvel_body_resid(SMODEL *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    
    int i;
    
    int DEBUG_NODE_ID = UNSET_INT;
    
    //if (perturb_node == UNSET_INT) DEBUG_NODE_ID = 2; // only perturb residual call
    
    int DEBUG_LOCAL = OFF;
    int DEBUG_PICKETS = OFF;
    if (DEBUG_NODE_ID != UNSET_INT) {
        for (i=0; i<mod->grid->elem3d[ie].nnodes; i++) {
            if (mod->grid->elem3d[ie].nodes[i] == DEBUG_NODE_ID - 1)  {DEBUG = ON; break;}
            else {DEBUG = OFF; DEBUG_LOCAL = OFF; DEBUG_PICKETS = OFF;}
        }
    }
#ifdef _DEBUG
    if (DEBUG == ON) DEBUG_LOCAL = ON;
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    time_t time1;  time(&time1);
#endif
    
    double t1 = 0.;
    
    SSW_3D *sw = mod->sw->d3;
    double dt = mod->dt;
    double tau_pg = mod->tau_pg;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GRID VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    SGRID *grid = mod->grid;
    SELEM_3D *elem3d = &(grid->elem3d[ie]);
    int nnodes = elem3d->nnodes;
    
    int isElementTetrahedron = TRUE;
    if (nnodes != NDONTET) isElementTetrahedron = FALSE;
    
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // velocity perturbation
    SVECT elem_vel[nnodes];  global_to_local_svect(sw->vel, elem_vel, elem3d->nodes, nnodes);
    if (perturb_var == PERTURB_W) {
        elem_vel[perturb_node].z += perturb_sign * perturbation;
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    double elem_dpl[nnodes];  global_to_local_dbl(sw->displacement, elem_dpl, elem3d->nodes, nnodes);
#ifdef _SEDIMENT
    // cjt :: add bed displacement to total displacement ... be carefull here, no need for db/dt now I think
    double elem_bed_dpl[nnodes];
    global_to_local_dbl(sw->bed_displacement, elem_bed_dpl, elem3d->nodes, nnodes);
    sarray_add_replace_dbl(elem_dpl,elem_bed_dpl,nnodes);
#endif
    
    SVECT elem_nds[nnodes];
    for (i=0; i<nnodes; i++) {
        elem_nds[i].x = grid->node[elem3d->nodes[i]].x;
        elem_nds[i].y = grid->node[elem3d->nodes[i]].y;
        elem_nds[i].z = grid->node[elem3d->nodes[i]].z + elem_dpl[i];
    }
    
    // calculate tetrehedral jacobians if needed
    SVECT grad_shp[nnodes];
    double grad_shp_x[nnodes], grad_shp_y[nnodes], grad_shp_z[nnodes];
    double elem_volume = 0.;
    SVECT elem_grad_u, elem_grad_v, elem_grad_w;
    if (isElementTetrahedron == TRUE) {
        elem_volume = get_tet_linear_djac_gradPhi2(NULL, elem_nds, grad_shp);
        grad_phi_dot_v(grad_shp, elem_vel, &elem_grad_u, &elem_grad_v, &elem_grad_w, nnodes);
        dumpVector(grad_shp, nnodes, grad_shp_x, grad_shp_y, grad_shp_z);
    } else {
        elem_volume = get_triprism_volume(elem_nds);
    }
    
    // dump velocities into doubles for function calls
    double u[nnodes], v[nnodes], w[nnodes];
    dumpVector(elem_vel, nnodes, u, v, w);
    
    // calculate elemental averages
    double elem_avg_strong_continuity = 0.;
    double elem_avg_grad_u_x = 0., elem_avg_grad_v_y = 0., elem_avg_grad_w_z = 0.;
    if (isElementTetrahedron == TRUE) {
        elem_avg_grad_u_x = elem_grad_u.x;
        elem_avg_grad_v_y = elem_grad_v.y;
        elem_avg_grad_w_z = elem_grad_w.z;
    } else {
        elem_avg_grad_u_x = integrate_triPrism_df(elem_nds, 1./elem_volume, u, 1);
        elem_avg_grad_v_y = integrate_triPrism_df(elem_nds, 1./elem_volume, v, 2);
        elem_avg_grad_w_z = integrate_triPrism_df(elem_nds, 1./elem_volume, w, 3);
    }
    elem_avg_strong_continuity = elem_avg_grad_u_x + elem_avg_grad_v_y + elem_avg_grad_w_z;
    
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
#if _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        printf("SW-WVEL BODY RESID :: ie: %d \t dt: %20.10f \t volume: %20.10f",ie,dt,elem_volume);
        if (perturb_var == PERTURB_W) {
            printf("\t perturbing W  || node: %d || perturbation: %20.10e\n",elem3d->nodes[perturb_node],perturb_sign*perturbation);
        }
        selem3d_printScreen(elem3d);
        printScreen_debug_vec("node locations: ",elem_nds, nnodes);
        printScreen_debug2_dbl("elem_displacement", elem_dpl, nnodes, elem3d->nodes);
        printScreen_debug_svect("elem_vel", elem_vel, nnodes, elem3d->nodes);
        
        if (isElementTetrahedron == TRUE) {
            printScreen_debug_svect("grad_shp", grad_shp, nnodes, elem3d->nodes);
        }
        printf("elem_avg_grad_u_x = %20.10e\n",elem_avg_grad_u_x);
        printf("elem_avg_grad_v_y = %20.10e\n",elem_avg_grad_v_y);
        printf("elem_avg_grad_w_z = %20.10e\n",elem_avg_grad_w_z);
        printf("elemental averages: continuity equation: %30.20e\n",elem_avg_strong_continuity);
    }
    if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                  FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    
    // initialize right hand side (to be returned)
    sarray_init_dbl(elem_rhs, nnodes);
    
    // temporary rhs
    double integral[nnodes], integral_X[nnodes], integral_Y[nnodes], integral_Z[nnodes];
    double rhs[nnodes]; sarray_init_dbl(rhs, nnodes);
    
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    CONTINUITY CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the continuity addition to the 3D SW elemental residual. \n
     * \note
     *********************************************************************************************/
    
    sarray_init_dbl(integral,nnodes);
    if (isElementTetrahedron == TRUE) {
        integrate_tetrahedron_gradPhi_dot_v(grad_shp, elem_volume, -dt, elem_vel, integral);
    } else {
        integrate_triPrism_gradPhi_dot_v(elem_nds, -dt, elem_vel, integral);
    }
    for (i=0; i<nnodes; i++) elem_rhs[i] += integral[i];
    
    
    // OLD WAY!
//    double constant = -0.25 * elem_volume * dt;
//    SVECT vel_sum; svect_init(&vel_sum);
//    vel_sum = svect_sum_array(elem_vel, NDONTET);
//    for (i=0; i<NDONTET; i++) {
//        integral[i] += constant * svect_dotp(vel_sum, grad_shp[i]);
//        elem_rhs[i] += constant * svect_dotp(vel_sum, grad_shp[i]);
//    }
    

#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        rhs_1dof("3D WVEL BODY || CONTINUITY: ",nnodes, ie, elem3d->nodes, integral);
        Is_DoubleArray_Inf_or_NaN(integral, nnodes ,__FILE__ ,__LINE__);
    }
#endif
    
    
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                      SUPG CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates the SUPG  addition to the elemental residual. \n
     *  \f{eqnarray*}{ \weakSwCsupgDDD{3d}{e}{i}{c} \f}
     *
     * \note CJT \:: use elementally averaged strong residuals to simplify the integrals
     * \note CJT \:: use elementally averaged velocities to also simplify the integrals
     *********************************************************************************************/
    
    sarray_init_dbl(integral,nnodes);
    if (isElementTetrahedron == TRUE) {
        double le = one_2 * (1./(fabs(grad_shp_z[0]) + fabs(grad_shp_z[1]) + fabs(grad_shp_z[2]) + fabs(grad_shp_z[3])));
        integrate_tetrahedron_fi(elem_volume, -dt * tau_pg * le * elem_avg_strong_continuity, grad_shp_z, integral);
        //if (DEBUG) printf("le: %20.10f \t tau_pg: %20.10f \n",le,tau_pg);
        
    } else {
        
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) { // make sure nodes are number correctly
            assert(fabs(elem_nds[0].x - elem_nds[3].x) < 1e-6); assert(fabs(elem_nds[0].y - elem_nds[3].y) < 1e-6);
            assert(fabs(elem_nds[1].x - elem_nds[4].x) < 1e-6); assert(fabs(elem_nds[1].y - elem_nds[4].y) < 1e-6);
            assert(fabs(elem_nds[2].x - elem_nds[5].x) < 1e-6); assert(fabs(elem_nds[2].y - elem_nds[5].y) < 1e-6);
        }
#endif
        // cjt :: for length, take average of 3 vertical sides
        double le = one_3 * (fabs(elem_nds[0].z - elem_nds[3].z) +
                             fabs(elem_nds[1].z - elem_nds[4].z) +
                             fabs(elem_nds[2].z - elem_nds[5].z));

        // direct integration (uses elementally averaged strong contininuity)
        integrate_triPrism_dphi(elem_nds, -dt * tau_pg * le * elem_avg_strong_continuity, integral, 3);
        
        // use quadrature and a linear continuity
//        int iqp, quad_order = 3; // quadrature order
//        SQUAD_PT *qp = NULL;
//        
//        SQUAD *quad = grid->quad_prism;
//        sarray_init_dbl(integral, nnodes);
//        for (iqp=0; iqp<quad[quad_order].n; iqp++) {
//            qp = &(quad[quad_order].pt[iqp]);
//            
//            // evalute triangular prism djac and cartesian shape function gradients at quadrature point
//            qp->djac = get_triprism_linear_djac_gradPhi(qp->xhat, qp->yhat, qp->zhat, elem_nds, qp->grad_shp);
//            
//            // evaluate continuity equation at quadrature points
//            double qp_grad_u_x = 0.;
//            double qp_grad_v_y = 0.;
//            double qp_grad_w_z = 0.;
//            for (i=0; i<NDONPRISM; i++) {
//                qp_grad_u_x += u[i] * qp->grad_shp[i].x;
//                qp_grad_v_y += v[i] * qp->grad_shp[i].y;
//                qp_grad_w_z += w[i] * qp->grad_shp[i].z;
//            }
//            double qp_Rcont = qp_grad_u_x + qp_grad_v_y + qp_grad_w_z;
//            
//            double t1 = -dt * le * qp->djac * qp->w;; //-dt * tau_pg * le * qp->djac * qp->w; // to match old trunk
//            for (i=0; i<NDONPRISM; i++) {
//                integral[i] += t1 * qp->grad_shp[i].z * qp_Rcont;
//            }
//        }
    }
    
    for (i=0; i<nnodes; i++) {
        mod->sw->d3->elem_rhs_supg_cont[i][ie] = integral[i]; // store for transport
        elem_rhs[i] += integral[i];
    }
    
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("3D WVEL BODY || SUPG: ",nnodes, ie, elem3d->nodes, integral);
            Is_DoubleArray_Inf_or_NaN(integral, nnodes ,__FILE__ ,__LINE__);
        }
#endif

    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                            EXTRA SW DA-CONTINUITY TERM CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Adds any extra 3D SW depth-averged continuity terms to the transport elemental residual. \n
     * note: CJT \:: already multiplied by dt
     *********************************************************************************************/
    //if (perturb_var != PERTURB_W) { // FOR NOW!
        for (i=0; i<nnodes; i++) {
            integral[i] = mod->sw->d3->elem_rhs_supg_dacont[i][ie];
            elem_rhs[i] += integral[i];
        }
        
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("3D WVEL BODY || DAC SUPG: ",nnodes, ie, elem3d->nodes, integral);
            Is_DoubleArray_Inf_or_NaN(rhs, nnodes ,__FILE__ ,__LINE__);
        }
#endif
    //}

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        rhs_1dof("3D WVEL TOTAL: ",nnodes, ie, elem3d->nodes, elem_rhs);
        Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    }
    
    time_t time2;  time(&time2);
    TIME_IN_WVEL_BODY_RESID += difftime(time2,time1);
#endif
    
    return;
}
