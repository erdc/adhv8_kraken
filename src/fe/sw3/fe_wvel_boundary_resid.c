/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes the boundary contributions for the 3D continuity equation.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] elem_rhs        the elemental residual array
 * @param[out] elem_rhs_noBed  the elemental residual array w/o bed surface contributions
 * @param[in]  mod             a pointer to the model struct
 * @param[in]  ie              the elemental id
 * @param[in]  perturbation     the Newton perturbation
 * @param[in]  perturb_node    the node to be perturbed
 * @param[in]  perturb_var     the variable to be perturbed
 * @param[in]  perturb_sign    the direction of Newton perturbation
 * @param[in]  PRESSURE        a flag to include pressure terms (TRUE/FALSE)
 * @param[in]  DEBUG           a debug flag
 *
 * \details Integrates the weak, discrete outflow boundary terms: \n
 * \f{eqnarray*}{
 *    R_{\,i,boundary}^{\,e} &=& dt * \bigg(
 *    \bcConv{sw}{e}{\phiddd{i}}{\velh} \,+
 *    \bcConv{b}{e}{\phiddd{i}}{\velh}  \,+
 *    \bcKinematic{\elev{}}{e}{\phiddd{i}}{\elev{h}}{\src{\elev{}}}
 *    \bigg)
 * \f}
 *
 * \note if done discetely consistent, the bed residual should be 0
 * \note if mask = 2, then the original equations are replaced with
 *         x_eq = 3D continuity
 *         y_eq = tangential momemtum equation
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_wvel_boundary_resid(SMODEL *mod, double *elem_rhs, double *elem_rhs_noBed, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    
    int DEBUG_LOCAL = OFF; 
    int DEBUG_PICKETS = OFF;
#ifdef _DEBUG
    if (DEBUG == ON) DEBUG_LOCAL = ON;
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    time_t time1;  time(&time1);
#endif
    
    int i;
    
    // aliases
    SSW_3D *sw = mod->sw->d3;
    double dt = mod->dt;
    double g = mod->gravity;
    STR_VALUE *str_values = mod->str_values;
    SFLAGS flags = mod->flag;
    double density = mod->density;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GRID VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    int isElementTriangle = FALSE;
    int isSidewall = FALSE;
    int isSurface = FALSE;
    int isBed = FALSE;
    
    SGRID *grid = mod->grid;
    SELEM_2D *elem2d = &(grid->elem2d[ie]);
    int nnodes = elem2d->nnodes;
    int string = elem2d->string;
    int hydro_bc_flag = str_values[string].flow.bc_flag;
    if (str_values[string].ol_flow.bc_flag == BCT_PRS_NEU) hydro_bc_flag = BCT_PRS_NEU;
    
    if (nnodes == NDONTRI) {isElementTriangle = TRUE;}
    else {isElementTriangle = FALSE;}
    
    if      (elem2d->bflag == 2) {isSidewall = TRUE;}
    else if (elem2d->bflag == 0) {isSurface = TRUE;}
    else if (elem2d->bflag == 1) {isBed = TRUE;}
    else {
        printf("bflag = %d\n",elem2d->bflag);
        tl_error(">> 2D Boundary face bflag is not defined.");
    }
    
#ifdef _DEBUG
    assert(nnodes == NDONQUAD || nnodes == NDONTRI);
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    SVECT elem_vel[nnodes];
    global_to_local_svect(sw->vel, elem_vel, elem2d->nodes, nnodes);
    if (perturb_var == PERTURB_W) {
        elem_vel[perturb_node].z += perturb_sign * perturbation;
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // local elevation displacements
    double elem_displacement[nnodes];
    global_to_local_dbl(sw->displacement, elem_displacement, elem2d->nodes, nnodes);
    
    // node positions
    SVECT elem_nds[nnodes];
    for (i=0; i<nnodes; i++) {
        elem_nds[i].x = grid->node[elem2d->nodes[i]].x;
        elem_nds[i].y = grid->node[elem2d->nodes[i]].y;
        elem_nds[i].z = grid->node[elem2d->nodes[i]].z + elem_displacement[i];
    }
#ifdef _DEBUG  // this must be true for analytic quadrilateral integrations to be correct
    if (isElementTriangle != TRUE) {
        assert(fabs(elem_nds[2].x-elem_nds[1].x) < 1e-6 && fabs(elem_nds[2].y-elem_nds[1].y) < 1e-6);
        assert(fabs(elem_nds[3].x-elem_nds[0].x) < 1e-6 && fabs(elem_nds[3].y-elem_nds[0].y) < 1e-6);
    }
#endif
    
    // get current and initial element areas
    SVECT elem_nds_fixed[nnodes]; // only used for quadrilaterals
    if (isElementTriangle == TRUE) {
        get_triangle_linear_djac_nrml_gradPhi(elem2d, NULL, elem_nds);
    } else {
        // note :: this assume that no displacements are hotstarted ...
        for (i=0; i<nnodes; i++) {
            elem_nds_fixed[i].x = grid->node[elem2d->nodes[i]].x;
            elem_nds_fixed[i].y = grid->node[elem2d->nodes[i]].y;
            elem_nds_fixed[i].z = grid->node[elem2d->nodes[i]].z;
        }
        elem2d->nrml = get_elem2d_normals(elem_nds);
    }
    
    // calculate normal velocities
    double elem_normal_vel[nnodes];
    for (i=0; i<nnodes; i++) {
        elem_normal_vel[i] = svect_dotp(elem_vel[i], elem2d->nrml);
    }
    
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        printf("SW-3D CONTINUITY BOUNDARY ELEM RESID :: ie: %d \t dt: %20.10f",ie,dt);
        if      (isSidewall == TRUE) {printf(" SIDEWALL FACE ");}
        else if (isSurface  == TRUE) {printf(" SURFACE FACE ");}
        else if (isBed      == TRUE) {printf(" BED FACE ");}
        if (perturb_var == PERTURB_W) {
            printf("\t perturbing W || node: %d || perturbation: %20.10e\n",elem2d->nodes[perturb_node],perturb_sign*perturbation);
        }
        if (isElementTriangle) {
            printf("djac2d: %30.20f \t djac2d3d: %30.20f \t djac2d3d_fixed: %30.20f \n",elem2d->djac,elem2d->djac3d,elem2d->djac3d_fixed);
        }
        printf("nodes: %d %d %d \n",elem2d->nodes[0],elem2d->nodes[1],elem2d->nodes[2]);
        printScreen_debug2_dbl("elem_displacement", elem_displacement, nnodes, elem2d->nodes);
        printScreen_debug_vec("node locations: ",elem_nds, nnodes);
        if (!isElementTriangle) printScreen_debug_vec("node locations_fixed: ",elem_nds_fixed, nnodes);
        printScreen_debug_svect("elem_vel", elem_vel, nnodes, elem2d->nodes);
        printScreen_debug2_dbl("elem_normal_vel", elem_normal_vel, nnodes, elem2d->nodes);
        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    }
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    double elem_rhs_bed[nnodes]; sarray_init_dbl(elem_rhs_bed, nnodes);
    sarray_init_dbl(elem_rhs, nnodes);
    sarray_init_dbl(elem_rhs_noBed, nnodes);
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   VELOCITY CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the convection boundary flux addition to the SW 2D elemental residual. \n
     * \note  There is an assumption here that the inflow or exit does not add or remove momentum.
     *********************************************************************************************/
    
    if (hydro_bc_flag == BCT_VEL_NEU || hydro_bc_flag == BCT_DIS_NEU) {

        // read normal velocity/discharge
        int isers = str_values[string].flow.isigma;
        double hydro_flux_normal = -sseries_get_value(isers, mod->series_head, 0); // (-) for flow into domain
        if (hydro_bc_flag == BCT_DIS_NEU) {
            hydro_flux_normal /= str_values[string].total_area;
        }
        
        if (isElementTriangle) {
            integrate_triangle_phi(elem2d->djac3d_fixed, dt * hydro_flux_normal, elem_rhs);
        } else {
            integrate_quadZ_phi(elem_nds_fixed, dt * hydro_flux_normal, elem_rhs);
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("3D CONTINUITY BCT_VEL_NEU/BCT_DIS_NEU: ",nnodes, ie, elem2d->nodes, elem_rhs);
            Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
        }
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    OUTFLOW CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the outflow addition to the elemental residual. \n
     * \note CJT \:: calculates implicit flow using implicitly calculated water elevation
     *********************************************************************************************/
    
    else if (hydro_bc_flag == BCT_OUTFLOW) {
        
        if (isElementTriangle) {
            integrate_triangle_phi_f(elem2d->djac3d, dt, elem_normal_vel, elem_rhs);
        } else {
            integrate_quadZ_phi_f(elem_nds, dt, elem_normal_vel, elem_rhs);
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("3D CONTINUITY BCT_OUTFLOW: ",nnodes, ie, elem2d->nodes, elem_rhs);
            Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
        }
#endif
        
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                            BED
     *--------------------------------------------------------------------------------------------
     * Calculates the bed flux contribution to the residual. \n
     * \note Should be mimimal if Kinematic BC is applied to the surface and this is consistent with DACONT
     * \note That said, we still need to add the flux to the residual, but we'll also store to track
     * \note When the sediment is included, this will need to include the bed velocity
     *********************************************************************************************/
    
    else if (hydro_bc_flag == BCT_BED) {
        
#ifdef _DEBUG
        assert(isElementTriangle == TRUE);
#endif
        
#ifdef _SEDIMENT
        // -----------------------------------------------------------------------------------------
        // depth-averaged continuity addition ------------------------------------------------------
        //  the bed surface change of the depth-averaged continuity equations here (note: dpl here is only from sediment)
        //        double elem_old_eta[NDONTRI];  // get wse at t(i+1/2)
        //        ELEM2D_GET_TPOSITION(elem_old_eta, elem_old_displacement, elem_older_displacement, mod->tau_temporal, dt, dt);
        //        double elem_eta[NDONTRI];  // get wse at t(i+3/2)
        //        ELEM2D_GET_TPOSITION(elem_eta, elem_displacement, elem_old_displacement, mod->tau_temporal, dt, dt);
        //
        //        double eta_integral[nnodes]; sarray_init_dbl(eta_integral, nnodes);
        //        double eta_old_integral[nnodes]; sarray_init_dbl(eta_old_integral, nnodes);
        //        integrate_triangle_phi(elem2d->djac2d, elem_eta, eta_integral);
        //        integrate_triangle_phi(elem2d->djac2d, elem_old_eta, eta_old_integral);
        //        for (i=0; i<nnodes; i++) {
        //            rhs[i] = eta_integral[i] - eta_old_integral[i];
        //        }
#endif
        integrate_triangle_phi_f(elem2d->djac3d, dt, elem_normal_vel, elem_rhs_bed);

#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("3D CONTINUITY BCT_BED: ",nnodes, ie, elem2d->nodes, elem_rhs_bed);
            Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
        }
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                       FREE SURFACE
     *--------------------------------------------------------------------------------------------
     * Calculates the bed source contribution to the residual. \n
     * \note When raining on the free surface, it was found that they had to account for momentum from
     *        rain or unrealistic w's emerged.
     * \note Assume pressure = 0 at surface
     *********************************************************************************************/
    else if (hydro_bc_flag == BCT_FRS || hydro_bc_flag == BCT_WATER_SOURCE) {

#ifdef _DEBUG
        assert(isElementTriangle == TRUE);
#endif
        
        double elem_old_displacement[nnodes], elem_older_displacement[nnodes];
        global_to_local_dbl(sw->old_displacement, elem_old_displacement, elem2d->nodes, elem2d->nnodes);
        global_to_local_dbl(sw->older_displacement, elem_older_displacement, elem2d->nodes, elem2d->nnodes);
#ifdef _SEDIMENT
        double elem_bed_displacement[nnodes], elem_old_bed_displacement[nnodes], elem_older_bed_displacement[nnodes];
        global_to_local_dbl(sw->bed_displacement, elem_bed_displacement, elem2d->nodes, elem2d->nnodes);
        global_to_local_dbl(sw->old_bed_displacement, elem_old_bed_displacement, elem2d->nodes, elem2d->nnodes);
        global_to_local_dbl(sw->older_bed_displacement, elem_older_bed_displacement, elem2d->nodes, elem2d->nnodes);
        // cjt :: add bed displacement to total displacement ... be carefull here, no need for db/dt now I think
        sarray_add_replace_dbl(elem_displacement, elem_bed_displacement, elem2d->nnodes);
        sarray_add_replace_dbl(elem_old_displacement, elem_old_bed_displacement, elem2d->nnodes);
        sarray_add_replace_dbl(elem_older_displacement, elem_older_bed_displacement, elem2d->nnodes);
#endif
        
        // read normal velocity/discharge
        int isers = str_values[string].flow.isigma;
        double hydro_flux_normal = -sseries_get_value(isers, mod->series_head, 0); // (-) for flow into domain
        if (hydro_bc_flag == BCT_DIS_NEU) {
            hydro_flux_normal /= (str_values[string].total_area);
        }
        
        double explicit_flow[nnodes];    sarray_init_dbl(explicit_flow,nnodes);
        integrate_triangle_phi(elem2d->djac, dt * hydro_flux_normal, explicit_flow);
        
        double elem_old_eta[NDONTRI];  // get wse at t(i+1/2)
        ELEM2D_GET_TPOSITION(elem_old_eta, elem_old_displacement, elem_older_displacement, mod->tau_temporal, dt, dt);
        double elem_eta[NDONTRI];      // get wse at t(i+3/2)
        ELEM2D_GET_TPOSITION(elem_eta, elem_displacement, elem_old_displacement, mod->tau_temporal, dt, dt);
        
        double eta_integral[nnodes];     sarray_init_dbl(eta_integral,nnodes);
        double eta_old_integral[nnodes]; sarray_init_dbl(eta_old_integral,nnodes);
        integrate_triangle_phi_f(elem2d->djac, 1., elem_eta, eta_integral);
        integrate_triangle_phi_f(elem2d->djac, 1., elem_old_eta, eta_old_integral);
        
        for (i=0; i<nnodes; i++) {
            elem_rhs[i] += (eta_integral[i] - eta_old_integral[i]) - explicit_flow[i];
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("3D CONTINUITY BCT_FRS: ",nnodes, ie, elem2d->nodes, elem_rhs);
            Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
        }
#endif

    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                             TAILWATER/PRESSURE CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the tailwater or pressure addition to the elemental residual. \n
     * \note
     *********************************************************************************************/
    
    else if (hydro_bc_flag == BCT_PRS_NEU) {
        
        if (isElementTriangle == TRUE) {
            integrate_triangle_phi_f(elem2d->djac3d, dt, elem_normal_vel, elem_rhs);
        } else {
            integrate_quadZ_phi_f(elem_nds, dt, elem_normal_vel, elem_rhs);
        }
        
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("3D CONTINUITY BCT_PRS_NEU: ",nnodes, ie, elem2d->nodes, elem_rhs);
            Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
        }
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                DIRICHLET BOUNDARY CONDITIONS
     *--------------------------------------------------------------------------------------------
     * Zeros matrix rows/columns when Dirichlet boundaries are used. \n
     * \note
     *********************************************************************************************/
//    for (i=0; i<nnodes; i++) {
//        int node_string = grid->node[ elem2d->nodes[i] ].string;
//        if (node_string > NORMAL) {
//            if (hydro_bc_flag == BCT_DIR || hydro_bc_flag == BCT_CEQ) {
//                elem_rhs[i] = 0;
//            }
//        }
//    }
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        rhs_1dof("3D CONTINUITY DIRICHLET BCS: ",nnodes, ie, elem2d->nodes, elem_rhs);
        Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    }
#endif
    
    
    // Here we add the bed flux residual to the overall residual and store the residual w/o bed flux for monitoring
    for (i=0; i<nnodes; i++) {
        elem_rhs_noBed[i] = elem_rhs[i];
        elem_rhs[i] += elem_rhs_bed[i];
    }
    
    
#ifdef _DEBUG
    time_t time2;  time(&time2);
    TIME_IN_WVEL_BOUNDARY_RESID += difftime(time2,time1);
#endif
}

//***************************************************************************//
//***************************************************************************//
//***************************************************************************//

