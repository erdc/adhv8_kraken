/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 2D diffusive wave elemental residual.
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
 *  \details Solves the following weak, discrete body terms of the diffusive wave equation: \n
 * \f[  R_{i,body}^e = \bodyTime{\,2d}{e}{\phidd{i}}{\depth{h}} \,-
 *                     \bodyConv{\,2d}{e}{\phidd{i}}{ (\depth{h} \, k \, \grad \elev{h}) } \,-
 *                     \bodySource{\,2d}{e}{\phidd{i}}{\srcGen{dw}{h}{}} \,+
 *                     \bodySUPG{\supg{i}{dw}{e}}
 * \f]
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_diffusive_body_resid(SMODEL *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    
    int DEBUG_LOCAL = OFF;
    
    /*
     if (ie==0 || ie==1 || ie==2) {
     DEBUG_LOCAL = ON;
     } else {
     DEBUG_LOCAL = OFF;
     }
     */
    
    int i;
    double t1, t2, t3, t4, t5;
    
    
    // aliases
    SSW_2D *sw2 = mod->sw->d2;
    SGRID *grid = mod->grid;
    SELEM_2D *elem2d = &(grid->elem2d[ie]); // be careful not to shallow copy here
    SVECT2D *grad_shp = elem2d->grad_shp;
    int nnodes = elem2d->nnodes;
    double alpha = mod->tau_temporal;
    double dt = mod->dt;
    double djac = elem2d->djac;
    double g = mod->gravity;
    double drying_lower_limit = mod->drying_lower_limit;
    int imat = elem2d->mat;
    STR_VALUE str_values = mod->str_values[elem2d->string];
    double roughness = str_values.fterms.manningsn;
    SQUAD *quad = grid->quad_rect; // for quadrilateral quadrature calculations
    
    int isTriangle = NO;  if (nnodes == NDONTRI) isTriangle = YES;
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        if (isTriangle == TRUE) assert(djac>SMALL);
        assert(alpha > -1E-6 && alpha < 1.0 + SMALL);
        assert(imat > -1);
        assert(perturb_var == UNSET_INT || perturb_var == PERTURB_H);
    }
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GRID VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    double node_z[nnodes];
    SNODE nodes[nnodes];
    for (i=0; i<nnodes; i++) {
        snode_copy(&(nodes[i]), grid->node[elem2d->nodes[i]]);
        node_z[i] = nodes[i].z;
    }
    
    SVECT elem_nds[nnodes];
    snode2svect(nodes, elem_nds, nnodes);
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    double elem_head[nnodes];
    global_to_local_dbl(sw2->head, elem_head, elem2d->nodes, nnodes);
    if (perturb_var == PERTURB_H) elem_head[perturb_node] += perturb_sign * perturbation;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    double elem_old_head[nnodes], elem_older_head[nnodes];
    global_to_local_dbl(sw2->old_head, elem_old_head, elem2d->nodes, nnodes);
    global_to_local_dbl(sw2->older_head, elem_older_head, elem2d->nodes, nnodes);
    
    double elem_avg_depth = 0., area = 0.;
    if (isTriangle == TRUE) {
        area = djac;
        elem_avg_depth = integrate_triangle_f(1.,1.,elem_head); // cjt :: djac = area cancel here
    } else {
        area = integrate_quadrilateral_area(elem_nds, 1.);
        elem_avg_depth = integrate_quadrilateral_f(elem_nds,1./area,elem_head);
    }
    
    // cjt :: going to be general here and assume velocities are not elementally constant.
    //        They are here because we are using the elementally averaged depth to calculated them.
    // gkc :: "vv" appears to be elementally averaged, but elem_vel[nodes] is not since depth at
    //        each node may be different.
    SVECT2D elem_vel[nnodes], vv;
    vv = getDiffusiveWaveVelocities(nnodes, elem_head, elem_avg_depth, node_z, grad_shp, roughness, mod->manning_units_constant, elem_vel);
    /* Gajanan gkc: WARNING: elem_vel is overridden below!!! Will set elem_vel[...] = vv */
    
    // for utilities and output
    double u[nnodes], v[nnodes];
    for (i=0; i<nnodes; i++) {
        //u[i] = elem_vel[i].x;  v[i] = elem_vel[i].y; /* Gajanan gkc: elem_vel calculation in getDiffusiveWaveVelocities possibly wrong. */
        if (elem_head[i] > 0.0){
            u[i] = vv.x;  v[i] = vv.y;
            elem_vel[i].x = vv.x;  elem_vel[i].y = vv.y; /* Gajanan gkc: elem_vel calculation in getDiffusiveWaveVelocities possibly wrong. */
        }
        else{
            u[i] = 0.0; v[i] = 0.0;
            elem_vel[i].x = 0.0;  elem_vel[i].y = 0.0; /* Gajanan gkc: This is likely unnecessary. */
        }
        
        /**********/
        //sw2->vel[elem2d->nodes[i]].x = vv.x; /* Most certainly wrong; elemental vel being assigned to nodes. */
        //sw2->vel[elem2d->nodes[i]].y = vv.y; /* Most certainly wrong; elemental vel being assigned to nodes. */
    }
    
    
    double elem_avg_u = 0., elem_avg_v = 0.;
    if (isTriangle == TRUE) { // cjt :: djac = area cancel here
        elem_avg_u = integrate_triangle_f(1.,1.,u);
        elem_avg_v = integrate_triangle_f(1.,1.,v);
    } else {
        elem_avg_u = integrate_quadrilateral_f(elem_nds,1./area,u);
        elem_avg_v = integrate_quadrilateral_f(elem_nds,1./area,v);
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        printf("DIFFUSIVE WAVE 2D ELEM RESID :: ie: %d \t dt: %20.10f \t area: %20.10f",ie,dt,area);
        if (perturb_var == PERTURB_H) {
            printf("\t perturbing H || node: %d || perturbation: %20.10e\n",nodes[perturb_node].id,perturb_sign*perturbation);
        }
        selem2d_printScreen(elem2d);
        
        printf("\n------------------ ELEMENTAL AVERAGES --------------------\n");
        printf("elem_avg_depth: %20.10e\n",elem_avg_depth);
        
        printf("\n------------------ FUNCTION GRADIENTS --------------------\n");
        
        printf("\n------------------- LINEAR VARIABLES ---------------------\n");
        printScreen_debug2_dbl("elem_head", elem_head, nnodes, elem2d->nodes);
        printScreen_debug2_dbl("elem_old_head", elem_old_head, nnodes, elem2d->nodes);
        printScreen_debug2_dbl("elem_older_head", elem_older_head, nnodes, elem2d->nodes);
        printScreen_debug_svec2d("elem_vel", elem_vel, nnodes, elem2d->nodes);
        
        printf("\n--------------------- OTHER STUFF ------------------------\n");
        printf("roughness: %20.10e \t mannings constant: %20.10e\n",roughness,mod->manning_units_constant);
        
    }
#endif
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    
    double rhs[nnodes];
    sarray_init_dbl(elem_rhs,nnodes);
    double vars[4]; // for passing doubles through wet-dry routines
    
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    TEMPORAL CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates the temporal addition to the elemental residual. \n
     * \f$ \bodyTime{\,2d}{e}{\phidd{i}}{\depth{h}} \f$
     *
     * \note GS \:: not sure if consistent or lumped mass should be used here
     $ \note CJT \:: returns below returns the elementally averaged temporal contribution
     *********************************************************************************************/
    
    double dhdt[nnodes];      // store for SUPG later
    double gwfactor[nnodes];      // store for DW-GW coupling. 0 for dry nodes, 1 for wet nodes
    for(i=0;i<nnodes;i++) gwfactor[i]=1.0;
    double head_32dt; // depth at t+(3/2)dt
    double head_12dt; // depth at t+(1/2)dt
    for (i=0; i<nnodes; i++) {
        head_32dt = get_second_order(elem_head[i], elem_old_head[i]);
        head_12dt = get_second_order(elem_old_head[i], elem_older_head[i]);
        dhdt[i] = (alpha * (head_32dt - head_12dt) + (1-alpha) * (elem_head[i]- elem_old_head[i]))/dt;
#ifdef _ADH_GROUNDWATER
#ifdef _DWGW_COUPLING
        if (mod->flag.GW_FLOW == ON){
            if (/*(head_32dt <= 0.0 || head_12dt <= 0.0) &&*/ (elem_head[i] <= 0.0 && elem_old_head[i] <= 0.0)){
                gwfactor[i] = 0.0;
            }
        }
#endif
#endif
    }
    
    sarray_init_dbl(rhs,nnodes);
    if (isTriangle == YES) {
//        DOF_3 temprhs[nnodes];
//        vars[2] = 0.0;
//        vars[3] = 0.0;
//        vars[0] = 1.0;
//
//        dof3_init_array(temprhs,nnodes);
//        vars[1] = (1. + alpha / 2.);
//        t1 = fe_sw2_wet_dry_wrapper(temprhs, elem_nds, elem_head, NULL, NULL, NULL, NULL, NULL, djac, ON, DEBUG, vars, fe_sw2_wd_continuity_temporal_triangle);
//        for (i=0; i<nnodes; i++) {
//            rhs[i] += vars[1] * temprhs[i].c_eq/dt;
//        }
//
//        dof3_init_array(temprhs,nnodes);
//        vars[1] = (-1.) * (1. + alpha);
//        t1 = fe_sw2_wet_dry_wrapper(temprhs, elem_nds, elem_head, NULL, NULL, NULL, NULL, NULL, djac, ON, DEBUG, vars, fe_sw2_wd_continuity_temporal_triangle);
//        for (i=0; i<nnodes; i++) {
//            rhs[i] += vars[1] * temprhs[i].c_eq/dt;
//        }
//
//        dof3_init_array(temprhs,nnodes);
//        vars[1] =  alpha / 2.;
//        t1 = fe_sw2_wet_dry_wrapper(temprhs, elem_nds, elem_head, NULL, NULL, NULL, NULL, NULL, djac, ON, DEBUG, vars, fe_sw2_wd_continuity_temporal_triangle);
//        for (i=0; i<nnodes; i++) {
//            rhs[i] += vars[1] * temprhs[i].c_eq/dt;
//        }

        //integrate_triangle_phi_f(djac, 1., dhdt, rhs);        // consistent mass matrix
        integrate_triangle_phi_f_lump(djac, 1., dhdt, rhs);   // lumped mass matrix
    } else {
        //integrate_quadrilateral_phi_f(elem_nds, 1., dhdt, rhs);       // consistent mass matrix
        integrate_quadrilateral_phi_f_lump(elem_nds, 1., dhdt, rhs);  // lumped mass matrix
    }
    for (i=0; i<nnodes; i++) elem_rhs[i] += gwfactor[i] * dt * rhs[i];
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        rhs_1dof("2D DIFFUSIVE WAVE || TEMPORAL: ",nnodes, ie, elem2d->nodes, rhs);
    }
    Is_DoubleArray_Inf_or_NaN(rhs, nnodes ,__FILE__ ,__LINE__);
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    ADVECTION CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the diffusion addition to the diffusive wave elemental residual. \n
     * \f$ \bodyConv{\,2d}{e}{\phidd{i}}{ (\depth{h} \, k \, \grad \elev{h}) } \f$
     *
     * \note GS  \:: Wrap in wet/dry integration
     * \note CJT \:: Includes SUPG contribution from convection
     * \note CJT \:: Until we include quadrilaterals in wetting and drying, this is only for triangles
     *********************************************************************************************/
    
    sarray_init_dbl(rhs,nnodes);
    if (isTriangle == YES) {
        DOF_3 dof3_rhs[nnodes]; dof3_init_array(dof3_rhs,nnodes);
        vars[0] = dt;
        vars[1] = mod->tau_pg;
        vars[2] = mod->manning_units_constant;
        vars[3] = roughness;
        t1 = fe_sw2_wet_dry_wrapper(dof3_rhs, elem_nds, elem_head, grad_shp, NULL, elem_vel, NULL, elem_head, djac, ON, DEBUG, vars, fe_diffusive_conv);
        for (i=0; i<nnodes; i++) {rhs[i] = dof3_rhs[i].c_eq;}
    } else {
        // cjt :: no wetting and drying here yet
        integrate_quadrilateral_gradPhi_dot_f_v(elem_nds, dt, elem_head, elem_vel, rhs); // cjt :: uses full depth
        integrate_quadrilateral_gradPhi_dot_v(elem_nds, dt * elem_avg_depth, elem_vel, rhs) ; // cjt :: uses elementally averaged depth
    }
    for (i=0; i<nnodes; i++) elem_rhs[i] += gwfactor[i] * rhs[i];

#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        printf("convection: v[0]: %20.10e\n", elem_head[0]);
        printf("convection: v[1]: %20.10e\n", elem_head[1]);
        printf("convection: v[2]: %20.10e\n", elem_head[2]);
        rhs_1dof("2D DIFFUSIVE WAVE || ADVECTION: ",nnodes, ie, elem2d->nodes, rhs);
    }
    Is_DoubleArray_Inf_or_NaN(rhs, nnodes ,__FILE__ ,__LINE__);
#endif
    
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                  SOURCE/SINK CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates the source (rain/evaporation/infiltration) addition to the Diffusive wave \n
     * elemental residual. \n
     * \f$ \bodySource{\,2d}{e}{\phidd{i}}{\srcGen{dw}{h}{}} \f$
     *
     * \note GS \:: For Diffusive Wave this is how I'm doing boundaries for now.
     *             If there is an inflow add it in as rain. This is not a problem because we
     *             don't bring flow in with momentum anyways.
     * \note GS \:: Wrapping this in wet/dry integration does not work.
     *********************************************************************************************/
    
    // water sourcing
    double elem_water_source = 0.;
    if (str_values.ol_flow.bc_flag == BCT_WATER_SOURCE) {
        int isers = str_values.ol_flow.isigma;
        elem_water_source = 1.0 * sseries_get_value(isers, mod->series_head, 0);
    }
    
    // infiltration
    double elem_infiltration = 0.;
    if (mod->green_ampt == TRUE && elem_avg_depth > 0) {
        elem_infiltration = mod->mat[imat].sw->hyd_conductivity * (mod->mat[imat].sw->psi + mod->mat[imat].sw->rooting_depth)/(mod->mat[imat].sw->rooting_depth);
    }
    
    // total water source
    double total_elem_source = elem_infiltration + elem_water_source;
    
    // now integrate and add to the element residual
    sarray_init_dbl(rhs,nnodes);
    if (isTriangle == YES) {
        integrate_triangle_phi(djac, total_elem_source, rhs);
    } else {
        integrate_quadrilateral_phi(elem_nds, total_elem_source, rhs);
    }
    for (i=0; i<nnodes; i++) {elem_rhs[i] -= dt * rhs[i];} /* No  gwfactor[i] here according to Kollett Maxwell paper. */
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        rhs_1dof("2D DIFFUSIVE WAVE || SOURCE: ",nnodes, ie, elem2d->nodes, rhs);
    }
    Is_DoubleArray_Inf_or_NaN(rhs, nnodes ,__FILE__ ,__LINE__);
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                      SUG CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the SUPG stabilization addition to the Diffusive wave elemental residual. \n
     * \f$
     * \supg{i}{dw}{e} = \supgEq{2d}{e}{\velb{h}}{\phidd{i}}{\resid{i}{}{dw}} =
     *                   \supgEq{2d}{e}{\velb{h}}{\phidd{i}}
     *                   {\lrpbb{\strongTrnsLinear{dw}{2d}{\depth{h}}{\velb{h}}}}
     * \f$
     *
     * \note CJT \:: The convection SUPG contribution is calculated in the wet/dry wrapper on Galerkin addition.  \n
     *               This is done for consistency, using the wet/dry portions of the element for both.
     * \note CJT \:: Elementally averaged values for temporal, convection and vel_x,y terms are used
     * \remark CJT \:: Where does the 0.7 come from?
     *********************************************************************************************/
    
    // calculate SUPG stabilization coefficient
    int tau_method_flag = 0; // not used
    int le_method_flag = 2; // area based
    double tau =  fe_get_supg_tau_sw(nnodes, elem_nds, g, elem_avg_depth, elem_avg_u, elem_avg_v,0., NULL, NULL, NULL, djac,mod->tau_pg, grid->ndim, tau_method_flag, le_method_flag);
    
    sarray_init_dbl(rhs,nnodes);
    if (isTriangle) {
        t1 =  integrate_triangle_f(1.,1.,dhdt);              // elementally averaged temporal term
        t2 =  total_elem_source;                             // elem averaged source term
        integrate_triangle_gradPhi_dot_v(grad_shp, djac, t1 - t2, elem_vel, rhs);
    } else {
        double elem_avg_temporal = integrate_quadrilateral_f(elem_nds, 1./area, dhdt);                 // elem averaged temporal term
        double elem_avg_source = total_elem_source;                                                   // elem averaged source term
        double elem_avg_convection = integrate_quadrilateral_grad_dot_h_v(elem_nds, 1./area, elem_vel, elem_head);  // elem averaged advection term
        double elem_avg_resid = elem_avg_temporal - elem_avg_source + elem_avg_convection;            // elementally averaged residual
        integrate_quadrilateral_gradPhi_dot_v(elem_nds, elem_avg_resid, elem_vel, rhs);
    }
//#ifdef _ADH_GROUNDWATER
//#ifdef _DWGWCOUPLING
//    if (mod->flag.GW_FLOW != ON){
//#endif
//#endif
        for (i=0; i<nnodes; i++) elem_rhs[i] += gwfactor[i] * tau * dt * rhs[i];
//#ifdef _ADH_GROUNDWATER
//#ifdef _DWGWCOUPLING
//    }
//#endif
//#endif
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        rhs_1dof("2D DIFFUSIVE WAVE || SUPG: ",nnodes, ie, elem2d->nodes, rhs);
    }
    Is_DoubleArray_Inf_or_NaN(rhs, nnodes ,__FILE__ ,__LINE__);
#endif
    
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    SET DIRICHLET BCS
     *--------------------------------------------------------------------------------------------
     * Sets the Dirichlet boundary conditions. \n
     *********************************************************************************************/
    int node_istring = UNSET_INT;
    for (i=0; i<nnodes; i++) {
        node_istring = grid->node[i].string;
        if (node_istring > NORMAL) {
            if (mod->str_values[node_istring].ol_flow.bc_flag == BCT_DIR || mod->str_values[node_istring].ol_flow.bc_flag == BCT_CEQ) {
                elem_rhs[i] = 0.0;
            }
        }
    }
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        rhs_1dof("2D BED LOAD || TOTAL: ",nnodes, ie, elem2d->nodes, elem_rhs);
    }
    Is_DoubleArray_Inf_or_NaN(rhs, nnodes ,__FILE__ ,__LINE__);
#endif
}































