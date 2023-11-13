/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_sw2_transport_elem_resid.c This file collections functions responsible for
 *          the 2D shallow water transport body contributions to the elemental residual.    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

// prototypes
void fe_sw2_bl_temporal(int nnodes, SVECT *elem_nds, double djac, double *elem_thickness, double *elem_c, double dt_factor, DOF_3 *rhs, char *string, int DEBUG, int DEBUG_LOCAL);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 2D bed load elemental residual.
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
 * @param[in]  pertubation   the Newton pertubation
 * @param[in]  perturb_node  the node to be perturbed
 * @param[in]  perturb_var   the variable to be perturbed
 * @param[in]  perturb_sign  the direction of Newton perturbation
 * @param[in]  DEBUG         a debug flag
 *
 * \note Bed load thickness = 1, so always wet
 *
 * \details Solves the following weak, discrete body terms of the diffusive wave equation: \n
 * \f[  R_{i,body}^e = \bodyTime{\,2d}{e}{\phidd{i}}{\cbl{h}{j}{bl}} \,-
 *                     \bodyConv{\,2d}{e}{\phidd{i}}{ (\velc{h}{bl} \, \cbl{h}{j}{bl}) } \,+
 *                     \bodyDiffusion{\,2d}{e}{\phidd{i}}{\diffTensor{bl}{h}{\diffTensorTRN}{j,}} \,- 
 *                     \bodySource{\,2d}{e}{\phidd{i}}{\srcGen{bl}{h}{j,}} \,+
 *                     \bodySUPG{\supg{i}{bl}{e}} 
 * \f]
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_bedload_body_resid(SMODEL *mod, double *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
#ifdef _SEDIMENT
    
    int i, DEBUG_LOCAL = OFF;
    
    // aliases
    SGRID *grid = mod->grid;
    SELEM_2D *elem2d = &(grid->elem2d[ie]);
    SVECT2D * grad_shp = elem2d->grad_shp;
    int nnodes = elem2d->nnodes;
    double dt = mod->dt;
    double gravity = mod->gravity;
    STR_VALUE *str_values = mod->str_values;
    int string = elem2d->string;
    SSED *sed = mod->sed;
    int ised = mod->ised;
    int imat = elem2d->mat;
    SQUAD *quad = grid->quad_rect; // for quadrilateral quadrature calculations
    double diffusion = mod->mat[imat].sed[ised].d_m;
    
    int i, DEBUG_LOCAL = OFF;
    double rhs[nnodes];          // local rhs for testing
    
    int isTriangle = NO;  if (nnodes == NDONTRI) isTriangle = YES;
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        assert(imat > 0);
        assert(perturb_var == UNSET_INT || perturb_var == C);
    }
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // ENSURE 3D CALL COMPATIBILITY
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
   
    double djac = 0;
    
    int node_bed_id[nnodes];
    if (grid->ndim == 3) {
        djac = elem2d->djac3d;   // 2d area of a triangular face of a 3d tet
        for (i; i<nnodes; i++) {
            node_bed_id[i] = grid->nodeID_3d_to_2d_bed[elem2d->nodes[i]];
        }
    } else {
        djac = elem2d->djac;
        for (i=0; i<nnodes; i++) {
            node_bed_id[i] = elem2d->nodes[i];
        }
    }
    
    double *depth = NULL;
    double *old_depth = NULL;
    double *older_depth = NULL;
    if (mod->flag.SW2_FLOW) {
        depth = mod->sw->d2->head;
        old_depth = mod->sw->d2->old_head;
        older_depth = mod->sw->d2->older_head;
    } else {
        depth = mod->sw->d3->depth;
        old_depth = mod->sw->d3->old_depth;
        older_depth = mod->sw->d3->old_depth;
    }
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        if (isTriangle == TRUE) assert(djac > SMALL);
    }
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GRID VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    SNODE nodes[nnodes];
    for (i=0; i<nnodes; i++) {
        snode_copy(&(nodes[i]), grid->node[elem2d->nodes[i]]);
    }
    
    SVECT elem_nds[nnodes];
    if (isTriangle == NO) snode2svect(nodes, elem_nds, nnodes);
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    double elem_c[nnodes];
    global_to_local_dbl(sed->bedload[ised].c, elem_c, node_bed_id, nnodes);
    if (perturb_var == PERTURB_C) elem_c[perturb_node] += perturb_sign * pertubation;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // depths
    double elem_head[nnodes], elem_old_head[nnodes], elem_older_head[nnodes];
    global_to_local_dbl(depth, elem_head, node_bed_id, nnodes);
    global_to_local_dbl(old_depth, elem_old_head, node_bed_id, nnodes);
    global_to_local_dbl(older_depth, elem_older_head, node_bed_id, nnodes);
    
    // concentrations
    double elem_old_c[nnodes], elem_older_c[nnodes];
    global_to_local_dbl(sed->bedload[ised].old_c, elem_old_c, node_bed_id, nnodes);
    global_to_local_dbl(sed->bedload[ised].older_c, elem_older_c, node_bed_id, nnodes);
    
    // thickness
    double elem_thickness[nnodes], elem_old_thickness[nnodes], elem_older_thickness[nnodes];
    global_to_local_dbl(sed->bedload[ised].thick, elem_thickness, node_bed_id, nnodes);
    global_to_local_dbl(sed->bedload[ised].old_thick, elem_old_thickness, node_bed_id, nnodes);
    global_to_local_dbl(sed->bedload[ised].older_thick, elem_older_thickness, node_bed_id, nnodes);
#ifdef _DEBUG
    for (i=0; i<nnodes; i++) {
        assert(fabs(elem_thickness[i] - 1) < SMALL && fabs(elem_old_thickness[i] - 1) < SMALL && fabs(elem_older_thickness[i] - 1) < SMALL);
    }
#endif
    
    // sources and sinks
    double elem_source[nnodes], elem_sink[nnodes];
    global_to_local_dbl(sed->bedload[ised].source, elem_source, node_bed_id, nnodes);
    global_to_local_dbl(sed->bedload[ised].sink, elem_sink, node_bed_id, nnodes);
    
    // bed velocity
    SVECT2D elem_bedload_vel[nnodes];
    global_to_local_svect2d(sed->bedload[ised].v, elem_bedload_vel, node_bed_id, nnodes);
    double u[nnodes], v[nnodes]; // for convienence
    for (i=0; i<nnodes; i++) {
        u[i] = elem_bedload_vel[i].x;  v[i] = elem_bedload_vel[i].y;
    }
    
    // bed shear
    double elem_bed_shear[nnodes];
    global_to_local_dbl(sed->bedload[ised].shear, elem_bed_shear, node_bed_id, nnodes);
    
    // calculate eddy viscosity tensor
    double diff_bl = 0.0, gcsh = 0.0;
    for (i = 0; i < nnodes; i++) {
        gcsh = 0.05 * 1.65 * g * mod->density * mod->sed->grain[ised].diameter;
        diff_bl += 5 * MAX(elem_head[i], 0.0) * sqrt(MAX(elem_bed_shear[i], gcsh) / mod->density);
    }
    diff_bl *= one_3;
    STENSOR2D d; stensor2d_init(&d);    // tensor material property
    double coef = diffusion;
    if (mod->mat[imat].sw->eev_mode == YES) {
        coef = mod->mat[imat].sw->eev_coef * diff_bl;
    }
    d.xx = coef;  d.yy = coef;  d.xy = 0.0;
    
    // elemental averages, function gradients, and 2d element areas
    double area = 0.;
    SVECT2D elem_vbed_avg,elem_grad_c, elem_grad_thick, elem_grad_u, elem_grad_v;
    double elem_t_avg = 0., elem_u_avg = 0., elem_v_avg = 0., elem_source_avg = 0., elem_sink_avg = 0., elem_strong_advection_avg = 0.;
    if (isTriangle == TRUE) {  // cjt :: djacs cancel here, since it is really area for a triangle
        area = djac;
        grad2d_f(grad_shp, elem_c, elem_grad_c, nnodes);
        grad2d_f(grad_shp, elem_thickness, elem_grad_thick, nnodes);
        grad2d_v(grad_shp, elem_bedload_vel, elem_grad_u, elem_grad_v, nnodes);
        
        elem_t_avg =      integrate_triangle_f(1.,1.,elem_thickness);
        elem_u_avg =      integrate_triangle_f(1.,1.,u);
        elem_v_avg =      integrate_triangle_f(1.,1.,v);
        elem_source_avg = integrate_triangle_f(1.,1.,elem_source);
        elem_sink_avg =   integrate_triangle_f(1.,1.,elem_sink);
        
        elem_vbed_avg.x = elem_u_avg; elem_vbed_avg.y = elem_v_avg;
        elem_strong_advection_avg = integrate_triangle_grad_dot_f_vbar(grad_shp, 1., elem_t_avg, elem_c, elem_vbed_avg);
    } else {
        area = integrate_quad_constant(elem_nds, 1.);
        elem_t_avg =       integrate_quadrilateral_f(elem_nds, 1./area, elem_thickness);
        elem_u_avg =       integrate_quadrilateral_f(elem_nds, 1./area, u);
        elem_v_avg =       integrate_quadrilateral_f(elem_nds, 1./area, v);
        elem_source_avg =  integrate_quadrilateral_f(elem_nds, 1./area, elem_source);
        elem_sink_avg =    integrate_quadrilateral_f(elem_nds, 1./area, elem_sink);
        elem_vbed_avg.x = elem_u_avg;
        elem_vbed_avg.y = elem_v_avg;
        elem_strong_advection_avg = integrate_quadrilateral_grad_dot_f_vbar(elem_nds, elem_t_avg/area, elem_c, elem_vbed_avg);
    }
    double elem_vmag_avg = svect2d_mag(elem_vbed_avg);
    
    // sets tau for SUPG contributions -- CJT :: why tau vs. tau2 ...
    double tau = 0., tau2 = 0.;
    double radius = sqrt(area);
    if (elem_vmag_avg > SMALL) {
        tau = (mod->tau_pg * radius / elem_vmag_avg);
    }
    if(svect2d_dotp(elem_vbed_avg, elem_vbed_avg) > SMALL) {
        tau2 = (mod->tau_pg * radius) / elem_vmag_avg);
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        printf("BED LOAD 2D ELEM RESID :: ie: %d \t dt: %20.10f \t area: %20.10f",ie,dt,area);
        if (perturb_var == PERTURB_C) {
            printf("\t perturbing C || node: %d || perturbation: %20.10e\n",nodes[perturb_node].id,perturb_sign*perturbation);
        }
        selem2d_printScreen(elem2d);
        
        printf("\n------------------ ELEMENTAL AVERAGES --------------------\n");
        printf("elem_t_avg: %20.10e\n",elem_t_avg);
        printf("elem_u_avg: %20.10e\n",elem_u_avg);
        printf("elem_v_avg: %20.10e\n",elem_v_avg);
        printf("elem_source_avg: %20.10e\n",elem_source_avg);
        printf("elem_sink_avg: %20.10e\n",elem_sink_avg);
        printf("elem_strong_advection_avg: %20.10e\n",elem_strong_advection_avg);
        printf("elem_vbed_avg: {%20.10e, %20.10e}\n",elem_vbed_avg.x, elem_vbed_avg.y);
        
        printf("\n------------------ FUNCTION GRADIENTS --------------------\n");
        printf("elem_grad_c: {%20.10e, %20.10e}\n",elem_grad_c.x,elem_grad_c.y);
        printf("elem_grad_thick: {%20.10e, %20.10e}\n",elem_grad_thick.x,elem_grad_thick.y);
        printf("elem_grad_bedload_u: {%20.10e, %20.10e}\n",elem_grad_u.x,elem_grad_u.y);
        printf("elem_grad_bedload_v: {%20.10e, %20.10e}\n",elem_grad_v.x,elem_grad_v.y);
        
        printf("\n------------------- LINEAR VARIABLES ---------------------\n");
        printScreen_debug2_dbl("elem_head", elem_head, nnodes, elem2d->node);
        printScreen_debug2_dbl("elem_old_head", elem_old_head, nnodes, elem2d->node);
        printScreen_debug2_dbl("elem_older_head", elem_older_head, nnodes, elem2d->node);
        printScreen_debug2_dbl("elem_c", elem_c, nnodes, elem2d->node);
        printScreen_debug2_dbl("elem_old_c", elem_old_c, nnodes, elem2d->node);
        printScreen_debug2_dbl("elem_older_c", elem_older_c, nnodes, elem2d->node);
        printScreen_debug2_dbl("elem_thickness", elem_thickness, nnodes, elem2d->node);
        printScreen_debug2_dbl("elem_old_thickness", elem_old_thickness, nnodes, elem2d->node);
        printScreen_debug2_dbl("elem_older_thickness", elem_older_thickness, nnodes, elem2d->node);
        printScreen_debug2_dbl("elem_source", elem_source, nnodes, elem2d->node);
        printScreen_debug2_dbl("elem_sink", elem_sink, nnodes, elem2d->node);
        printScreen_debug2_dbl("elem_bed_shear", elem_bed_shear, nnodes, elem2d->node);
        printScreen_debug_svect2d("elem_bedload_vel", elem_bedload_vel, nnodes, elem2d->node);
        
        printf("\n--------------------- OTHER STUFF ------------------------\n");
        printf("eddy viscosity tensor: "); stensor2dai_printScreen(d_total);
        printf("supg :: tau: %10.5e \t tau2: %10.5e\n",tau,tau2);

    }
#endif
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    sarray_dbl_init(elem_rhs, nnodes);
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   TEMPORAL CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the temporal addition to the bed load elemental residual. \n
     * \f$ R_{i,body}^e =  \bodyTime{\,2d}{e}{\phidd{i}}{\ctrns{h}{j}{bl} \,\thick{h}} \f$
     * \note GLB \:: 11-13 remove pg time term: it is unstable for some reason
     * \note CJT \:: only uses 1rst order time discretization for now
     *********************************************************************************************/
    
    // ++ t(n+1) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_sw2_bl_temporal(nnodes, elem_nds, djac, elem_thickness, elem_c, +1., elem_rhs, "2D TRANSPORT || TEMPORAL T(N+1): ", DEBUG, DEBUG_LOCAL);
    
    // ++ t(n) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_sw2_bl_temporal(nnodes, elem_nds, djac, elem_old_thickness, elem_old_c, -1., elem_rhs, "2D TRANSPORT || TEMPORAL T(N): ", DEBUG, DEBUG_LOCAL);

    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   ADVECTION CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the advection addition to the bed load elemental residual. \n
     * \f$ R_{i,body}^e =  \bodyConv{\,2d}{e}{\phidd{i}}{ (\velc{h}{bl} \, {\ctrns{h}{j}{bl} \,\thick{h}}) } \f$
     *********************************************************************************************/

    // GALERKIN
    sarray_dbl_init(rhs,nnodes);
    if (isTriangle == YES) {
        integrate_triangle_gradPhi_dot_f_g_v(grad_shp, djac, 1., elem_c, elem_thickness, elem_bedload_vel, rhs);
    } else {
        integrate_quadrilateral_gradPhi_dot_f_g_v(elem_nds, 1., elem_c, elem_thickness, elem_bedload_vel, rhs);
    }
    for (i=0; i<nnodes; i++) {elem_rhs[i] -= dt * rhs[i];}
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        rhs_1dof("2D BED LOAD || ADVECTION: ",nnodes, ie, elem2d->nodes, rhs);
    }
#endif

    // SUPG
    sarray_dbl_init(rhs,nnodes);
    if (isTriangle == TRUE) {
        integrate_triangle_gradPhi_dot_vbar(grad_shp, djac, elem_strong_advection_avg, elem_vbed_avg, rhs);
    } else {
        integrate_quadrilateral_gradPhi_dot_vbar(elem_nds, elem_strong_advection_avg, elem_vbed_avg, rhs);
    }
    for (i=0; i<nnodes; i++) elem_rhs[i] += tau2 * dt * rhs[i];
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        rhs_1dof("2D BED LOAD || SUPG ADVECTION: ",nnodes, ie, elem2d->nodes, rhs);
    }
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                 DIFFUSION CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the diffusion addition to the bed load elemental residual. \n
     * \f$ R_{i,body}^e =  \bodyDiffusion{\,2d}{e}{\phidd{i}}{\diffTensor{bl}{h}{\diffTensorTRN}{j,}} \f$
     * \note CJT \:: uses elementally averaged thickness
     * \note CJT \:: must use quadrature for quadrilateral integration
     *********************************************************************************************/
    
    // GALERKIN
    sarray_dbl_init(rhs,nnodes);
    if (isTriangle == TRUE) {
        SVECT2D diffusive_flux; VT_2D_TENS2D_VECT2D_PROD(diffusive_flux, d, elem_grad_c);
        integrate_triangle_gradPhi_dot_v(grad_shp, djac, thickness_avg, diffusive_flux, rhs);
    } else {
        integrate_quadrilateral_gradPhi_dot_Df(elem_nds, quad, thickness_avg, d, elem_c, rhs);
    }
    for (i=0; i<nnodes; i++) elem_rhs[i] += dt * rhs[i];
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        rhs_1dof("2D BED LOAD || DIFFUSION: ",nnodes, ie, elem2d->nodes, rhs);
    }
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                 SOURCE/SINK CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the advection addition to the bed load elemental residual. \n
     * \f$ R_{i,body}^e =  \bodySource{\,2d}{e}{\phidd{i}}{\srcGen{bl}{h}{j,}} \f$
     * \note GLB \:: 0213 for sediment, use average elemental loading
     *********************************************************************************************/
    
    double *elem_total_source[nnodes];
    double elem_total_source_avg = elem_source_avg + elem_sink_avg;
    for (i=0; i<nnodes; i++) elem_total_source[i] = elem_total_source_avg;

    // GALERKIN
    sarray_dbl_init(rhs,nnodes);
    if (isTriangle == TRUE) {
        integrate_triangle_phi_f_lump(djac, 1., elem_total_source, rhs);
    } else {
        integrate_quadrilateral_phi_f_lump(elem_nds, 1., elem_total_source, rhs);
    }
    for (i=0; i<nnodes; i++) {elem_rhs[i] -= dt * rhs[i];}
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        rhs_1dof("2D BED LOAD || SOURCE/SINK: ",nnodes, ie, elem2d->nodes, rhs);
    }
#endif
    
    // SUPG
    sarray_dbl_init(rhs,nnodes);
    if (isTriangle == TRUE) {
        integrate_triangle_gradPhi_dot_vbar(grad_shp, djac, elem_total_source_avg, elem_vbed_avg, rhs);
    } else {
        integrate_quadrilateral_gradPhi_dot_vbar(elem_nds, elem_total_source_avg, elem_vbed_avg, rhs);
    }
    t1 = tau * dt;
    for (i=0; i<nnodes; i++) {elem_rhs[i] -= t1 * rhs[i];}
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        rhs_1dof("2D BED LOAD || SUPG SOURCE/SINK: ",nnodes, ie, elem2d->nodes, rhs);
    }
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    SET DIRICHLET BCS
     *--------------------------------------------------------------------------------------------
     * Sets the Dirichlet boundary conditions. \n
     *********************************************************************************************/
    int node_istring = UNSET_INT;
    for (i=0; i<nnodes; i++) {
        node_istring = grid->nodes[node_bed_id[i]].string;
        if (node_istring > NORMAL) {
            if (mod->str_values[node_istring].sed[mod->ised].bc_flag == BCT_DIR || mod->str_values[node_istring].sed[mod->ised].bc_flag == BCT_CEQ) {
                elem_rhs[i] = 0.0;
            }
        }
    }

#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        rhs_1dof("2D BED LOAD || TOTAL: ",nnodes, ie, elem2d->nodes, elem_rhs);
    }
#endif
    
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 2D transport body temporal contributions to the shallow water equations.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_sw2_bl_temporal(int nnodes, SVECT *elem_nds, double djac, double *elem_thickness, double *elem_c, double dt_factor, DOF_3 *rhs, char *string, int DEBUG, int DEBUG_LOCAL) {
    
    int i;
    
    sarray_dbl_init(rhs,nnodes);
    
    // CJT :: group these two variables.  Should be fine, especially since thick = 1
    double ch_group[nnodes];
    for (i=0; i<nnodes; i++) ch_group[i] = elem_c[i] * elem_thickness[i];
    
    if (isTriangle == TRUE) {
        integrate_triangle_phi_f_lump(djac, 1., ch_group, rhs);
    } else {
        integrate_quadrilateral_phi_f_lump(elem_nds, 1., ch_group, rhs);
    }
    for (i=0; i<nnodes; i++) {elem_rhs[i] += dt_factor * rhs[i]};
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        rhs_1dof(string, nnodes, ie, elem2d->nodes, rhs);
    }
#endif
}

