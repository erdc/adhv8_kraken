/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_2d_transport_body_resid.c This file collections functions responsible for
 *          the 2D shallow water transport body contributions to the elemental residual.    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

// prototypes
void fe_2d_trns_temporal(SELEM_2D *elem2d, SVECT *elem_nds, double djac, double *elem_head, double *elem_c, double dt_factor, double *elem_rhs, char *string, int DEBUG_LOCAL);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 2D general transport body contributions to the elemental residual.
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
 * \details Solves the following weak, discrete body terms of the diffusive wave equation: \n
 * \f[  R_{i,body}^{\,e} = \bodyTime{\,2d}{e}{\phidd{i}}{\, \ctrns{h}{j}{trns} \, \depth{h}} \,-
 *                     \bodyConv{\,2d}{e}{\phidd{i}}{ (\velc{h}{trns} \, \ctrns{h}{j}{trns} \, \depth{h}) } \,+
 *                     \bodyDiffusion{\,2d}{e}{\phidd{i}}{\diffTensor{trns}{h}{\diffTensorTRN}{j,}} \,-
 *                     \bodySource{\,2d}{e}{\phidd{i}}{\srcGen{trns}{h}{j,}} \,+
 *                     \bodySUPG{\supg{i}{trns}{e}}
 * \f]
 *
 * \note CJT \:: Old velocities used here are NOT consitituent velocities, they are fluid.  This is weird to me.
 *               To fix this, we would need to store old mfcf and cvf factors though.
 * \note CJT \:: wd_flag = 0   // fully wet element
 * \note CJT \::wd_flag = 1,2 // partial dry element
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_2d_transport_body_resid(SMODEL *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    
    
    int i,j;
    
    double t1, t2, t3, t4, t5;
    
    
    int DEBUG_LOCAL = OFF;
    int DEBUG_PICKETS = OFF;
    int DEBUG_NODE_ID = UNSET_INT;
    int DEBUG_NODE_PE = OFF;
#ifdef _DEBUG
    if (DEBUG_NODE_ID != UNSET_INT) {
        for (i=0; i<mod->grid->elem2d[ie].nnodes; i++) {
            if (mod->grid->elem2d[ie].nodes[i] == DEBUG_NODE_ID - 1 && mod->grid->smpi->myid == DEBUG_NODE_PE) {
                DEBUG = ON;
                break;
            } else {
                DEBUG = OFF;
                DEBUG_LOCAL = OFF;
                DEBUG_PICKETS = OFF;
            }
        }
    }
    if (DEBUG == ON) DEBUG_LOCAL = ON;
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    time_t time1;  time(&time1);
#endif
    
    // aliases
    SSW_2D *sw2 = mod->sw->d2;
    double alpha = mod->tau_temporal;
    double dt = mod->dt;
    double g = mod->gravity;
    double drying_lower_limit = mod->drying_lower_limit;
    int imat = mod->grid->elem2d[ie].mat;
    STR_VALUE str_values = mod->str_values[mod->grid->elem2d[ie].string];
    SMAT_TRN mat_trans;
    int index;
    double *c = NULL, *c_old = NULL, *c_older = NULL, *source = NULL, *sink = NULL, *mfcf = NULL;
    SVECT2D *vcf = NULL;
    
    if (mod->is_sediment_running) {
#ifdef _SEDIMENT
        SSUSLOAD *susload = NULL;
        index = mod->ised;
        susload = &(mod->sed->susload[index]);
        mat_trans = mod->mat[imat].sed[index];
        c = susload->c;
        c_old = susload->c_old;
        c_older = susload->c_older;
        source = susload->source;
        sink = susload->sink;
        mfcf = susload->mfcf;
        vcf = susload->vcf;
#endif
    } else {
        SCON *con = NULL;
        index = mod->itrns;
        con = &(mod->con[index]);
        mat_trans = mod->mat[imat].trn[index];
        c = con->concentration;
        c_old = con->old_concentration;
        c_older = con->older_concentration;
        source = con->source;
        sink = con->sink;
        mfcf = con->mfcf;
        vcf = con->vcf;
    }
    
#ifdef _DEBUG
        assert(perturb_node < mod->grid->elem2d[ie].nnodes);
        assert(alpha > -1E-6 && alpha <= 1.0);
        assert(imat >= 0);
        assert(perturb_var == PERTURB_NONE || perturb_var == PERTURB_C);
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GRID VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    SGRID *grid = mod->grid;
    int wd_flag = grid->wd_flag[ie];
    SQUAD *quad = grid->quad_rect; // for quadrilateral quadrature calculations
    SELEM_2D *elem2d = &(grid->elem2d[ie]); // be careful not to shallow copy here
    int nnodes = elem2d->nnodes;
    
    double djac = 0.;
    SVECT2D *grad_shp = NULL;
    int isElementTriangle = NO;
    if (nnodes == NDONTRI) {
        isElementTriangle = TRUE;
        djac = elem2d->djac;
        grad_shp = elem2d->grad_shp;
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            assert(djac > SMALL);
        }
#endif
    }
    
    SVECT elem_nds[nnodes];
    SNODE elem_nodes[nnodes];
    for (i=0; i<nnodes; i++) {
        elem_nds[i].x = grid->node[elem2d->nodes[i]].x;
        elem_nds[i].y = grid->node[elem2d->nodes[i]].y;
        elem_nds[i].z = grid->node[elem2d->nodes[i]].z;
        snode_copy(&(elem_nodes[i]), grid->node[elem2d->nodes[i]]);
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    double elem_c[nnodes];
    global_to_local_dbl(c, elem_c, elem2d->nodes, nnodes);
    if (perturb_var == PERTURB_C) elem_c[perturb_node] += perturb_sign * perturbation;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // constituent concentrations
    double elem_c_old[nnodes], elem_c_older[nnodes];
    global_to_local_dbl(c_old, elem_c_old, elem2d->nodes, nnodes);
    global_to_local_dbl(c_older, elem_c_older, elem2d->nodes, nnodes);
    
    // water depths
    double elem_head[nnodes], elem_old_head[nnodes], elem_older_head[nnodes];
    global_to_local_dbl(sw2->head, elem_head, elem2d->nodes, nnodes);
    global_to_local_dbl(sw2->old_head, elem_old_head, elem2d->nodes, nnodes);
    global_to_local_dbl(sw2->older_head, elem_older_head, elem2d->nodes, nnodes);
    
    // water velocities
    SVECT2D elem_vel[nnodes], elem_old_vel[nnodes];
    global_to_local_svect2d(sw2->old_vel, elem_old_vel, elem2d->nodes, nnodes);
    global_to_local_svect2d(sw2->vel, elem_vel, elem2d->nodes, nnodes);
    
    // sediment velocity
    SVECT2D elem_vel_con[nnodes];
    for (i=0; i<nnodes; i++) {
        j = elem2d->nodes[i];
        elem_vel_con[i].x = (elem_vel[i].x + vcf[j].x) * mfcf[j];
        elem_vel_con[i].y = (elem_vel[i].y + vcf[j].y) * mfcf[j];
    }
    double elem_u_con[nnodes], elem_v_con[nnodes];
    dumpVector2D(elem_vel_con, nnodes, elem_u_con, elem_v_con);
    
    // sources/sinks
    double elem_source[nnodes], elem_sink[nnodes];
    global_to_local_dbl(source, elem_source, elem2d->nodes, nnodes);
    global_to_local_dbl(sink, elem_sink, elem2d->nodes, nnodes);
    
    // non-wet/dry elemental averages, area and function/vector gradients
    double area = 0., elem_h_avg = 0., elem_h_old_avg = 0., elem_source_avg = 0., elem_sink_avg = 0., elem_c_avg = 0.;
    double elem_u_con_avg = 0., elem_v_con_avg = 0.;
    SVECT2D elem_grad_u_con, elem_grad_v_con, elem_grad_c, elem_vel_con_avg;
    svect2d_init(&elem_grad_u_con); svect2d_init(&elem_grad_v_con); svect2d_init(&elem_grad_c);
    if (isElementTriangle == TRUE) {
        area = djac;
        
        // elemental averages - djacs cancel here
        elem_c_avg =      integrate_triangle_f(1.,1.,elem_c);
        elem_h_avg =      integrate_triangle_f(1.,1.,elem_head);
        elem_h_old_avg =  integrate_triangle_f(1.,1.,elem_old_head);
        elem_u_con_avg =  integrate_triangle_f(1.,1.,elem_u_con);
        elem_v_con_avg =  integrate_triangle_f(1.,1.,elem_v_con);
        elem_source_avg = integrate_triangle_f(1.,1.,elem_source);
        elem_sink_avg =   integrate_triangle_f(1.,1.,elem_sink);
        
        
        // function gradients
        grad2d_phi_f(grad_shp, elem_c, &elem_grad_c, nnodes);
        grad2d_phi_dot_v(grad_shp, elem_vel_con, &elem_grad_u_con, &elem_grad_v_con, nnodes);
        
    } else {
        area = integrate_quadrilateral_area(elem_nds, 1.);
        
        // elemental averages
        elem_c_avg =      integrate_quadrilateral_f(elem_nds,1./area,elem_c);
        elem_h_avg =      integrate_quadrilateral_f(elem_nds,1./area,elem_head);
        elem_h_old_avg =  integrate_quadrilateral_f(elem_nds,1./area,elem_old_head);
        elem_u_con_avg =  integrate_quadrilateral_f(elem_nds,1./area,elem_u_con);
        elem_v_con_avg =  integrate_quadrilateral_f(elem_nds,1./area,elem_v_con);
        elem_source_avg = integrate_quadrilateral_f(elem_nds,1./area,elem_source);
        elem_sink_avg =   integrate_quadrilateral_f(elem_nds,1./area,elem_sink);
    }
    elem_vel_con_avg.x = elem_u_con_avg;
    elem_vel_con_avg.y = elem_v_con_avg;
    double elem_vel_con_avg_mag = svect2d_mag(elem_vel_con_avg);
    
    // wet/dry factors and elemental averages
    double factor = 1., factor_old = 1., factor_older = 1.;
    double elem_h_wd_avg = elem_h_avg, elem_c_wd_avg = elem_c_avg;
    SVECT2D elem_vel_con_wd_avg = elem_vel_con_avg;
    if (isElementTriangle == TRUE && wd_flag != 0) {
        fe_sw2_wd_average(elem_nds, elem_head, elem_vel_con, elem_head, djac, &elem_h_wd_avg, &elem_vel_con_wd_avg);
        //fe_sw2_wd_average(elem_nds, elem_old_head, elem_vel_old, elem_old_head, djac, &elem_h_old_wd_avg, &elem_old_vel_wd_avg);
        fe_sw2_wd_average(elem_nds, elem_head, NULL, elem_c, djac, &elem_c_wd_avg, NULL);
        
        factor = fe_sw2_wet_dry_factor(elem_nds, elem_head, djac);
        factor_old = fe_sw2_wet_dry_factor(elem_nds, elem_old_head, djac);
        factor_older = fe_sw2_wet_dry_factor(elem_nds, elem_older_head, djac);
    }
    double elem_vel_con_wd_avg_mag = svect2d_mag(elem_vel_con_wd_avg);
    //SVECT2D elem_old_vel_wd_avg_mag = svect2d_mag(elem_old_vel_wd_avg);
    
    // friction variables
    double h_fric = MAX(elem_h_avg, SMALL);
    double roughness = fe_sw2_get_roughness(mod, h_fric, elem_vel_con_wd_avg_mag, elem2d->string, UNUSED);
    
    // loads the diffusion material coefficients :: GLB dxx and dyy have to be zero because of the way we have defined anisotripic diffusion
    STENSOR2D d; stensor2d_init(&d);
    STENSOR2D dvor; stensor2d_init(&dvor);
    STENSOR2D d_total; stensor2d_init(&d_total);
    double ev_st = 0., ev_tr = 0.;
    d.xx = 0.0;
    d.yy = 0.0;
    d.xy = mat_trans.d_m;
    
    if (mod->mat[imat].sw->EEVF == TRUE) {
        // use estimated eddy viscosity algorithm
        if (isElementTriangle == TRUE) {
            fe_sw2_get_EEVF(mod->mat[imat].sw->eev_mode, mod->mat[imat].sw->eev_coef, g, drying_lower_limit, h_fric, roughness, area, elem_grad_u_con, elem_grad_v_con, elem_vel_con_avg_mag, wd_flag, &ev_st, &ev_tr);
        } else {
            // use elementally averaged gradients here
            SVECT2D elem_grad_u_con_avg = integrate_quadrilateral_gradF(elem_nds, 1./area, elem_u_con);
            SVECT2D elem_grad_v_con_avg = integrate_quadrilateral_gradF(elem_nds, 1./area, elem_v_con);
            
            fe_sw2_get_EEVF(mod->mat[imat].sw->eev_mode, mod->mat[imat].sw->eev_coef, g, drying_lower_limit, h_fric, roughness, area, elem_grad_u_con_avg, elem_grad_v_con_avg, elem_vel_con_avg_mag, wd_flag, &ev_st, &ev_tr);
        }
        d.xx = ev_st;
        d.xy = ev_tr;
        d.yy = ev_st;
    }
    
    // non-wetting and drying velocities and unit normals
    SVECT2D vdir;
    vdir.x = elem_u_con_avg / MAX(elem_vel_con_avg_mag, 1.E-8);
    vdir.y = elem_v_con_avg / MAX(elem_vel_con_avg_mag, 1.E-8);
    
    if (mod->flag.VORTICITY == ON && isElementTriangle == TRUE) {
        SVECT2D vor_dir;
        vor_dir.x = -vdir.y;
        vor_dir.y =  vdir.x;
        double vor_c = one_3 * (  mod->con[mod->vorticity_id].old_concentration[elem2d->nodes[0]]
                                + mod->con[mod->vorticity_id].old_concentration[elem2d->nodes[1]]
                                + mod->con[mod->vorticity_id].old_concentration[elem2d->nodes[2]]);
        double vor_nbv = 6. * vor_c * elem_h_wd_avg / sqrt (6. * mod->con[mod->vorticity_id].property[1]);
        vor_nbv = fabs(vor_nbv);
        vor_nbv = MIN(vor_nbv,(2.0 * elem_vel_con_wd_avg_mag));
        dvor.xx = 0.5 * vor_nbv * elem_h_wd_avg;
        dvor.yy = 0.5 * vor_nbv * elem_h_wd_avg;
        dvor.xy = 0.0;
    }
    
    // total diffusion diffusion
    d_total = d; // shallow copy ok here
    //if (mod->vorticity_id > UNSET_INT) {
    //    void stensor2d_add_replace(&d_total, dvor);
    //}
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        printf("TRANSPORT 2D ELEM RESID :: ie: %d \t dt: %20.10f \t area: %20.10f",ie,dt,area);
        if (perturb_var == PERTURB_C) {
            printf("\t perturbing C || node: %d || perturbation: %20.10e\n",elem2d->nodes[perturb_node],perturb_sign*perturbation);
        }
        selem2d_printScreen(elem2d);
        
        printf("\n------------------ ELEMENTAL AVERAGES --------------------\n");
        printf("elem_c_avg: %20.10e  elem_c_wd_avg: %20.10e \n",elem_c_avg,elem_c_wd_avg);
        printf("elem_h_avg: %20.10e  elem_h_wd_avg: %20.10e \n",elem_h_avg,elem_h_wd_avg);
        //printf("elem_h_old_avg: %20.10e  elem_h_wd_avg: %20.10e \n",elem_h_old_avg,elem_h_old_wd_avg);
        printf("elem_u_con_avg: %20.10e  elem_u_con_wd_avg: %20.10e \n",elem_vel_con_avg.x,elem_vel_con_wd_avg.x);
        printf("elem_v_con_avg: %20.10e  elem_v_con_wd_avg: %20.10e \n",elem_vel_con_avg.y,elem_vel_con_wd_avg.y);
        //printf("elem_old_u_avg: %20.10e  elem_old_u_wd_avg: %20.10e \n",elem_old_vel_avg.x,elem_old_vel_wd_avg.x);
        //printf("elem_old_v_avg: %20.10e  elem_old_v_wd_avg: %20.10e \n",elem_old_vel_avg.y,elem_old_vel_wd_avg.y);
        printf("elem_source_avg: %20.10e \t elem_sink_avg: %20.10e \n",elem_source_avg, elem_sink_avg);
        printf("elem_vel_con_avg_mag: %20.10e \t elem_vel_con_wd_avg_mag: %20.10e \n",elem_vel_con_avg_mag,elem_vel_con_wd_avg_mag);
        //printf("elem_old_vel_avg_mag: %20.10e \t elem_old_vel_wd_avg_mag: %20.10e \n",elem_old_vel_avg_mag,elem_old_vel_wd_avg_mag);
        
        printf("\n------------------ FUNCTION GRADIENTS --------------------\n");
        printf("elem_grad_c: {%20.10e, %20.10e}\n",elem_grad_c.x,elem_grad_c.y);
        printf("elem_grad_vel.x: {%20.10e, %20.10e}\n",elem_grad_u_con.x,elem_grad_u_con.y);
        printf("elem_grad_vel.y: {%20.10e, %20.10e}\n",elem_grad_v_con.x,elem_grad_v_con.y);
        
        printf("\n------------------- LINEAR VARIABLES ---------------------\n");
        printScreen_debug_vec("node locations: ",elem_nds, nnodes);
        printScreen_debug2_dbl("elem_head", elem_head, nnodes, elem2d->nodes);
        printScreen_debug2_dbl("elem_old_head", elem_old_head, nnodes, elem2d->nodes);
        printScreen_debug2_dbl("elem_older_head", elem_older_head, nnodes, elem2d->nodes);
        printScreen_debug2_dbl("elem_c", elem_c, nnodes, elem2d->nodes);
        printScreen_debug2_dbl("elem_c_old", elem_c_old, nnodes, elem2d->nodes);
        printScreen_debug2_dbl("elem_c_older", elem_c_older, nnodes, elem2d->nodes);
        printScreen_debug_svec2d("elem_vel_con", elem_vel_con, nnodes, elem2d->nodes);
        printScreen_debug_svec2d("elem_vel", elem_vel, nnodes, elem2d->nodes);
        printScreen_debug_svec2d("elem_old_vel", elem_old_vel, nnodes, elem2d->nodes);
        printScreen_debug2_dbl("elem_source", elem_source, nnodes, elem2d->nodes);
        printScreen_debug2_dbl("elem_sink", elem_sink, nnodes, elem2d->nodes);
        
        printf("\n--------------------- OTHER STUFF ------------------------\n");
        printf("wet-dry factors: current: %20.10f \t old: %20.10f \t older: %20.10f\n",factor,factor_old,factor_older);
        //printf("eddy viscosity tensor: "); stensor2dai_printScreen(d_total);
    }
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                  FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    
    double rhs[nnodes];
    DOF_3 dof3_rhs[nnodes]; // for wet/dry wrapper interface
    sarray_init_dbl(elem_rhs, nnodes);
    sarray_init_dbl(rhs, nnodes);
    double vars[2]; // for passing doubles into wet-dry algorithms
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   TEMPORAL CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the temporal addition to the bed load elemental residual. \n
     * \f$ R_{i,body}^{\,e} =  \bodyTime{\,2d}{e}{\phidd{i}}{\ctrns{h}{j}{} \,\depth{h}} \f$
     * \note
     *********************************************************************************************/
    
    // ++ t(n+1) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_2d_trns_temporal(elem2d, elem_nds, djac, elem_head, elem_c, (1. + alpha / 2.), elem_rhs, "2D TRANSPORT || TEMPORAL T(N+1): ", DEBUG_LOCAL);
#ifdef _DEBUG
    Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
#endif
    
    // ++ t(n) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_2d_trns_temporal(elem2d, elem_nds, djac, elem_old_head, elem_c_old, (-1.) * (1. + alpha), elem_rhs, "2D TRANSPORT || TEMPORAL T(N): ", DEBUG_LOCAL);
#ifdef _DEBUG
    Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
#endif
    
    // ++ t(n-1) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_2d_trns_temporal(elem2d, elem_nds, djac, elem_older_head, elem_c_older, (alpha / 2.), elem_rhs, "2D TRANSPORT || TEMPORAL T(N-1): ", DEBUG_LOCAL);
#ifdef _DEBUG
    Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                     ADVECTION CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the advection addition to the 2d transport elemental residual. \n
     * \f$ R_{i,body}^{\,e} = \bodyConv{\,2d}{e}{\phidd{i}}{ (\velc{h}{trns} \, \ctrns{h}{j}{trns} \, \depth{h}) } \f$
     * \note Triangles wrap in wet/dry integration
     * \note Quadrilaterals use full depth
     *********************************************************************************************/
    
    sarray_init_dbl(rhs,nnodes);
    if (isElementTriangle == TRUE) {
        dof3_init_array(dof3_rhs,nnodes);
        vars[0] = dt;
        vars[1] = -1;
        t1 = fe_sw2_wet_dry_wrapper(dof3_rhs, elem_nds, elem_head, grad_shp, NULL, elem_vel_con, NULL, elem_c, djac, ON, DEBUG, vars, fe_2d_transport_wd_convection_triangle);
        for (i=0; i<nnodes; i++) {rhs[i] = dof3_rhs[i].c_eq;}
    } else {
        integrate_quadrilateral_gradPhi_dot_f_g_v(elem_nds, -dt, elem_c, elem_head, elem_vel_con, rhs);
    }
    for (i=0; i<nnodes; i++) {elem_rhs[i] += rhs[i];}
#ifdef _DEBUG
    Is_DoubleArray_Inf_or_NaN(rhs ,nnodes ,__FILE__ ,__LINE__);
    if (DEBUG_LOCAL == ON) {
        rhs_1dof("2D TRANSPORT || ADVECTION: ",nnodes, ie, elem2d->nodes, rhs);
    }
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                 DIFFUSION CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the diffusion addition to the 2d transport elemental residual. \n
     * \f$ R_{i,body}^{\,e} =  \bodyDiffusion{\,2d}{e}{\phidd{i}}{\diffTensor{}{h}{\diffTensorTRN}{j,}} \f$
     * \note CJT \:: uses average depth
     * \note SUPG does not contribute here, since first order elements are used
     *********************************************************************************************/
    
    SVECT2D diffusive_flux; svect2d_init(&diffusive_flux);
    double factor_diff = MAX(factor,.01);
    double depth =  MAX(elem_h_wd_avg , SMALL);
    double depth_diff = MAX(depth, drying_lower_limit * .01);
    
    sarray_init_dbl(rhs,nnodes);
    if (isElementTriangle == TRUE) {
        VT_2D_TENS2D_VECT2D_PROD_AI(diffusive_flux, d_total, elem_grad_c, vdir);
        integrate_triangle_gradPhi_dot_vbar(grad_shp, dt * depth_diff * factor_diff, djac, diffusive_flux, rhs);
    } else {
        integrate_quadrilateral_gradPhi_dot_Df(elem_nds, quad, dt * depth_diff * factor_diff, d_total, elem_c, rhs);
    }
    for (i=0; i<nnodes; i++) {elem_rhs[i] += rhs[i];}
#ifdef _DEBUG
    Is_DoubleArray_Inf_or_NaN(rhs ,nnodes ,__FILE__ ,__LINE__);
    if (DEBUG_LOCAL == ON) {
        printf("djac: %20.10f dt: %20.10f factor_diff: %20.10f depth: %20.10f depth_diff: %20.10f\n",djac,dt,factor_diff,depth,depth_diff);
        printf("diffusion tensor: xx: %20.10f yy: %20.10f xy: %20.10f \n",d_total.xx, d_total.yy, d_total.xy);
        printf("diffusive_flux: %20.10f %20.10f\n",diffusive_flux.x,diffusive_flux.y);
        printf("elem_grad_c: %20.10f %20.10f\n",elem_grad_c.x,elem_grad_c.y);
        printf("vdir: %20.10f %20.10f\n",vdir.x,vdir.y);
        rhs_1dof("2D BED LOAD || DIFFUSION: ",nnodes, ie, elem2d->nodes, rhs);
        
    }
#endif
    
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                 SOURCE/SINK CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the advection addition to the 2d transport  elemental residual. \n
     * \f$ R_{i,body}^{\,e} =  \bodySource{\,2d}{e}{\phidd{i}}{\srcGen{}{h}{j,}} \f$
     * \note GLB \:: 0213 for sediment, use average elemental loading
     *********************************************************************************************/
    
    double elem_avg_total_source = 0.;
    if (mod->is_sediment_running) {
        
        // combine sources and sinks
        elem_avg_total_source = elem_source_avg + elem_sink_avg;
        
        // GALERKIN
        sarray_init_dbl(rhs,nnodes);
        if (isElementTriangle == TRUE) {
            integrate_triangle_phi(djac, -dt * factor * elem_avg_total_source, rhs);
        } else {
            integrate_quadrilateral_phi(elem_nds, -dt * factor * elem_avg_total_source, rhs);
        }
        for (i=0; i<nnodes; i++) {elem_rhs[i] += rhs[i];}
#ifdef _DEBUG
        Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("2D TRANSPORT || SOURCE/SINK: ",nnodes, ie, elem2d->nodes, rhs);
            Is_DoubleArray_Inf_or_NaN(rhs ,nnodes ,__FILE__ ,__LINE__);
        }
#endif
    }
    
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                     VORTICITY CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the vorticity addition to the 2d transport elemental residual. \n
     *
     * \note CJT \:: vorticity source is treated as constant over the element
     * \note CJT \:: vorticity sink is lumped
     *********************************************************************************************/
    
    double xdestruct = 0., ydestruct = 0., elem_avg_vorticity_source = 0.;
    if (mod->is_sediment_running == FALSE && index == mod->vorticity_id) {
        
        // 1 indicates that transport is calling
        SVECT2D prod_destruct = tl_bendway_correction(grad_shp, elem_old_vel, elem_old_head, elem_c, mod->con[mod->vorticity_id].property[1],
                                                      mod->con[mod->vorticity_id].property[2], drying_lower_limit, roughness, mod->density, 1, ie);
        xdestruct += prod_destruct.x;
        ydestruct += prod_destruct.y;
        if (elem_head[0] < drying_lower_limit || elem_head[1] < drying_lower_limit || elem_head[2] < drying_lower_limit) {
            xdestruct = 0.;
            ydestruct = -0.1 * sqrt(g * h_fric)/h_fric;
        }
        
        double elem_decay[nnodes];
        for (i=0; i<nnodes; i++) {elem_decay[i] = -ydestruct * elem_c[i];}
        
        // GALERKIN SOURCE
        sarray_init_dbl(rhs,nnodes);
        t1 = -xdestruct * elem_h_wd_avg * dt * factor;
        if (isElementTriangle == TRUE) {
            integrate_triangle_phi(djac, t1, rhs);
        } else {
            integrate_quadrilateral_phi(elem_nds, t1, rhs);
        }
        for (i=0; i<nnodes; i++) {elem_rhs[i] += rhs[i];}
#ifdef _DEBUG
        Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("2D TRANSPORT || VORTICITY SOURCE: ",nnodes, ie, elem2d->nodes, rhs);
        }
#endif
        
        // GALERKIN SINK
        sarray_init_dbl(rhs,nnodes);
        t1 = elem_h_wd_avg * dt * factor;
        if (isElementTriangle == TRUE) {
            integrate_triangle_phi_f_lump(djac, t1, elem_decay, rhs);
        } else {
            integrate_quadrilateral_phi_f_lump(elem_nds, t1, elem_decay, rhs);
        }
        
        for (i=0; i<nnodes; i++) {elem_rhs[i] += rhs[i];}
#ifdef _DEBUG
        Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("2D TRANSPORT || VORTICITY SINK: ",nnodes, ie, elem2d->nodes, rhs);
        }
#endif
        
        // calculate elemental average values for GLS contributions and combine for one integration call
        double elem_avg_xdestruction = -xdestruct; // since it is elementally constant
        double elem_avg_decay = 0.;
        if (isElementTriangle == TRUE) {
            elem_avg_decay = integrate_triangle_f(1.,1.,elem_decay); // cjt :: djacs cancel here
        } else {
            elem_avg_decay = integrate_quadrilateral_f(elem_nds, 1./area, elem_decay);
        }
        elem_avg_vorticity_source = elem_h_wd_avg * (elem_avg_xdestruction + elem_avg_decay);
        
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    SUPG CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the streamline upwind Petrov-Galerkin addition to the 2D transport elemental residual. \n
     * \note CJT \:: temporal term has been omitted for some reason
     *********************************************************************************************/
    
    // the SUPG tau parameter
    double grad_shp_x[nnodes], grad_shp_y[nnodes];
    for (i=0; i<nnodes; i++) {
        grad_shp_x[i] = grad_shp[i].x;
        grad_shp_y[i] = grad_shp[i].y;
    }
    int tau_method_flag = 2; // OLD ADH WAY (probably best to use method 1)
    int le_method_flag = 2; // method for calculating elemental length :: 2 uses area
    double total_diffusion = 0; /// not used in tau_method_flag = 2
    double tau = fe_get_supg_tau(nnodes, elem_nds, total_diffusion, elem_u_con_avg, elem_v_con_avg, 0.,
                                 grad_shp_x, grad_shp_y, NULL, area, mod->tau_pg, 2, tau_method_flag, le_method_flag);
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        printf("tau: %20.10f\n",tau);
    }
#endif
    
    // get elementally averaged strong transport residual
    double elem_avg_temporal = 0., elem_avg_convection = 0., elem_avg_source = 0.;
    if (isElementTriangle == TRUE) {
        elem_avg_temporal = 0.;
        elem_avg_convection = elem_h_wd_avg * svect2d_dotp(elem_vel_con_wd_avg,elem_grad_c);;
        elem_avg_source = elem_avg_total_source + elem_avg_vorticity_source;
        
    } else {
        elem_avg_temporal = 0.;
        elem_avg_convection = integrate_quadrilateral_grad_dot_h_g_v(elem_nds, 1./area, elem_vel_con, elem_c, elem_head);
        elem_avg_source = elem_avg_total_source + elem_avg_vorticity_source;
    }
    double elem_avg_strong_transport_resid = elem_avg_temporal + elem_avg_convection + factor * elem_avg_source;
    
    sarray_init_dbl(rhs,nnodes);
    if (isElementTriangle == TRUE) {
        integrate_triangle_gradPhi_dot_vbar(grad_shp, tau * dt * elem_avg_strong_transport_resid, djac, elem_vel_con_wd_avg, rhs);
    } else {
        integrate_quadrilateral_gradPhi_dot_vbar(elem_nds, tau * dt * elem_avg_strong_transport_resid, elem_vel_con_wd_avg, rhs);
    }
    for (i=0; i<nnodes; i++) {elem_rhs[i] += rhs[i];}
#ifdef _DEBUG
    Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
    if (DEBUG_LOCAL == ON) {
        rhs_1dof("2D TRANSPORT SUPG: ",nnodes, ie, elem2d->nodes, rhs);
    }
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                               EXTRA SW CONTINUITY TERM CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Adds any extra 2D SW terms (SUPG, shock capturing, etc) to the transport elemental residual. \n
     * note: \:: should this use fully wet or wet-dry elementally averaged c?
     *********************************************************************************************/
    
    if (debug.no_hydro == OFF) {
        sarray_init_dbl(rhs,nnodes);
        for (i=0; i<nnodes; i++) {
            rhs[i] = elem_c_wd_avg * mod->sw->d2->elem_rhs_dacont_extra_terms[i][ie];
            elem_rhs[i] += rhs[i];
        }
#ifdef _DEBUG
        Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
        if (DEBUG_LOCAL == ON) {
            printf("elem_c_wd_avg: %20.10f\n",elem_c_wd_avg);
            rhs_1dof("2D TRANSPORT || SW CONTINUITY: ",nnodes, ie, elem2d->nodes, rhs);
        }
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    SET DIRICHLET BCS
     *--------------------------------------------------------------------------------------------
     * Sets the Dirichlet boundary conditions. \n
     *********************************************************************************************/
    int istring = UNSET_INT;
    INPUT_DATA *str_value;
    for (i=0; i<nnodes; i++) {
        istring = grid->node[elem2d->nodes[i]].string;
        if (mod->is_sediment_running) {
            str_value = &(mod->str_values[istring].sed[mod->ised]);
        } else {
            str_value = &(mod->str_values[istring].trans[mod->itrns]);
        }
        if (istring > NORMAL) {
            if (str_value->bc_flag == BCT_DIR || str_value->bc_flag == BCT_CEQ) {
                elem_rhs[i] = 0.0;
            }
        }
    }
#ifdef _DEBUG
    Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
    if (DEBUG_LOCAL == ON) {
        rhs_1dof("2D TRANSPORT || TOTAL RHS AFER DB: ",nnodes, ie, elem2d->nodes, elem_rhs);
        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    }
    
    time_t time2;  time(&time2);
    TIME_IN_2D_TRANSPORT_BODY_RESID += difftime(time2,time1);
#endif
    
    return;
    
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

void fe_2d_trns_temporal(SELEM_2D *elem2d, SVECT *elem_nds, double djac, double *elem_head, double *elem_c, double dt_factor, double *elem_rhs, char *string, int DEBUG_LOCAL) {
    
    int i;
    int nnodes = elem2d->nnodes;
    
    double rhs[nnodes]; sarray_init_dbl(rhs, nnodes);
    DOF_3 dof3_rhs[nnodes]; dof3_init_array(dof3_rhs,NDONTRI);
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        assert(elem2d != NULL);
        assert(elem_nds != NULL);
        assert(djac > 0);
        assert(nnodes > 0 && nnodes < NDONQUAD);
        tl_check_all_pickets(__FILE__,__LINE__);
    }
#endif
    
    double vars[2];
    vars[0] = 1; // dt
    vars[1] = 1; // c
    
    if (nnodes == NDONTRI) {
        double t1 = fe_sw2_wet_dry_wrapper(dof3_rhs, elem_nds, elem_head, NULL, NULL, NULL, NULL, elem_c, djac, ON, DEBUG_LOCAL, vars, fe_2d_transport_wd_temporal_triangle);
        for (i=0; i<nnodes; i++) {rhs[i] = dof3_rhs[i].c_eq;}
    } else {
        integrate_quadrilateral_phi_f_g(elem_nds, 1., elem_c, elem_head, rhs);
    }
    for (i=0; i<nnodes; i++) {elem_rhs[i] += dt_factor * rhs[i];}
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        rhs_1dof(string, nnodes, elem2d->id, elem2d->nodes, rhs);
        Is_DoubleArray_Inf_or_NaN(rhs,nnodes,__FILE__ ,__LINE__);
    }
#endif
}
