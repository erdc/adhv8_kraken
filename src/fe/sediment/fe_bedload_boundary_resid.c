/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns 1D bed load elemental residual.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] elem_rhs      the 1D elemental residual array
 * @param[in]  mod           a pointer to the model struct
 * @param[in]  ie            the elemental id
 * @param[in]  pertubation   the Newton pertubation
 * @param[in]  perturb_node  the node to be perturbed
 * @param[in]  perturb_var   the variable to be perturbed
 * @param[in]  perturb_sign  the direction of Newton perturbation
 * @param[in]  DEBUG         a debug flag
 *
 *  \details Integrates the weak, discrete outflow boundary terms: \n
 *  \f{eqnarray*}{
 *   R_{i,boundary}^e &=& dt * \bigg( \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \ctbjh)} \bigg)
 *  \f}
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_bedload_boundary_resid(SMODEL *mod, double *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
#ifdef _SEDIMENT
    
    int i, DEBUG_LOCAL = OFF;

    // aliases
    SSED *sed = mod->sed;
    int ised = mod->ised;
    SGRID *grid = mod->grid;
    SELEM_1D *elem1d = &(grid->elem1d[ie]); // be careful not to shallow copy here
    int nnodes = elem1d->nnodes;
    double dt = mod->dt;
    double djac = elem1d->djac;
    double gravity = mod->gravity;
    STR_VALUE *str_values = mod->str_values;
    int string = elem1d->string;
    int nnodes2d = grid->elem2d[elem1d->elem2d].nnodes;  // the total number of nodes on the 2D element
    
#ifdef _DEBUG
    assert(nnodes == NDONSEG); // only 1D element is a line segment
    assert(djac > 0); // only 1D element is a line segment
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // ENSURE 3D CALL COMPATIBILITY
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    int node_bed_id[nnodes2d];
    if (grid->ndim == 3) {
        for (i; i<nnodes; i++) {
            node_bed_id[i] = grid->nodeID_3d_to_2d_bed[elem1d->nodes[i]];
        }
    } else {
        for (i=0; i<nnodes; i++) {
            node_bed_id[i] = elem->elem1d->nodes[i];
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

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    double elem_c[nnodes];
    global_to_local_dbl(sed->bedload[ised].c, elem_c, node_bed_id, nnodes);
    if (perturb_var == PERTURB_C) elem_c += perturb_sign * pertubation;
    
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
    double elem_c_avg = one_2 * (elem->c[0] + elem->c[1]);
    double elem_c_old_avg = one_2 * (elem->old_c[0] + elem->old_c[1]);
    
    // thickness
    double elem_thickness[nnodes], elem_old_thickness[nnodes], elem_older_thickness[nnodes];
    global_to_local_dbl(sed->bedload[ised].thick, elem_thickness, node_bed_id, nnodes);
    global_to_local_dbl(sed->bedload[ised].old_thick, elem_old_thickness, node_bed_id, nnodes);
    global_to_local_dbl(sed->bedload[ised].older_thick, elem_older_thickness, node_bed_id, nnodes);
    double elem_t_avg = one_2 * (elem->thickness[0] + elem->thickness[1]);
    
    // sources and sinks
    double elem_source[nnodes], elem_sink[nnodes];
    global_to_local_dbl(sed->bedload[ised].source, elem_source, node_bed_id, nnodes);
    global_to_local_dbl(sed->bedload[ised].sink, elem_sink, node_bed_id, nnodes);
    double elem_source_avg = one_2 * (elem->source[0] + elem_source[1]);
    double elem_sink_avg = one_2 * (elem->sink[0] + elem->sink[1]);
    
    // bed velocity
    SVECT2D elem_bedload_vel[nnodes];
    global_to_local_svect2d(sed->bedload[ised].v, elem_bedload_vel, node_bed_id, nnodes);
    
    // bed velocity normal to 1d element
    double elem_bedvel_nrml[nnodes];
    elem_bedvel_nrml[0] = VT_2D_DOT(elem_bedload_vel[0], elem1d->nrml);
    elem_bedvel_nrml[1] = VT_2D_DOT(elem_bedload_vel[1], elem1d->nrml);

    // bed shear (not currently used?)
    double elem_bed_shear[nnodes];
    global_to_local_dbl(sed->bedload[ised].shear, elem_bed_shear, node_bed_id, nnodes);
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        printf("BED LOAD 1D ELEM RESID :: ie: %d \t dt: %20.10f \t djac: %20.10f \n",ie,dt,djac);
        if (perturb_var == PERTURB_C) {
            printf("\t perturbing C || node: %d || perturbation: %20.10e\n",nodes[perturb_node].id,perturb_sign*perturbation);
        }
        printf("elemental averages:  c: %10.5e \t c_old: %10.5e \t thickness: %10.5e \t source: %10.5e \t sink: %10.5e \n",
               elem_c_avg, elem_c_old_avg, elem_t_avg, elem_source_avg, elem_sink_avg);
        printf("\n--------------------------------------------------------- \n");
        printScreen_debug2_dbl("elem_head", elem_head, nnodes, elem1d->node);
        printScreen_debug2_dbl("elem_old_head", elem_old_head, nnodes, elem1d->node);
        printScreen_debug2_dbl("elem_older_head", elem_older_head, nnodes, elem1d->node);
        printScreen_debug2_dbl("elem_c", elem_c, nnodes, elem1d->node);
        printScreen_debug2_dbl("elem_old_c", elem_old_c, nnodes, elem1d->node);
        printScreen_debug2_dbl("elem_older_c", elem_older_c, nnodes, elem1d->node);
        printScreen_debug2_dbl("elem_thickness", elem_thickness, nnodes, elem1d->node);
        printScreen_debug2_dbl("elem_old_thickness", elem_old_thickness, nnodes, elem1d->node);
        printScreen_debug2_dbl("elem_older_thickness", elem_older_thickness, nnodes, elem1d->node);
        printScreen_debug2_dbl("elem_source", elem_source, nnodes, elem1d->node);
        printScreen_debug2_dbl("elem_sink", elem_sink, nnodes, elem1d->node);
        printScreen_debug2_dbl("elem_bed_shear", elem_bed_shear, nnodes, elem1d->node);
        printScreen_debug2_dbl("elem_bedvel_nrml", elem_bedvel_nrml, nnodes, elem1d->node);
        printScreen_debug_svect2d("elem_bedload_vel", elem_bedload_vel, nnodes, elem1d->node);
    }
#endif
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    sarray_dbl_init(elem_rhs, grid->elem2d[elem1d->elem2d].nnodes);
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                 NATURAL BOUNDARY CONDITIONS
     *                                   VELOCITY CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the velocity/discharge addition to the elemental residual. \n
     * \note  If the hydro flux is into the model (-), then the constituent flux must be assigned - it is set to 0 if not.
     * \note  If the hydro flux is out of the model, then the assigned constituent flux is ignored and calculated implicitly.
     *********************************************************************************************/
    if(str_values[string].ol_flow.bc_flag == BCT_VEL_NEU) {

        int isers = str_values[string].ol_flow.isigma;
        double water_flux = -1.0 * sseries_get_value(isers, mod->series_head, 0);

        if(water_flux < -SMALL) {
            double con_user = 0.;
            if(str_values[string].sed[ised].bc_flag == BCT_NEU) {
                isers = str_values[string].sed[ised].isigma;
                con_user = 0.; //sseries_get_value(isers, mod->series_head, 0) / con.property[0];
                //con_user *= 1.E-6;
            }
            integrate_line_phi(djac, con_user * water_flux * dt, elem_rhs);
        } else {
            integrate_line_phi_f(djac, water_flux * dt, elem_c, elem_rhs)
        }
#ifdef DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) {
            rhs_1dof("1D BED LOAD || BCT_VEL_NEU: ", nnodes, ie, elem1d->nodes, elem_rhs);
        }
#endif

    }

    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    OUTFLOW CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the outflow addition to the elemental residual. \n
     * \note CJT \:: calculates implicit flow
     *********************************************************************************************/
    else if(str_values[string].ol_flow.bc_flag == BCT_OUTFLOW) {
        integrate_line_phi_f_g_h(djac, dt, elem_c[0], elem_c[1], elem_thickness[0], elem_thickness[1], elem_bedvel_nrml[0], elem_bedvel_nrml[1], elem_rhs);
#ifdef DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) {
            rhs_1dof("1D BED LOAD || OUTFLOW: ", nnodes, ie, elem1d->nodes, elem_rhs);
        }
#endif
    }

    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                             TAILWATER/PRESSURE CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the tailwater or pressure addition to the elemental residual. \n
     * \note  If the implicit hydro flux is into the model (-), then the constituent flux must be assigned - it is set to 0 if not.
     * \note  If the implicit hydro flux is out of the model, then the assigned constituent flux is ignored and calculated implicitly.
     *********************************************************************************************/
    else if (str_values[string].ol_flow.bc_flag == BCT_PRS_NEU || str_values[string].ol_flow.bc_flag == BCT_DIS_NEU) {
        double flux = one_2 * (elem_bedvel_nrml[0] * elem_head[0] + elem_bedvel_nrml[1] * elem_head[1]);
        if(flux >= 0.) {
             integrate_line_phi_f_g_h(djac, dt, elem_c[0], elem_c[1], elem_thickness[0], elem_thickness[1], elem_bedvel_nrml[0], elem_bedvel_nrml[1], elem_rhs);
        } else {
            /*! flow is into the model, so we need to have a concentration defined */
            double con_user =0.;
            if(str_values[string].sed[ised].bc_flag == BCT_NEU) {
                int isers = str_values[string].sed[ised].isigma;
                con_user = 0.; //sseries_get_value(isers, mod->series_head, 0) / con.property[0];
                //con_user *= 1.E-6;
            }
            integrate_line_phi_f_g(djac, con_user * dt, elem_thickness, elem_bedvel_nrml, elem_rhs);
        }
#ifdef DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) {
            rhs_1dof("1D BED LOAD || BCT_PRS_NEU: ", nnodes, ie, elem1d->nodes, elem_rhs);
        }
#endif
    }

    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                     SET DIRICHLET BCS
     *-------------------------------------------------------------------------------------------*/
    for (i=0; i<nnodes; i++) {
        if (string > NORMAL) {
            if(str_values[string].sed[ised].bc_flag == BCT_DIR || str_values[string].sed[ised].bc_flag == BCT_CEQ) {
                elem_rhs[i] = 0;
            }
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif
}
