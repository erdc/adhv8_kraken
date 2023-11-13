/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes the boundary contributions for the 3D depth-averaged continuity and velocity system solve.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] elem_rhs            (SDOF *)   the elemental residual array
 * @param[in]  mod                 (SMODEL *) a pointer to the model struct
 * @param[in]  ie                  (int)      the elemental id
 * @param[in]  pertubation         (double)   the Newton pertubation
 * @param[in]  perturb_node        (int)      the local node ID to be perturbed
 * @param[in]  perturb_var         (int)      an integer ID for the variable to be perturbed
 * @param[in]  perturb_sign        (int)      +/- sign for the direction of Newton perturbation
 * @param[in]  node_in_column_flag (int *)    flags nodes that share a column with the perturbed node
 * @param[in]  DEBUG               (int)      a debug flag
 *
 * \details Integrates the weak, discrete outflow boundary terms: \n
 * \f{eqnarray*}{
 *    \residDA{I}{}{c} = & \sum\limits_{i \in C(I)}
 *      \bigg[  \sum\limits_e \bigg[\intesw{\phiddd{i} ( \velh \cdot \nrml)}  \bigg] \bigg] + \,
 *              \sum\limits_e \bigg[\intes{\phiddd{i} \bigg( \deriv{\elev{}}{t} - S_{\elev{}} \bigg) n_z\,} + \,
 *                           \inteb{\phiddd{i} \bigg( \deriv{\bed{}}{t} - S_{\bed{}} \bigg) n_z\,} \bigg] \\
 *     \resid{i}{}{mx}  = & \sum\limits_{e}
 *                          \bigg[ -\bcStress{\eta}{e}{\phidd{i}}{\lrpb{\tau_{winds,x}^{\,h} + \tau_{waves,x}^{\,h}}} \,+
 *                                  \bcConv{}{e}{\phidd{i}}{(\velr^{\,h} \, u^{\,h})} \, +
 *                                  \bcPressure{}{e}{\phidd{i}}{P^{\,h}}{x} \bigg] \\
 *
 *     \resid{i}{}{my}  = & \sum\limits_{e}
 *                          \bigg[ -\bcStress{\eta}{e}{\phidd{i}}{\lrpb{\tau_{winds,y}^{\,h} + \tau_{waves,y}^{\,h}}} \,+
 *                                  \bcConv{}{e}{\phidd{i}}{(\velr^{\,h} \, v^{\,h})} \, +
 *                                  \bcPressure{}{e}{\phidd{i}}{P^{\,h}}{y} \bigg]
 * \f}
 *
 * \note if mask = 2, then the original equations are replaced with
 *         x_eq = 3D continuity
 *         y_eq = tangential momemtum equation
 * \note displacement *cannot* be hotstarted here!
 * \note CJT \:: Do NOT return until the end of this routine, as DB conditions are applied at the end.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

static int printFieldWidth = 20;
static int printPrecision  = 10;

static int DEBUG_LOCAL = OFF;
static int DEBUG_PICKETS = OFF;
static int DEBUG_NODE_ID = OFF;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_hvel_boundary_resid(SMODEL *mod, DOF_3 *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int *node_in_column_flag, int DEBUG) {
    
    int i;
    
#ifdef _DEBUG
    if (DEBUG_NODE_ID != UNSET_INT) {
        for (i=0; i<mod->grid->elem2d[ie].nnodes; i++) {
            if (mod->grid->elem2d[ie].nodes[i] == DEBUG_NODE_ID - 1)  {DEBUG = ON; break;} else {DEBUG = OFF; DEBUG_LOCAL = OFF; DEBUG_PICKETS = OFF;}
        }
    }
    if (DEBUG == ON) DEBUG_LOCAL = ON;
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    time_t time1;  time(&time1);
#endif
    
    int PRS_FLAG = ON, nd1, nd2;
    
    // aliases
    SSW_3D *sw = mod->sw->d3;
    double dt = mod->dt;
    double g = mod->gravity;
    STR_VALUE *str_values = mod->str_values;
    SFLAGS flags = mod->flag;
    double density = mod->density;
    
    assert(mod->nsys == 3);
    assert(mod->nsys_sq == 9);
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GRID VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    SGRID *grid = mod->grid;
    SELEM_2D *elem2d = &(grid->elem2d[ie]);
    int string = elem2d->string;
    int imat = elem2d->mat;
    SMAT_SW mat = *mod->mat[imat].sw;
    int hydro_bc_flag = str_values[string].flow.bc_flag;
    if (str_values[string].ol_flow.bc_flag == BCT_PRS_NEU) hydro_bc_flag = BCT_PRS_NEU;
    
    int isElementTriangle = FALSE;
    int isSidewall = FALSE;
    int isSurface = FALSE;
    int isBed = FALSE;
    if      (elem2d->bflag == 2) {isSidewall = TRUE;}
    else if (elem2d->bflag == 0) {isSurface = TRUE;}
    else if (elem2d->bflag == 1) {isBed = TRUE;}
    else {
        printf("bflag = %d\n",elem2d->bflag);
        tl_error(">> 2D Boundary face bflag is not defined.");
    }
    
    int **nds_on_edge;
    int nnodes = UNSET_INT, nnodes_quad = UNSET_INT;
    if (elem2d->nnodes == NDONTRI) {
        isElementTriangle = TRUE;
        nnodes = NDONTRI;
        nnodes_quad = 2*NDONTRI; // for 2D triangle
        nds_on_edge = grid->nd_on_TriEdge;
    } else {
        isElementTriangle = FALSE;
        nnodes = NDONQUAD;
        nnodes_quad = 2*NDONQUAD; // for 2D prism
        nds_on_edge = grid->nd_on_QuadEdge;
    }
    
    // get local 2d nodes
    SNODE nodes[nnodes];
    for (i=0; i<nnodes; i++) {
        snode_copy(&(nodes[i]), grid->node[elem2d->nodes[i]]);
    }
    // get bottom and surface node
    int iseg;
    int elem_surface_nodeID[nnodes];
    int elem_bottom_nodeID[nnodes];
    for(i=0; i<nnodes; i++) {
        iseg = find_vertical_segment(grid, elem2d->nodes[i], grid->vertical_hash);
        elem_surface_nodeID[i] = grid->vertical_list[iseg]->id;
        elem_bottom_nodeID[i] = (grid->vertical_list[iseg]->prev)->id;
    }
    
    // boundary condition masks
    int mask[mod->nsys * nnodes], istart;
    for (i=0; i<nnodes; i++) {
        istart = mod->nsys * i;
        mask[istart]   = mod->bc_mask[elem2d->nodes[i] * mod->nsys];      // x_eq
        mask[istart+1] = mod->bc_mask[elem2d->nodes[i] * mod->nsys + 1];  // y_eq
        mask[istart+2] = mod->bc_mask[elem2d->nodes[i] * mod->nsys + 2];  // c_eq
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    int perturbed_surf_node = UNSET_INT; // the surface node of the column containing the perturbed node
    double elem_displacement[nnodes];
    global_to_local_dbl(sw->displacement, elem_displacement, elem2d->nodes, nnodes);
    if (perturb_var == PERTURB_DPL) {
        double elem_dpl_perturb[nnodes];
        global_to_local_dbl(sw->dpl_perturbation, elem_dpl_perturb, elem2d->nodes, nnodes);
        
        for (i=0; i<nnodes; i++) {
            if (node_in_column_flag[i] == 1) {
                elem_displacement[i] += perturb_sign * elem_dpl_perturb[i];
            }
        }
        
        // get the perturbed surface node
        for(i=0; i<nnodes; i++) {
            if (node_in_column_flag[i] == 1) { // use this column surface perturbation
                perturbed_surf_node = elem_surface_nodeID[i];
                break;
            }
        }
    }
    
    SVECT elem_vel[nnodes];
    global_to_local_svect(sw->vel, elem_vel, elem2d->nodes, nnodes);
    if (perturb_var == PERTURB_U) {
        PRS_FLAG = ON; //OFF; // CJT ::: TEMPORARY!
        elem_vel[perturb_node].x += perturb_sign * pertubation;
    } else if (perturb_var == PERTURB_V) {
        PRS_FLAG = ON; //OFF;
        elem_vel[perturb_node].y += perturb_sign * pertubation;
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    
    // node vertex and midpoint locations
    SVECT elem_nds[nnodes_quad], elem_nds_fixed[nnodes_quad];
    elem_get_midpt_locations2(nodes, elem_displacement, elem_nds, nnodes, nnodes_quad, nds_on_edge);
    
    // get vector node storage (for convienence)
    if (isElementTriangle == TRUE) {
        get_triangle_linear_djac_nrml_gradPhi(elem2d, NULL, elem_nds);
    } else {
#ifdef _DEBUG // this must be true for analytic quadrilateral integrations to be correct
        assert(fabs(elem_nds[2].x-elem_nds[1].x) < 1e-6 && fabs(elem_nds[2].y-elem_nds[1].y) < 1e-6);
        assert(fabs(elem_nds[3].x-elem_nds[0].x) < 1e-6 && fabs(elem_nds[3].y-elem_nds[0].y) < 1e-6);
#endif
        double elem_dpl_fixed[4] = {0., 0., 0., 0.}; // needs to be initial dpl, 0 should be ok as long as dpl is not hotstarted
        elem_get_midpt_locations2(nodes, elem_dpl_fixed, elem_nds_fixed, nnodes, nnodes_quad, nds_on_edge);
        elem2d->nrml = get_elem2d_normals(elem_nds);
    }
    double nx = elem2d->nrml.x, ny = elem2d->nrml.y, nz = elem2d->nrml.z;
    
    // calculate normal velocities
    double elem_normal_vel[nnodes];
    for (i=0; i<nnodes; i++) {
        elem_normal_vel[i] = svect_dotp(elem_vel[i], elem2d->nrml);
    }
    
    // for convienence
    double u[nnodes], v[nnodes], w[nnodes];
    for (i=0; i<nnodes; i++) {
        u[i] = elem_vel[i].x;
        v[i] = elem_vel[i].y;
        w[i] = elem_vel[i].z;
    }
    
    // 2D element tangent vector
    SVECT2D elem_tanvec[nnodes];
    global_to_local_svect2d(sw->tanvec, elem_tanvec, elem2d->nodes, elem2d->nnodes);
    
    // local pressures
    double elem_pressure[nnodes_quad], elem_prs_plus[nnodes_quad], elem_prs_minus[nnodes_quad];
    sarray_init_dbl(elem_pressure, nnodes_quad);
    
    if (PRS_FLAG == ON) {
        
        int valueID1, valueID2;
        double ret_val = 0., prs_value[5] = {0., 0., 0., 0., 0.};
        
        // -- at vertex
        global_to_local_dbl(sw->prs, elem_pressure, elem2d->nodes, nnodes);
        
        // -- at midpoints
        for(i=nnodes; i<nnodes_quad; i++) {
            nd1 = nds_on_edge[i - nnodes][0];
            nd2 = nds_on_edge[i - nnodes][1];
            ret_val = get_value_midpt_list(grid, grid->midpt_list, elem2d->nodes[nd1], elem2d->nodes[nd2], prs_value);
            elem_pressure[i] = prs_value[0];
        }
        
        if (perturb_var == PERTURB_DPL) {
            
            // perturb vertex nodes in perturbed column
            if (perturb_sign == 1) { // + perturb
                global_to_local_dbl(sw->prs_plus,  elem_prs_plus,  elem2d->nodes, nnodes);
                for(i=0; i<nnodes; i++) {
                    if(node_in_column_flag[i] == 1) {
                        elem_pressure[i] = elem_prs_plus[i];
                    }
                }
                valueID1 = 1; valueID2 = 2;
            } else if (perturb_sign == -1) { // - perturb
                global_to_local_dbl(sw->prs_minus, elem_prs_minus, elem2d->nodes, nnodes);
                for(i=0; i<nnodes; i++) {
                    if(node_in_column_flag[i] == 1) {
                        elem_pressure[i] = elem_prs_minus[i];
                    }
                }
                valueID1 = 3; valueID2 = 4;
            } else {tl_error("perturbation sign is bad.");}
            
            // peturb midpoint nodes affected by the perturbed column
            for(i=nnodes; i<nnodes_quad; i++) {
                nd1 = nds_on_edge[i - nnodes][0];
                nd2 = nds_on_edge[i - nnodes][1];
                ret_val = get_value_midpt_list(grid, grid->midpt_list, elem2d->nodes[nd1], elem2d->nodes[nd2], prs_value);
                
                if(elem_surface_nodeID[nd1] < elem_surface_nodeID[nd2]) {
                    if(node_in_column_flag[nd1] == 1) {
                        elem_pressure[i] = prs_value[valueID1];
                    } else if (node_in_column_flag[nd2] == 1) {
                        elem_pressure[i] = prs_value[valueID2];
                    } else {
                        elem_pressure[i] = prs_value[0];
                    }
                } else {
                    if(node_in_column_flag[nd1] == 1) {
                        elem_pressure[i] = prs_value[valueID2];
                    } else if(node_in_column_flag[nd2] == 1) {
                        elem_pressure[i] = prs_value[valueID1];
                    } else {
                        elem_pressure[i] = prs_value[0];
                    }
                }
            }
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        printf("SW-3D HVEL BOUNDARY ELEM RESID :: ie: %d \t dt: %20.10f",ie,dt);
        if      (isSidewall == TRUE) {printf(" SIDEWALL FACE ");}
        else if (isSurface  == TRUE) {printf(" SURFACE FACE ");}
        else if (isBed      == TRUE) {printf(" BED FACE ");}
        if      (perturb_var == PERTURB_U)   {printf("\t perturbing U || node: %d\n",nodes[perturb_node].id);}
        else if (perturb_var == PERTURB_V)   {printf("\t perturbing V || node: %d\n",nodes[perturb_node].id);}
        else if (perturb_var == PERTURB_DPL) {printf("\t perturbing DPL || node: %d\n",nodes[perturb_node].id);}
        if (isElementTriangle) {
            printf("nodes: %d %d %d \n",elem2d->nodes[0],elem2d->nodes[1],elem2d->nodes[2]);
            printf("djac2d: %30.20f \t djac2d3d: %30.20f \t djac2d3d_fixed: %30.20f \n",elem2d->djac,elem2d->djac3d,elem2d->djac3d_fixed);
        }
        printf("\n--------------------------------------------------------- \n");
        printf("Node Locations With Pressure Values --------------- \n");
        printf("nx = %20.10f; ny = %20.10f; nz = %20.10f;\n",elem2d->nrml.x,elem2d->nrml.y,elem2d->nrml.z);
        for (i=0; i<nnodes; i++) {
            if (elem_pressure[i] < -1e-6) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> Pressure is negative on element vertex node.");
            }
            printf("node_id = %d; x%d = %*.*e; y%d = %*.*e; z%d = %*.*e; p%d = %*.*e;\n",elem2d->nodes[i]+1,i+1,printFieldWidth,printPrecision,elem_nds[i].x,i+1,printFieldWidth,printPrecision,elem_nds[i].y,i+1,printFieldWidth,printPrecision,elem_nds[i].z,i+1,printFieldWidth,printPrecision,elem_pressure[i]);
        }
        for (i=nnodes; i<nnodes_quad; i++) {
            nd1 = nds_on_edge[i - nnodes][0];
            nd2 = nds_on_edge[i - nnodes][1];
            if (elem_pressure[i] < -1e-6) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> Pressure is negative on element midpoint node.");
            }
            printf("x%d = %*.*e; y%d = %*.*e; z%d = %*.*e; p%d = %*.*e; # mdpt of element nodes [%d,%d] \n",i+1,printFieldWidth,printPrecision,elem_nds[i].x,i+1,printFieldWidth,printPrecision,elem_nds[i].y,i+1,printFieldWidth,printPrecision,elem_nds[i].z,i+1,printFieldWidth,printPrecision,elem_pressure[i],nd1+1,nd2+1);
        }
        printf("--------------------------------------------------------- \n");
        printf("--------------------------------------------------------- \n");
        for (i=0; i<nnodes * mod->nsys; i++) printf("mask[%d]: %d\n",i,mask[i]);
        for (i=0; i<nnodes; i++) printf("tanvec[%d]: %*.*e \t %*.*e\n",i,printFieldWidth,printPrecision,elem_tanvec[i].x,printFieldWidth,printPrecision,elem_tanvec[i].y);
        printf("normals: %*.*e %*.*e %*.*e\n",printFieldWidth,printPrecision,nx,printFieldWidth,printPrecision,ny,printFieldWidth,printPrecision,nz);
        printScreen_debug2_dbl("elem_displacement", elem_displacement, nnodes, elem2d->nodes);
        printScreen_debug2_dbl("elem_normal_vel", elem_normal_vel, nnodes, elem2d->nodes);
        printScreen_debug_svect("elem_vel", elem_vel, nnodes, elem2d->nodes);
        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    }
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    dof3_init_array(elem_rhs, elem2d->nnodes);
    
    double prs_integral[nnodes]; sarray_init_dbl(prs_integral, nnodes);
    double dprs_integral[nnodes_quad]; sarray_init_dbl(dprs_integral, nnodes_quad);
    double elem_dprs[nnodes_quad]; sarray_init_dbl(elem_dprs, nnodes_quad);
    double implicit_flow[nnodes]; sarray_init_dbl(implicit_flow, nnodes);
    double explicit_flow[nnodes]; sarray_init_dbl(explicit_flow, nnodes);
    double flux_x[nnodes]; sarray_init_dbl(flux_x, nnodes);
    double flux_y[nnodes]; sarray_init_dbl(flux_y, nnodes);
    double friction_integral_x[nnodes]; sarray_init_dbl(friction_integral_x, nnodes);
    double friction_integral_y[nnodes]; sarray_init_dbl(friction_integral_y, nnodes);
    
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
        
        // -----------------------------------------------------------------------------------------
        // depth-averaged continuity addition ------------------------------------------------------
        sarray_init_dbl(explicit_flow, nnodes);
        sarray_init_dbl(implicit_flow, nnodes);
        if (isElementTriangle) {
            integrate_triangle_phi_f(elem2d->djac3d, dt, elem_normal_vel, implicit_flow);
            integrate_triangle_phi(elem2d->djac3d_fixed, dt * hydro_flux_normal, explicit_flow);
        } else {
            integrate_quadZ_phi_f(elem_nds, dt, elem_normal_vel, implicit_flow);
            integrate_quadZ_phi(elem_nds_fixed, dt * hydro_flux_normal, explicit_flow);
        }
#ifdef _DEBUG
        if (DEBUG) {
            for (i=0; i<nnodes; i++) {
                printf("implicit flow[%d]: %20.10f\n",i,implicit_flow[i]);
            }
            for (i=0; i<nnodes; i++) {
                printf("explicit flow[%d]: %20.10f\n",i,explicit_flow[i]);
            }
        }
#endif
        for (i=0; i<nnodes; i++) {
            if (mask[i * mod->nsys] == YES) {
                elem_rhs[i].c_eq = implicit_flow[i];
            } else {
                elem_rhs[i].c_eq = explicit_flow[i];
            }
        }
        
        // -----------------------------------------------------------------------------------------
        // depth-averaged momentum addition --------------------------------------------------------
        sarray_init_dbl(prs_integral, nnodes);
        if (isElementTriangle) {
            integrate_triangle_phi_fquad(elem2d->djac3d, dt * g, elem_pressure, prs_integral); // quadratic
        } else {
            integrate_quadZ_phi_f(elem_nds, dt * g, elem_pressure, prs_integral); // linear
        }
        double drag = 0.;
        for (i=0; i<nnodes; i++) {
            if (mask[i * mod->nsys] != 2) {
                elem_rhs[i].x_eq = prs_integral[i]*nx + friction_integral_x[i]*drag;
                elem_rhs[i].y_eq = prs_integral[i]*ny + friction_integral_y[i]*drag;
            } else {
                elem_rhs[i].x_eq = implicit_flow[i] - explicit_flow[i];
                elem_rhs[i].y_eq = prs_integral[i] * (nx * elem_tanvec[i].x + ny * elem_tanvec[i].y);
                elem_rhs[i].y_eq += drag * (friction_integral_x[i] * elem_tanvec[i].x + friction_integral_y[i] * elem_tanvec[i].y);
            }
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_3dof("3D HVEL BCT_VEL_NEU/BCT_DIS_NEU: ",nnodes, ie, elem2d->nodes, elem_rhs);
            Is_DoubleArray_Inf_or_NaN(implicit_flow,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(explicit_flow,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(prs_integral,nnodes,__FILE__ ,__LINE__);
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
        
        // -----------------------------------------------------------------------------------------
        // depth-averaged continuity addition ------------------------------------------------------
        sarray_init_dbl(implicit_flow, nnodes);
        if (isElementTriangle) {
            integrate_triangle_phi_f(elem2d->djac3d, dt, elem_normal_vel, implicit_flow);
        } else {
            integrate_quadZ_phi_f(elem_nds, dt, elem_normal_vel, implicit_flow);
        }
        for (i=0; i<nnodes; i++) {elem_rhs[i].c_eq += implicit_flow[i];}
        
        // -----------------------------------------------------------------------------------------
        // depth-averaged momentum addition --------------------------------------------------------
        sarray_init_dbl(flux_x, nnodes);
        sarray_init_dbl(flux_y, nnodes);
        sarray_init_dbl(prs_integral, nnodes);
        if (isElementTriangle) {
            integrate_triangle_phi_f_g(elem2d->djac3d, dt, elem_normal_vel, u, flux_x);
            integrate_triangle_phi_f_g(elem2d->djac3d, dt, elem_normal_vel, v, flux_y);
            integrate_triangle_phi_fquad(elem2d->djac3d, dt, elem_pressure, prs_integral); // quadratic
        } else {
            integrate_quadZ_phi_f_g(elem_nds, dt, elem_normal_vel, u, flux_x);
            integrate_quadZ_phi_f_g(elem_nds, dt, elem_normal_vel, v, flux_y);
            integrate_quadZ_phi_f(elem_nds, dt, elem_pressure, prs_integral); // linear
        }
        for (i=0; i<nnodes; i++) {
            if (mask[i * mod->nsys] != 2) {
                elem_rhs[i].x_eq += flux_x[i] + prs_integral[i] * g * nx;
                elem_rhs[i].y_eq += flux_y[i] + prs_integral[i] * g * ny;
            } else {
                elem_rhs[i].y_eq += elem_tanvec[i].x * flux_x[i] + elem_tanvec[i].y * flux_y[i];
                elem_rhs[i].y_eq += prs_integral[i] * g * (nx * elem_tanvec[i].x + ny * elem_tanvec[i].y);
            }
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_3dof("3D HVEL BCT_OUTFLOW: ",nnodes, ie, elem2d->nodes, elem_rhs);
            Is_DoubleArray_Inf_or_NaN(implicit_flow,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(flux_x,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(flux_y,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(prs_integral,nnodes,__FILE__ ,__LINE__);
        }
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                        BED FLUX
     *--------------------------------------------------------------------------------------------
     * Calculates the bed flux contribution to the residual. \n
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
        //            elem_elem_rhs[i].c_eq = eta_integral[i] - eta_old_integral[i];
        //        }
#endif
        
        
        // -----------------------------------------------------------------------------------------
        // depth-averaged momentum addition --------------------------------------------------------
        
        double avg_depth = one_3*(elem_pressure[0] + elem_pressure[1] + elem_pressure[2]);  // since elem_pressure - P/(g*rho) and P = rho * g * h
        double drag = fe_sw3_get_roughness(mod, avg_depth, string);
        double normal_discharge = 0.; // for now
        double avg_u = sarray_avg_dbl(u, nnodes);
        double avg_v = sarray_avg_dbl(v, nnodes);
        double avg_w = sarray_avg_dbl(w, nnodes);
        double avg_vel_mag = sqrt(avg_u*avg_u + avg_v*avg_v + avg_w*avg_w);
        
        sarray_init_dbl(friction_integral_x, nnodes);
        sarray_init_dbl(friction_integral_y, nnodes);
        sarray_init_dbl(prs_integral, nnodes);
        integrate_triangle_phi_f(elem2d->djac3d, dt * avg_vel_mag, u, friction_integral_x);
        integrate_triangle_phi_f(elem2d->djac3d, dt * avg_vel_mag, v, friction_integral_y);
        integrate_triangle_phi_fquad(elem2d->djac3d, dt * g, elem_pressure, prs_integral); // quadratic // note: if prism prs is linear, should be linear here too
        
        for (i=0; i<nnodes; i++) {
            if (mask[i * mod->nsys] != 2) {
                elem_rhs[i].x_eq = prs_integral[i]*nx + friction_integral_x[i]*drag;
                elem_rhs[i].y_eq = prs_integral[i]*ny + friction_integral_y[i]*drag;
            }
            else {
                elem_rhs[i].y_eq =  prs_integral[i] * (nx * elem_tanvec[i].x + ny * elem_tanvec[i].y);
                elem_rhs[i].y_eq += drag * (friction_integral_x[i] * elem_tanvec[i].x + friction_integral_y[i] * elem_tanvec[i].y);
            }
            
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_3dof("3D HVEL BCT_BED: ",nnodes, ie, elem2d->nodes, elem_rhs);
            Is_DoubleArray_Inf_or_NaN(implicit_flow,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(friction_integral_x,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(friction_integral_y,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(prs_integral,nnodes,__FILE__ ,__LINE__);
        }
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                          FREE SURFACE
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
        SVECT2D elem_surface_stress;
        
        double elem_old_displacement[nnodes], elem_older_displacement[nnodes];
        global_to_local_dbl(sw->old_displacement, elem_old_displacement, elem2d->nodes, elem2d->nnodes);
        global_to_local_dbl(sw->older_displacement, elem_older_displacement, elem2d->nodes, elem2d->nnodes);
        
#ifdef _SEDIMENT
        double elem_bed_displacement[nnodes], elem_old_bed_displacement[nnodes],elem_older_bed_displacement[nnodes];
        global_to_local_dbl(sw->bed_displacement, elem_bed_displacement, elem2d->nodes, elem2d->nnodes);
        global_to_local_dbl(sw->old_bed_displacement, elem_old_bed_displacement, elem2d->nodes, elem2d->nnodes);
        global_to_local_dbl(sw->older_bed_displacement, elem_older_bed_displacement, elem2d->nodes, elem2d->nnodes);
        // cjt :: add bed displacement to total displacement ... be carefull here, no need for db/dt now I think
        sarray_add_replace_dbl(elem_displacement, elem_bed_displacement, elem2d->nnodes);
        sarray_add_replace_dbl(elem_old_displacement, elem_old_bed_displacement, elem2d->nnodes);
        sarray_add_replace_dbl(elem_older_displacement, elem_older_bed_displacement, elem2d->nnodes);
#endif
        
        // -----------------------------------------------------------------------------------------
        // depth-averaged continuity addition ------------------------------------------------------
        
        // read normal velocity/discharge
        int isers = str_values[string].flow.isigma;
        double hydro_flux_normal = -sseries_get_value(isers, mod->series_head, 0); // (-) for flow into domain
        if (hydro_bc_flag == BCT_DIS_NEU) {
            hydro_flux_normal /= str_values[string].total_area;
        }
        sarray_init_dbl(explicit_flow, nnodes);
        integrate_triangle_phi(elem2d->djac3d_fixed, dt * hydro_flux_normal, explicit_flow);
        
        double elem_old_eta[NDONTRI];  // get wse at t(i+1/2)
        ELEM2D_GET_TPOSITION(elem_old_eta, elem_old_displacement, elem_older_displacement, mod->tau_temporal, dt, dt);
        double elem_eta[NDONTRI];      // get wse at t(i+3/2)
        ELEM2D_GET_TPOSITION(elem_eta, elem_displacement, elem_old_displacement, mod->tau_temporal, dt, dt);
        
        double eta_integral[nnodes]; sarray_init_dbl(eta_integral, nnodes);
        double eta_old_integral[nnodes]; sarray_init_dbl(eta_old_integral, nnodes);
        integrate_triangle_phi_f(elem2d->djac, 1., elem_eta, eta_integral);
        integrate_triangle_phi_f(elem2d->djac, 1., elem_old_eta, eta_old_integral);
        
        for (i=0; i<nnodes; i++) {
            elem_rhs[i].c_eq = (eta_integral[i] - eta_old_integral[i]) - explicit_flow[i];
        }
        
        // -----------------------------------------------------------------------------------------
        // depth-averaged momentum addition --------------------------------------------------------
        if (flags.WIND == ON || flags.WAVE == ON) {
            
            SVECT2D elem_wave_stress, elem_wind_stress;
            svect2d_init(&elem_wave_stress);
            svect2d_init(&elem_wind_stress);
            
            int grid2d_nodeIDs[MAX_NNODES_ON_ELEM2D]; // map from 3d node ids to 2d grid node ids
            for(i=0; i<elem2d->nnodes; i++) {
                grid2d_nodeIDs[i] = grid->nodeID_3d_to_2d_sur[ elem2d->nodes[i] ];
            }
            
            if (flags.WIND) {
                int nws = OFF;
#ifdef _WINDLIB
                nws = mod->windlib->nws;
#endif
                double avg_depth = one_3*(elem_pressure[0] + elem_pressure[1] + elem_pressure[2]); // since elem_pressure - P/(g*rho) and P = rho * g * h
                elem_wind_stress = swind_elem2d_get_local_stress(3, sw->winds, elem2d, grid2d_nodeIDs, avg_depth, g, density, mat.windatt, mat.wind_flag, flags.WIND_LIBRARY, nws);
            }
            if (flags.WAVE) {
                elem_wave_stress = swave_elem2d_get_local_stress(3, flags.CSTORM_WSID, sw->waves, *elem2d, grid2d_nodeIDs);
            }
            
            
            
            elem_surface_stress.x = elem_wind_stress.x + elem_wave_stress.x;
            elem_surface_stress.y = elem_wind_stress.y + elem_wave_stress.y;
            
            sarray_init_dbl(flux_x, nnodes);
            sarray_init_dbl(flux_y, nnodes);
            integrate_triangle_phi(elem2d->djac3d, dt * elem_surface_stress.x, flux_x);
            integrate_triangle_phi(elem2d->djac3d, dt * elem_surface_stress.y, flux_y);
            
            for (i = 0; i < nnodes; i++) {
                if (mask[i * mod->nsys] != 2) {
                    elem_rhs[i].x_eq = -flux_x[i];
                    elem_rhs[i].y_eq = -flux_y[i];
                } else {
                    elem_rhs[i].y_eq = -(flux_x[i] * elem_tanvec[i].x + flux_y[i] * elem_tanvec[i].y);
                }
            }
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            printf("elem_surface_stress: %30.20e \t %30.20e\n",elem_surface_stress.x,elem_surface_stress.y);
            rhs_3dof("3D HVEL BCT_FRS/SOURCE: ",nnodes, ie, elem2d->nodes, elem_rhs);
            Is_DoubleArray_Inf_or_NaN(eta_integral,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(eta_old_integral,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(flux_x,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(flux_y,nnodes,__FILE__ ,__LINE__);
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
        
        // -----------------------------------------------------------------------------------------
        // depth-averaged continuity addition ------------------------------------------------------
        sarray_init_dbl(implicit_flow, nnodes);
        if (isElementTriangle == TRUE) {
            integrate_triangle_phi_f(elem2d->djac3d, dt, elem_normal_vel, implicit_flow);
        } else {
            integrate_quadZ_phi_f(elem_nds, dt, elem_normal_vel, implicit_flow);
        }
        for (i=0; i<nnodes; i++) {elem_rhs[i].c_eq = implicit_flow[i];}
        
        // -----------------------------------------------------------------------------------------
        // depth-averaged momentum addition --------------------------------------------------------
        
        // Use boundary elevation to force boundary pressure
        if (PRS_FLAG == ON) {
            
            // get user water elevation
            int isers = str_values[string].ol_flow.isigma;
            double user_elevation = sseries_get_value(isers, mod->series_head, 0);
            
            // get implicit water elevation
            double elem_surface_dpl[nnodes];
            global_to_local_dbl(sw->displacement, elem_surface_dpl, elem_surface_nodeID, elem2d->nnodes);
            if (perturb_var == PERTURB_DPL) {
                for(i=0; i<elem2d->nnodes; i++) {
                    if (node_in_column_flag[i] == 1) {
                        elem_surface_dpl[i] += perturb_sign * sw->dpl_perturbation[elem_surface_nodeID[i]] ; // used in 2d PRS_NEU
                    }
                }
            }
            double current_elevation[nnodes];
            for (i=0; i<nnodes; i++) {
                current_elevation[i] = mod->grid->node[ elem_surface_nodeID[i] ].z + elem_surface_dpl[i];
            }
            
            // baroclinic density
            double elem_surface_density[nnodes];
            if (flags.BAROCLINIC == 1) {
                for(i=0; i<nnodes; i++) {
                    elem_surface_density[i] = sw->density[elem_surface_nodeID[i]];
                }
            }
            // normalize density
            if (mod->flag.BAROCLINIC > 0) {
                sarray_scale_replace_dbl( elem_surface_density, 1./mod->density, nnodes);
            } else {
                sarray_init_value_dbl(elem_surface_density, nnodes, 1.);
            }
            
            // now calculate a new pressure integral which equals the different of the user and current elevations
            double elem_dprs[nnodes_quad];
            for (i = 0; i<nnodes; i++) {
                elem_dprs[i] = elem_surface_density[i] * (user_elevation - current_elevation[i]);
            }
            
            for (i=nnodes; i<nnodes_quad; i++) {
                nd1 = nds_on_edge[i-nnodes][0];
                nd2 = nds_on_edge[i-nnodes][1];
                elem_dprs[i] = one_2 * (elem_dprs[nd1] + elem_dprs[nd2]);
            }
            
            // now calculated pressure(current elevation) + dpressure(boundary condition)
            sarray_add_replace_dbl(elem_pressure, elem_dprs, nnodes_quad);
            
            sarray_init_dbl(prs_integral, nnodes);
            if (isElementTriangle == TRUE) {
                integrate_triangle_phi_fquad(elem2d->djac3d, dt * g, elem_pressure, prs_integral); // quadratic
            } else {
                integrate_quadZ_phi_f(elem_nds, dt * g, elem_pressure, prs_integral); // linear
            }
            
            for (i=0; i<nnodes; i++) {
                if (mask[i * mod->nsys] != 2) {
                    elem_rhs[i].x_eq = prs_integral[i] * nx;
                    elem_rhs[i].y_eq = prs_integral[i] * ny;
                }
                else {
                    elem_rhs[i].y_eq = prs_integral[i] * (elem_tanvec[i].x * nx + elem_tanvec[i].y * ny);
                }
            }
        }
        
        // Now, we only let water momentum flux *out* of the domain, otherwise, no penetration (cjt :: should we not sum_norm_vel > 0 c_eq too then??)
        double sum_normal_vel = sarray_sum_dbl(elem_normal_vel,nnodes);
        if (sum_normal_vel > 0.) { // out of the domain
            
            // calculate normal relative velocities
            SVECT elem_rel_vel[nnodes];
            double normal_rvel[nnodes], elem_old_displacement[nnodes], elem_older_displacement[nnodes];
            global_to_local_dbl(sw->old_displacement, elem_old_displacement, elem2d->nodes, elem2d->nnodes);
            global_to_local_dbl(sw->older_displacement, elem_older_displacement, elem2d->nodes, elem2d->nnodes);
            elem_get_relative_velocities(elem_vel, elem_displacement, elem_old_displacement, elem_older_displacement, dt, nnodes, elem_rel_vel);
            for (i=0; i<nnodes; i++) {
                normal_rvel[i] = svect_dotp(elem_rel_vel[i], elem2d->nrml);
            }
            
            sarray_init_dbl(flux_x, nnodes);
            sarray_init_dbl(flux_y, nnodes);
            if (isElementTriangle == TRUE) {
                integrate_triangle_phi_f_g(elem2d->djac3d, dt, normal_rvel, u, flux_x);
                integrate_triangle_phi_f_g(elem2d->djac3d, dt, normal_rvel, v, flux_y);
            } else {
                integrate_quadZ_phi_f_g(elem_nds, dt, normal_rvel, u, flux_x);
                integrate_quadZ_phi_f_g(elem_nds, dt, normal_rvel, v, flux_y);
            }
            
            for (i=0; i<nnodes; i++) {
                if (mask[i * mod->nsys] != 2) {
                    elem_rhs[i].x_eq += flux_x[i];
                    elem_rhs[i].y_eq += flux_y[i];
                }
                else {
                    elem_rhs[i].y_eq += (flux_x[i] * elem_tanvec[i].x + flux_y[i] * elem_tanvec[i].y);
                }
            }
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_3dof("3D HVEL BCT_PRS_NEU: ",nnodes, ie, elem2d->nodes, elem_rhs);
            Is_DoubleArray_Inf_or_NaN(implicit_flow,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(prs_integral,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(flux_x,nnodes,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(flux_y,nnodes,__FILE__ ,__LINE__);
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
    //        if (mask[i * mod->nsys] == YES) {
    //            elem_rhs[i].x_eq = 0.;
    //        }
    //        if (mask[i * mod->nsys + 1] == YES) {
    //            elem_rhs[i].y_eq = 0.;
    //        }
    //        if (hydro_bc_flag == BCT_FRS) {
    //            if (mask[i * mod->nsys + 2] == YES) {
    //                elem_rhs[i].c_eq = 0.;
    //            }
    //        }
    //    }
    //#ifdef _DEBUG
    //    if (DEBUG_LOCAL == ON) {
    //        rhs_3dof("3D CONTINUITY DIRICHLET BCS: ",nnodes, ie, elem2d->nodes, elem_rhs);
    //        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    //    }
    //#endif
    
    
#ifdef _DEBUG
    time_t time2;  time(&time2);
    TIME_IN_HVEL_BOUNDARY_RESID += difftime(time2,time1);
#endif
    
}

//***************************************************************************//
//***************************************************************************//
//***************************************************************************//

