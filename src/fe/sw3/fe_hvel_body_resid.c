/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_hvel_body_resid.c This file collections functions responsible for
 *         the 3D shallow water h,u,v equation contributions to the elemental residual.     */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

DOF_3 elem_rhs_temporal[NDONPRISM]; // to keep track of total temporal term contributions
void fe_hvel_temporal(int ie, int nnodes, SELEM_3D *elem3d, SVECT *elem_nds, double djac, double *u, double *v, double dt_factor, DOF_3 *elem_rhs, char *string, int DEBUG_LOCAL);


static int printFieldWidth = 30;
static int printPrecision  = 20;

static int DEBUG_LOCAL = OFF;
static int DEBUG_PICKETS = OFF;
static int DEBUG_NODE_ID = OFF;
static int DEBUG_NODE_PE = OFF;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 3D shallow water h,u,v equation contributions to the elemental residual.
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
 * \note cjt \:: all residuals are multiplied by dt
 * \note cjt \:: I have found that assuming an average continuity equation for prisms prevents convergence.
 *                 So, it is assumed linear here.  An average momentum equation seems fine though, for now.
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_hvel_body_resid(SMODEL *mod, DOF_3 *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int *node_in_column_flag, int DEBUG) {
    
    int i;
    
#ifdef _DEBUG
    if (DEBUG_NODE_ID != UNSET_INT) {
        for (i=0; i<mod->grid->elem3d[ie].nnodes; i++) {
            if (mod->grid->elem3d[ie].nodes[i] == DEBUG_NODE_ID - 1 && mod->grid->smpi->myid == DEBUG_NODE_PE) {
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
    
    int PRESSURE_FLAG = ON;
    int j, nd1, nd2, perturbed_surf_node = UNSET_INT;
    double t1 = 0., constant = 0.;
    
    SSW_3D *sw = mod->sw->d3;
    double g = mod->gravity;
    double dt = mod->dt;
    double tau_pg = mod->tau_pg;
    STR_VALUE *str_values = mod->str_values;
    double alpha = mod->tau_temporal;
    double dtTerm1 = (1. + alpha / 2.);
    double dtTerm2 = (-1.) * (1. + alpha);
    double dtTerm3 = (alpha / 2.) ;
    
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        assert(alpha > -1E-6 && alpha < 1.0001);
        assert(perturb_var == PERTURB_NONE || perturb_var == PERTURB_U || perturb_var == PERTURB_V || perturb_var == PERTURB_DPL);
    }
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GRID VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    SGRID *grid = mod->grid;
    SQUAD *quad = grid->quad_prism;
    SELEM_3D *elem3d = &(grid->elem3d[ie]);
    int nnodes = elem3d->nnodes;
    int nnodes_quad = elem3d->nnodes_quad;
    int icol = elem3d->icol;
    int imat = elem3d->mat;
    SMAT_SW mat = *mod->mat[imat].sw;
    SELEM_2D *elem2d_sur = &(grid->elem2d[grid->elem2d_sur[icol]]);
    SELEM_2D *elem2d_bed = &(grid->elem2d[grid->elem2d_bed[icol]]);
#ifdef _DEBUG
    if(nnodes <= 0 || nnodes > MAX_NNODES_ON_ELEM3D) {
        selem3d_printScreen(elem3d);
        printf("nnodes: %d\n",nnodes);
        tl_check_all_pickets(__FILE__,__LINE__);
        tl_error("ERROR: nnodes");
    }
    if(nnodes_quad <= 0 || nnodes_quad > MAX_NNODES_ON_ELEM3D_QUAD) {
        selem3d_printScreen(elem3d);
        printf("nnodes_quad: %d\n",nnodes_quad);
        tl_check_all_pickets(__FILE__,__LINE__);
        tl_error("ERROR: nnodes_quad");
    }
    assert(elem2d_sur->nnodes == NDONTRI);
    assert(elem2d_bed->nnodes == NDONTRI);
    assert(icol<grid->nelems2d && icol > -1);
    assert(grid->elem2d_sur[icol]<grid->nelems2d && grid->elem2d_sur[icol]>-1);
    assert(grid->elem2d_bed[icol]<grid->nelems2d && grid->elem2d_bed[icol]>-1);
#endif
    
    int iseg;
    SNODE elem_nodes[nnodes];
    int elem_surface_nodeID[nnodes], elem_bottom_nodeID[nnodes];
    for(i=0; i<nnodes; i++) {
        iseg = find_vertical_segment(grid, elem3d->nodes[i], grid->vertical_hash);
        elem_surface_nodeID[i] = grid->vertical_list[iseg]->id;            // surface node ID
        elem_bottom_nodeID[i] = (grid->vertical_list[iseg]->prev)->id;     // bed node ID
#ifdef _DEBUG
        assert(elem_surface_nodeID[i] < grid->nnodes && elem_surface_nodeID[i] > -1);
        assert(elem_bottom_nodeID[i] < grid->nnodes && elem_bottom_nodeID[i] > -1);
#endif
        snode_copy(&(elem_nodes[i]), grid->node[elem3d->nodes[i]]);
    }
    
    int isElementTetrahedron = TRUE;
    if (nnodes != NDONTET) isElementTetrahedron = FALSE;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // displacement perturbation
    double elem_dpl[nnodes];  global_to_local_dbl(sw->displacement, elem_dpl, elem3d->nodes, nnodes);
    if (perturb_var == PERTURB_DPL) {
        double elem_dpl_perturb[MAX_NNODES_ON_ELEM3D];
        global_to_local_dbl(sw->dpl_perturbation, elem_dpl_perturb, elem3d->nodes, nnodes);
        
        for(i=0; i<nnodes; i++) {
            if (node_in_column_flag[i] == 1) {
                elem_dpl[i] += perturb_sign * elem_dpl_perturb[i];
            }
        }
        
        // get the perturbed surface node
        for(i=0; i<nnodes; i++) {
            if (node_in_column_flag[i] == 1) { // use this column surface perturbation
                perturbed_surf_node = elem_surface_nodeID[i]; // can just use perturb_node for this
                break;
            }
        }
    }
    
    // velocity perturbation
    SVECT elem_vel[nnodes];  global_to_local_svect(sw->vel, elem_vel, elem3d->nodes, nnodes);
    if (perturb_var == PERTURB_U) {
        PRESSURE_FLAG = OFF;
        elem_vel[perturb_node].x += perturb_sign * perturbation;
    } else if (perturb_var == PERTURB_V) {
        PRESSURE_FLAG = OFF;
        elem_vel[perturb_node].y += perturb_sign * perturbation;
    } else if (perturb_var == PERTURB_W) {
        PRESSURE_FLAG = OFF;
        elem_vel[perturb_node].z += perturb_sign * perturbation;
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // calculate pressures
    double elem_pressure[nnodes_quad];
    sarray_init_dbl(elem_pressure,nnodes_quad);
    
    if (PRESSURE_FLAG == ON) {
        
        int valueID1, valueID2;
        double ret_val = 0.;
        double prs_value[5] = {0., 0., 0., 0., 0.};
        
        // -- at vertex
        global_to_local_dbl(sw->prs, elem_pressure, elem3d->nodes, nnodes);
        
        // -- at midpoints
        for(i=nnodes; i<nnodes_quad; i++) {
            nd1 = elem3d->edges[i - nnodes][0];
            nd2 = elem3d->edges[i - nnodes][1];
            ret_val = get_value_midpt_list(grid, grid->midpt_list, elem3d->nodes[nd1], elem3d->nodes[nd2], prs_value);
            elem_pressure[i] = prs_value[0];
        }
        
        if (perturb_var == PERTURB_DPL) {
            // perturb vertex nodes in perturbed column
            double elem_prs_plus[nnodes], elem_prs_minus[nnodes];
            if (perturb_sign == 1) { // + perturb
                global_to_local_dbl(sw->prs_plus,  elem_prs_plus,  elem3d->nodes, nnodes);
                for(i=0; i<nnodes; i++) {
                    if(node_in_column_flag[i] == 1) {
                        elem_pressure[i] = elem_prs_plus[i];
                    }
                }
                valueID1 = 1; valueID2 = 2;
            } else if (perturb_sign == -1) { // - perturb
                global_to_local_dbl(sw->prs_minus, elem_prs_minus, elem3d->nodes, nnodes);
                for(i=0; i<nnodes; i++) {
                    if(node_in_column_flag[i] == 1) {
                        elem_pressure[i] = elem_prs_minus[i];
                    }
                }
                valueID1 = 3; valueID2 = 4;
            } else {tl_error("perturbation sign is bad.");}
            
            // peturb midpoint nodes affected by the perturbed column
            for(i=nnodes; i<nnodes_quad; i++) {
                nd1 = elem3d->edges[i - nnodes][0];
                nd2 = elem3d->edges[i - nnodes][1];
                ret_val = get_value_midpt_list(grid, grid->midpt_list, elem3d->nodes[nd1], elem3d->nodes[nd2], prs_value);
                
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
    
    // displacements
    double elem_dpl_old[nnodes], elem_dpl_older[nnodes];
    global_to_local_dbl(sw->old_displacement,   elem_dpl_old,   elem3d->nodes, nnodes);
    global_to_local_dbl(sw->older_displacement, elem_dpl_older, elem3d->nodes, nnodes);
#ifdef _SEDIMENT
    // cjt :: add bed displacement to total displacement ... be carefull here, no need for db/dt now I think
    double elem_bed_dpl[nnodes], elem_bed_dpl_old[nnodes], elem_bed_dpl_older[nnodes];
    global_to_local_dbl(sw->bed_displacement,       elem_bed_dpl,       elem3d->nodes, nnodes);
    global_to_local_dbl(sw->old_bed_displacement,   elem_bed_dpl_old,   elem3d->nodes, nnodes);
    global_to_local_dbl(sw->older_bed_displacement, elem_bed_dpl_older, elem3d->nodes, nnodes);
    
    sarray_add_replace_dbl(elem_dpl_old,   elem_bed_dpl_old,   nnodes);
    sarray_add_replace_dbl(elem_dpl_older, elem_bed_dpl_older, nnodes);
#endif
    
    // velocities
    SVECT elem_vel_old[nnodes], elem_vel_older[nnodes];
    global_to_local_svect(sw->old_vel, elem_vel_old, elem3d->nodes, nnodes);
    global_to_local_svect(sw->older_vel, elem_vel_older, elem3d->nodes, nnodes);
    
    // grid and relative velocities
    SVECT elem_rel_vel[nnodes], elem_grid_vel[nnodes];
    svect_init_array(elem_grid_vel, nnodes);
    for (i = 0; i < nnodes; i++) elem_grid_vel[i].z = (elem_dpl[i] - elem_dpl_old[i])/dt;
    svect_subtract_array2(elem_rel_vel, elem_vel, elem_grid_vel, nnodes);
    
    // nodal vectors with displacements
    SVECT elem_nds[nnodes_quad], elem_nds_old[nnodes_quad], elem_nds_older[nnodes_quad];
    elem_get_midpt_locations2(elem_nodes, elem_dpl,       elem_nds,       nnodes, nnodes_quad, elem3d->edges);
    elem_get_midpt_locations2(elem_nodes, elem_dpl_old,   elem_nds_old,   nnodes, nnodes_quad, elem3d->edges);
    elem_get_midpt_locations2(elem_nodes, elem_dpl_older, elem_nds_older, nnodes, nnodes_quad, elem3d->edges);
    
    // baroclinic calculations
    double elem_density[nnodes]; sarray_init_value_dbl(elem_density, nnodes, 1.);
    if (mod->flag.BAROCLINIC == 1) {
        global_to_local_dbl(sw->density, elem_density, elem3d->nodes, nnodes);
    }
    // normalize density (only used in turbulence)
    if (mod->flag.BAROCLINIC > 0) {
        ELEM3D_SCALE(elem_density, 1. / mod->density);
    }
    
    // calculate tetrehedral jacobians if needed
    SVECT grad_shp[nnodes], grad_shp_old[nnodes], grad_shp_older[nnodes];
    double elem_djac = 0., elem_djac_old = 0., elem_djac_older = 0.;
    double elem_volume = 0., elem_volume_old = 0., elem_volume_older = 0.;
    SVECT elem_grad_u, elem_grad_v, elem_grad_w, elem_grad_p, elem_grad_den;
    if (isElementTetrahedron == TRUE) {
        elem_djac = get_tet_linear_djac_gradPhi2(NULL, elem_nds, grad_shp);
        elem_djac_old = get_tet_linear_djac_gradPhi2(NULL, elem_nds_old, grad_shp_old);
        elem_djac_older = get_tet_linear_djac_gradPhi2(NULL, elem_nds_older, grad_shp_older);
        
        grad_phi_dot_v(grad_shp, elem_vel, &elem_grad_u, &elem_grad_v, &elem_grad_w, nnodes);
        grad_phi_f(grad_shp, elem_pressure, &elem_grad_p, nnodes); // should be a quadratic gradient
        grad_phi_f(grad_shp, elem_density, &elem_grad_den, nnodes);
        
        elem_volume = elem_djac;
        elem_volume_old = elem_djac_old;
        elem_volume_older = elem_djac_older;
    } else {
        elem_volume = get_triprism_volume(elem_nds);
        elem_volume_old = get_triprism_volume(elem_nds_old);
        elem_volume_older = get_triprism_volume(elem_nds_older);
    }
    
    // calculate 2d surface normals and jacobians
    // note :: bed dpl does not affect 2d surface element
    // note :: this is always on triangles
    double elem2d_sur_dpl[NDONTRI];
    global_to_local_dbl(sw->displacement, elem2d_sur_dpl, elem2d_sur->nodes, elem2d_sur->nnodes);
    if (perturb_var == PERTURB_DPL) {
        for(i=0; i<elem2d_sur->nnodes; i++) {
            if (elem2d_sur->nodes[i] == perturbed_surf_node) {
                elem2d_sur_dpl[i] += perturb_sign * sw->dpl_perturbation[elem2d_sur->nodes[i]];
            }
        }
    }
    SVECT elem_nds_sur[NDONTRI];
    for (i=0; i<NDONTRI; i++) {
        elem_nds_sur[i].x = grid->node[ elem2d_sur->nodes[i] ].x;
        elem_nds_sur[i].y = grid->node[ elem2d_sur->nodes[i] ].y;
        elem_nds_sur[i].z = grid->node[ elem2d_sur->nodes[i] ].z + elem2d_sur_dpl[i];
    }
    get_triangle_linear_djac_nrml_gradPhi(elem2d_sur, NULL, elem_nds_sur);
    
    // average column depth on 3 nodes of surface triangle
    double elem_avg_depth = tl_find_avg_column_depth(grid, elem2d_sur->nodes, sw->displacement);
    if (perturb_var == PERTURB_DPL) {
        elem_avg_depth += perturb_sign * one_3 * sw->dpl_perturbation[perturbed_surf_node];
    }
#ifdef _SEDIMENT // if sediment is displacing, add to depth
    double elem_bed_dpl[MAX_NNODES_ON_ELEM2D];
    global_to_local_dbl(sw->bed_displacement, elem_bed_dpl, elem2d_bed->nodes, elem2d_bed->nnodes);
    elem_avg_depth += one_3 * (elem_bed_dpl[0] + elem_bed_dpl[1] + elem_bed_dpl[2]);
#endif
    
#ifdef _DEBUG
    if (elem_avg_depth <= 0) {
        DEBUG_LOCAL = ON;
        printf("\n");
        printf("surface node 1 id: %d surface node z: %*.*e surface node dpl: %*.*e\n",elem2d_sur->nodes[0],
               printFieldWidth,printPrecision,elem_nds_sur[0].z,
               printFieldWidth,printPrecision,sw->displacement[elem2d_sur->nodes[0]]);
        printf("surface node 2 id: %d surface node z: %*.*e surface node dpl: %*.*e\n",elem2d_sur->nodes[1],
               printFieldWidth,printPrecision,elem_nds_sur[1].z,
               printFieldWidth,printPrecision,sw->displacement[elem2d_sur->nodes[1]]);
        printf("surface node 3 id: %d surface node z: %*.*e surface node dpl: %*.*e\n",elem2d_sur->nodes[2],
               printFieldWidth,printPrecision,elem_nds_sur[2].z,
               printFieldWidth,printPrecision,sw->displacement[elem2d_sur->nodes[2]]);
    }
#endif 
    
    // depth averaged velocities
    SVECT2D elem_depth_avg_vel[NDONTRI], elem_avg_davel;
    for (i=0; i<elem2d_sur->nnodes; i++) {
        elem_depth_avg_vel[i].x = sw->depth_avg_vel[ grid->nodeID_3d_to_2d_sur[elem2d_sur->nodes[i]] ].x;
        elem_depth_avg_vel[i].y = sw->depth_avg_vel[ grid->nodeID_3d_to_2d_sur[elem2d_sur->nodes[i]] ].y;
    }
    elem_avg_davel = svect2d_average_array(elem_depth_avg_vel, elem2d_sur->nnodes); // surface (triangle) elementally averaged depth-averaged velocity
    
    // coriolis force
    double rads = DEGREE2RAD * mat.coriolis;
    double elem_angular_speed = 2. * EARTH_ROTAT * sin(rads);
    
    // dump velocities into doubles for function calls
    double u[nnodes], v[nnodes], w[nnodes], u_old[nnodes], v_old[nnodes], w_old[nnodes], u_older[nnodes], v_older[nnodes], w_older[nnodes];
    double ur[nnodes], vr[nnodes], wr[nnodes];
    dumpVector(elem_vel, nnodes, u, v, w);
    dumpVector(elem_vel_old,  nnodes, u_old, v_old, w_old);
    dumpVector(elem_vel_older,  nnodes, u_older, v_older, w_older);
    dumpVector(elem_rel_vel, nnodes, ur, vr, wr);
    
    double grad_shp_x[nnodes], grad_shp_y[nnodes], grad_shp_z[nnodes];
    if (isElementTetrahedron == TRUE) {
        dumpVector(grad_shp, nnodes, grad_shp_x, grad_shp_y, grad_shp_z);
    }
    
    // OLD ADH WAY :: **************************************************
    double elem_djac_12dt = 0., elem_djac_32dt = 0.;
    double elem_dpl_12dt[nnodes], elem_dpl_32dt[nnodes];
    SVECT elem_vel_12dt[nnodes], elem_vel_32dt[nnodes], elem_nds_12dt[nnodes_quad], elem_nds_32dt[nnodes_quad];
    elem_get_tposition(elem_dpl_12dt, elem_dpl_old, elem_dpl_older, mod->tau_temporal, dt, dt,nnodes);
    elem_get_tposition(elem_dpl_32dt, elem_dpl, elem_dpl_old, mod->tau_temporal, dt, dt,nnodes);
    elem_get_tposition_vect(elem_vel_12dt, elem_vel_old, elem_vel_older, mod->tau_temporal, dt, dt, nnodes);
    elem_get_tposition_vect(elem_vel_32dt, elem_vel, elem_vel_old, mod->tau_temporal, dt, dt, nnodes);
    elem_get_midpt_locations2(elem_nodes, elem_dpl_12dt, elem_nds_12dt, nnodes, nnodes_quad, elem3d->edges);
    elem_get_midpt_locations2(elem_nodes, elem_dpl_32dt, elem_nds_32dt, nnodes, nnodes_quad, elem3d->edges);
    double u_12dt[nnodes], v_12dt[nnodes], w_12dt[nnodes], u_32dt[nnodes], v_32dt[nnodes], w_32dt[nnodes];
    dumpVector(elem_vel_12dt, nnodes, u_12dt, v_12dt, w_12dt);
    dumpVector(elem_vel_32dt, nnodes, u_32dt, v_32dt, w_32dt);
    // *******************************************************************
    
    // calculate elemental averages
    SVECT2D elem_avg_strong_mom;
    SVECT elem_avg_vel, elem_avg_vel_old, elem_avg_vel_older, elem_avg_vel_rel;
    SVECT elem_avg_grad_density,elem_avg_grad_u, elem_avg_grad_v, elem_avg_grad_w;
    double elem_avg_prs, elem_avg_density, elem_avg_dudt, elem_avg_dvdt;
    double elem_avg_strong_continuity, elem_avg_u, elem_avg_v, elem_avg_w, elem_avg_ur, elem_avg_vr, elem_avg_wr;
    double elem_avg_u_old, elem_avg_v_old, elem_avg_u_older, elem_avg_v_older;
    if (isElementTetrahedron == TRUE) {
        
        /* Gajanan gkc bug fix - adding below to enable turbulence calculations and remove nan errors due to undefined elem_avg_grad_u/v/w */
        elem_avg_grad_u.x = elem_grad_u.x;
        elem_avg_grad_u.y = elem_grad_u.y;
        elem_avg_grad_u.z = elem_grad_u.z;
        
        elem_avg_grad_v.x = elem_grad_v.x;
        elem_avg_grad_v.y = elem_grad_v.y;
        elem_avg_grad_v.z = elem_grad_v.z;
        
        elem_avg_grad_w.x = elem_grad_w.x;
        elem_avg_grad_w.y = elem_grad_w.y;
        elem_avg_grad_w.z = elem_grad_w.z;
        /* Gajanan gkc bug fix - adding above - Please check this, Corey! */
        
        elem_avg_grad_density.x = elem_grad_den.x;
        elem_avg_grad_density.y = elem_grad_den.y;
        elem_avg_grad_density.z = elem_grad_den.z;
        
        // quadratic pressure average
        double sum1 = elem_pressure[0] + elem_pressure[1] + elem_pressure[2] + elem_pressure[3];
        double sum2 = elem_pressure[4] + elem_pressure[5] + elem_pressure[6] + elem_pressure[7] + elem_pressure[8] + elem_pressure[9];
        elem_avg_prs = one_20 * (-sum1 + 4. * sum2);
        
        elem_avg_vel_old = svect_average_array(elem_vel_old, nnodes);
        elem_avg_vel_older = svect_average_array(elem_vel_older, nnodes);
        elem_avg_vel = svect_average_array(elem_vel, nnodes);
        elem_avg_vel_rel = svect_average_array(elem_rel_vel, nnodes);
        ELEM3D_LOCAL_AVG(elem_density, elem_avg_density);
        
        elem_avg_dudt = (dtTerm1 * elem_avg_vel.x + dtTerm2 * elem_avg_vel_old.x + dtTerm3 * elem_avg_vel_older.x)/dt;
        elem_avg_dvdt = (dtTerm1 * elem_avg_vel.y + dtTerm2 * elem_avg_vel_old.y + dtTerm3 * elem_avg_vel_older.y)/dt;
        
        // note :: gradients and djac are constants across tet element, so avg is just value
        elem_avg_strong_continuity = elem_grad_u.x + elem_grad_v.y + elem_grad_w.z;
        elem_avg_strong_mom.x = elem_avg_dudt + svect_dotp(elem_avg_vel_rel, elem_grad_u) - elem_avg_vel.y * elem_angular_speed;
        elem_avg_strong_mom.y = elem_avg_dvdt + svect_dotp(elem_avg_vel_rel, elem_grad_v) + elem_avg_vel.x * elem_angular_speed;
        
        // note : grad_p should really be linear here, since p is quadratic
        if (PRESSURE_FLAG == ON) {
            elem_avg_strong_mom.x += g * elem_grad_p.x;
            elem_avg_strong_mom.y += g * elem_grad_p.y;
        }
    } else {
        
        elem_avg_density = integrate_triPrism_f(elem_nds, 1./elem_volume, elem_density);
        elem_avg_prs = integrate_triPrism_f(elem_nds, 1./elem_volume, elem_pressure); // linear pressure average
        elem_avg_u =   integrate_triPrism_f(elem_nds, 1./elem_volume, u);
        elem_avg_v =   integrate_triPrism_f(elem_nds, 1./elem_volume, v);
        elem_avg_w =   integrate_triPrism_f(elem_nds, 1./elem_volume, w);
        elem_avg_ur =  elem_avg_u; // no lateral grid movement
        elem_avg_vr =  elem_avg_v; // no lateral grid movement
        elem_avg_wr =      integrate_triPrism_f(elem_nds, 1./elem_volume, wr);
        elem_avg_u_old =   integrate_triPrism_f(elem_nds_old, 1./elem_volume_old, u_old);
        elem_avg_v_old =   integrate_triPrism_f(elem_nds_old, 1./elem_volume_old, v_old);
        elem_avg_u_older = integrate_triPrism_f(elem_nds_older, 1./elem_volume_older, u_older);
        elem_avg_v_older = integrate_triPrism_f(elem_nds_older, 1./elem_volume_older, v_older);
        elem_avg_vel.x = elem_avg_u;
        elem_avg_vel.y = elem_avg_v;
        elem_avg_vel.z = elem_avg_w;
        elem_avg_vel_rel.x = elem_avg_ur;
        elem_avg_vel_rel.y = elem_avg_vr;
        elem_avg_vel_rel.z = elem_avg_wr;
        elem_avg_grad_density = integrate_triPrism_df_full(elem_nds, 1./elem_volume, elem_density);
        elem_avg_grad_u = integrate_triPrism_df_full(elem_nds, 1./elem_volume, u);
        elem_avg_grad_v = integrate_triPrism_df_full(elem_nds, 1./elem_volume, v);
        elem_avg_grad_w = integrate_triPrism_df_full(elem_nds, 1./elem_volume, w);
        
        elem_avg_dudt = (dtTerm1 * elem_avg_u + dtTerm2 * elem_avg_u_old + dtTerm3 * elem_avg_u_older)/dt;
        elem_avg_dvdt = (dtTerm1 * elem_avg_v + dtTerm2 * elem_avg_v_old + dtTerm3 * elem_avg_v_older)/dt;
        
        elem_avg_strong_continuity = elem_avg_grad_u.x + elem_avg_grad_v.y + elem_avg_grad_w.z;
        elem_avg_strong_mom.x = elem_avg_dudt + svect_dotp(elem_avg_vel_rel, elem_avg_grad_u) - elem_avg_v * elem_angular_speed;
        elem_avg_strong_mom.y = elem_avg_dvdt + svect_dotp(elem_avg_vel_rel, elem_avg_grad_v) + elem_avg_u * elem_angular_speed;
        
        // OLD WAY :: seems to work better, at least for corilois verification
        double elem_avg_u_32dt =   integrate_triPrism_f(elem_nds, 1./elem_volume, u_32dt);
        double elem_avg_v_32dt =   integrate_triPrism_f(elem_nds, 1./elem_volume, v_32dt);
        double elem_avg_u_12dt =   integrate_triPrism_f(elem_nds, 1./elem_volume, u_12dt);
        double elem_avg_v_12dt =   integrate_triPrism_f(elem_nds, 1./elem_volume, v_12dt);
        elem_avg_dudt = (elem_avg_u_32dt - elem_avg_u_12dt)/dt;
        elem_avg_dvdt = (elem_avg_v_32dt - elem_avg_v_12dt)/dt;
        elem_avg_strong_mom.x = (elem_avg_u_32dt - elem_avg_u_12dt)/dt + (elem_avg_ur*elem_avg_grad_u.x + elem_avg_vr*elem_avg_grad_u.y + elem_avg_wr*elem_avg_grad_u.z) - elem_avg_v * elem_angular_speed;
        elem_avg_strong_mom.y = (elem_avg_v_32dt - elem_avg_v_12dt)/dt + (elem_avg_ur*elem_avg_grad_v.x + elem_avg_vr*elem_avg_grad_v.y + elem_avg_wr*elem_avg_grad_v.z) + elem_avg_u * elem_angular_speed;
        
        // note : grad_p should really be linear here, since p is quadratic
        if (PRESSURE_FLAG == ON) {
            elem_avg_strong_mom.x += integrate_triPrism_df(elem_nds, g/elem_volume, elem_pressure, 1);
            elem_avg_strong_mom.y += integrate_triPrism_df(elem_nds, g/elem_volume, elem_pressure, 2);
        }
        
        //printf("temporal ::  %*.*e \n",printFieldWidth,printPrecision, (elem_avg_u_32dt - elem_avg_u_12dt)/dt);
        //printf("advection :: %*.*e \n",printFieldWidth,printPrecision, elem_avg_ur*elem_avg_grad_u.x + elem_avg_vr*elem_avg_grad_u.y + elem_avg_wr*elem_avg_grad_u.z);
        //printf("coriolis ::  %*.*e \n",printFieldWidth,printPrecision, -elem_avg_v * elem_angular_speed);
        //printf("pressure ::  %*.*e \n",printFieldWidth,printPrecision, integrate_triPrism_df(elem_nds, g/elem_volume, elem_pressure, 1));
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        printf("SW-HVEL BODY RESID :: ie: %d \t dt: %*.*f \t volume: %*.*f \t volume_old: %*.*f \t volume_older: %*.*f",
               ie,
               printFieldWidth,printPrecision,dt,
               printFieldWidth,printPrecision,elem_volume,
               printFieldWidth,printPrecision,elem_volume_old,
               printFieldWidth,printPrecision,elem_volume_older);
        if (perturb_var == PERTURB_DPL) {
            printf("\t perturbing DPL  || node: %d || perturbation: %-*.*e",elem3d->nodes[perturb_node],printFieldWidth,printPrecision,perturb_sign*perturbation);
        } else if (perturb_var == PERTURB_U) {
            printf("\t perturbing U    || node: %d || perturbation: %-*.*e",elem3d->nodes[perturb_node],printFieldWidth,printPrecision,perturb_sign*perturbation);
        } else if (perturb_var == PERTURB_V) {
            printf("\t perturbing V    || node: %d || perturbation: %-*.*e",elem3d->nodes[perturb_node],printFieldWidth,printPrecision,perturb_sign*perturbation);
        }
        printf("\n");
        //selem3d_printScreen(elem3d);
        for (i=0; i<nnodes_quad; i++) printf("nds %d :: x: %-*.*f y: %-*.*f z: %-*.*f prs: %-*.*f\n",
                                             i,
                                             printFieldWidth,printPrecision,elem_nds[i].x,
                                             printFieldWidth,printPrecision,elem_nds[i].y,
                                             printFieldWidth,printPrecision,elem_nds[i].z,
                                             printFieldWidth,printPrecision,elem_pressure[i] );
        printScreen_debug_vec("node locations old: ",elem_nds_old, nnodes);
        printScreen_debug_vec("node locations older: ",elem_nds_older, nnodes);
        printScreen_debug2_dbl("elem_displacement", elem_dpl, nnodes, elem3d->nodes);
        printScreen_debug2_dbl("elem_displacement_old", elem_dpl_old, nnodes, elem3d->nodes);
        printScreen_debug2_dbl("elem_displacement_older", elem_dpl_older, nnodes, elem3d->nodes);
        printScreen_debug_svect("elem_vel", elem_vel, nnodes, elem3d->nodes);
        printScreen_debug_svect("elem_vel_old", elem_vel_old, nnodes, elem3d->nodes);
        printScreen_debug_svect("elem_vel_older", elem_vel_older, nnodes, elem3d->nodes);
        printScreen_debug_svect("elem_grid_vel", elem_grid_vel, nnodes, elem3d->nodes);
        printScreen_debug_svect("elem_rel_vel", elem_rel_vel, nnodes, elem3d->nodes);
        printScreen_debug2_dbl("elem_density", elem_density, nnodes, elem3d->nodes);
        
        if (isElementTetrahedron == TRUE) {
            printScreen_debug_svect("grad_shp", grad_shp, nnodes, elem3d->nodes);
            printf("elem_avg_prs: %-*.*e\n",printFieldWidth,printPrecision,elem_avg_prs);
            printf("elem_avg_density: %-*.*e\n",printFieldWidth,printPrecision,elem_avg_density);
            printf("elem_avg_vel:     x = %-*.*e \t y = %-*.*e \t z = %-*.*e\n",
                   printFieldWidth,printPrecision,elem_avg_vel.x,
                   printFieldWidth,printPrecision,elem_avg_vel.y,
                   printFieldWidth,printPrecision,elem_avg_vel.z);
            printf("elem_avg_vel_rel: x = %-*.*e \t y = %-*.*e \t z = %-*.*e\n",
                   printFieldWidth,printPrecision,elem_avg_vel_rel.x,
                   printFieldWidth,printPrecision,elem_avg_vel_rel.y,
                   printFieldWidth,printPrecision,elem_avg_vel_rel.z);
        } else {
            printf("elem_avg_prs: %-*.*e\n",printFieldWidth,printPrecision,elem_avg_prs);
            printf("elem_avg_density: %-*.*e\n",printFieldWidth,printPrecision,elem_avg_density);
            printf("elem_avg_vel:     x = %-*.*e \t y = %-*.*e \t z = %-*.*e\n",
                   printFieldWidth,printPrecision,elem_avg_vel.x,
                   printFieldWidth,printPrecision,elem_avg_vel.y,
                   printFieldWidth,printPrecision,elem_avg_vel.z);
            printf("elem_avg_vel_rel: x = %-*.*e \t y = %-*.*e \t z = %-*.*e\n",
                   printFieldWidth,printPrecision,elem_avg_vel_rel.x,
                   printFieldWidth,printPrecision,elem_avg_vel_rel.y,
                   printFieldWidth,printPrecision,elem_avg_vel_rel.z);
            printf("elem_avg_grad_u:  x = %-*.*e \t y = %-*.*e \t z = %-*.*e\n",
                   printFieldWidth,printPrecision,elem_avg_grad_u.x,
                   printFieldWidth,printPrecision,elem_avg_grad_u.y,
                   printFieldWidth,printPrecision,elem_avg_grad_u.z);
            printf("elem_avg_grad_v:  x = %-*.*e \t y = %-*.*e \t z = %-*.*e\n",
                   printFieldWidth,printPrecision,elem_avg_grad_v.x,
                   printFieldWidth,printPrecision,elem_avg_grad_v.y,
                   printFieldWidth,printPrecision,elem_avg_grad_v.z);
            printf("elem_avg_grad_w:  x = %-*.*e \t y = %-*.*e \t z = %-*.*e\n",
                   printFieldWidth,printPrecision,elem_avg_grad_w.x,
                   printFieldWidth,printPrecision,elem_avg_grad_w.y,
                   printFieldWidth,printPrecision,elem_avg_grad_w.z);
        }
        printf("elem_avg_depth: %-*.*e\n",printFieldWidth,printPrecision,elem_avg_depth);
        printf("elem_avg_davel: %-*.*e \t %-*.*e\n",printFieldWidth,printPrecision,elem_avg_davel.x,printFieldWidth,printPrecision,elem_avg_davel.y);
        printf("elemental averages: continuity equation: %-*.*e || momentum equation: %-*.*e %-*.*e\n",
               printFieldWidth,printPrecision,elem_avg_strong_continuity,
               printFieldWidth,printPrecision,elem_avg_strong_mom.x,
               printFieldWidth,printPrecision,elem_avg_strong_mom.y);
        
        assert(elem_avg_depth > 0);
    }
    if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                  FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    
    // initialize right hand side (to be returned)
    dof3_init_array(elem_rhs, nnodes);
    dof3_init_array(elem_rhs_temporal, NDONPRISM);
    
    // temporary rhs
    double integral[nnodes], integral_X[nnodes], integral_Y[nnodes], integral_Z[nnodes];
    DOF_3 rhs[nnodes]; dof3_init_array(rhs, nnodes);
    
    // zero SUPG da-cont terms for wvel addition
    for (i=0; i<nnodes; i++) mod->sw->d3->elem_rhs_supg_dacont[i][ie] = 0.;
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                              DEPTH-AVERAGED CONTINUITY CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the depth-averaged continuity addition to the 3D SW elemental residual. \n
     * \note
     *********************************************************************************************/
    
    sarray_init_dbl(integral,nnodes);
    if (isElementTetrahedron == TRUE) {
        integrate_tetrahedron_gradPhi_dot_v(grad_shp, elem_djac, -dt, elem_vel, integral);
    } else {
        integrate_triPrism_gradPhi_dot_v(elem_nds, -dt, elem_vel, integral);
        
        //        // old prism code way
        //        integrate_triPrism_dphi_f(elem_nds, dt, u, integral_X, 1);
        //        integrate_triPrism_dphi_f(elem_nds, dt, v, integral_Y, 2);
        //        integrate_triPrism_dphi_f(elem_nds, dt, w, integral_Z, 3);
        //        for (i=0; i<nnodes; i++) integral[i] = integral_X[i] + integral_Y[i] + integral_Z[i];
    }
    for (i=0; i<nnodes; i++) elem_rhs[i].c_eq += integral[i];
    
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON || dof3_debug_continue(elem_rhs, nnodes)) {
        rhs_3dof("3D HVEL BODY || DA-C: ",nnodes, ie, elem3d->nodes, elem_rhs);
    }
    dof3_debug(rhs, nnodes, __FILE__, __LINE__);
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   MOMENTUM TEMPORAL CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the temporal addition to the 3D SW momentum elemental residual. \n
     * \note
     *********************************************************************************************/
    
    // ++ t(n+1) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_hvel_temporal(ie, nnodes, elem3d, elem_nds, elem_djac, u, v, dtTerm1, elem_rhs, "3D HVEL || TEMPORAL T(N+1) ", DEBUG_LOCAL);
    
    // ++ t(n) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_hvel_temporal(ie, nnodes, elem3d, elem_nds_old, elem_djac_old, u_old, v_old, dtTerm2, elem_rhs, "3D HVEL || TEMPORAL T(N) ", DEBUG_LOCAL);
    
    // ++ t(n-1) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_hvel_temporal(ie, nnodes, elem3d, elem_nds_older, elem_djac_older, u_older, v_older, dtTerm3, elem_rhs, "3D HVEL || TEMPORAL T(N-1) ", DEBUG_LOCAL);
    
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        rhs_3dof("3D SW BODY || TOTAL TEMPORAL RHS: ",nnodes, ie, elem3d->nodes, elem_rhs_temporal);
    }
    dof3_debug(rhs, nnodes, __FILE__, __LINE__);
#endif
    
    
    
    // OLD ADH WAY
    
    //    /* ------------------------------------------------------------------------------------------*/
    //    /* SECOND ORDER TERM 1 :: alpha * (3/2 u^{n+1) - 1/2 u^{n})                                  */
    //    /* ------------------------------------------------------------------------------------------*/
    //
    //    if (mod->tau_temporal > 1e-6) {
    //
    //        constant = mod->tau_temporal;
    //        if (isElementTetrahedron == TRUE) {
    //            //integrate_tetrahedron_phi_f(constant * elem_djac_32dt, u_32dt, integral_X);
    //            //integrate_tetrahedron_phi_f(constant * elem_djac_32dt, v_32dt, integral_Y);
    //        }
    //        else {
    //            double integral_X[nnodes];  sarray_init_dbl(integral_X,nnodes);
    //            double integral_Y[nnodes];  sarray_init_dbl(integral_Y,nnodes);
    //            integrate_triPrism_phi_f(elem_nds_32dt, constant, u_32dt, integral_X);
    //            integrate_triPrism_phi_f(elem_nds_32dt, constant, v_32dt, integral_Y);
    //        }
    //        for (i=0; i<nnodes; i++) {
    //            elem_rhs[i].x_eq += integral_X[i];
    //            elem_rhs[i].y_eq += integral_Y[i];
    //        }
    //#ifdef _DEBUG
    //        if (DEBUG_LOCAL == ON || dof3_debug_continue(elem_rhs_temporal, nnodes)) {
    //            dof3_init_array(rhs, nnodes);
    //            for (i=0; i<nnodes; i++) {
    //                rhs[i].x_eq = integral_X[i];
    //                rhs[i].y_eq = integral_Y[i];
    //            }
    //            rhs_3dof("temporal 1", nnodes, ie, elem3d->nodes, rhs);
    //        }
    //#endif
    //
    //        /* ------------------------------------------------------------------------------------------*/
    //        /* SECOND ORDER TERM 2 ::  -alpha * (3/2 u^{n} - 1/2 u^{n-1})                                */
    //        /* ------------------------------------------------------------------------------------------*/
    //
    //        constant = -mod->tau_temporal;
    //        if (isElementTetrahedron == TRUE) {
    //            //integrate_tetrahedron_phi_f(constant * elem_djac_12dt, u_12dt, integral_X);
    //            //integrate_tetrahedron_phi_f(constant * elem_djac_12dt, v_12dt, integral_Y);
    //        }
    //        else {
    //            double integral_X[nnodes];  sarray_init_dbl(integral_X,nnodes);
    //            double integral_Y[nnodes];  sarray_init_dbl(integral_Y,nnodes);
    //            integrate_triPrism_phi_f(elem_nds_12dt, constant, u_12dt, integral_X);
    //            integrate_triPrism_phi_f(elem_nds_12dt, constant, v_12dt, integral_Y);
    //        }
    //        for (i=0; i<nnodes; i++) {
    //            elem_rhs[i].x_eq += integral_X[i];
    //            elem_rhs[i].y_eq += integral_Y[i];
    //        }
    //#ifdef _DEBUG
    //        if (DEBUG_LOCAL == ON || dof3_debug_continue(elem_rhs_temporal, nnodes)) {
    //            dof3_init_array(rhs, nnodes);
    //            for (i=0; i<nnodes; i++) {
    //                rhs[i].x_eq = integral_X[i];
    //                rhs[i].y_eq = integral_Y[i];
    //            }
    //            rhs_3dof("temporal 2", nnodes, ie, elem3d->nodes, rhs);
    //        }
    //#endif
    //    }
    //
    //    /* ------------------------------------------------------------------------------------------*/
    //    /* FIRST ORDER TERM 1 ::  (1-alpha) * u^{n+1)                                                */
    //    /* ------------------------------------------------------------------------------------------*/
    //
    //    sarray_init_dbl(integral_X,nnodes);
    //    sarray_init_dbl(integral_Y,nnodes);
    //    constant = 1. - mod->tau_temporal;
    //    if (isElementTetrahedron == TRUE) {
    //        //integrate_tetrahedron_phi_f(constant * elem_djac, u, integral_X);
    //        //integrate_tetrahedron_phi_f(constant * elem_djac, v, integral_Y);
    //    }
    //    else {
    //        integrate_triPrism_phi_f(elem_nds, constant, u, integral_X);
    //        integrate_triPrism_phi_f(elem_nds, constant, v, integral_Y);
    //    }
    //    for (i=0; i<nnodes; i++) {
    //        elem_rhs[i].x_eq += integral_X[i];
    //        elem_rhs[i].y_eq += integral_Y[i];
    //    }
    //#ifdef _DEBUG
    //if (DEBUG_LOCAL == ON || dof3_debug_continue(elem_rhs_temporal, nnodes)) {
    //        dof3_init_array(rhs, nnodes);
    //        for (i=0; i<nnodes; i++) {
    //            rhs[i].x_eq = integral_X[i];
    //            rhs[i].y_eq = integral_Y[i];
    //        }
    //        rhs_3dof("temporal 3", nnodes, ie, elem3d->nodes, rhs);
    //    }
    //#endif
    //
    //    /* ------------------------------------------------------------------------------------------*/
    //    /* FIRST ORDER TERM 2 ::  -(1-alpha) * u^{n}                                                 */
    //    /* ------------------------------------------------------------------------------------------*/
    //
    //    sarray_init_dbl(integral_X,nnodes);
    //    sarray_init_dbl(integral_Y,nnodes);
    //    constant = -(1. - mod->tau_temporal);
    //    if (isElementTetrahedron == TRUE) {
    //        //integrate_tetrahedron_phi_f(constant * elem_djac_old, u_old, integral_X);
    //        //integrate_tetrahedron_phi_f(constant * elem_djac_old, v_old, integral_Y);
    //    }
    //    else {
    //        integrate_triPrism_phi_f(elem_nds_old, constant, u_old, integral_X);
    //        integrate_triPrism_phi_f(elem_nds_old, constant, v_old, integral_Y);
    //    }
    //    for (i=0; i<nnodes; i++) {
    //        elem_rhs[i].x_eq += integral_X[i];
    //        elem_rhs[i].y_eq += integral_Y[i];
    //    }
    //#ifdef _DEBUG
    //if (DEBUG_LOCAL == ON || dof3_debug_continue(elem_rhs_temporal, nnodes)) {
    //        dof3_init_array(rhs, nnodes);
    //        for (i=0; i<nnodes; i++) {
    //            rhs[i].x_eq = integral_X[i];
    //            rhs[i].y_eq = integral_Y[i];
    //        }
    //        rhs_3dof("temporal 4", nnodes, ie, elem3d->nodes, rhs);
    //    }
    //#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    ADVECTION CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the advection addition to the 3D SW momentum elemental residual. \n
     * \note
     *********************************************************************************************/
    
#ifdef _DEBUG
    if (debug.no_advection == OFF) {
#endif
        sarray_init_dbl(integral_X,nnodes);
        sarray_init_dbl(integral_Y,nnodes);
        if (isElementTetrahedron) {
            integrate_tetrahedron_gradPhi_dot_f_v(grad_shp, elem_djac, -dt, u, elem_rel_vel, integral_X);
            integrate_tetrahedron_gradPhi_dot_f_v(grad_shp, elem_djac, -dt, v, elem_rel_vel, integral_Y);
        } else {
            
            integrate_triPrism_gradPhi_dot_f_v(elem_nds, -dt, u, elem_rel_vel, integral_X);
            integrate_triPrism_gradPhi_dot_f_v(elem_nds, -dt, v, elem_rel_vel, integral_Y);
            
            //            // old prism code way
            //            double xx[nnodes], yy[nnodes], zz[nnodes];
            //            sarray_init_dbl(xx,nnodes);
            //            sarray_init_dbl(yy,nnodes);
            //            sarray_init_dbl(zz,nnodes);
            //            integrate_triPrism_dphi_f_g(elem_nds, -dt, ur, u, xx, 1);
            //            integrate_triPrism_dphi_f_g(elem_nds, -dt, vr, u, yy, 2);
            //            integrate_triPrism_dphi_f_g(elem_nds, -dt, wr, u, zz, 3);
            //            for (i=0; i<nnodes; i++) integral_X[i] = xx[i] + yy[i] + zz[i];
            //
            //            sarray_init_dbl(xx,nnodes);
            //            sarray_init_dbl(yy,nnodes);
            //            sarray_init_dbl(zz,nnodes);
            //            integrate_triPrism_dphi_f_g(elem_nds, -dt, ur, v, xx, 1);
            //            integrate_triPrism_dphi_f_g(elem_nds, -dt, vr, v, yy, 2);
            //            integrate_triPrism_dphi_f_g(elem_nds, -dt, wr, v, zz, 3);
            //            for (i=0; i<nnodes; i++) integral_Y[i] = xx[i] + yy[i] + zz[i];
        }
        for (i=0; i<nnodes; i++) {
            elem_rhs[i].x_eq += integral_X[i];
            elem_rhs[i].y_eq += integral_Y[i];
        }
#ifdef _DEBUG
    }
    if (DEBUG_LOCAL == ON || dof3_debug_continue(elem_rhs, nnodes)) {
        for (i=0; i<nnodes; i++) {rhs[i].c_eq = 0.; rhs[i].x_eq = integral_X[i]; rhs[i].y_eq = integral_Y[i];}
        rhs_3dof("3D HVEL BODY || ADVECTION: ",nnodes, ie, elem3d->nodes, rhs);
    }
    dof3_debug(rhs, nnodes, __FILE__, __LINE__);
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    PRESSURE CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the pressure addition to the 3D SW momentum elemental residual. \n
     * \note
     *********************************************************************************************/
    
    if (PRESSURE_FLAG) {
        
        sarray_init_dbl(integral_X,nnodes);
        sarray_init_dbl(integral_Y,nnodes);
        if (isElementTetrahedron) {
            // constant pressure
            integrate_tetrahedron_fi(elem_djac, elem_avg_prs, grad_shp_x, integral_X);
            integrate_tetrahedron_fi(elem_djac, elem_avg_prs, grad_shp_y, integral_Y);
        } else {
            // constant pressure || note ::  is not stable for some reason ...
            //integrate_triPrism_dphi(elem_nds, -dt * g * elem_avg_prs, integral_X, 1);
            //integrate_triPrism_dphi(elem_nds, -dt * g * elem_avg_prs, integral_Y, 2);
            
            // linear pressure || note :: works, but seems to cap out early on NTL convergence
            //integrate_triPrism_dphi_f(elem_nds, -dt * g, elem_pressure, integral_X, 1);
            //integrate_triPrism_dphi_f(elem_nds, -dt * g, elem_pressure, integral_Y, 2);
            
            // linear pressure using quadrature ||note ::  slow, but the convergence rates seem better than analytic, probably due to precision
            //integrate_triPrism_dphi_f_quadrature(elem_nds, grid->quad_prism, -dt * g, elem_pressure, integral_X, integral_Y, 3);
            integrate_triPrism_dphi_f_quadrature(elem_nds, grid->quad_prism, 1., elem_pressure, integral_X, integral_Y, 3);
        }
        for (i=0; i<nnodes; i++) {
            elem_rhs[i].x_eq += -dt * g * integral_X[i];
            elem_rhs[i].y_eq += -dt * g * integral_Y[i];
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON || dof3_debug_continue(elem_rhs, nnodes)) {
            for (i=0; i<nnodes; i++) {
                rhs[i].c_eq = 0.; rhs[i].x_eq = -dt * g * integral_X[i]; rhs[i].y_eq = -dt * g * integral_Y[i];
            }
            rhs_3dof("3D HVEL BODY || PRESSURE: ",nnodes, ie, elem3d->nodes, rhs);
        }
        dof3_debug(rhs, nnodes, __FILE__, __LINE__);
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    DIFFUSION CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the diffusion addition to the 3D SW momentum elemental residual. \n
     * \note
     *********************************************************************************************/
    
    /* Compute requirements for the turbulence models */
    /* Mode 0 = Smagorinkski */
    /* Mode 1 = MY 2 */
    /* Mode 2 = MY 2.5 */
    /* Mode 3 = MY 3 */
    /* Mode 4 = K-e  GSavant */
    /* Modes are open ended, we should be able to include other models using the transport capabilities and Operator-Splitting */
    double evxz = 0., evyz = 0., tur_viscosity = 0., dist_above_bed = 0., v_mag_actual = 0., avg_den = 0.,avg_factor = 0.;
    double exx = 0., exy = 0., exz = 0., eyx = 0., eyy = 0., eyz = 0.;
    double djac_avg = 0.;
    if (mat.turbulence_model_xy == 0) {
        tur_viscosity = tur_smag(elem_avg_grad_u.x, elem_avg_grad_u.y, elem_avg_grad_v.x, elem_avg_grad_v.y, elem_volume, mat.smag_coeff);
    }
    if (mat.turbulence_model_xy == 1) {
        double v_mag_actual = VECT3D_MAG_SAFE(elem_avg_vel);
        tur_viscosity = tur_ws(v_mag_actual, elem_avg_depth, mat.smag_coeff);
    }
    if (mat.turbulence_model_z == 1) {
        double dist_above_bed = elem_avg_depth - elem_avg_prs/g; // this is the distance of element above bed
        evxz = tur_MY_2(g, elem_avg_depth, dist_above_bed, elem_avg_grad_u.z, elem_avg_grad_v.z, elem_avg_grad_density.z, elem_avg_density, mat.supression_func, 0);
        evyz = tur_MY_2(g, elem_avg_depth, dist_above_bed, elem_avg_grad_u.z, elem_avg_grad_v.z, elem_avg_grad_density.z, elem_avg_density, mat.supression_func, 1);
        evxz /= elem_avg_density;
        evyz /= elem_avg_density;
    }
    
    /* Store elemental eddy viscosity */
    grid->hyd_eddy[ie] = evxz + one_3 * (mat.ev.xz + mat.ev.zz + mat.ev.yz);
    
    double multiply = 1.;
    if (elem_avg_grad_density.z > 0) {
        multiply = 10.;
    }
    exx = multiply * (mat.ev.xx + mod->viscosity + tur_viscosity);
    exy = multiply * (mat.ev.xy + mod->viscosity + tur_viscosity);
    exz = multiply * (mat.ev.xz + mod->viscosity + evxz);
    eyx = multiply * (mat.ev.xy + mod->viscosity + tur_viscosity);
    eyy = multiply * (mat.ev.yy + mod->viscosity + tur_viscosity);
    eyz = multiply * (mat.ev.yz + mod->viscosity + evyz);
    
    sarray_init_dbl(integral_X,nnodes);
    sarray_init_dbl(integral_Y,nnodes);
    if (isElementTetrahedron) {
        SVECT diffusive_flux_x, diffusive_flux_y;
        diffusive_flux_x.x = exx * elem_grad_u.x;  diffusive_flux_x.y = exy * elem_grad_u.y;  diffusive_flux_x.z = exz * elem_grad_u.z;
        diffusive_flux_y.x = eyx * elem_grad_v.x;  diffusive_flux_y.y = eyy * elem_grad_v.y;  diffusive_flux_y.z = eyz * elem_grad_v.z;
#ifdef _DEBUG
        if (DEBUG_LOCAL) {
            printf("diffusive_flux_x = %*.*e \t %*.*e \t %*.*e\n",
                   printFieldWidth,printPrecision,diffusive_flux_x.x,
                   printFieldWidth,printPrecision,diffusive_flux_x.y,
                   printFieldWidth,printPrecision,diffusive_flux_x.z);
            printf("diffusive_flux_y = %*.*e \t %*.*e \t %*.*e\n",
                   printFieldWidth,printPrecision,diffusive_flux_y.x,
                   printFieldWidth,printPrecision,diffusive_flux_y.y,
                   printFieldWidth,printPrecision,diffusive_flux_y.z);
        }
#endif
        integrate_tetrahedron_gradPhi_dot_vcon(grad_shp, elem_djac, dt, diffusive_flux_x, integral_X);
        integrate_tetrahedron_gradPhi_dot_vcon(grad_shp, elem_djac, dt, diffusive_flux_y, integral_Y);
    }
    else {
        
        SVECT ex, ey;
        ex.x = exx; ex.y = exy; ex.z = exz;
        ey.x = eyx; ey.y = eyy; ey.z = eyz;
        integrate_triPrism_gradPhi_dot_DgradV(elem_nds, grid->quad_prism, dt, ex, ey, u, v, integral_X, integral_Y);
    }
    for (i=0; i<nnodes; i++) {
        elem_rhs[i].x_eq += integral_X[i];
        elem_rhs[i].y_eq += integral_Y[i];
    }
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON || dof3_debug_continue(elem_rhs, nnodes)) {
        printf("exx = %*.*e; exy = %*.*e; exz = %*.*e\n",printFieldWidth,printPrecision,exx,printFieldWidth,printPrecision,exy,printFieldWidth,printPrecision,exz);
        printf("eyx = %*.*e; eyy = %*.*e; eyz = %*.*e\n",printFieldWidth,printPrecision,eyx,printFieldWidth,printPrecision,eyy,printFieldWidth,printPrecision,eyz);
        for (i=0; i<nnodes; i++) {
            rhs[i].c_eq = 0.; rhs[i].x_eq = integral_X[i]; rhs[i].y_eq = integral_Y[i];
        }
        rhs_3dof("3D HVEL BODY || DIFFUSION: ",nnodes, ie, elem3d->nodes, rhs);
    }
    dof3_debug(rhs, nnodes, __FILE__, __LINE__);
#endif
    
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    CORIOLIS CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the coriolis addition to the 3D SW momentum elemental residual. \n
     * \note
     *********************************************************************************************/
    
    if (isElementTetrahedron == TRUE) {
        integrate_tetrahedron_phi_f(elem_djac, -dt*elem_angular_speed, v, integral_X);
        integrate_tetrahedron_phi_f(elem_djac,  dt*elem_angular_speed, u, integral_Y);
    } else {
        integrate_triPrism_phi_f(elem_nds, -dt*elem_angular_speed, v, integral_X);
        integrate_triPrism_phi_f(elem_nds,  dt*elem_angular_speed, u, integral_Y);
    }
    for (i=0; i<nnodes; i++) {
        elem_rhs[i].x_eq += integral_X[i];
        elem_rhs[i].y_eq += integral_Y[i];
    }
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON || dof3_debug_continue(elem_rhs, nnodes)) {
        for (i=0; i<nnodes; i++) {
            rhs[i].c_eq = 0.; rhs[i].x_eq = integral_X[i]; rhs[i].y_eq = integral_Y[i];
        }
        rhs_3dof("3D HVEL BODY || CORIOLIS: ",nnodes, ie, elem3d->nodes, rhs);
    }
    dof3_debug(rhs, nnodes, __FILE__, __LINE__);
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
    
    DOF_3 rhs_mom_supg[nnodes]; dof3_init_array(rhs_mom_supg, nnodes);
    DOF_3 rhs_dacont_supg[nnodes]; dof3_init_array(rhs_dacont_supg, nnodes);
    
    // get dimensionalization coefficients
    double tau3d = fe_get_supg_tau_sw(nnodes, elem_nds, g, elem_avg_depth, elem_avg_vel_rel.x, elem_avg_vel_rel.y, elem_avg_vel_rel.z, grad_shp_x, grad_shp_y, grad_shp_z, elem_volume, tau_pg, 3, UNSET_INT, UNSET_INT);
    double tau3d_C  = tau3d * g * elem_avg_depth;
    double tau2d = fe_get_supg_tau_sw(nnodes, elem_nds, g, elem_avg_depth, elem_avg_davel.x, elem_avg_davel.y, 0.,NULL, NULL, NULL, elem2d_sur->djac, tau_pg, 2, UNSET_INT, 2);

    
    if (isElementTetrahedron == TRUE) {
        
        double upwind[nnodes];
        
        // depth-averaged continuity supg
        sarray_init_dbl(integral,nnodes);
        for(i=0; i<nnodes; i++){upwind[i] = elem_avg_davel.x * grad_shp_x[i] + elem_avg_davel.y * grad_shp_y[i];}
        integrate_tetrahedron_fi(elem_djac, tau2d * dt * elem_avg_strong_continuity, upwind, integral);
        for (i=0; i<nnodes; i++) {rhs_dacont_supg[i].c_eq = integral[i];}
        
        // depth-averaged continuity supg contribution to momentum
        sarray_init_dbl(integral_X,nnodes);  sarray_init_dbl(integral_Y,nnodes);
        integrate_tetrahedron_fi(elem_djac, tau3d_C * dt * elem_avg_strong_continuity, grad_shp_x, integral_X);
        integrate_tetrahedron_fi(elem_djac, tau3d_C * dt * elem_avg_strong_continuity, grad_shp_y, integral_Y);
        for (i=0; i<nnodes; i++) {
            rhs_dacont_supg[i].x_eq = integral_X[i];
            rhs_dacont_supg[i].y_eq = integral_Y[i];
        }
        
        // momentum supg
        sarray_init_dbl(integral_X,nnodes);
        sarray_init_dbl(integral_Y,nnodes);
        for(i=0; i<nnodes; i++){upwind[i] = svect_dotp(elem_avg_vel_rel, grad_shp[i]);}
        integrate_tetrahedron_fi(elem_djac, tau3d * dt * elem_avg_strong_mom.x, upwind, integral_X);
        integrate_tetrahedron_fi(elem_djac, tau3d * dt * elem_avg_strong_mom.y, upwind, integral_Y);
        for (i=0; i<nnodes; i++) {
            rhs_mom_supg[i].x_eq = integral_X[i];
            rhs_mom_supg[i].y_eq = integral_Y[i];
        }
        
        // momentum supg contribution to depth-averaged continuity
        sarray_init_dbl(integral,nnodes);
        for(i=0; i<nnodes; i++){upwind[i] = grad_shp_x[i] * elem_avg_strong_mom.x + grad_shp_y[i] * elem_avg_strong_mom.y;}
        integrate_tetrahedron_fi(elem_djac, tau2d * dt, upwind, integral);
        for (i=0; i<nnodes; i++) {rhs_mom_supg[i].c_eq = integral[i];}
        
    }  else {
        
        //        // use direct integration
        //        sarray_init_dbl(integral_X,nnodes); sarray_init_dbl(integral_Y,nnodes); sarray_init_dbl(integral_Z,nnodes);
        //        integrate_triPrism_dphi(elem_nds, 1., integral_X, 1);
        //        integrate_triPrism_dphi(elem_nds, 1., integral_Y, 2);
        //        integrate_triPrism_dphi(elem_nds, 1., integral_Z, 3);
        //
        //        // depth-averaged continuity supg
        //        t1 = tau2d * dt * elem_avg_strong_continuity;
        //        for (i=0; i<nnodes; i++) {rhs_dacont_supg[i].c_eq += t1 * (elem_avg_davel.x * integral_X[i] + elem_avg_davel.y * integral_Y[i]);}
        //
        //        // depth-averaged continuity supg contribution to momentum
        //        t1 = tau3d_C * dt * elem_avg_strong_continuity;
        //        for (i=0; i<nnodes; i++){
        //            rhs_dacont_supg[i].x_eq += t1 * integral_X[i];
        //            rhs_dacont_supg[i].y_eq += t1 * integral_Y[i];
        //        }
        //
        //        // momentum supg
        //        t1 = tau3d * dt;
        //        for (i=0; i<nnodes; i++){
        //            rhs_mom_supg[i].x_eq += t1 * elem_avg_strong_mom.x * (elem_avg_vel_rel.x * integral_X[i] + elem_avg_vel_rel.y * integral_Y[i] + elem_avg_vel_rel.z * integral_Z[i]);
        //            rhs_mom_supg[i].y_eq += t1 * elem_avg_strong_mom.y * (elem_avg_vel_rel.x * integral_X[i] + elem_avg_vel_rel.y * integral_Y[i] + elem_avg_vel_rel.z * integral_Z[i]);
        //        }
        //
        //        // momentum supg contribution to depth-averaged continuity
        //        t1 = tau2d * dt;
        //        for (i=0; i<nnodes; i++) {rhs_mom_supg[i].c_eq += t1 * (elem_avg_strong_mom.x * integral_X[i] + elem_avg_strong_mom.y * integral_Y[i]);}
        
        // use quadrature
        double constant_con, constant_conMoM;
        int iqp, quad_order = 2; // quadrature order
        SQUAD_PT *qp = NULL;
        
        for (iqp=0; iqp<quad[quad_order].n; iqp++) {
            qp = &(quad[quad_order].pt[iqp]);
            
            // evalute triangular prism djac and cartesian shape function gradients at quadrature point
            qp->djac = get_triprism_linear_djac_gradPhi(qp->xhat, qp->yhat, qp->zhat, elem_nds, qp->grad_shp);
            
            // evaluate continuity equation at quadrature points
            double qp_grad_u_x = 0.;
            double qp_grad_v_y = 0.;
            double qp_grad_w_z = 0.;
            for (i=0; i<NDONPRISM; i++) {
                qp_grad_u_x += u[i] * qp->grad_shp[i].x;
                qp_grad_v_y += v[i] * qp->grad_shp[i].y;
                qp_grad_w_z += w[i] * qp->grad_shp[i].z;
            }
            double qp_Rcont = qp_grad_u_x + qp_grad_v_y + qp_grad_w_z;
            constant_con    = tau2d    * qp_Rcont;
            constant_conMoM = tau3d_C  * qp_Rcont;
            
            t1 = dt * qp->djac * qp->w; // 1rst-order integrand + gradients = 2nd-order
            for (i=0; i<nnodes; i++) {
                
                // depth-averaged continuity supg
                rhs_dacont_supg[i].c_eq += t1 * constant_con * (elem_avg_davel.x * qp->grad_shp[i].x + elem_avg_davel.y * qp->grad_shp[i].y);
                
                // depth-averaged continuity supg contribution to momentum (these significantly effect convergence)
                rhs_dacont_supg[i].x_eq += t1 * constant_conMoM  * qp->grad_shp[i].x;
                rhs_dacont_supg[i].y_eq += t1 * constant_conMoM  * qp->grad_shp[i].y;
                
                // momentum supg
                //printf("constant_x: %*.*e \t constant_y: %*.*e\n",tau3d * elem_avg_strong_mom.x ,tau3d * elem_avg_strong_mom.y);
                double dum = tau3d * elem_avg_strong_mom.x ;
                rhs_mom_supg[i].x_eq += t1 * dum * svect_dotp(elem_avg_vel_rel, qp->grad_shp[i]);
                dum = tau3d * elem_avg_strong_mom.y;
                rhs_mom_supg[i].y_eq += t1 * dum * svect_dotp(elem_avg_vel_rel, qp->grad_shp[i]);
                
                // momentum supg contribution to depth-averaged continuity || this guys seems to need more stabilization than tets for some reason ...
                rhs_mom_supg[i].c_eq += t1 * tau2d * (qp->grad_shp[i].x * elem_avg_strong_mom.x + qp->grad_shp[i].y * elem_avg_strong_mom.y);
                
            }
        }
    }
    
    for (i=0; i<nnodes; i++) {
        mod->sw->d3->elem_rhs_supg_dacont[i][ie] = rhs_dacont_supg[i].c_eq + rhs_mom_supg[i].c_eq;
        elem_rhs[i].c_eq += mod->sw->d3->elem_rhs_supg_dacont[i][ie];
        elem_rhs[i].x_eq += rhs_dacont_supg[i].x_eq + rhs_mom_supg[i].x_eq;
        elem_rhs[i].y_eq += rhs_dacont_supg[i].y_eq + rhs_mom_supg[i].y_eq;
    }
    
#ifdef _DEBUG
    if (DEBUG_LOCAL || dof3_debug_continue(rhs_dacont_supg, nnodes) || dof3_debug_continue(rhs_mom_supg, nnodes)) {
        printf("SUPG coefficients || tau2d: %*.*e tau3d: %*.*e tau3d_c: %*.*e\n",
               printFieldWidth,printPrecision,tau2d,
               printFieldWidth,printPrecision,tau3d,
               printFieldWidth,printPrecision,tau3d_C);
        //printf("elem_avg_depth: %*.*e \t elem_avg_vel_rel.x: %*.*e \t  elem_avg_vel_rel.y: %*.*e \t  elem_avg_vel_rel.z: %*.*e \t deph_avg_vel: %*.*e %*.*e\n",
        //       elem_avg_depth, elem_avg_vel_rel.x, elem_avg_vel_rel.y, elem_avg_vel_rel.z, elem_avg_davel.x, elem_avg_davel.y);
        rhs_3dof("3D HVEL BODY || DAC SUPG: ", nnodes, ie, elem3d->nodes, rhs_dacont_supg);
        rhs_3dof("3D HVEL BODY || MOM SUPG: ", nnodes, ie, elem3d->nodes, rhs_mom_supg);
        //        dof3_init_array(rhs, nnodes);
        //        for (i=0; i<nnodes; i++) {
        //            rhs[i].c_eq += rhs_mom_supg[i].c_eq;
        //            rhs[i].x_eq += rhs_dacont_supg[i].x_eq + rhs_mom_supg[i].x_eq;
        //            rhs[i].y_eq += rhs_dacont_supg[i].y_eq + rhs_mom_supg[i].y_eq;
        //        }
        //        rhs_3dof("3D HVEL BODY || OLD ADH MOMENTUM: ", nnodes, ie, elem3d->nodes, rhs);
        //
        //        dof3_init_array(rhs, nnodes);
        //        for (i=0; i<nnodes; i++) {
        //            rhs[i].c_eq += rhs_dacont_supg[i].c_eq;
        //        }
        //        rhs_3dof("3D HVEL BODY || OLD ADH DA-CONT: ", nnodes, ie, elem3d->nodes, rhs);
    }
    dof3_debug(rhs_dacont_supg, nnodes, __FILE__, __LINE__);
    dof3_debug(rhs_mom_supg, nnodes, __FILE__, __LINE__);
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
        rhs_3dof("3D SW BODY || TOTAL RHS: ",nnodes, ie, elem3d->nodes, elem_rhs);
        printf("\n");
        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    }
    dof3_debug(elem_rhs, nnodes, __FILE__, __LINE__);
    
    time_t time2;  time(&time2);
    TIME_IN_HVEL_BODY_RESID += difftime(time2,time1);
#endif
    
    //exit(-1);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 3D HVEL body temporal contributions to the shallow water equations.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_hvel_temporal(int ie, int nnodes, SELEM_3D *elem3d, SVECT *elem_nds, double djac, double *u, double *v, double dt_factor, DOF_3 *elem_rhs, char *string, int DEBUG_LOCAL) {
    int i;
    double integral_X[nnodes];  sarray_init_dbl(integral_X,nnodes);
    double integral_Y[nnodes];  sarray_init_dbl(integral_Y,nnodes);
    if (nnodes == NDONTET) {
        integrate_tetrahedron_phi_f(djac, dt_factor, u, integral_X);
        integrate_tetrahedron_phi_f(djac, dt_factor, v, integral_Y);
    }
    else {
        integrate_triPrism_phi_f(elem_nds, dt_factor, u, integral_X);
        integrate_triPrism_phi_f(elem_nds, dt_factor, v, integral_Y);
    }
    for (i=0; i<nnodes; i++) {
        elem_rhs[i].x_eq += integral_X[i];
        elem_rhs[i].y_eq += integral_Y[i];
    }
#ifdef _DEBUG
    for (i=0; i<nnodes; i++) {
        elem_rhs_temporal[i].c_eq += 0.0;
        elem_rhs_temporal[i].x_eq += integral_X[i];
        elem_rhs_temporal[i].y_eq += integral_Y[i];
    }
    if (DEBUG_LOCAL == ON || dof3_debug_continue(elem_rhs_temporal, nnodes)) {
        rhs_3dof(string, nnodes, ie, elem3d->nodes, elem_rhs_temporal);
    }
    dof3_debug(elem_rhs_temporal, nnodes, __FILE__, __LINE__);
#endif
}


















