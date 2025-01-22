#include "adh.h"
// PROTOTYPES *******************************************************************************//
//*******************************************************************************************//
double fe_3m2d_dry_wet_factornew(SVECT2D * x, double s, double *h, double djac, int *num_dry);
#define SCALE_FOR_MASTER_TET_VOLUME (1.0/6.0) // tetrahedral djac to volume scale
#define SCALE_FOR_MASTER_TRI_VOLUME (1.0/2.0) // tetrahedral djac to volume scale
//*******************************************************************************************//
//*******************************************************************************************//

static double water_sources = 0.;
static double total_time_mass_flux = 0.;
static void global_to_local_svect(SVECT *global, SVECT *local, int *nodes, int nnodes) ;
static void global_to_local_dbl(double *global, double *local, int *nodes, int nnodes);
static void global_to_local_svect2d(SVECT2D *global, SVECT2D *local, int *nodes, int nnodes);
static double get_elem3d_volume(SVECT *node, int nnodes);
static SVECT get_elem2d_normals(SVECT *nd);
static void get_triangle_linear_djac_nrml_gradPhi(SELEM_2D *elem2d, SNODE *nd_SNODE, SVECT *nd_SVECT);
static double get_triprism_volume(SVECT *node);
//*******************************************************************************************//
//*******************************************************************************************//
//*******************************************************************************************//
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     stores local 3D vector values from the global arrays
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] the local/elemental values
 *  @param[in] global the global values
 *  @param[in] nodes an array if integer global node IDs
 *  @param[in] nnodes the number of nodes on the element
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void global_to_local_svect(SVECT *global, SVECT *local, int *nodes, int nnodes) {
    int i=0;
    for (i=0; i<nnodes; i++) {
        local[i] = global[nodes[i]];
    }
}
inline void global_to_local_dbl(double *global, double *local, int *nodes, int nnodes) {
    int i=0;
    for (i=0; i<nnodes; i++) {
        local[i] = global[nodes[i]];
    }
}
inline void global_to_local_svect2d(SVECT2D *global, SVECT2D *local, int *nodes, int nnodes) {
    int i=0;
    for (i=0; i<nnodes; i++) {
        local[i] = global[nodes[i]];
    }
}
inline double get_triprism_volume(SVECT *node) {
    double d1 = node[1].x - node[2].x, d2 = node[0].x - node[2].x, d3 = node[0].x - node[1].x;
    double prism_term = d1*node[0].y - d2*node[1].y + d3*node[2].y;
    double t1 = node[0].z + node[1].z + node[2].z;
    t1 -= (node[3].z + node[4].z + node[5].z);
    return (one_6 * prism_term * t1);
}
// returns a 3D element volume
inline double get_elem3d_volume(SVECT *node, int nnodes) {
    if (nnodes == NDONTET) {
        double d1 = node[1].x - node[2].x, d2 = node[0].x - node[2].x, d3 = node[0].x - node[1].x;
        double d4 = node[2].x - node[3].x, d5 = node[1].x - node[3].x, d6 = node[0].x - node[3].x;
        return (one_6 * ((d4*node[1].y - d5*node[2].y + d1*node[3].y)*node[0].z -
                         (d4*node[0].y - d6*node[2].y + d2*node[3].y)*node[1].z +
                         (d5*node[0].y - d6*node[1].y + d3*node[3].y)*node[2].z -
                         (d1*node[0].y - d2*node[1].y + d3*node[2].y)*node[3].z));
    } else if (nnodes == NDONPRISM) {
        return get_triprism_volume(node);
    } else {
        tl_error("bad call to elem2d_get_volume");
    }
    return (-1);
}

// returns a normal vector to the 2D element plane
inline SVECT get_elem2d_normals(SVECT *nd) {
    SVECT side1, side2;         /* two sides - the cross product is the normal */
    double magnitude = 0.;      /* the magnitude of the normal */
    SVECT nrml;
    svect_init(&nrml);
    svect_init(&side1);
    svect_init(&side2);
    
    /* computes two side vectors */
    side1.x = nd[1].x - nd[0].x;
    side1.y = nd[1].y - nd[0].y;
    side1.z = nd[1].z - nd[0].z;
    side2.x = nd[2].x - nd[0].x;
    side2.y = nd[2].y - nd[0].y;
    side2.z = nd[2].z - nd[0].z;
    
    /* computes their cross product */
    nrml.x = side1.y * side2.z - side1.z * side2.y;
    nrml.y = side1.z * side2.x - side1.x * side2.z;
    nrml.z = side1.x * side2.y - side1.y * side2.x;
    magnitude = svect_mag(nrml);
    
    /* normalizes the normal */
    svect_scale_replace(&nrml, 1./magnitude);
    
    return nrml;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// these are all constants on a triangle and can be calculated once and stored if grid points do not move
// note :: nd_SNODE is a global node vector
inline void get_triangle_linear_djac_nrml_gradPhi(SELEM_2D *elem2d, SNODE *nd_SNODE, SVECT *nd_SVECT) {
    
    SVECT side1, side2;     /* two sides - the cross product is the normal */
    double magnitude;       /* the magnitude of the normal */
    
#ifdef _DEBUG
    if (nd_SVECT == NULL) assert(nd_SNODE != NULL);
    if (nd_SNODE == NULL) assert(nd_SVECT != NULL);
#endif
    
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    if (nd_SNODE != NULL) { /* Gajanan gkc : Changing the to allow passing mod->grid->node as nd_SNODE in the function. */
        x1 = nd_SNODE[elem2d->nodes[0]].x; x2 = nd_SNODE[elem2d->nodes[1]].x; x3 = nd_SNODE[elem2d->nodes[2]].x;
        y1 = nd_SNODE[elem2d->nodes[0]].y; y2 = nd_SNODE[elem2d->nodes[1]].y; y3 = nd_SNODE[elem2d->nodes[2]].y;
        z1 = nd_SNODE[elem2d->nodes[0]].z; z2 = nd_SNODE[elem2d->nodes[1]].z; z3 = nd_SNODE[elem2d->nodes[2]].z;
    } else {
        x1 = nd_SVECT[0].x; x2 = nd_SVECT[1].x; x3 = nd_SVECT[2].x;
        y1 = nd_SVECT[0].y; y2 = nd_SVECT[1].y; y3 = nd_SVECT[2].y;
        z1 = nd_SVECT[0].z; z2 = nd_SVECT[1].z; z3 = nd_SVECT[2].z;
    }
    
    //printf("x1: %20.10f y1: %20.10f z1: %20.10f\n",x1,y1,z1);
    //printf("x2: %20.10f y2: %20.10f z2: %20.10f\n",x2,y2,z2);
    //printf("x3: %20.10f y3: %20.10f z3: %20.10f\n",x3,y3,z3);
    
    /* initialize the variables */
    svect_init(&elem2d->nrml);
    svect_init(&side1);
    svect_init(&side2);
    
    /* computes two side vectors */
    side1.x = x2 - x1;
    side1.y = y2 - y1;
    side1.z = z2 - z1;
    side2.x = x3 - x1;
    side2.y = y3 - y1;
    side2.z = z3 - z1;
    
    /* computes 2d element un-normalized normal and normal mag */
    elem2d->nrml.x = side1.y * side2.z - side1.z * side2.y;
    elem2d->nrml.y = side1.z * side2.x - side1.x * side2.z;
    elem2d->nrml.z = side1.x * side2.y - side1.y * side2.x;
    magnitude = svect_mag(elem2d->nrml);
    
    /* normalizes the normal */
    svect_scale_replace(&(elem2d->nrml), 1./magnitude);
    
    /* sets the jacobian */
    elem2d->djac3d = SCALE_FOR_MASTER_TRI_VOLUME * magnitude;
    assert(fabs(elem2d->djac3d) > SMALL);
    
    /* if 3d jacobian for original element has not been set yet, do it */
    if (fabs(elem2d->djac3d_fixed) < 1e-7) {
        elem2d->djac3d_fixed = elem2d->djac3d;
    }
    
    /* computes the 2d jacobian */
    elem2d->djac = (side1.x * side2.y) - (side2.x * side1.y);
    
    /* sets the 2d jacobian and the 2d shape functions */
    // in 3d, this jacobian is (-)!!!!! causes problem in bedload 

    if(fabs(elem2d->djac) > SMALL) {
        elem2d->grad_shp[0].x = (side1.y - side2.y);
        elem2d->grad_shp[0].y = (side2.x - side1.x);
        elem2d->grad_shp[1].x =  side2.y;
        elem2d->grad_shp[1].y = -side2.x;
        elem2d->grad_shp[2].x = -side1.y;
        elem2d->grad_shp[2].y =  side1.x;
        svect2d_scale_replace_array(elem2d->grad_shp, (1./elem2d->djac), elem2d->nnodes);
    } else {
        elem2d->grad_shp[0].x = UNSET_FLT;
        elem2d->grad_shp[0].y = UNSET_FLT;
        elem2d->grad_shp[1].x = UNSET_FLT;
        elem2d->grad_shp[1].y = UNSET_FLT;
        elem2d->grad_shp[2].x = UNSET_FLT;
        elem2d->grad_shp[2].y = UNSET_FLT;
    }
    
    /* scales the 2d jacobian to make it the area */
    elem2d->djac *= SCALE_FOR_MASTER_TRI_VOLUME;
}

// computes the water mass over a 2d element
// integral_0^1 integral_0^(1-x) (h0 (1-x-y)+h1 x+h2 y) dy dx = 1/6 (h0+h1+h2)
// cjt :: just use implicit flow of AdH for boundaries!!
double tl_find_grid_mass_elem2d(double density, STR_VALUE *str, SSERIES *series_head, double *depth, SGRID *grid, SFLAGS flags) {
    
    int i, j, ie = UNSET_INT, isers = UNSET_INT;
    double factor = 0.;
    SVECT2D coord[NDONQUAD], elem_vel[NDONQUAD];
    SVECT elem_nds[NDONQUAD];
    double elem_head[NDONQUAD];
    double elem_hneg[NDPRFC], factor1 = 0., factor2 = 0.;
    int num_dry, nelem_dry[4] = {0,0,0,0}, tot_elem_dry[4] = {0,0,0,0};
    
    double grid_mass = 0.;      /* the return variable */
    for (ie = 0; ie < grid->nelems2d; ie++) {
        //if(grid->elem2d[ie].my_pe > 1){
            
            for (i=0; i<grid->elem2d[ie].nnodes; i++) {
                elem_head[i] = depth[grid->elem2d[ie].nodes[i]];
                coord[i].x = grid->node[grid->elem2d[ie].nodes[i]].x;
                coord[i].y = grid->node[grid->elem2d[ie].nodes[i]].y;
                elem_nds[i].x = coord[i].x;
                elem_nds[i].y = coord[i].y;
                elem_nds[i].z = grid->node[grid->elem2d[ie].nodes[i]].z;
            }
            
            // ADD BODY MASS
            if (grid->elem2d[ie].nnodes == NDONTRI) {
                
                // CJT :: triangle mass collection includes this crazy stuff for wetting and drying
                // CJT :: instead of integrating and possible having negative mass, this actually adds the
                //        negative mass as a positive total
                for (i = 0; i < NDONTRI; i++) {
                    elem_hneg[i] = 0.0;
                }
                factor = fe_3m2d_dry_wet_factornew(coord, 0., elem_head, grid->elem2d[ie].djac, &num_dry);
                
                nelem_dry[num_dry]++;
                switch (num_dry) {
                    case 1:
                        for (i = 0; i < NDONTRI; i++) {
                            if (elem_head[i] < 0.0) {
                                elem_hneg[i] = elem_head[i];
                            }
                        }
                        factor1 = 1.0;
                        factor2 = -1.0 * (1.0 - factor);
                        break;
                    case 2:
                        for (i = 0; i < NDONTRI; i++) {
                            if (elem_head[i] < 0.0) {
                                elem_hneg[i] = elem_head[i];
                                elem_head[i] = 0.0;
                            }
                        }
                        factor1 = factor;
                        factor2 = 0.0;
                        break;
                    default:
                        factor1 = factor;
                        factor2 = 0.0;
                        break;
                }
                grid_mass += integrate_triangle_f(grid->elem2d[ie].djac, factor1, elem_head) +
                integrate_triangle_f(grid->elem2d[ie].djac, factor2, elem_hneg);
                
            } else {
                // quadrilateral
                grid_mass += integrate_quadrilateral_f(elem_nds, 1., elem_head);
            }
        //}
    }
#ifdef _MESSG
    grid_mass = messg_dsum(grid_mass, grid->smpi->ADH_COMM);
#endif
    
    return grid_mass * density;
    
}

//*******************************************************************************************//
//*******************************************************************************************//
//*******************************************************************************************//
// calculates the error from initial and current grid mass
double tl_find_grid_mass_error_elem2d(double density, double *depth, SVECT2D *vel, SGRID *grid, SFLAGS flags, double initial_grid_mass, SSERIES *series_head, STR_VALUE *str, double dt, double *total_time_mass_flux_T) {

    int ie = UNSET_INT, isers = UNSET_INT, i=0;
    double elem_water_source = 0., body_source = 0., factor = 1., total_mass_edge_flux = 0., elem_avg_depth = 0., constant = 0.;
    double elem_head[NDONQUAD]; sarray_init_dbl(elem_head,NDONQUAD);
    SVECT elem_nds[NDONQUAD]; svect_init_array(elem_nds,NDONQUAD);
    SVECT2D elem_vel[NDONQUAD]; svect2d_init_array(elem_vel,NDONQUAD);
    double normal_vel[NDONQUAD]; sarray_init_dbl(normal_vel,NDONQUAD);
    
    // get total grid mass
    double grid_mass = tl_find_grid_mass_elem2d(density, str, series_head, depth, grid, flags); // cjt :: already gathered
    // calculate total mass from body sources/sinks
    //for (ie = 0; ie < grid->nelems2d; ie++) {
//        if (str[grid->elem2d[ie].string].ol_flow.bc_flag == BCT_WATER_SOURCE) {
//            if (grid->elem2d[ie].nnodes == NDONTRI) {
//                for (i=0; i<grid->elem2d[ie].nnodes; i++) {
//                    elem_head[i] = depth[grid->elem2d[ie].nodes[i]];
//                    elem_nds[i].x = grid->node[grid->elem2d[ie].nodes[i]].x;
//                    elem_nds[i].y = grid->node[grid->elem2d[ie].nodes[i]].y;
//                    elem_nds[i].z = grid->node[grid->elem2d[ie].nodes[i]].z;
//                }
//                factor = fe_sw2_wet_dry_factor(elem_nds, elem_head, grid->elem2d[ie].djac);
//            }
//            isers = str[grid->elem2d[ie].string].ol_flow.isigma;
//            elem_water_source = sseries_get_value(isers, series_head, 0);
//            body_source += elem_water_source * grid->elem2d[ie].djac * dt * density * factor; //(1./2.) * 2 * dt;
//            i++;
//        }
//    }
    
    // cjt :: implicit flow of mass through external boundaries
    for (ie = 0; ie < grid->nelems1d; ie++) {
        global_to_local_svect2d(vel, elem_vel, grid->elem1d[ie].nodes, grid->elem1d[ie].nnodes);
        for (i = 0; i<grid->elem1d[ie].nnodes; i++) {
            elem_head[i] = depth[grid->elem1d[ie].nodes[i]];
            normal_vel[i] = svect2d_dotp(elem_vel[i], grid->elem1d[ie].nrml);
        }
        
        // cjt :: (-) normal flux is into the domain here, so this turns into a positive mass flux
        constant = density * dt * grid->elem1d->djac;
        total_mass_edge_flux -= integrate_line_h_f(grid->elem1d->djac, density * dt, normal_vel, elem_head);
    }
    
    double total_mass_flux = 0.;
#ifdef _MESSG
    total_mass_flux = messg_dsum(body_source + total_mass_edge_flux, grid->smpi->ADH_COMM);
#else
    total_mass_flux  = body_source + total_mass_edge_flux;
#endif
    *total_time_mass_flux_T += total_mass_flux; // add this time to previous time-step calls
    double grid_mass_error = fabs(grid_mass - *total_time_mass_flux_T - initial_grid_mass);
#ifdef _MESSG
    if(grid->smpi->myid==0) printf("\ninitial mass: %20.10e  mass flux out/into domain: %20.10e current_mass: %20.10e mass_error: %20.10e  mass_rel_error: %20.10e\n",initial_grid_mass, *total_time_mass_flux_T, grid_mass, grid_mass_error,100*(grid_mass_error)/(initial_grid_mass));
#endif   
    return grid_mass_error;
}

//*******************************************************************************************//
//*******************************************************************************************//
//*******************************************************************************************//
// computes the water mass over a 3d tet element
double tl_find_grid_mass_elem3d(double density, SGRID *grid, double *displacement) {
    
    int ie = UNSET_INT, i, j;
    double grid_mass = 0.;
    
    SVECT elem_nds[NDONPRISM];
    for (ie = 0; ie < grid->nelems3d; ie++) {
        //if(grid->elem3d[ie].my_pe > 1){ // has to be a resident element
            for (j=0; j<grid->elem3d[ie].nnodes; j++) {
                elem_nds[j].x = grid->node[grid->elem3d[ie].nodes[j]].x;
                elem_nds[j].y = grid->node[grid->elem3d[ie].nodes[j]].y;
                elem_nds[j].z = grid->node[grid->elem3d[ie].nodes[j]].z + displacement[grid->elem3d[ie].nodes[j]];
            }
            grid_mass += density * get_elem3d_volume(elem_nds, grid->elem3d[ie].nnodes);
        //}
    }
#ifdef _MESSG
    grid_mass = messg_dsum(grid_mass, grid->smpi->ADH_COMM);
#endif
    return grid_mass;
}

//*******************************************************************************************//
//*******************************************************************************************//
//*******************************************************************************************//
// calculates the error from initial and current grid mass
// Int[ rho * dV(t + dt)] (current grid volume) + Sum_{ndt}(Int[rho * (v(t+dt) dot n) dA] * dt)  (total boundary flux over time) - Int[ rho * dV(t=0)] (initial grid volume)

// note :: not for barolinic simulations yet
// note :: I tried to use dpl(t+3/2) and dpl(t+1/2), but results were much worse for some reason
double tl_find_3d_grid_mass_error(STR_VALUE *str, SSERIES *series_head, double initial_grid_mass, double density, SGRID *grid, SVECT *vel, double *displacement, double *old_displacement, double *older_displacement, double dt, double *new_grid_mass, double *total_time_mass_flux_T) {
    
    int i, j, string = UNSET_INT, ie = UNSET_INT, isers = UNSET_INT;
    double flux = 0., area2d3d = 0.;
    SVECT elem_nds[NDONPRISM];
    SELEM_2D *elem2d;
    
    // get total grid mass
    *new_grid_mass = tl_find_grid_mass_elem3d(density, grid, displacement);
    
    // get mass leaving/entering domain through boundaries
    double total_mass_flux = 0.;
    double elem_dpl[NDONPRISM], elem_old_dpl[NDONPRISM];
    SVECT elem_vel[NDONPRISM], elem_rel_vel[NDONPRISM];
    for (ie = 0; ie < grid->nelems2d; ie++) {
        
        // the kinematic SW bc and ALE grid on surface and bed means no flow there
        //if (grid->elem2d[ie].bflag == 0 || grid->elem2d[ie].bflag == 1) continue;
        
        //if(grid->elem2d[ie].my_pe > 1) { // must be a resident element
            
            elem2d = &grid->elem2d[ie];
            string = elem2d->string;
            
            global_to_local_svect(vel, elem_vel, elem2d->nodes, elem2d->nnodes);
            global_to_local_dbl(displacement, elem_dpl, elem2d->nodes, elem2d->nnodes);
            global_to_local_dbl(old_displacement, elem_old_dpl, elem2d->nodes, elem2d->nnodes);
            
            for (j=0; j<elem2d->nnodes; j++) {
                elem_nds[j].x = grid->node[elem2d->nodes[j]].x;
                elem_nds[j].y = grid->node[elem2d->nodes[j]].y;
                elem_nds[j].z = grid->node[elem2d->nodes[j]].z + elem_dpl[j];
            }
            
            if (elem2d->nnodes == NDONTRI) {
                get_triangle_linear_djac_nrml_gradPhi(elem2d, NULL, elem_nds);
            } else {
                elem2d->nrml = get_elem2d_normals(elem_nds);
            }
            
            double normal_rel_vel[elem2d->nnodes];
            for (i = 0; i<elem2d->nnodes; i++) {
                elem_rel_vel[i].x = elem_vel[i].x;
                elem_rel_vel[i].y = elem_vel[i].y;
                elem_rel_vel[i].z = elem_vel[i].z - (elem_dpl[i] - elem_old_dpl[i])/dt;
                normal_rel_vel[i] = svect_dotp(elem_rel_vel[i], elem2d->nrml);
            }
            
            // cjt :: (-) normal flux is into the domain here, so this turns into a positive mass flux
            if (elem2d->nnodes == NDONTRI) {
                total_mass_flux += -integrate_triangle_f(elem2d->djac3d, density * dt, normal_rel_vel);
            } else {
                total_mass_flux += -integrate_quadZ_f(elem_nds, density * dt, normal_rel_vel);
            }
        //}
    }
    
#ifdef _MESSG
    total_mass_flux = messg_dsum(total_mass_flux, grid->smpi->ADH_COMM);
#endif
    
    *total_time_mass_flux_T += total_mass_flux; // add this time to previous time-step calls
    double grid_mass_error = (*new_grid_mass - *total_time_mass_flux_T) - initial_grid_mass;
    
    return grid_mass_error;
}

//*******************************************************************************************//
//*******************************************************************************************//
//*******************************************************************************************//

double fe_3m2d_dry_wet_factornew(SVECT2D *x, double s, double *h, double djac, int *num_dry) {
    
    int i, is_dry, total_vtx, passcode, only_vtx;
    double elem_djac, factor, diff_h[3];
    
    /*here's where the magic happens - passcode is a number that will be used
     in both its binary and decimal representations to flag which vertices are dry
     total_vtx is the number of dry vertices (either 0, 1, 2, or 3) */
    
    passcode = 0;
    total_vtx = 0;
    
    /*checking to see if the vertices are dry or not
     if a vertex is dry, increment total_vtx and set the corresponding
     bit of the binary representation of passcode to '1' by using the binary
     shift operator */
    
    for(i = 0; i < 3; i++)
    {
        is_dry = (h[i] < 0.0) ? 1 : 0;
        total_vtx += is_dry;
        passcode += (is_dry << i);
    }
    
    /*if total_vtx = 3, we're all dry and no integration needs to take place
     if total vtx = 0, we're all wet and the entire element may be integrated as normal
     otherwise (total_vtx=1 or 2), we can prepare to create a new element by performing
     mathematical magic using the decimal representation of passcode to find the 'unique' element
     i.e., that which is the only one of its kind (dry or wet) in the original element */
    
    if(total_vtx == 3) {
        *num_dry = 3;
        return 0.0;
    } else if(total_vtx == 0) {
        *num_dry = 0;
        return 1.0;
    }
    
    if(total_vtx != 1)
        passcode = (7 - passcode);
    
    only_vtx = passcode / 2;
    
    {
        double xsi, new_h[3];
        SVECT2D new_x[3];
        /*loop through the vertices and find only_vtx then calculate position of the 2 new vertecies */
        for(i = 0; i < 3; i++)
        {
            if(i != only_vtx)
            {
                xsi = h[i] / (h[i] - h[only_vtx]);
                FIXED_POS(x[i].x, x[only_vtx].x, xsi, new_x[i].x);
                FIXED_POS(x[i].y, x[only_vtx].y, xsi, new_x[i].y);
                new_h[i] = 0;
            }
            else
            {
                new_x[i] = x[i];
                new_h[i] = h[i];
            }
        }
        
        elem_djac = djac;
        TRI_AREA(new_x[0], new_x[1], new_x[2], djac);
        /* it is possible to have a tiny portion of the element that is dry */
        /* if it is extremely small we can have problems when the gradients */
        /* are calculated.  In reality when the dry portion is really small */
        /* then the element should be considered all wet and so we return   */
        /* immediately with factor = 1.0                                    */
        if(djac / elem_djac < SMALL) {
            *num_dry = (i=total_vtx%2)==1 ? 0 : 3;
            return (double)(i);
        }
        factor = (total_vtx >> 1) - (total_vtx % 2);
    }
    *num_dry = total_vtx;
    factor *= djac / elem_djac;
    return factor + ((factor < 0) ? 1.0 : 0.0);
}
