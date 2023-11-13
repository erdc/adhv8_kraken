/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// max local node allocated static veriables for debugging purposes

static int printFieldWidth = 13;
static int printPrecision  = 9;

static int DEBUG_LOCAL = OFF;
static int DEBUG_PICKETS = OFF;
static int DEBUG_NODE_ID = UNSET_INT;
static int DEBUG_EXIT_ON_NAN = ON;

static SELEM_2D *elem2d = NULL;
static int nnodes = UNSET_INT, nnodes_quad = UNSET_INT, string = UNSET_INT;

static double dt=0., gravity=0., elem_area=0., density = 0., elem_avg_density = 0.;
static int isElementTriangle = FALSE;
static int isSidewall = FALSE;
static int isSurface = FALSE;
static int isBed = FALSE;
static SNODE elem_nodes[NDONQUAD];

static SVECT2D elem_tanvec[NDONQUAD];
static SVECT elem_nds[NDONQUAD], elem_nds_fixed[NDONQUAD], elem_vel[NDONQUAD];
static double elem_normal_vel[NDONQUAD], elem_u[NDONQUAD], elem_v[NDONQUAD], elem_w[NDONQUAD];
static double elem_dpl[NDONQUAD], elem_prs[NDONQUAD], elem_rho[NDONQUAD];
static double nx = 0., ny = 0., nz = 0.;

static double elem_dpl_fixed[4] = {0., 0., 0., 0.}; // needs to be initial dpl, 0 should be ok as long as dpl is not hotstarted

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// file prototypes
int fe_ns3_boundary_test_for_NAN_or_INF();
void fe_ns3_boundary_debug(double perturbation, int perturb_node, int perturb_var, int perturb_sign);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes the boundary residual contributions for the 3D navier stokes model.
 *  \author    Charlie Berger, Ph.D.
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
 *  \details Integrates the weak, discrete outflow boundary terms: \n
 *  \f{eqnarray*}{
 *   R_{i,boundary}^e &=& dt * \bigg( \bcConv{}{e}{\phidd{i}}{(\velh_t \, \cjh)} \bigg)
 *                     =  dt * \bigg( \bcConv{}{e}{\phidd{i}}{(q_n^{hydro} \, \cjh)} \bigg)
 *  \f}
 *
 * \note CJT \:: Assumes that no displacements are hotstarted!
 * \note CJT \:: In the original AdH NS, when Charlie perturbed pressure, he pertubed displacement on a moving grid
 *               intead (for static grid, he perturbed pressure).  I found this odd and am not doing it here.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_ns3_boundary_resid(SMODEL *mod, DOF_4 *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    
    int i,print_debug_info = OFF,INCLUDE_PRESSURE = ON;
    
    // DO NOT TAMPER WITH +++++++++++++++++++++++++++
    // for time-dependent debugging
    if (mod->t_prev > DEBUG_TIME) DEBUG_LOCAL = ON;
    
    // for node-dependent debugging
    for (i=0; i<mod->grid->elem2d[ie].nnodes; i++) {
        if (mod->grid->elem2d[ie].nodes[i] == DEBUG_NODE_ID) {
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
    dt = mod->dt;
    gravity = mod->gravity;
    STR_VALUE *str_values = mod->str_values;
    SFLAGS flags = mod->flag;
    density = mod->density;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GRID VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    SGRID *grid = mod->grid;
    int **nds_on_edge = NULL;
    
    elem2d = &(grid->elem2d[ie]);
    string = elem2d->string;
    isElementTriangle = FALSE;
    nnodes = 0; nnodes_quad = 0.;
    if (elem2d->nnodes == NDONTRI) {
        isElementTriangle = TRUE;
        nnodes = NDONTRI;
        nnodes_quad = 2*NDONTRI; // for 2D triangle
        nds_on_edge = grid->nd_on_TriEdge;
    } else {
        nnodes = NDONQUAD;
        nnodes_quad = 2*NDONQUAD; // for 2D prism
        nds_on_edge = grid->nd_on_QuadEdge;
    }
    
    isSidewall = FALSE;
    isSurface = FALSE;
    isBed = FALSE;
    if      (elem2d->bflag == 2) {isSidewall = TRUE;}
    else if (elem2d->bflag == 0) {isSurface = TRUE;}
    else if (elem2d->bflag == 1) {isBed = TRUE;}
    else {
        printf("bflag = %d\n",elem2d->bflag);
        tl_error(">> 2D Boundary face bflag is not defined.");
    }
    
    int gnodeID_perturbed = elem2d->nodes[perturb_node];
    int hydro_bc_flag = str_values[string].flow.bc_flag;
    
    for(i=0; i<nnodes;i++) {
        snode_copy(&(elem_nodes[i]), grid->node[elem2d->nodes[i]]);
    }
    
    // boundary condition masks
    int mask[mod->nsys * nnodes], istart;
    for (i=0; i<nnodes; i++) {
        istart = mod->nsys * i;
        mask[istart]   = mod->bc_mask[elem2d->nodes[i] * mod->nsys];      // x_eq
        mask[istart+1] = mod->bc_mask[elem2d->nodes[i] * mod->nsys + 1];  // y_eq
        mask[istart+2] = mod->bc_mask[elem2d->nodes[i] * mod->nsys + 2];  // z_eq
        mask[istart+3] = mod->bc_mask[elem2d->nodes[i] * mod->nsys + 3];  // c_eq
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // pressure perturbation
    global_to_local_dbl(ns->prs, elem_prs, elem2d->nodes, nnodes);
    if (perturb_var == PERTURB_P) {elem_prs[perturb_node] += perturb_sign * perturbation;}
    
    // displacement perturbation (for moving grid nodes only)
    global_to_local_dbl(ns->displacement, elem_dpl, elem2d->nodes, nnodes);
    if (perturb_var == PERTURB_D) {elem_dpl[perturb_node] += perturb_sign * perturbation;}
    if (mod->flag.MG == OFF) {assert(fabs(elem_dpl[0]) < SMALL && fabs(elem_dpl[1]) < SMALL && fabs(elem_dpl[2]) < SMALL );}
    
    // velocity perturbation
    global_to_local_svect(ns->vel, elem_vel, elem2d->nodes, nnodes);
    //for (i=0; i<nnodes; i++) printf("elem_vel[%d]: %20.10f %20.10f %20.10f\n",i,elem_vel[i].x,elem_vel[i].y,elem_vel[i].z);
    
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
    
    // moving grid displacements only on surface and if turned on (add w perturbation)
    int gnodeID = UNSET_INT;
    if (mod->flag.MG == ON) {
        for (i=0; i<nnodes; i++) {
            gnodeID = elem2d->nodes[i];
            
            // add displacement perturbation from W
            if (grid->node[gnodeID].bflag == 0 && perturb_var == PERTURB_W) elem_dpl[i] += (perturb_sign * perturbation * dt);
            
            // adjust pressures in nodes below surface if W or D is perturbed
            if (grid->node[perturb_node].bflag == 0 && grid->node[gnodeID].bflag != 0) {
                if (perturb_var == PERTURB_D) { // perturb pressure below surface
                    if (fabs(elem_nodes[i].x - elem_nodes[perturb_node].x) < SMALL && fabs(elem_nodes[i].y - elem_nodes[perturb_node].y) < SMALL) {
                        elem_prs[i] += mod->density * mod->gravity * (perturb_sign * perturbation);
                    }
                }
                if (perturb_var == PERTURB_W) { // perturb pressure below surfae
                    if (fabs(elem_nodes[i].x - elem_nodes[perturb_node].x) < SMALL && fabs(elem_nodes[i].y - elem_nodes[perturb_node].y) < SMALL) {
                        elem_prs[i] += mod->density * mod->gravity * (perturb_sign * perturbation * dt);
                    }
                }
            }
        }
    }
    
    // node vertex and midpoint locations
    elem_get_midpt_locations2(elem_nodes, elem_dpl, elem_nds, nnodes, nnodes_quad, nds_on_edge);
    
    // density
    global_to_local_dbl(ns->density, elem_rho, elem2d->nodes, nnodes);
    
    // get phi gradients, etc.
    if (isElementTriangle == TRUE) {
        
        get_triangle_linear_djac_nrml_gradPhi(elem2d, NULL, elem_nds); // update normals and 2d/3d jacobians
        elem_area = elem2d->djac3d;
        elem_avg_density = one_3 * (elem_rho[0] + elem_rho[1] + elem_rho[2]);
        
    } else {
        
        elem_area = integrate_quadrilateral_area(elem_nds, 1.);
        elem_get_midpt_locations2(elem_nodes, elem_dpl_fixed, elem_nds_fixed, nnodes, nnodes_quad, nds_on_edge);
        elem2d->nrml = get_elem2d_normals(elem_nds);
        
    }
    nx = elem2d->nrml.x;
    ny = elem2d->nrml.y;
    nz = elem2d->nrml.z;
    
    // calculate normal velocities
    for (i=0; i<nnodes; i++) {
        elem_normal_vel[i] = svect_dotp(elem_vel[i], elem2d->nrml);
    }
    
    // for convienence
    for (i=0; i<nnodes; i++) {
        //printf("elem_vel[%d]: %20.10f %20.10f %20.10f\n",i,elem_vel[i].x,elem_vel[i].y,elem_vel[i].z);
        elem_u[i] = elem_vel[i].x;
        elem_v[i] = elem_vel[i].y;
        elem_w[i] = elem_vel[i].z;
    }
    
    // 2D element tangent vector
    global_to_local_svect2d(ns->tanvec, elem_tanvec, elem2d->nodes, elem2d->nnodes);
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    print_debug_info = OFF;
    
    assert(ie == elem2d->id);
    assert(nnodes == NDONQUAD || nnodes == NDONTRI);
    assert(perturb_var == PERTURB_NONE || perturb_var == PERTURB_P || perturb_var == PERTURB_U ||
           perturb_var == PERTURB_V || perturb_var == PERTURB_W || perturb_var == PERTURB_D);
    
    if (isElementTriangle == TRUE) {
        if (elem2d->djac3d < SMALL) {print_debug_info = ON; DEBUG_EXIT_ON_NAN = ON;}
        if (elem2d->djac3d_fixed < SMALL) {print_debug_info = ON; DEBUG_EXIT_ON_NAN = ON;}
    } else {
        // this order must be true for analytic quadrilateral integrations to be correct
        assert(fabs(elem_nds[2].x-elem_nds[1].x) < 1e-6 && fabs(elem_nds[2].y-elem_nds[1].y) < 1e-6);
        assert(fabs(elem_nds[3].x-elem_nds[0].x) < 1e-6 && fabs(elem_nds[3].y-elem_nds[0].y) < 1e-6);
    }
    if (elem_area < SMALL) {print_debug_info = ON; DEBUG_EXIT_ON_NAN = ON;}
    
    
    for (i=0; i<nnodes; i++) {
        if (mod->flag.MG == OFF) {
            if (fabs(elem_dpl[i]) > SMALL || fabs(elem_dpl_fixed[i]) > SMALL) {
                printf("ERROR :: elem_dpl[i] :: %20.10f \t elem_dpl_fixed[i] :: %20.10f\n",elem_dpl[i],elem_dpl_fixed[i]);
                tl_error("elem_dpl should be 0 for a fixed grid!\n");
            }
        }
    }
    
    // if there is a NAN or INF, print arrays and possibly exit
    if (print_debug_info == OFF) print_debug_info = fe_ns3_boundary_test_for_NAN_or_INF();
    if (print_debug_info == ON) {
        tl_check_all_pickets(__FILE__,__LINE__);
        fe_ns3_boundary_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
    
    if (DEBUG_LOCAL == ON) fe_ns3_boundary_debug(perturbation, perturb_node, perturb_var, perturb_sign);
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    dof4_init_array(elem_rhs, nnodes);
    
    double prs_integral[nnodes]; sarray_init_dbl(prs_integral, nnodes);
    double implicit_flow[nnodes]; sarray_init_dbl(implicit_flow, nnodes);
    double explicit_flow[nnodes]; sarray_init_dbl(explicit_flow, nnodes);
    double flux_x[nnodes]; sarray_init_dbl(flux_x, nnodes);
    double flux_y[nnodes]; sarray_init_dbl(flux_y, nnodes);
    double flux_z[nnodes]; sarray_init_dbl(flux_z, nnodes);
    double friction_integral_x[nnodes]; sarray_init_dbl(friction_integral_x, nnodes);
    double friction_integral_y[nnodes]; sarray_init_dbl(friction_integral_y, nnodes);
    double friction_integral_z[nnodes]; sarray_init_dbl(friction_integral_z, nnodes);
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   VELOCITY CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the convection boundary flux addition to the 2D elemental residual. \n
     * \note
     *********************************************************************************************/
    
    if (hydro_bc_flag == BCT_VEL_NEU || hydro_bc_flag == BCT_DIS_NEU) {
        
        // read normal velocity/discharge
        int isers = str_values[string].flow.isigma;
        double hydro_flux_normal = -sseries_get_value(isers, mod->series_head, 0); // (-) for flow into domain
        if (hydro_bc_flag == BCT_DIS_NEU) {
            hydro_flux_normal /= str_values[string].total_area;
        }
        assert(hydro_flux_normal < NOT_QUITE_SMALL); // use this BC only for inflow
        
        // -----------------------------------------------------------------------------------------
        // continuity addition ---------------------------------------------------------------------
        if (isElementTriangle) {
            integrate_triangle_phi_f(elem2d->djac3d, 1., elem_normal_vel, implicit_flow);
            integrate_triangle_phi(elem2d->djac3d_fixed,hydro_flux_normal, explicit_flow);
        } else {
            if (isSidewall == TRUE) {
                integrate_quadZ_phi_f(elem_nds, 1., elem_normal_vel, implicit_flow);
                integrate_quadZ_phi(elem_nds_fixed, hydro_flux_normal, explicit_flow);
            } else {
                integrate_quadrilateral_phi_f(elem_nds, 1., elem_normal_vel, implicit_flow);
                integrate_quadrilateral_phi(elem_nds_fixed, hydro_flux_normal, explicit_flow);
            }
        }
        for (i=0; i<nnodes; i++) {
            if (mask[i * mod->nsys] == YES) {
                // if there is DB on the momentum, then let mass flow implicitly
                elem_rhs[i].c_eq = implicit_flow[i]; // * (1./mod->density);
                //tag();
                //exit(-1);
            } else {
                elem_rhs[i].c_eq = explicit_flow[i]; // * (1./mod->density);
            }
        }
        
        // -----------------------------------------------------------------------------------------
        // momentum addition -----------------------------------------------------------------------
        double avg_u=0., avg_v=0., avg_w=0., avg_rho=0., avg_vel_mag = 0.; // elemental averages
        if (isElementTriangle) {
            avg_u = integrate_triangle_f(1., 1., elem_u); // elem_area's cancel
            avg_v = integrate_triangle_f(1., 1., elem_v); // elem_area's cancel
            avg_w = integrate_triangle_f(1., 1., elem_w); // elem_area's cancel
            avg_rho = integrate_triangle_f(1., 1., elem_rho); // elem_area's cancel
            avg_vel_mag = sqrt(avg_u*avg_u + avg_v*avg_v + avg_w*avg_w);
            
            if (INCLUDE_PRESSURE == ON) {
                integrate_triangle_phi_f(elem2d->djac3d, 1., elem_prs, prs_integral); // linear
                ////integrate_triangle_phi_fquad(elem2d->djac3d, 1., elem_prs, prs_integral); // quadratic
            }
            integrate_triangle_phi_f(elem2d->djac3d, avg_vel_mag, elem_u, friction_integral_x); // 1/2 rho f u |u|  :: friction term
            integrate_triangle_phi_f(elem2d->djac3d, avg_vel_mag, elem_v, friction_integral_y); // 1/2 rho f v |u|  :: friction term
            integrate_triangle_phi_f(elem2d->djac3d, avg_vel_mag, elem_w, friction_integral_z); // 1/2 rho f w |u|  :: friction term
        } else {
            if (isSidewall == TRUE) {
                avg_u = integrate_quadZ_f(elem_nds, 1./elem_area, elem_u);
                avg_v = integrate_quadZ_f(elem_nds, 1./elem_area, elem_v);
                avg_w = integrate_quadZ_f(elem_nds, 1./elem_area, elem_w);
                avg_rho = integrate_quadZ_f(elem_nds, 1./elem_area, elem_rho);
                avg_vel_mag = sqrt(avg_u*avg_u + avg_v*avg_v + avg_w*avg_w);
                
                integrate_quadZ_phi_f(elem_nds, 1., elem_prs, prs_integral); // linear
                integrate_quadZ_phi_f(elem_nds, avg_vel_mag, elem_u, friction_integral_x);
                integrate_quadZ_phi_f(elem_nds, avg_vel_mag, elem_v, friction_integral_y);
                integrate_quadZ_phi_f(elem_nds, avg_vel_mag, elem_w, friction_integral_z);
            } else {
                avg_u = integrate_quadrilateral_f(elem_nds, 1./elem_area, elem_u);
                avg_v = integrate_quadrilateral_f(elem_nds, 1./elem_area, elem_v);
                avg_w = integrate_quadrilateral_f(elem_nds, 1./elem_area, elem_w);
                avg_rho = integrate_quadrilateral_f(elem_nds, 1./elem_area, elem_rho);
                avg_vel_mag = sqrt(avg_u*avg_u + avg_v*avg_v + avg_w*avg_w);
                
                integrate_quadrilateral_phi_f(elem_nds, 1., elem_prs, prs_integral); // linear
                integrate_quadrilateral_phi_f(elem_nds, avg_vel_mag, elem_u, friction_integral_x);
                integrate_quadrilateral_phi_f(elem_nds, avg_vel_mag, elem_v, friction_integral_y);
                integrate_quadrilateral_phi_f(elem_nds, avg_vel_mag, elem_w, friction_integral_z);
            }
        }
        assert(avg_rho > 0);
        
        double drag = 0.0; //one_2 * avg_rho * mod->str_values[string].roughness;  // 1/2 rho f u |u|  :: friction term
        for (i=0; i<nnodes; i++) {
            if (ROTATE == OFF || mask[i * mod->nsys] != 2) {
                
                
                //if (elem_nodes[i].bflag == 1 && nz > 0) {nz = -nz;}
                //if (elem_nodes[i].id==5) printf("BOUNDARY PRESSURE :: integral_X: %20.10f integral_Y: %20.10f integral_Z: %20.10f\n",
                //                                prs_integral[i]*nx*(dt/elem_avg_density),
                //                                prs_integral[i]*ny*(dt/elem_avg_density),
                //                                prs_integral[i]*nz*(dt/elem_avg_density));
                
                elem_rhs[i].x_eq = prs_integral[i]*nx * (1./mod->density) + drag*friction_integral_x[i] * (1./mod->density) ;
                elem_rhs[i].y_eq = prs_integral[i]*ny * (1./mod->density) + drag*friction_integral_y[i] * (1./mod->density) ;
                elem_rhs[i].z_eq = prs_integral[i]*nz * (1./mod->density) + drag*friction_integral_z[i] * (1./mod->density) ;
            } else {
                //exit(-1);
                
                //elem_rhs[i].x_eq = implicit_flow[i] - explicit_flow[i];
                //elem_rhs[i].y_eq = prs_integral[i] * (nx * elem_tanvec[i].x + ny * elem_tanvec[i].y);
                //elem_rhs[i].y_eq += drag * (friction_integral_x[i] * elem_tanvec[i].x + friction_integral_y[i] * elem_tanvec[i].y);
                //elem_rhs[i].z_eq = prs_integral[i]*nz + drag*friction_integral_z[i];
                
                //elem_rhs[i].x_eq  = implicit_flow[i] - explicit_flow[i];
                //elem_rhs[i].y_eq  = (dt/elem_avg_density) * prs_integral[i] * (nx * elem_tanvec[i].x + ny * elem_tanvec[i].y);
                //elem_rhs[i].y_eq += (dt/elem_avg_density) * drag * (friction_integral_x[i] * elem_tanvec[i].x + friction_integral_y[i] * elem_tanvec[i].y);
                //elem_rhs[i].z_eq  = prs_integral[i]*nz * (dt/elem_avg_density)  + drag*friction_integral_z[i] * (dt/elem_avg_density) ;
                
                elem_rhs[i].x_eq  = (1./mod->density) * (implicit_flow[i] - explicit_flow[i]);
                elem_rhs[i].y_eq  = (1./mod->density) * prs_integral[i] * (nx * elem_tanvec[i].x + ny * elem_tanvec[i].y);
                elem_rhs[i].y_eq += (1./mod->density) * drag * (friction_integral_x[i] * elem_tanvec[i].x + friction_integral_y[i] * elem_tanvec[i].y);
                elem_rhs[i].z_eq  = prs_integral[i]*nz * (1./mod->density) + drag*friction_integral_z[i] * (1./mod->density) ;
            }
        }
        
        // -----------------------------------------------------------------------------------------
        // debug -----------------------------------------------------------------------------------
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_4dof("3D NS BCT_VEL_NEU/BCT_DIS_NEU: ",nnodes, ie, elem2d->nodes, elem_rhs, printFieldWidth, printPrecision);
        }
        if (dof4_debug_continue(elem_rhs,nnodes) == 1) {
            tl_check_all_pickets(__FILE__,__LINE__);
            rhs_4dof("3D NS BCT_VEL_NEU/BCT_DIS_NEU: ", nnodes, elem2d->id, elem2d->nodes, elem_rhs, printFieldWidth, printPrecision);
            printf("\n");
            fe_ns3_boundary_debug(perturbation, perturb_node, perturb_var, perturb_sign);
            if (DEBUG_EXIT_ON_NAN) exit(-1);
        }
#endif
        
        /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         *                                    OUTFLOW CONTRIBUTION
         *--------------------------------------------------------------------------------------------
         * Calculates the outflow addition to the elemental residual. \n
         * \note CJT \:: calculates implicit flow using implicitly calculated water elevation
         *********************************************************************************************/
        
    } else if (hydro_bc_flag == BCT_OUTFLOW) {
        
        // -----------------------------------------------------------------------------------------
        // continuity addition ---------------------------------------------------------------------
        
        if (isElementTriangle) {
            integrate_triangle_phi_f(elem2d->djac3d, 1., elem_normal_vel, implicit_flow);
        } else {
            if (isSidewall == TRUE) {
                integrate_quadZ_phi_f(elem_nds, 1., elem_normal_vel, implicit_flow);
            } else {
                integrate_quadrilateral_phi_f(elem_nds, 1., elem_normal_vel, implicit_flow);
            }
        }
        for (i=0; i<nnodes; i++) {elem_rhs[i].c_eq = implicit_flow[i];} // * (1./mod->density);}
        
        // -----------------------------------------------------------------------------------------
        // momentum addition -----------------------------------------------------------------------
        
        if (isElementTriangle) {
            //if (INCLUDE_PRESSURE == ON) {
            integrate_triangle_phi_f(elem2d->djac3d, 1., elem_prs, prs_integral); // linear
            //integrate_triangle_phi_fquad(elem2d->djac3d, dt, elem_prs, prs_integral); // quadratic
            
            //}
            
            integrate_triangle_phi_f_g(elem2d->djac3d, 1., elem_normal_vel, elem_u, flux_x); // comes from advection I.B.P.
            integrate_triangle_phi_f_g(elem2d->djac3d, 1., elem_normal_vel, elem_v, flux_y); // comes from advection I.B.P.
            integrate_triangle_phi_f_g(elem2d->djac3d, 1., elem_normal_vel, elem_w, flux_z); // comes from advection I.B.P.
        } else {
            if (isSidewall == TRUE) {
                integrate_quadZ_phi_f_g(elem_nds, 1., elem_normal_vel, elem_u, flux_x);
                integrate_quadZ_phi_f_g(elem_nds, 1., elem_normal_vel, elem_v, flux_y);
                integrate_quadZ_phi_f_g(elem_nds, 1., elem_normal_vel, elem_w, flux_z);
                integrate_quadZ_phi_f(elem_nds, 1., elem_prs, prs_integral); // linear
                
            } else {
                integrate_quadrilateral_phi_f_g(elem_nds, 1., elem_normal_vel, elem_u, flux_x);
                integrate_quadrilateral_phi_f_g(elem_nds, 1., elem_normal_vel, elem_v, flux_y);
                integrate_quadrilateral_phi_f_g(elem_nds, 1., elem_normal_vel, elem_w, flux_z);
                integrate_quadrilateral_phi_f(elem_nds, 1., elem_prs, prs_integral); // linear
                
            }
        }
        for (i=0; i<nnodes; i++) {
            if (ROTATE == OFF || mask[i * mod->nsys] != 2) {
                //elem_rhs[i].x_eq = elem_avg_density * flux_x[i] + prs_integral[i] * nx;
                //elem_rhs[i].y_eq = elem_avg_density * flux_y[i] + prs_integral[i] * ny;
                //elem_rhs[i].z_eq = elem_avg_density * flux_z[i] + prs_integral[i] * nz;
                
                elem_rhs[i].x_eq = flux_x[i] + prs_integral[i] * nx * (1./mod->density);
                elem_rhs[i].y_eq = flux_y[i] + prs_integral[i] * ny * (1./mod->density);
                elem_rhs[i].z_eq = flux_z[i] + prs_integral[i] * nz * (1./mod->density);
            } else {
                //elem_rhs[i].y_eq += elem_tanvec[i].x * flux_x[i] + elem_tanvec[i].y * flux_y[i];
                //elem_rhs[i].y_eq += prs_integral[i] * (nx * elem_tanvec[i].x + ny * elem_tanvec[i].y);
                //elem_rhs[i].z_eq += flux_z[i] + prs_integral[i] * nz;
                
                elem_rhs[i].y_eq  = elem_tanvec[i].x * flux_x[i] + elem_tanvec[i].y * flux_y[i];
                elem_rhs[i].y_eq += (1./mod->density) * prs_integral[i] * (nx * elem_tanvec[i].x + ny * elem_tanvec[i].y);
                elem_rhs[i].z_eq  = (1./mod->density) * prs_integral[i] * nz + flux_z[i] ;
                
                //exit(-1);
            }
        }
        
        // -----------------------------------------------------------------------------------------
        // debug -----------------------------------------------------------------------------------
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_4dof("3D NS BCT_OUTFLOW ",nnodes, ie, elem2d->nodes, elem_rhs, printFieldWidth, printPrecision);
        }
        if (dof4_debug_continue(elem_rhs,nnodes) == 1) {
            tl_check_all_pickets(__FILE__,__LINE__);
            rhs_4dof("3D NS OUTFLOW: ", nnodes, elem2d->id, elem2d->nodes, elem_rhs, printFieldWidth, printPrecision);
            printf("\n");
            fe_ns3_boundary_debug(perturbation, perturb_node, perturb_var, perturb_sign);
            if (DEBUG_EXIT_ON_NAN) exit(-1);
        }
#endif
    }
    
#ifdef _DEBUG
    time_t time2;  time(&time2);
    TIME_IN_NS3_BOUNDARY_RESID += difftime(time2,time1);
#endif
    
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Checks 3D NS boundary variables for NANS or INFS.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int fe_ns3_boundary_test_for_NAN_or_INF() {
    
    int print_debug_info = OFF;
    
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_u, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_u has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_v, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_v has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_w, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_w has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_dpl, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_dpl has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_nds, nnodes,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_nds has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_nds_fixed, nnodes,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_nds_fixed has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_prs, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_prs has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_normal_vel, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_normal_velocity has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_area,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_area has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_Double_Inf_or_NaN_noExit(nx,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("nx has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(ny,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("ny has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(nz,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("nz has either a NAN or INF, printing debug info \n\n");
    
    //print_debug_info = Is_Double_Inf_or_NaN_noExit(elem2d->djac,__FILE__ ,__LINE__);
    //if (print_debug_info != NO) printf("elem2d->djac has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem2d->djac3d,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem2d->djac3d has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem2d->djac3d_fixed,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem2d->djac3d_fixed has either a NAN or INF, printing debug info \n\n");
    
    
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem2d->djac,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem2d->djac has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem2d->djac,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem2d->djac has either a NAN or INF, printing debug info \n\n");
    
    if (print_debug_info == NO) print_debug_info = OFF;
    
    return print_debug_info;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Prints to screen DEBUG info for the 3D NS boundary residual.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_ns3_boundary_debug(double perturbation, int perturb_node, int perturb_var, int perturb_sign) {
    
    printf("3D TRANSPORT BOUNDARY ELEM RESID :: ie: %d \t dt: %*.*f",elem2d->id+1,printFieldWidth,printPrecision,dt);
    
    if      (isSidewall == TRUE) {printf(" SIDEWALL FACE ");}
    else if (isSurface  == TRUE) {printf(" SURFACE FACE ");}
    else if (isBed      == TRUE) {printf(" BED FACE ");}
    
    if (perturb_var == PERTURB_P) {
        printf("\t perturbing P || node: %d || perturbation: %-*.*e\n",elem2d->nodes[perturb_node],
               printFieldWidth,printPrecision,perturb_sign*perturbation);
    } else if (perturb_var == PERTURB_U) {
        printf("\t perturbing U || node: %d || perturbation: %-*.*e\n",elem2d->nodes[perturb_node],
               printFieldWidth,printPrecision,perturb_sign*perturbation);
    } else if (perturb_var == PERTURB_V) {
        printf("\t perturbing v || node: %d || perturbation: %-*.*e\n",elem2d->nodes[perturb_node],
               printFieldWidth,printPrecision,perturb_sign*perturbation);
    }else if (perturb_var == PERTURB_W) {
        printf("\t perturbing W || node: %d || perturbation: %-*.*e\n",elem2d->nodes[perturb_node],
               printFieldWidth,printPrecision,perturb_sign*perturbation);
    }else if (perturb_var == PERTURB_D) {
        printf("\t perturbing D || node: %d || perturbation: %-*.*e\n",elem2d->nodes[perturb_node],
               printFieldWidth,printPrecision,perturb_sign*perturbation);
    }else if (perturb_var == PERTURB_NONE) {
        printf("\t perturbing nothing\n");
    }
    
    
    if (isElementTriangle) {
        printf("djac2d: %-*.*f \t djac2d3d: %-*.*f \t djac2d3d_fixed: %-*.*f \n",
               printFieldWidth,printPrecision,elem2d->djac,
               printFieldWidth,printPrecision,elem2d->djac3d,
               printFieldWidth,printPrecision,elem2d->djac3d_fixed);
    }
    
    printf("elem_area: %-*.*f \t elem_normal: {%-*.*f, %-*.*f, %-*.*f}\n",
           printFieldWidth,printPrecision,elem_area,
           printFieldWidth,printPrecision,nx,
           printFieldWidth,printPrecision,ny,
           printFieldWidth,printPrecision,nz);
    
    printScreen_debug_vec("node locations: ",elem_nds, elem2d->nnodes);
    printScreen_debug_vec("node locations_fixed: ",elem_nds_fixed, elem2d->nnodes);
    printScreen_debug_svec2d("elem_tanvec", elem_tanvec, elem2d->nnodes, elem2d->nodes);
    printScreen_debug2_dbl("elem_prs", elem_prs, elem2d->nnodes, elem2d->nodes);
    printScreen_debug2_dbl("elem_dpl", elem_dpl, elem2d->nnodes, elem2d->nodes);
    printScreen_debug2_dbl("elem_normal_vel", elem_normal_vel, elem2d->nnodes, elem2d->nodes);
    printScreen_debug2_dbl("elem_u", elem_u, elem2d->nnodes, elem2d->nodes);
    printScreen_debug2_dbl("elem_v", elem_v, elem2d->nnodes, elem2d->nodes);
    printScreen_debug2_dbl("elem_w", elem_w, elem2d->nnodes, elem2d->nodes);
    printScreen_debug2_dbl("elem_rho", elem_rho, elem2d->nnodes, elem2d->nodes);
}
