/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_3d_transport_boundary_resid.c This file collections functions responsible for
 *          the 3D shallow water transport boundary contributions to the elemental residual.*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// max local node allocated static veriables for debugging purposes

static int DEBUG_LOCAL = OFF;
static int DEBUG_PICKETS = OFF;
static int DEBUG_NODE_ID = UNSET_INT;
static int DEBUG_EXIT_ON_NAN = ON;

static SELEM_2D *elem2d = NULL;
static int nnodes = UNSET_INT;

static double dt=0., gravity=0., c_inv=SMALL;
static int isElementTriangle = FALSE;
static int isSidewall = FALSE;
static int isSurface = FALSE;
static int isBed = FALSE;

static double elem_c[NDONQUAD], elem_c_old[NDONQUAD], elem_c_older[NDONQUAD], elem_dpl[NDONQUAD];
static SVECT elem_nds[NDONQUAD], elem_vel[NDONQUAD], elem_nds_fixed[NDONQUAD];
static double elem_source[NDONQUAD], elem_sink[NDONQUAD], elem_total_source[NDONQUAD];
static double elem_normal_vel[NDONQUAD];
static double elem_avg_nrml_vel = 0;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// file prototypes
int fe_3d_transport_boundary_test_for_NAN_or_INF();
void fe_3d_transport_boundary_debug(double perturbation, int perturb_node, int perturb_var, int perturb_sign);
double djac_quadz(double x1,double x2,double y1,double y2,double z1,double z2,double z3,double z4, double xhat) {
    return (one_8*sqrt(x1*x1 - 2*x1*x2 + x2*x2 + y1*y1 - 2*y1*y2 + y2*y2)*((xhat - 1)*z1 - (xhat + 1)*z2 + (xhat + 1)*z3 - (xhat - 1)*z4));
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes the boundary residual contributions for the 3D transport model.
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
 *  \details Integrates the weak, discrete outflow boundary terms: \n
 *  \f{eqnarray*}{
 *   R_{i,boundary}^e &=& dt * \bigg( \bcConv{}{e}{\phidd{i}}{(\velh_t \, \cjh)} \bigg)
 *                     =  dt * \bigg( \bcConv{}{e}{\phidd{i}}{(q_n^{hydro} \, \cjh)} \bigg)
 *  \f}
 *
 * \note CJT \:: Assumes that no displacements are hotstarted!
 * \note CJT \:: I think we need to rotate the boundaries here like SW
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_3d_transport_boundary_resid(SMODEL *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    
    int i;
    
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
    
    // determine which hydro is running
    double *dpl=NULL, *dpl_old=NULL, *dpl_older=NULL;
    SVECT *vel=NULL, *vel_old=NULL, *vel_older=NULL;
    if (mod->flag.SW_FLOW == ON) {
        dpl = mod->sw->d3->displacement;
        dpl_old = mod->sw->d3->old_displacement;
        dpl_older = mod->sw->d3->older_displacement;
        vel = mod->sw->d3->vel;
        vel_old = mod->sw->d3->old_vel;
        vel_older = mod->sw->d3->older_vel;
    } else if (mod->flag.NS_FLOW == ON) {
        dpl = mod->ns->d3->displacement;
        dpl_old = mod->ns->d3->old_displacement;
        dpl_older = mod->ns->d3->older_displacement;
        vel = mod->ns->d3->vel;
        vel_old = mod->ns->d3->old_vel;
        vel_older = mod->ns->d3->older_vel;
    }
    
    dt = mod->dt;
    gravity = mod->gravity;
    
    double *c=NULL, *old_c=NULL, *older_c=NULL, *source=NULL, *sink=NULL;
    int string = mod->grid->elem2d[ie].string;
    
    // determine which transport is running
    INPUT_DATA *str_value_hydro = &(mod->str_values[string].flow);
    INPUT_DATA *str_value_con = NULL;
    int trns_index = UNSET_INT;
    c_inv = SMALL;
    if (mod->is_sediment_running) {
#ifdef _SEDIMENT
        trns_index = mod->ised;
        c = mod->sed->susload[trns_index].c;
        old_c = mod->sed->susload[trns_index].old_c;
        older_c = mod->sed->susload[trns_index].older_c;
        source = mod->sed->susload[trns_index].source;
        sink = mod->sed->susload[trns_index].sink;
        c_inv = 1.E-6 / mod->sed->grain[trns_index].reference_c;
        str_value_con = &(mod->str_values[mod->grid->elem2d[ie].string].sed[trns_index]);
#endif
    } else {
        trns_index = mod->itrns;
        c = mod->con[trns_index].concentration;
        old_c = mod->con[trns_index].old_concentration;
        older_c = mod->con[trns_index].older_concentration;
        source = mod->con[trns_index].source;
        sink = mod->con[trns_index].sink;
        c_inv = 1. / mod->con[trns_index].property[0];
        str_value_con = &(mod->str_values[mod->grid->elem2d[ie].string].trans[trns_index]);
    }
    int hydro_bc_flag = str_value_hydro->bc_flag;
    int con_bc_flag   = str_value_con->bc_flag;
    if (mod->str_values[string].ol_flow.bc_flag == BCT_PRS_NEU) hydro_bc_flag = BCT_PRS_NEU;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GRID VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    SGRID *grid = mod->grid;
    elem2d = &(grid->elem2d[ie]);
    nnodes = elem2d->nnodes;
    
    isElementTriangle = FALSE;
    isSidewall = FALSE;
    isSurface = FALSE;
    isBed = FALSE;
    if (nnodes == NDONTRI) {isElementTriangle = TRUE;}
    else {isElementTriangle = FALSE;}
    
    if      (elem2d->bflag == 2) {isSidewall = TRUE;}
    else if (elem2d->bflag == 0) {isSurface = TRUE;}
    else if (elem2d->bflag == 1) {isBed = TRUE;}
    else {
        printf("bflag = %d\n",elem2d->bflag);
        tl_error(">> 2D Boundary face bflag is not defined.");
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    sarray_init_dbl(elem_c,NDONQUAD);
    global_to_local_dbl(c, elem_c, elem2d->nodes, nnodes);
    if (perturb_var == PERTURB_C) elem_c[perturb_node] += perturb_sign * perturbation;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // local elevation displacements
    sarray_init_dbl(elem_dpl,NDONQUAD);
    global_to_local_dbl(dpl, elem_dpl, elem2d->nodes, nnodes);
    
    // node positions
    svect_init_array(elem_nds,NDONQUAD);
    for (i=0; i<nnodes; i++) {
        elem_nds[i].x = grid->node[elem2d->nodes[i]].x;
        elem_nds[i].y = grid->node[elem2d->nodes[i]].y;
        elem_nds[i].z = grid->node[elem2d->nodes[i]].z + elem_dpl[i];
    }
    
    // current and initial element areas
    if (isElementTriangle == TRUE) {
        get_triangle_linear_djac_nrml_gradPhi(elem2d, NULL, elem_nds);
    }
    
    // note :: this assume that no displacements are hotstarted ...
    svect_init_array(elem_nds_fixed,NDONQUAD);
    for (i=0; i<nnodes; i++) {
            elem_nds_fixed[i].x = grid->node[elem2d->nodes[i]].x;
            elem_nds_fixed[i].y = grid->node[elem2d->nodes[i]].y;
            elem_nds_fixed[i].z = grid->node[elem2d->nodes[i]].z;
    }
    elem2d->nrml = get_elem2d_normals(elem_nds);
    
    // sources/sinks
    sarray_init_dbl(source,NDONQUAD);
    sarray_init_dbl(sink,NDONQUAD);
    sarray_init_dbl(elem_total_source,NDONQUAD);
    global_to_local_dbl(source, elem_source, elem2d->nodes, nnodes);
    global_to_local_dbl(sink, elem_sink, elem2d->nodes, nnodes);
    for (i=0; i<nnodes; i++) {
        elem_total_source[i] = elem_sink[i] + elem_source[i];
    }
    
    // water velocities
    svect_init_array(elem_vel,NDONQUAD);
    global_to_local_svect(vel, elem_vel, elem2d->nodes, nnodes);

    
    // if this is a clay sediment, add clay settling velocity to vz
#ifdef _SEDIMENT
    if (mod->sed->grain[trns_index].type == CLA) {
        for (i=0; i<nnodes; i++) {
            elem_vel[i].z -= mod->sed->grain[trns_index].clay.settling_velocity;
        }
    }
#endif
    
    // calculate normal velocities and elementally averaged normal velocity
    elem_avg_nrml_vel = 0.;
    sarray_init_dbl(elem_normal_vel,NDONQUAD);
    for (i = 0; i<nnodes; i++) {
        elem_normal_vel[i] = svect_dotp(elem_vel[i], elem2d->nrml);
    }
    elem_avg_nrml_vel = sarray_avg_dbl(elem_normal_vel, nnodes);
    
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    assert(ie == elem2d->id);
    assert(nnodes == NDONQUAD || nnodes == NDONTRI);
    assert(perturb_var == PERTURB_NONE || perturb_var == PERTURB_C);
    
    if (isElementTriangle == TRUE) {
        assert(elem2d->djac3d > 0);
        assert(elem2d->djac3d_fixed > 0);
    } else {
        // this order must be true for analytic quadrilateral integrations to be correct
        assert(fabs(elem_nds[2].x-elem_nds[1].x) < 1e-6 && fabs(elem_nds[2].y-elem_nds[1].y) < 1e-6);
        assert(fabs(elem_nds[3].x-elem_nds[0].x) < 1e-6 && fabs(elem_nds[3].y-elem_nds[0].y) < 1e-6);
    }
    
    // if there is a NAN or INF, print arrays and possibly exit
    int print_debug_info = fe_3d_transport_boundary_test_for_NAN_or_INF();
    if (print_debug_info != NO) {
        tl_check_all_pickets(__FILE__,__LINE__);
        fe_3d_transport_boundary_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
    
    if (DEBUG_LOCAL == ON) fe_3d_transport_boundary_debug(perturbation, perturb_node, perturb_var, perturb_sign);
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    
#endif
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    sarray_init_dbl(elem_rhs, nnodes);
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   VELOCITY CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the convection boundary flux addition to the 2D elemental residual. \n
     * \note  If the hydro flux is into the model (-), then the constituent flux must be assigned - it is set to 0 if not.
     * \note  If the hydro flux is out of the model, then the assigned constituent flux is ignored and calculated implicitly.
     * \note  There is an assumption here that the inflow or exit does not add or remove momentum.
     *********************************************************************************************/
    
    if (hydro_bc_flag == BCT_VEL_NEU || hydro_bc_flag == BCT_DIS_NEU) {
        
        int isers = str_value_hydro->isigma;
        double hydro_flux_normal = -1.0 * sseries_get_value(isers, mod->series_head, 0);
        
        if (hydro_bc_flag == BCT_DIS_NEU) {
            hydro_flux_normal /= mod->str_values[string].total_area;
        }
        
        if (hydro_flux_normal < -SMALL) {
            // get user input consitutent flux if given, otherwise set to 0
            if(con_bc_flag > NORMAL) {
                int isers = str_value_con->isigma;
                double user_c = sseries_get_value(isers, mod->series_head, 0) * c_inv;
                //printf("user_c: %20.10f \t c_inv: %20.10f\n",user_c,c_inv);
                if (isElementTriangle) {
                    integrate_triangle_phi(elem2d->djac3d_fixed, dt * hydro_flux_normal * user_c, elem_rhs);
                } else {
                    //integrate_quadZ_phi(elem_nds_fixed, dt * hydro_flux_normal * user_c, elem_rhs);
                    
                    // use quadrature
                    SQUAD *quad = grid->quad_rect;
                    int iqp, quad_order = 2; // quadrature order
                    SQUAD_PT *qp = NULL;
                    
                    for (iqp=0; iqp<quad[quad_order].n; iqp++) {
                        qp = &(quad[quad_order].pt[iqp]);
                        qp->djac = djac_quadz(elem_nds[0].x,elem_nds[1].x,elem_nds[0].y,elem_nds[1].y,elem_nds[0].z,elem_nds[1].z,elem_nds[2].z,elem_nds[3].z,qp->xhat);
                        double t1 = dt * hydro_flux_normal * qp->djac * qp->w * user_c;
                        for (i=0; i<nnodes; i++) {
                            elem_rhs[i] += t1 * qp->lshape[i];
                        }
                    }
                    
                }
            }
        } else {
            if (isElementTriangle) {
                integrate_triangle_phi_f(elem2d->djac3d_fixed, dt * hydro_flux_normal, elem_c, elem_rhs);
            } else {
                //integrate_quadZ_phi_f(elem_nds_fixed, dt * hydro_flux_normal, elem_c, elem_rhs);
                
                SQUAD *quad = grid->quad_rect;
                int iqp, quad_order = 3; // quadrature order
                SQUAD_PT *qp = NULL;
                
                // use quadrature
                for (iqp=0; iqp<quad[quad_order].n; iqp++) {
                    qp = &(quad[quad_order].pt[iqp]);
                    qp->djac = djac_quadz(elem_nds[0].x,elem_nds[1].x,elem_nds[0].y,elem_nds[1].y,elem_nds[0].z,elem_nds[1].z,elem_nds[2].z,elem_nds[3].z,qp->xhat);
                    double qp_c = SQUAD_get_function(qp, elem_c, nnodes);
                    double t1 = dt * hydro_flux_normal * qp->djac * qp->w * qp_c;
                    for (i=0; i<nnodes; i++) {
                        elem_rhs[i] += t1 * qp->lshape[i];
                    }
                }
            }
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("3D TRANSPORT BCT_VEL_NEU/BCT_DIS_NEU: ",nnodes, ie, elem2d->nodes, elem_rhs);
        }
        if (Is_DoubleArray_Inf_or_NaN_noExit(elem_rhs,nnodes,__FILE__ ,__LINE__) != NO) {
            tl_check_all_pickets(__FILE__,__LINE__);
            rhs_1dof("3D TRANSPORT BCT_VEL_NEU/BCT_DIS_NEU: ", nnodes, elem2d->id, elem2d->nodes, elem_rhs);
            printf("\n");
            fe_3d_transport_boundary_debug(perturbation, perturb_node, perturb_var, perturb_sign);
            if (DEBUG_EXIT_ON_NAN) exit(-1);
        }
#endif
    }
    
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    OUTFLOW CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the outflow addition to the elemental residual. \n
     * \note CJT \:: calculates implicit flow using implicitly calculated water elevation
     *********************************************************************************************/
    
    else if(hydro_bc_flag == BCT_OUTFLOW) {
        if (isElementTriangle) {
            integrate_triangle_phi_f_g(elem2d->djac3d_fixed, dt, elem_normal_vel, elem_c, elem_rhs);
        } else {
            integrate_quadZ_phi_f_g(elem_nds_fixed, dt, elem_normal_vel, elem_c, elem_rhs);
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("3D TRANSPORT BCT_OUTFLOW: ",nnodes, ie, elem2d->nodes, elem_rhs);
        }
        if (Is_DoubleArray_Inf_or_NaN_noExit(elem_rhs,nnodes,__FILE__ ,__LINE__) != NO) {
            tl_check_all_pickets(__FILE__,__LINE__);
            rhs_1dof("3D TRANSPORT BCT_OUTFLOW: ", nnodes, elem2d->id, elem2d->nodes, elem_rhs);
            printf("\n");
            fe_3d_transport_boundary_debug(perturbation, perturb_node, perturb_var, perturb_sign);
            if (DEBUG_EXIT_ON_NAN) exit(-1);
        }
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    FREE SURFACE SOURCING
     *--------------------------------------------------------------------------------------------
     * Calculates the bed source contribution to the residual. \n
     * \note  If the implicit hydro flux is into the model (-), then the constituent flux must be assigned - it is set to 0 if not.
     * \note  If the implicit hydro flux is out of the model (+), then the assigned constituent flux is ignored and calculated implicitly.
     *********************************************************************************************/
    
    else if(hydro_bc_flag == BCT_FRS) {
        // get water flux
        int isers = str_value_hydro->isigma;
        double hydro_flux_normal = sseries_get_value(isers, mod->series_head, 0);
        if (hydro_flux_normal > 0) {
            // water is leaving the domain through the surface (evaporation, etc)
        } else {
            // water is fluxed into the domain through the surface
            if (con_bc_flag > NORMAL) {
                int isers = str_value_con->isigma;
                double user_c = sseries_get_value(isers, mod->series_head, 0) * c_inv;
                if (isElementTriangle) {
                    integrate_triangle_phi(elem2d->djac3d_fixed, dt * hydro_flux_normal * user_c, elem_rhs);
                } else {
                    integrate_quadZ_phi(elem_nds_fixed, dt * hydro_flux_normal * user_c, elem_rhs);
                }
            }
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("3D TRANSPORT BCT_FRS: ",nnodes, ie, elem2d->nodes, elem_rhs);
        }
        if (Is_DoubleArray_Inf_or_NaN_noExit(elem_rhs,nnodes,__FILE__ ,__LINE__) != NO) {
            tl_check_all_pickets(__FILE__,__LINE__);
            rhs_1dof("3D TRANSPORT BCT_FRS: ", nnodes, elem2d->id, elem2d->nodes, elem_rhs);
            printf("\n");
            fe_3d_transport_boundary_debug(perturbation, perturb_node, perturb_var, perturb_sign);
            if (DEBUG_EXIT_ON_NAN) exit(-1);
        }
#endif
    }
    
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                        BED SOURCING
     *--------------------------------------------------------------------------------------------
     * Calculates the bed source contribution to the residual. \n
     * \note This comes from Sedlib.  Currently this uses lumping for triangles.
     *********************************************************************************************/
    
    else if (hydro_bc_flag == BCT_BED) {
        if (mod->is_sediment_running) {
            if (isElementTriangle) {
                integrate_triangle_phi_f_lump(elem2d->djac3d_fixed, dt, elem_total_source, elem_rhs);
            } else {
                integrate_quadZ_phi_f(elem_nds_fixed, dt, elem_total_source, elem_rhs);
            }
#ifdef _DEBUG
            if (DEBUG_LOCAL == ON) {
                rhs_1dof("3D TRANSPORT BCT_BED: ",nnodes, ie, elem2d->nodes, elem_rhs);
            }
            if (Is_DoubleArray_Inf_or_NaN_noExit(elem_rhs,nnodes,__FILE__ ,__LINE__) != NO) {
                tl_check_all_pickets(__FILE__,__LINE__);
                rhs_1dof("3D TRANSPORT BCT_BED: ", nnodes, elem2d->id, elem2d->nodes, elem_rhs);
                printf("\n");
                fe_3d_transport_boundary_debug(perturbation, perturb_node, perturb_var, perturb_sign);
                if (DEBUG_EXIT_ON_NAN) exit(-1);
            }
#endif
        }
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                             TAILWATER/PRESSURE CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the tailwater or pressure addition to the elemental residual. \n
     * \note  If the implicit hydro flux is into the model (-), then the constituent flux must be assigned - it is set to 0 if not.
     * \note  If the implicit hydro flux is out of the model, then the assigned constituent flux is ignored and calculated implicitly.
     *********************************************************************************************/
    
    else if (hydro_bc_flag == BCT_PRS_NEU) {
        if (elem_avg_nrml_vel > 0) { // flow out of the domain
            if (isElementTriangle) {
                integrate_triangle_phi_f_g(elem2d->djac3d, dt, elem_normal_vel, elem_c, elem_rhs);
            } else {
                //integrate_quadZ_phi_f_g(elem_nds, dt, elem_normal_vel, elem_c, elem_rhs);
                
                // use quadrature
                SQUAD *quad = grid->quad_rect;
                int iqp, quad_order = 4; // quadrature order
                SQUAD_PT *qp = NULL;
                
                for (iqp=0; iqp<quad[quad_order].n; iqp++) {
                    qp = &(quad[quad_order].pt[iqp]);
                    qp->djac = djac_quadz(elem_nds[0].x,elem_nds[1].x,elem_nds[0].y,elem_nds[1].y,elem_nds[0].z,elem_nds[1].z,elem_nds[2].z,elem_nds[3].z,qp->xhat);
                    double qp_c =       SQUAD_get_function(qp, elem_c, nnodes);
                    double qp_norm_v =       SQUAD_get_function(qp, elem_normal_vel, nnodes);
                    double t1 = dt * qp->djac * qp->w * qp_c * qp_norm_v;
                    for (i=0; i<nnodes; i++) {
                        elem_rhs[i] += t1 * qp->lshape[i];
                    }
                }
            }
        } else {                     // flow into the domain
            if (con_bc_flag > NORMAL) {
                int isers = str_value_con->isigma;
                double user_c = sseries_get_value(isers, mod->series_head, 0) * c_inv;
                if (isElementTriangle) {
                    integrate_triangle_phi_f(elem2d->djac3d_fixed, dt * user_c, elem_normal_vel, elem_rhs);
                } else {
                    //integrate_quadZ_phi_f(elem_nds_fixed, dt * user_c, elem_normal_vel, elem_rhs);
                    
                    // use quadrature
                    SQUAD *quad = grid->quad_rect;
                    int iqp, quad_order = 3; // quadrature order
                    SQUAD_PT *qp = NULL;
                    
                    for (iqp=0; iqp<quad[quad_order].n; iqp++) {
                        qp = &(quad[quad_order].pt[iqp]);
                        qp->djac = djac_quadz(elem_nds[0].x,elem_nds[1].x,elem_nds[0].y,elem_nds[1].y,elem_nds[0].z,elem_nds[1].z,elem_nds[2].z,elem_nds[3].z,qp->xhat);
                        double qp_norm_v =       SQUAD_get_function(qp, elem_normal_vel, nnodes);
                        double t1 = dt * qp->djac * qp->w * user_c * qp_norm_v;
                        for (i=0; i<nnodes; i++) {
                            elem_rhs[i] += t1 * qp->lshape[i];
                        }
                    }
                }
            }
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("3D TRANSPORT BCT_PRS_NEU: ",nnodes, ie, elem2d->nodes, elem_rhs);
        }
        if (Is_DoubleArray_Inf_or_NaN_noExit(elem_rhs,nnodes,__FILE__ ,__LINE__) != NO) {
            tl_check_all_pickets(__FILE__,__LINE__);
            rhs_1dof("3D TRANSPORT BCT_PRS_NEU: ", nnodes, elem2d->id, elem2d->nodes, elem_rhs);
            printf("\n");
            fe_3d_transport_boundary_debug(perturbation, perturb_node, perturb_var, perturb_sign);
            if (DEBUG_EXIT_ON_NAN) exit(-1);
        }
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                DIRICHLET BOUNDARY CONDITIONS
     *--------------------------------------------------------------------------------------------
     * Zeros matrix rows/columns when Dirichlet boundaries are used. \n
     * \note
     *********************************************************************************************/
    for (i=0; i<nnodes; i++) {
        int node_string = grid->node[ elem2d->nodes[i] ].string;
        if (node_string > NORMAL) {
            if (mod->is_sediment_running) {
                str_value_con = &(mod->str_values[node_string].sed[trns_index]);
            } else {
                str_value_con = &(mod->str_values[node_string].trans[trns_index]);
            }
            if (str_value_con->bc_flag == BCT_DIR || str_value_con->bc_flag == BCT_CEQ) {
                elem_rhs[i] = 0;
            }
        }
    }
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        rhs_1dof("3D TRANSPORT DIRICHLET BCS: ",nnodes, ie, elem2d->nodes, elem_rhs);
        Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    }
    if (Is_DoubleArray_Inf_or_NaN_noExit(elem_rhs,nnodes,__FILE__ ,__LINE__) != NO) {
        tl_check_all_pickets(__FILE__,__LINE__);
        rhs_1dof("3D TRANSPORT DIRICHLET BCS: ", nnodes, elem2d->id, elem2d->nodes, elem_rhs);
        printf("\n");
        fe_3d_transport_boundary_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
#endif
    
    
#ifdef _DEBUG
    time_t time2;  time(&time2);
    TIME_IN_3D_TRANSPORT_BOUNDARY_RESID += difftime(time2,time1);
#endif
    
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Checks 3D transport boundary variables for NANS or INFS.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int fe_3d_transport_boundary_test_for_NAN_or_INF() {
    
    int print_debug_info = NO;
    
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_c, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_c has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_c_old, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_c_old has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_c_older, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_c_older has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_dpl, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_dpl has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_vel, nnodes, __FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_vel has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_nds, nnodes,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_nds has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_nds_fixed, nnodes,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_nds_fixed has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_source, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_source has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_sink, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_sink has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_total_source, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_total_source has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_nrml_vel,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_nrml_vel has either a NAN or INF, printing debug info \n\n");
    
    return print_debug_info;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Prints to screen DEBUG info for the 3D transport boundary residual.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_3d_transport_boundary_debug(double perturbation, int perturb_node, int perturb_var, int perturb_sign) {
    
    printf("3D TRANSPORT BOUNDARY ELEM RESID :: ie: %d \t dt: %20.10f",elem2d->id,dt);
    if      (isSidewall == TRUE) {printf(" SIDEWALL FACE ");}
    else if (isSurface  == TRUE) {printf(" SURFACE FACE ");}
    else if (isBed      == TRUE) {printf(" BED FACE ");}
    if (perturb_var == PERTURB_C) {
        printf("\t perturbing C || node: %d || perturbation: %20.10e\n",elem2d->nodes[perturb_node],perturb_sign*perturbation);
    }
    printf("elem_avg_nrml_vel: %20.10e\n",elem_avg_nrml_vel);
    if (isElementTriangle) {
        printf("djac2d: %30.20f \t djac2d3d: %30.20f \t djac2d3d_fixed: %30.20f \n",elem2d->djac,elem2d->djac3d,elem2d->djac3d_fixed);
    }
    printf("nodes: %d %d %d \n",elem2d->nodes[0],elem2d->nodes[1],elem2d->nodes[2]);
    printScreen_debug_vec("node locations: ",elem_nds, elem2d->nnodes);
    printScreen_debug_vec("node locations_fixed: ",elem_nds_fixed, elem2d->nnodes);
    printScreen_debug2_dbl("elem_c", elem_c, elem2d->nnodes, elem2d->nodes);
    printScreen_debug2_dbl("elem_source", elem_source, elem2d->nnodes, elem2d->nodes);
    printScreen_debug2_dbl("elem_sink", elem_sink, elem2d->nnodes, elem2d->nodes);
    printScreen_debug2_dbl("elem_dpl", elem_dpl, elem2d->nnodes, elem2d->nodes);
    printScreen_debug2_dbl("elem_normal_vel", elem_normal_vel, elem2d->nnodes, elem2d->nodes);
    printScreen_debug_svect("elem_vel", elem_vel, elem2d->nnodes, elem2d->nodes);
}
