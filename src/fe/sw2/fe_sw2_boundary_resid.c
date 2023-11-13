
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_sw2_1d_elem_resid.c This file collections functions responsible for
 *          the 2D shallow water boundary contributions to the elemental residual.          */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static int printFieldWidth = 20;
static int printPrecision  = 10;

static int DEBUG_LOCAL = OFF;
static int DEBUG_PICKETS = OFF;
static int DEBUG_NODE_ID = UNSET_INT;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_sw2_1d_explicit_flow    (DOF_3 *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, double flux, double h_avg, double *u, double *v, int DEBUG, int DEBUG_LOCAL);
void fe_sw2_1d_implicit_flow    (DOF_3 *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, SVECT2D *v, double *h, int DEBUG, int DEBUG_LOCAL);
void fe_sw2_1d_wall_friction    (DOF_3 *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, double *u, double *v, double *resistance, int DEBUG, int DEBUG_LOCAL);
void fe_sw2_1d_pressure         (DOF_3 *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, double *h, double g, int DEBUG, int DEBUG_LOCAL);
void fe_sw2_1d_density_pressure (DOF_3 *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, double *h, double *d, double g, int DEBUG, int DEBUG_LOCAL);
void fe_sw2_1d_hvel             (DOF_3 *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, SVECT2D *v, SVECT2D user_velocity, double *h, int DEBUG, int DEBUG_LOCAL);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes the 1D element residual for the SW-2D model.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  \note CJT \:: Uses the SW2D model struct
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
 *   \residDA{i}{}{c}   &=& dt * \intbe{\phidd{i} \, q_n}{e} = dt * \intbe{\phidd{i} \, (\velb{h} \cdot \nrml) \, \depth{h}}{e} \\
 *   \residDA{i}{}{mx}  &=& dt * \bigg[
 *                                \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \ub{h} \, \depth{h})}
 *                             +  \bcFriction{}{e}{\phidd{i}}{\ub{h}}{\velb{h}}
 *                             +  \bcPressureDD{}{e}{\phidd{i}}{(\depth{h})}{x}
 *                             +  \bcPressureDD{}{e}{\phidd{i} \,\rho^h}{(\depth{h})}{x}
 *                               \bigg] \\
 *   \residDA{i}{}{my}  &=& dt * \bigg[
 *                                \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \ub{h} \, \depth{h})}
 *                             +  \bcFriction{}{e}{\phidd{i}}{\vb{h}}{\velb{h}}
 *                             +  \bcPressureDD{}{e}{\phidd{i}}{(\depth{h})}{y}
 *                             +  \bcPressureDD{}{e}{\phidd{i} \,\rho^h}{(\depth{h})}{y}
 *                               \bigg]
 *  \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_boundary_resid(SMODEL *mod, DOF_3 *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    
    int i;
    
#ifdef _DEBUG
    if (DEBUG_NODE_ID != UNSET_INT) {
        for (i=0; i<mod->grid->elem1d[ie].nnodes; i++) {
            if (mod->grid->elem1d[ie].nodes[i] == DEBUG_NODE_ID - 1)  {DEBUG = ON; break;} else {DEBUG = OFF; DEBUG_LOCAL = OFF; DEBUG_PICKETS = OFF;}
        }
    }
    if (DEBUG == ON) DEBUG_LOCAL = ON;
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    time_t time1;  time(&time1);
#endif
    
    
    // aliases
    SSW_2D *sw2 = mod->sw->d2;
    SELEM_1D elem1d = mod->grid->elem1d[ie];
    int string = elem1d.string;
    int nnodes = elem1d.nnodes;
    SVECT2D nrml = elem1d.nrml;
    double dt = mod->dt;
    double djac = elem1d.djac;
    double g = mod->gravity;
    STR_VALUE *str_values = mod->str_values;
    SWEIR_C *weir = mod->weir;
    SFLAP_C *flap = mod->flap;
    SSLUICE_C *sluice = mod->sluice;
    
    int PRESSURE_FLAG = ON;
    
#ifdef _DEBUG
    assert(nnodes == NDONSEG); // only 1D element is a line segment
    assert(djac > 0); // only 1D element is a line segment
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    if (perturb_var == PERTURB_NONE) perturb_sign = 0;
    
    double elem_head[nnodes];
    global_to_local_dbl(sw2->head, elem_head, elem1d.nodes, nnodes);
    if (perturb_var == PERTURB_H) elem_head[perturb_node] += perturb_sign * perturbation;
    
    
    SVECT2D elem_vel[nnodes];
    global_to_local_svect2d(sw2->vel, elem_vel, elem1d.nodes, nnodes);
    if (perturb_var == PERTURB_U) {
        elem_vel[perturb_node].x += perturb_sign * perturbation;
        PRESSURE_FLAG = OFF;
    }
    if (perturb_var == PERTURB_V) {
        elem_vel[perturb_node].y += perturb_sign * perturbation;
        PRESSURE_FLAG = OFF;
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    SVECT2D elem_old_vel[nnodes], elem_older_vel[nnodes];
    global_to_local_svect2d(sw2->old_vel, elem_old_vel, elem1d.nodes, nnodes);
    global_to_local_svect2d(sw2->older_vel, elem_older_vel, elem1d.nodes, nnodes);
    double u[NDONSEG], v[NDONSEG];
    for (i=0; i<nnodes; i++) { u[i] = elem_vel[i].x; v[i] = elem_vel[i].y; }
    
    double elem_old_head[nnodes], elem_older_head[nnodes];
    global_to_local_dbl(sw2->old_head, elem_old_head, elem1d.nodes, nnodes);
    global_to_local_dbl(sw2->older_head, elem_older_head, elem1d.nodes, nnodes);
    
    double elem_density[nnodes]; sarray_init_dbl(elem_density,nnodes);
    if (mod->flag.BAROCLINIC == 1) {
        global_to_local_dbl(sw2->density, elem_density, elem1d.nodes, nnodes);
        for (i = 0; i < nnodes; i++) {
            if (elem_head[i] < 0) {
                elem_density[i] = 0;
            }
        }
    }
    double u_tan = one_2 * (elem_vel[0].x + elem_vel[1].x) * nrml.y - one_2 * (elem_vel[0].y + elem_vel[1].y) * nrml.x;
    double mag_utan = fabs(u_tan);
    double h_avg = one_2 * (elem_head[0] + elem_head[1]); if (h_avg < SMALL) h_avg = SMALL;
    double u_avg = one_2 * (u[0] + u[1]);
    double v_avg = one_2 * (v[0] + v[1]);
    
    double roughness[nnodes];
    roughness[0] = fe_sw2_get_roughness(mod, elem_head[0], 0.0, string, UNUSED);
    roughness[1] = fe_sw2_get_roughness(mod, elem_head[1], 0.0, string, UNUSED);
    
    double resistance[nnodes];
    resistance[0] = one_2 * roughness[0] * mag_utan;
    resistance[1] = one_2 * roughness[1] * mag_utan;
    
    double elem_node_z[nnodes];
    elem_node_z[0] = mod->grid->node[elem1d.nodes[0]].z + sw2->bed_displacement[elem1d.nodes[0]];
    elem_node_z[1] = mod->grid->node[elem1d.nodes[1]].z + sw2->bed_displacement[elem1d.nodes[1]];
    
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        printf("SW-2D 1D ELEM RESID :: ie: %d \t dt: %20.10f \t djac: %20.10f nrml: {%10.5e, %10.5e}\n",ie,dt,djac,nrml.x,nrml.y);
        if (perturb_var == PERTURB_H) {
            printf("\t perturbing H || node: %d || perturbation: %20.10e\n",elem1d.nodes[perturb_node],perturb_sign*perturbation);
        } else if (perturb_var == PERTURB_U) {
            printf("\t perturbing U || node: %d || perturbation: %20.10e\n",elem1d.nodes[perturb_node],perturb_sign*perturbation);
        } else if (perturb_var == PERTURB_V) {
            printf("\t perturbing V || node: %d || perturbation: %20.10e\n",elem1d.nodes[perturb_node],perturb_sign*perturbation);
        }
        printf("h_avg: %20.10e \t u_avg: %20.10e \t v_avg: %20.10e \t u_tan: %20.10e\n",h_avg, u_avg, v_avg, u_tan);
        printf("\n--------------------------------------------------------- \n");
        printScreen_debug2_dbl("roughness", roughness, nnodes, elem1d.nodes);
        printScreen_debug2_dbl("elem_density", elem_density, nnodes, elem1d.nodes);
        printScreen_debug2_dbl("elem_head", elem_head, nnodes, elem1d.nodes);
        printScreen_debug2_dbl("elem_old_head", elem_old_head, nnodes, elem1d.nodes);
        printScreen_debug2_dbl("elem_older_head", elem_older_head, nnodes, elem1d.nodes);
        printScreen_debug_svec2d("elem_vel", elem_vel, nnodes, elem1d.nodes);
        printScreen_debug_svec2d("elem_old_vel", elem_old_vel, nnodes, elem1d.nodes);
        printScreen_debug_svec2d("elem_older_vel", elem_older_vel, nnodes, elem1d.nodes);
    }
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    dof3_init_array(elem_rhs, nnodes);
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   VELOCITY CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the convection boundary flux addition to the SW 2D elemental residual. \n
     * \note
     *********************************************************************************************/
    if (str_values[string].ol_flow.bc_flag == BCT_VEL_NEU) {
        int isers = str_values[string].ol_flow.isigma;
        double flux = -1.0 * sseries_get_value(isers, mod->series_head, 0);
        fe_sw2_1d_explicit_flow(elem_rhs, elem1d, ie, nrml, djac, dt, flux, h_avg, u, v, DEBUG, DEBUG_LOCAL);
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_VEL_NEU :: flux: %20.10e \n",flux);
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    OUTFLOW CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the outflow addition to the elemental residual. \n
     * \note CJT \:: calculates implicit flow using implicitly calculated water elevation
     *********************************************************************************************/
    else if (str_values[string].ol_flow.bc_flag == BCT_OUTFLOW) {
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_OUTFLOW \n");
#endif
        fe_sw2_1d_implicit_flow(elem_rhs, elem1d, ie, nrml, djac, dt, elem_vel, elem_head, DEBUG, DEBUG_LOCAL);
    }
    
#ifdef _ADH_STRUCTURES
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   WEIR CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates a weir addition to the elemental residual. \n
     *
     *********************************************************************************************/
    else if(str_values[string].ol_flow.bc_flag == BCT_WEIRU) {
        
        int iwn = str_values[string].weir_num;
        double flux = weir[iwn].fluxu;
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_VEL_NEU :: flux: %20.10e \n",flux);
#endif
        fe_sw2_1d_explicit_flow(elem_rhs, elem1d, ie, nrml, djac, dt, flux, h_avg, u, v, DEBUG, DEBUG_LOCAL);
        fe_sw2_1d_wall_friction(elem_rhs, elem1d, ie, nrml, djac, dt, u, v, resistance, DEBUG, DEBUG_LOCAL);
    }
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   WEIR CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates a wier addition to the elemental residual. \n
     *
     *********************************************************************************************/
    else if(str_values[string].ol_flow.bc_flag == BCT_WEIRD) {
        int iwn = str_values[string].weir_num;
        double flux = weir[iwn].fluxd;
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_VEL_NEU :: flux: %20.10e \n",flux);
#endif
        fe_sw2_1d_explicit_flow(elem_rhs, elem1d, ie, nrml, djac, dt, flux, h_avg, u, v, DEBUG, DEBUG_LOCAL);
        fe_sw2_1d_wall_friction(elem_rhs, elem1d, ie, nrml, djac, dt, u, v, resistance, DEBUG, DEBUG_LOCAL);
    }
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   FLAPU CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates a flag gate addition to the elemental residual. \n
     *
     *********************************************************************************************/
    else if(str_values[string].ol_flow.bc_flag == BCT_FLAPU) {
        int iwn = str_values[string].flap_num;
        double flux = flap[iwn].fluxu;
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_VEL_NEU :: flux: %20.10e \n",flux);
#endif
        fe_sw2_1d_explicit_flow(elem_rhs, elem1d, ie, nrml, djac, dt, flux, h_avg, u, v, DEBUG, DEBUG_LOCAL);
        fe_sw2_1d_wall_friction(elem_rhs, elem1d, ie, nrml, djac, dt, u, v, resistance, DEBUG, DEBUG_LOCAL);
    }
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   FLAPD CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates a flap gate addition to the elemental residual. \n
     *
     *********************************************************************************************/
    else if(str_values[string].ol_flow.bc_flag == BCT_FLAPD) {
        int iwn = str_values[string].flap_num;
        double flux = flap[iwn].fluxd;
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_VEL_NEU :: flux: %20.10e \n",flux);
#endif
        fe_sw2_1d_explicit_flow(elem_rhs, elem1d, ie, nrml, djac, dt, flux, h_avg, u, v, DEBUG, DEBUG_LOCAL);
        fe_sw2_1d_wall_friction(elem_rhs, elem1d, ie, nrml, djac, dt, u, v, resistance, DEBUG, DEBUG_LOCAL);
    }
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   SLUICEU CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates a sluice addition to the elemental residual. \n
     *
     *********************************************************************************************/
    else if(str_values[string].ol_flow.bc_flag == BCT_SLUICEU) {
        int iwn = str_values[string].sluice_num;
        double flux = sluice[iwn].fluxu;
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_VEL_NEU :: flux: %20.10e \n",flux);
#endif
        fe_sw2_1d_explicit_flow(elem_rhs, elem1d, ie, nrml, djac, dt, flux, h_avg, u, v, DEBUG, DEBUG_LOCAL);
        fe_sw2_1d_wall_friction(elem_rhs, elem1d, ie, nrml, djac, dt, u, v, resistance, DEBUG, DEBUG_LOCAL);
    }
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   SLUICED CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates a sluice addition to the elemental residual. \n
     *
     *********************************************************************************************/
    else if(str_values[string].ol_flow.bc_flag == BCT_SLUICED) {
        int iwn = str_values[string].sluice_num;
        double flux = sluice[iwn].fluxd;
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_VEL_NEU :: flux: %20.10e \n",flux);
#endif
        fe_sw2_1d_explicit_flow(elem_rhs, elem1d, ie, nrml, djac, dt, flux, h_avg, u, v, DEBUG, DEBUG_LOCAL);
        fe_sw2_1d_wall_friction(elem_rhs, elem1d, ie, nrml, djac, dt, u, v, resistance, DEBUG, DEBUG_LOCAL);
    }
#endif // _ADH_STRUCTURES
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                 DISCHARGE CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the natural discharge addition to the elemental residual. \n
     * \note
     *********************************************************************************************/
    else if (str_values[string].ol_flow.bc_flag == BCT_DIS_NEU) {
        if (elem_head[0] < SMALL || elem_head[1] < SMALL)  return;
        double discharge[2] = { 0., 0. };
        double troughness[2] = { 0., 0. };
        double fact = 0., convey = 0., three_2 = 3.0 / 2.0;
        
        convey = str_values[string].conveyance;
        troughness[0] = fe_sw2_get_roughness(mod, elem_head[0], mag_utan, string, UNUSED);
        troughness[1] = fe_sw2_get_roughness(mod, elem_head[1], mag_utan, string, UNUSED);
        
        int isers = str_values[string].ol_flow.isigma;
        fact = -1.0 * sseries_get_value(isers, mod->series_head,0) * sqrt(2. * g / troughness[0]) / convey;
        discharge[0] = fact * pow(elem_head[0], three_2);
        fact = -1.0 * sseries_get_value(isers, mod->series_head,0) * sqrt(2. * g / troughness[1]) / convey;
        discharge[1] = fact * pow(elem_head[0], three_2);

        double elem_avg_discharge = one_2 * (discharge[0] + discharge[1]);
        fe_sw2_1d_explicit_flow(elem_rhs, elem1d, ie, nrml, djac, dt, elem_avg_discharge, h_avg, u, v, DEBUG, DEBUG_LOCAL);
        fe_sw2_1d_wall_friction(elem_rhs, elem1d, ie, nrml, djac, dt, u, v, resistance, DEBUG, DEBUG_LOCAL);
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_DIS_NEU :: discharge: %20.10e, %20.10e \n",discharge[0],discharge[1]);
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                             TAILWATER/PRESSURE CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the tailwater or pressure addition to the elemental residual. \n
     * \note CJT \:: calculates implicit flow using user provides water elevation
     *********************************************************************************************/
    else if (str_values[string].ol_flow.bc_flag == BCT_PRS_NEU) {
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_PRS_NEU \n");
#endif
        
        double tail_elev = 0.;
        double elem_depth[2] = { 0., 0. };
        
        int isers = str_values[string].ol_flow.isigma;
        tail_elev = sseries_get_value(isers, mod->series_head,0);
        
        elem_depth[0] = tail_elev - elem_node_z[0];
        elem_depth[1] = tail_elev - elem_node_z[0];
        if (elem_depth[0] < SMALL || elem_depth[1] < SMALL || elem_head[0] < SMALL || elem_head[1] < SMALL) return;
        
        fe_sw2_1d_implicit_flow(elem_rhs, elem1d, ie, nrml, djac, dt, elem_vel, elem_head, DEBUG, DEBUG_LOCAL);
        
        if (PRESSURE_FLAG == ON) {
            fe_sw2_1d_pressure(elem_rhs, elem1d, ie, nrml, djac, dt, elem_depth, g, DEBUG, DEBUG_LOCAL);
            fe_sw2_1d_pressure(elem_rhs, elem1d, ie, nrml, djac, dt, elem_head, -g, DEBUG, DEBUG_LOCAL);
            
            if (mod->flag.BAROCLINIC != OFF) { // density impacts
                fe_sw2_1d_density_pressure(elem_rhs, elem1d, ie, nrml, djac, dt, elem_depth, elem_density,  g, DEBUG, DEBUG_LOCAL);
                fe_sw2_1d_density_pressure(elem_rhs, elem1d, ie, nrml, djac, dt, elem_head,  elem_density, -g, DEBUG, DEBUG_LOCAL);
            }
            
        }
    }
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                      ELEVATION + VELOCITIES BOUNDARY CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates applied tailwater + velocity natural boundaries to the elemental residual. \n
     * \note CJT \::
     *********************************************************************************************/
    else if (str_values[string].ol_flow.bc_flag == BCT_OVH_NEU) {
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_OVH_NEU \n");
#endif
        
        // get user elevation
        int isers = str_values[string].ol_flow.isigma;
        double tail_elev = sseries_get_value(isers, mod->series_head,0);
        double elem_depth[2] = { 0., 0. };
        elem_depth[0] = tail_elev - elem_node_z[0];
        elem_depth[1] = tail_elev - elem_node_z[1];
        if (elem_depth[0] < SMALL || elem_depth[1] < SMALL) {
            tl_error("assigned tailwater elevation is lower than the bed.");
        }
        
        // get user velocities
        SVECT2D vel_def;
        isers = str_values[string].ol_flow.ivx;
        vel_def.x = sseries_get_value(isers, mod->series_head,0);
        isers = str_values[string].ol_flow.ivy;
        vel_def.y = sseries_get_value(isers, mod->series_head,0);
        
        double user_nrml_vel = vel_def.x * nrml.x + vel_def.y * nrml.y;
        if (user_nrml_vel < -SMALL) { // assigned flow is into the model
            fe_sw2_1d_hvel(elem_rhs, elem1d, ie, nrml, djac, dt, elem_vel, vel_def, elem_head, DEBUG, DEBUG_LOCAL);
        } else {             // assigned flow is out of the model
            fe_sw2_1d_implicit_flow(elem_rhs, elem1d, ie, nrml, djac, dt, elem_vel, elem_head, DEBUG, DEBUG_LOCAL);
        }
        
        if (PRESSURE_FLAG == ON) {
            fe_sw2_1d_pressure(elem_rhs, elem1d, ie, nrml, djac, dt, elem_depth, g, DEBUG, DEBUG_LOCAL);
            fe_sw2_1d_pressure(elem_rhs, elem1d, ie, nrml, djac, dt, elem_head, -g, DEBUG, DEBUG_LOCAL);
            if (mod->flag.BAROCLINIC != OFF) { // density impacts
                fe_sw2_1d_density_pressure(elem_rhs, elem1d, ie, nrml, djac, dt, elem_depth, elem_density,  g, DEBUG, DEBUG_LOCAL);
                fe_sw2_1d_density_pressure(elem_rhs, elem1d, ie, nrml, djac, dt, elem_head,  elem_density, -g, DEBUG, DEBUG_LOCAL);
            }
        }
    }
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                DIRICHLET VELOCITY CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates Dirichlet velocity contributions to the elemental residual. \n
     * \note
     *********************************************************************************************/
    else if (str_values[string].ol_flow.bc_flag == BCT_DB_VEL) {
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_DB_VEL \n");
#endif
        fe_sw2_1d_implicit_flow(elem_rhs, elem1d, ie, nrml, djac, dt, elem_vel, elem_head, DEBUG, DEBUG_LOCAL);
    }
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                          MONOLITHIC COUPLING INTERFACE CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates monolithich coupling contributions to the elemental residual. \n
     * \note
     *********************************************************************************************/
    else if (str_values[string].ol_flow.bc_flag == BCT_HYBRID_INTERNAL) {
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_HYBRID_INTERNAL \n");
#endif
        if (PRESSURE_FLAG == ON) {
            // the following line is added to eliminate the integration done in 2d in fe_sw2_boundpressure()
            fe_sw2_1d_pressure(elem_rhs, elem1d, ie, nrml, djac, dt, elem_head, -g, DEBUG, DEBUG_LOCAL);
            if (mod->flag.BAROCLINIC != OFF) { // density impacts
                fe_sw2_1d_density_pressure(elem_rhs, elem1d, ie, nrml, djac, dt, elem_head,  elem_density, -g, DEBUG, DEBUG_LOCAL);
            }
        }
    }
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    else {
        printf("Error in boundary type specification\n The boundary type is %d on string %d\n", str_values[string].ol_flow.bc_flag, string + 1);
        tl_error("The overland flow 1D elem. boundary conditions didn't match available types.");
    }
    
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                DIRICHLET BOUNDARY CONDITIONS
     *--------------------------------------------------------------------------------------------
     * Zeros matrix rows/columns when Dirichlet boundaries are used. \n
     * \note
     *********************************************************************************************/
    for (i=0; i<nnodes; i++) {
        string = mod->grid->node[ elem1d.nodes[i] ].string;
        if (string > NORMAL) {
            if (str_values[string].ol_flow.bc_flag == BCT_PRS_DIR)
                elem_rhs[i].c_eq = 0.0;
            else if (str_values[string].ol_flow.bc_flag == BCT_VEL_DIR) {
                elem_rhs[i].x_eq = 0.0;
                elem_rhs[i].y_eq = 0.0;
            }
            else if (str_values[string].ol_flow.bc_flag == BCT_VEL_PRS_DIR) {
                elem_rhs[i].c_eq = 0.0;
                elem_rhs[i].x_eq = 0.0;
                elem_rhs[i].y_eq = 0.0;
            }
        }
    }
    
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        rhs_3dof("1D SW BOUNDARY || TOTAL RHS AFER DB",nnodes, ie, elem1d.nodes, elem_rhs);
        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    }
    
    time_t time2;  time(&time2);
    TIME_IN_2D_SW_BOUNDARY_RESID += difftime(time2,time1);
#endif
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 1D element explicit flow contributions to 2D-SW equations
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elem_rhs    the 1D elemental residual array
 * @param[in]  ie            the 1D element ID
 * @param[in]  nrml          the 1D element normal
 * @param[in]  djac          the 1D element jacobian
 * @param[in]  dt            the current model time-step
 * @param[in]  flux          a user-supplied flux
 * @param[in]  h_avg         the elementally averaged depth
 * @param[in]  u             x-velocity
 * @param[in]  v             y-velocity
 *
 * \f{eqnarray*}{  \residDA{i}{}{c}   &=& dt * \intbe{\phidd{i} \, q_n}{e}    \f}
 *
 * For the momentum contributions: \n
 * -- if the flux is into the domain, the normal flux is split into x,y velocities, so that
 *  \f{eqnarray*}{
 *     u^* = (q_n / \overline{\depth{h}}) \, n_x \\
 *     v^* = (q_n / \overline{\depth{h}}) \, n_y
 *  \f}
 * and added to the residual as
 *  \f{eqnarray*}{
 *   \residDA{i}{}{mx}  &=& dt * \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \ub{h} \, \depth{h})} = dt * \intbe{\phidd{i} \, q_n \, u^*}{e} \\
 *   \residDA{i}{}{my}  &=& dt * \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \vb{h} \, \depth{h})} = dt * \intbe{\phidd{i} \, q_n \, v^*}{e}
 *  \f}
 * -- if the flux is out of the domain, an implicity u,v is used so that
 *  \f{eqnarray*}{
 *   \residDA{i}{}{mx}  &=& dt * \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \ub{h} \, \depth{h})} = dt * \intbe{\phidd{i} \, q_n \, u}{e} \\
 *   \residDA{i}{}{my}  &=& dt 8 \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \vb{h} \, \depth{h})} = dt * \intbe{\phidd{i} \, q_n \, v}{e}
 *  \f}
 *
 *  \details  Forms the 1d element matrix for the flux terms in the shallow water equations. \n
 *            This is the unit discharge out of the boundary. Positive is out, Negative is in.
 * NOTE: CJT \:: all units here must be equivalent to Int [d(uh)/dt * dxdy] = m^4/s^2
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_1d_explicit_flow(DOF_3 *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, double flux, double h_avg, double *u, double *v, int DEBUG, int DEBUG_LOCAL) {
    
    int i;
    double rhs_c_eq[NDONSEG], rhs_x_eq[NDONSEG], rhs_y_eq[NDONSEG]; // local rhs
    sarray_init_dbl(rhs_c_eq, NDONSEG);
    sarray_init_dbl(rhs_x_eq, NDONSEG);
    sarray_init_dbl(rhs_y_eq, NDONSEG);
    
    integrate_line_phi(djac, 1., rhs_c_eq);
    
    if (flux < 0) {    // flow is into the model
        double u_star = 0., v_star = 0.;
        if (h_avg > 0) {
            u_star = flux * nrml.x / h_avg; // average directional velocity
            v_star = flux * nrml.y / h_avg; // average directional velocity
        }
        integrate_line_phi(djac, u_star, rhs_x_eq);
        integrate_line_phi(djac, v_star, rhs_y_eq);
        
    } else {           // flow is out of the model
        integrate_line_phi_f(djac, 1., u, rhs_x_eq);
        integrate_line_phi_f(djac, 1., v, rhs_y_eq);
    }
    
    double t1 = dt * flux;
    for (i=0; i<NDONSEG; i++) {
        elem_rhs[i].c_eq += t1 * rhs_c_eq[i];
        elem_rhs[i].x_eq += t1 * rhs_x_eq[i];
        elem_rhs[i].y_eq += t1 * rhs_y_eq[i];
    }
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        DOF_3 rhs[NDONSEG];
        for (i=0; i<NDONSEG; i++) {
            rhs[i].c_eq = t1 * rhs_c_eq[i];
            rhs[i].x_eq = t1 * rhs_x_eq[i];
            rhs[i].y_eq = t1 * rhs_y_eq[i];
        }
        rhs_3dof("1D SW EXPLICIT FLOW: ",NDONSEG, ie, elem1d.nodes, rhs);
    }
#endif
}



/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 1D element implicit flow contributions to 2D-SW equations
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elem_rhs    the 1D elemental residual array
 * @param[in]  ie            the 1D element ID
 * @param[in]  nrml          the 1D element normal
 * @param[in]  djac          the 1D element jacobian
 * @param[in]  dt            the current model time-step
 * @param[in]  v             the hydrodynamic velocity
 * @param[in]  h             the water depth
 *
 * \f{eqnarray*}{  \residDA{i}{}{c}   &=& dt * \intbe{\phidd{i} \, q_n}{e} \\
 *                 \residDA{i}{}{mx}  &=& dt * \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \ub{h} \, \depth{h})} \\
 *                 \residDA{i}{}{my}  &=& dt * \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \vb{h} \, \depth{h})}
 * \f}
 *
 * \note CJT\::  Gaurav *thinks the "r" variable is a prox for v so he can eliminate eddies forming on the boundary
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_1d_implicit_flow(DOF_3 *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, SVECT2D *v, double *h, int DEBUG, int DEBUG_LOCAL) {
    
    int i;
    double rhs_c_eq[NDONSEG], rhs_x_eq[NDONSEG], rhs_y_eq[NDONSEG]; // local rhs
    sarray_init_dbl(rhs_c_eq, NDONSEG);
    sarray_init_dbl(rhs_x_eq, NDONSEG);
    sarray_init_dbl(rhs_y_eq, NDONSEG);
    
    double u1 =  v[0].x, u2 = v[1].x;
    double v1 =  v[0].y, v2 = v[1].y;
    double h1 =  h[0],   h2 = h[1];
    double ur1 = v[0].x, ur2 = v[1].x;
    double vr1 = v[0].y, vr2 = v[1].y;
    
    
    double dir_vect[NDONSEG];
    dir_vect[0] = v[0].x * nrml.x + v[0].y * nrml.y;
    dir_vect[1] = v[1].x * nrml.x + v[1].y * nrml.y;
    if (dir_vect[0] < SMALL) {
        ur1 = 0.;
        vr1 = 0.;
    }
    if (dir_vect[1] < SMALL) {
        ur2 = 0.;
        vr2 = 0.;
    }
    
    integrate_line_phi_h_v_dot_n(djac, 1., h, v, nrml, rhs_c_eq); // phi_i * u * h * nx + phi_i * v * h * ny
    integrate_line_phi_f_g_h(djac, nrml.x, ur1, ur2, u1, u2, h1, h2, rhs_x_eq); // phi_i * ur * u * h * nx
    integrate_line_phi_f_g_h(djac, nrml.y, vr1, vr2, v1, v2, h1, h2, rhs_y_eq); // phi_i * vr * v * h * ny
    integrate_line_phi_f_g_h(djac, nrml.y, ur1, ur2, v1, v2, h1, h2, rhs_x_eq); // phi_i * ur * v * h * ny
    integrate_line_phi_f_g_h(djac, nrml.x, vr1, vr2, u1, u2, h1, h2, rhs_y_eq); // phi_i * vr * u * h * nx
    
    for (i=0; i<NDONSEG; i++) {
        elem_rhs[i].c_eq += dt * rhs_c_eq[i];
        elem_rhs[i].x_eq += dt * rhs_x_eq[i];
        elem_rhs[i].y_eq += dt * rhs_y_eq[i];
    }
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        DOF_3 rhs[NDONSEG];
        for (i=0; i<NDONSEG; i++) {
            rhs[i].c_eq = dt * rhs_c_eq[i];
            rhs[i].x_eq = dt * rhs_x_eq[i];
            rhs[i].y_eq = dt * rhs_y_eq[i];
        }
        rhs_3dof("1D SW IMPLICIT FLOW: ",NDONSEG, ie, elem1d.nodes, rhs);
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 1D element wall friction contributions to 2D-SW equations
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elem_rhs    the 1D elemental residual array
 * @param[in]  ie            the 1D element ID
 * @param[in]  nrml          the 1D element normal
 * @param[in]  djac          the 1D element jacobian
 * @param[in]  dt            the current model time-step
 * @param[in]  u             x-velocity
 * @param[in]  v             y-velocity
 * @param[in]  resistance    resistance
 *
 *  \details Integrates the weak, discrete outflow boundary terms: \n
 *  \f{eqnarray*}{
 *   \residDA{i}{}{mx}  &=& dt * \bcFriction{}{e}{\phidd{i}}{\ub{h}}{\velb{h}} \\
 *   \residDA{i}{}{my}  &=& dt * \bcFriction{}{e}{\phidd{i}}{\vb{h}}{\velb{h}}
 *  \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_1d_wall_friction(DOF_3 *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, double *u, double *v, double *resistance, int DEBUG, int DEBUG_LOCAL) {
    
    int i;
    double rhs_x_eq[NDONSEG], rhs_y_eq[NDONSEG]; // local rhs
    sarray_init_dbl(rhs_x_eq, NDONSEG);
    sarray_init_dbl(rhs_y_eq, NDONSEG);
    
    integrate_line_phi_f(djac, fabs(nrml.y), u, rhs_x_eq);
    integrate_line_phi_f(djac, fabs(nrml.x), v, rhs_y_eq);
    
    for (i=0; i<NDONSEG; i++) {
        elem_rhs[i].x_eq += dt * resistance[i] * rhs_x_eq[i];
        elem_rhs[i].y_eq += dt * resistance[i] * rhs_y_eq[i];
    }
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        DOF_3 rhs[NDONSEG];
        for (i=0; i<NDONSEG; i++) {
            rhs[i].c_eq = 0;
            rhs[i].x_eq = dt * resistance[i] * rhs_x_eq[i];
            rhs[i].y_eq = dt * resistance[i] * rhs_y_eq[i];
        }
        rhs_3dof("1D SW WALL FRICTION: ",NDONSEG, ie, elem1d.nodes, rhs);
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Returns the 1D element pressure contributions to 2D-SW equations
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout] elem_rhs    the 1D elemental residual array
 * @param[in]  ie            the 1D element ID
 * @param[in]  nrml          the 1D element normal
 * @param[in]  djac          the 1D element jacobian
 * @param[in]  dt            the current model time-step
 * @param[in]  h             the water depth
 * @param[in]  g             gravity
 *
 * \f{eqnarray*}{
 *                 \residDA{i}{}{mx}  &=& dt * \bcPressureDD{}{e}{\phidd{i}}{(\depth{h})}{x} \\
 *                 \residDA{i}{}{my}  &=& dt * \bcPressureDD{}{e}{\phidd{i}}{(\depth{h})}{y}
 * \f}
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_1d_pressure(DOF_3 *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, double *h, double g, int DEBUG, int DEBUG_LOCAL) {
    
    int i;
    double rhs[NDONSEG] = {0., 0.};
    
    
    integrate_line_phi_h_h(djac, 1., h[0], h[1], rhs);
    
    double t1 = dt * one_2 * g;
    for (i=0; i<NDONSEG; i++) {
        elem_rhs[i].x_eq += t1 * nrml.x * rhs[i];
        elem_rhs[i].y_eq += t1 * nrml.y * rhs[i];
    }
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        DOF_3 rhs2[NDONSEG];
        for (i=0; i<NDONSEG; i++) {
            rhs2[i].c_eq = 0;
            rhs2[i].x_eq = t1 * nrml.x * rhs[i];
            rhs2[i].y_eq = t1 * nrml.y * rhs[i];
        }
        rhs_3dof("1D SW PRESSURE: ",NDONSEG, ie, elem1d.nodes, rhs2);
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Returns the 1D element density pressure contributions to 2D-SW equations
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Gary Brown, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout] elem_rhs    the 1D elemental residual array
 * @param[in]  ie            the 1D element ID
 * @param[in]  nrml          the 1D element normal
 * @param[in]  djac          the 1D element jacobian
 * @param[in]  dt            the current model time-step
 * @param[in]  h             the water depth
 * @param[in]  d             water density
 * @param[in]  g             gravity
 *
 * \f{eqnarray*}{
 *                 \residDA{i}{}{mx}  &=& dt * \bcPressureDD{}{e}{\phidd{i}}{(\rho^h \, \depth{h})}{x} \\
 *                 \residDA{i}{}{my}  &=& dt * \bcPressureDD{}{e}{\phidd{i}}{(\rho^h \, \depth{h})}{y}
 * \f}
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_1d_density_pressure(DOF_3 *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, double *h, double *d, double g, int DEBUG, int DEBUG_LOCAL) {
    
    int i;
    double rhs[NDONSEG] = {0., 0.};
    
    integrate_line_phi_h_h_g(djac, 1., h[0], h[1], d[0], d[1], rhs);
    
    double t1 = dt * one_2 * g;
    for (i=0; i<NDONSEG; i++) {
        elem_rhs[i].x_eq += t1 * nrml.x * rhs[i];
        elem_rhs[i].y_eq += t1 * nrml.y * rhs[i];
    }
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        DOF_3 rhs2[NDONSEG];
        for (i=0; i<NDONSEG; i++) {
            rhs2[i].c_eq = 0;
            rhs2[i].x_eq = t1 * nrml.x * rhs[i];
            rhs2[i].y_eq = t1 * nrml.y * rhs[i];
        }
        rhs_3dof("1D SW DENSITY PRESSURE: ",NDONSEG, ie, elem1d.nodes, rhs2);
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 1D element explicit flow and pressure contributions to 2D-SW equations
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elem_rhs    the 1D elemental residual array
 * @param[in]  ie            the 1D element ID
 * @param[in]  nrml          the 1D element normal
 * @param[in]  djac          the 1D element jacobian
 * @param[in]  dt            the current model time-step
 * @param[in]  v             the water velocity
 * @param[in]  user_velocity the water velocity prescribed by user
 * @param[in]  h             the water depth
 *
 * \f{eqnarray*}{  \residDA{i}{}{c}   &=& dt * \intbe{\phidd{i} \, (\velb{h} \cdot \nrml) \, \depth{h}}{e} \\
 *                 \residDA{i}{}{mx}  &=& dt * \intbe{\phidd{i} \, q_n \, \ub{}}{e}\\
 *                 \residDA{i}{}{my}  &=& dt * \intbe{\phidd{i} \, q_n \, \vb{}}{e}
 * \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_1d_hvel(DOF_3 *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, SVECT2D *v, SVECT2D user_velocity, double *h, int DEBUG, int DEBUG_LOCAL){
    
    int i;
    double rhs[NDONSEG] = {0., 0.};
    integrate_line_phi_h_v_dot_n(djac, 1., h, v, nrml, rhs); // phi_i * u * h * nx + phi_i * v * h * ny
    
    for (i=0; i<NDONSEG; i++) {
        elem_rhs[i].c_eq += dt * rhs[i];
        elem_rhs[i].x_eq += dt * rhs[i] * user_velocity.x;
        elem_rhs[i].y_eq += dt * rhs[i] * user_velocity.y;
    }
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        DOF_3 rhs2[NDONSEG];
        for (i=0; i<NDONSEG; i++) {
            rhs2[i].c_eq = dt * rhs[i];
            rhs2[i].x_eq = dt * rhs[i] * user_velocity.x;
            rhs2[i].y_eq = dt * rhs[i] * user_velocity.y;
        }
        rhs_3dof("1D SW HVEL: ",NDONSEG, ie, elem1d.nodes, rhs2);
    }
#endif
}







































