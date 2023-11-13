/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes the boundary residual contributions for the SW-2D transport model.
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
 * @param[in]  PRESSURE      a flag to include pressure terms (TRUE/FALSE)
 * @param[in]  DEBUG         a debug flag
 *
 *  \details Integrates the weak, discrete outflow boundary terms: \n
 *  \f{eqnarray*}{
 *   \overline{R}_{i,boundary}^e &=& dt * \bigg( \bcConv{}{e}{\phidd{i}}{(\velcDA{h}{t} \, \cjhDA \, \depth{h})} \bigg)
 *                     =  dt * \bigg( \bcConv{}{e}{\phidd{i}}{(q_n^{hydro} \, \cjhDA)} \bigg)
 *  \f}
 * 
 * \note Requires water velocity and depth
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_2d_transport_boundary_resid(SMODEL *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    
    int DEBUG_LOCAL = OFF;
    int DEBUG_PICKETS = OFF;
#ifdef _DEBUG
    if (DEBUG == ON) DEBUG_LOCAL = ON;
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    time_t time1;  time(&time1);
#endif
    
    int i;
    
    // aliases
    SSW_2D *sw = mod->sw->d2;
    double dt = mod->dt;
    double gravity = mod->gravity;
    
    int index = UNSET_INT;
    double c_inv = SMALL;
    double *c, *old_c, *older_c;
    int string = mod->grid->elem1d[ie].string;
    INPUT_DATA *str_value_hydro = &(mod->str_values[string].ol_flow);
    INPUT_DATA *str_value_con;
    if (mod->is_sediment_running) {
#ifdef _SEDIMENT
        index = mod->ised;
        c = mod->sed->susload[index].c;
        old_c = mod->sed->susload[index].old_c;
        older_c = mod->sed->susload[index].older_c;
        c_inv = 1.E-6 / mod->sed->grain[index].reference_c;
        str_value_con = &(mod->str_values[string].sed[index]);
#endif
    } else {
        index = mod->itrns;
        c = mod->con[index].concentration;
        old_c = mod->con[index].old_concentration;
        older_c = mod->con[index].older_concentration;
        c_inv = 1. / mod->con[index].property[0];
        str_value_con = &(mod->str_values[string].trans[index]);
    }
    int con_bc_flag = str_value_con->bc_flag;
    int hydro_bc_flag = str_value_hydro->bc_flag;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GRID VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    SGRID *grid = mod->grid;
    SELEM_1D *elem1d = &(grid->elem1d[ie]);
    SVECT2D nrml = elem1d->nrml;
    double djac = elem1d->djac;
    int nnodes = elem1d->nnodes;

#ifdef _DEBUG
    assert(nnodes == NDONSEG); // only 1D element is a line segment
    assert(djac > 0); // only 1D element is a line segment
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    double elem_c[nnodes];
    global_to_local_dbl(c, elem_c, elem1d->nodes, nnodes);
    if (perturb_var == PERTURB_C) elem_c[perturb_node] += perturb_sign * perturbation;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // constituent concentrations
    double elem_old_c[nnodes], elem_older_c[nnodes];
    global_to_local_dbl(old_c, elem_old_c, elem1d->nodes, nnodes);
    global_to_local_dbl(older_c, elem_older_c, elem1d->nodes, nnodes);
    
    // water depth
    double elem_head[nnodes], elem_old_head[nnodes], elem_older_head[nnodes];
    global_to_local_dbl(sw->head, elem_head, elem1d->nodes, nnodes);
    global_to_local_dbl(sw->old_head, elem_old_head, elem1d->nodes, nnodes);
    global_to_local_dbl(sw->older_head, elem_older_head, elem1d->nodes, nnodes);
 
    // water velocities
    SVECT2D elem_vel[nnodes], elem_old_vel[nnodes], elem_older_vel[nnodes];
    global_to_local_svect2d(sw->vel, elem_vel, elem1d->nodes, nnodes);
    global_to_local_svect2d(sw->old_vel, elem_old_vel, elem1d->nodes, nnodes);
    global_to_local_svect2d(sw->older_vel, elem_older_vel, elem1d->nodes, nnodes);
    double u[NDONSEG], v[NDONSEG];
    for (i=0; i<nnodes; i++) { u[i] = elem_vel[i].x; v[i] = elem_vel[i].y; }
    
    // sediment velocities
    SVECT2D elem_vel_sed[nnodes];
    if (mod->is_sediment_running) {
#ifdef _SEDIMENT
        int gnode;
        for (i=0; i<nnodes; i++) {
            gnode = elem1d->nodes[i];
            elem_vel_sed[i].x = (elem_vel[i].x + mod->sed->susload[index].vcf[gnode].x) * mod->sed->susload[index].mfcf[gnode];
            elem_vel_sed[i].y = (elem_vel[i].y + mod->sed->susload[index].vcf[gnode].y) * mod->sed->susload[index].mfcf[gnode];
        }
#endif
    } else {
        for (i=0; i<nnodes; i++) {
            elem_vel_sed[i].x = elem_vel[i].x;
            elem_vel_sed[i].y = elem_vel[i].y;
        }
    }
    
    // sediment velocity normal to the 1d element
    double elem_vel_sed_nrml[nnodes];
    elem_vel_sed_nrml[0] = svect2d_dotp(elem_vel_sed[0], nrml);
    elem_vel_sed_nrml[1] = svect2d_dotp(elem_vel_sed[1], nrml);
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        printf("SW-2D BOUNDARY ELEM RESID :: ie: %d \t dt: %20.10f \t djac: %20.10f nrml: {%10.5e, %10.5e}\n",ie,dt,djac,nrml.x,nrml.y);
        if (perturb_var == PERTURB_C) {
            printf("\t perturbing C || node: %d || perturbation: %20.10e\n",elem1d->nodes[perturb_node],perturb_sign*perturbation);
        }
        printScreen_debug2_dbl("elem_c", elem_c, nnodes, elem1d->nodes);
        printScreen_debug2_dbl("elem_c", elem_vel_sed_nrml, nnodes, elem1d->nodes);
        printScreen_debug2_dbl("elem_old_c", elem_old_c, nnodes, elem1d->nodes);
        printScreen_debug2_dbl("elem_older_c", elem_older_c, nnodes, elem1d->nodes);
        printScreen_debug2_dbl("elem_head", elem_head, nnodes, elem1d->nodes);
        printScreen_debug2_dbl("elem_old_head", elem_old_head, nnodes, elem1d->nodes);
        printScreen_debug2_dbl("elem_older_head", elem_older_head, nnodes, elem1d->nodes);
        printScreen_debug_svec2d("elem_sed_vel", elem_vel, nnodes, elem1d->nodes);
        printScreen_debug_svec2d("elem_vel", elem_vel, nnodes, elem1d->nodes);
        printScreen_debug_svec2d("elem_old_vel", elem_old_vel, nnodes, elem1d->nodes);
        printScreen_debug_svec2d("elem_older_vel", elem_older_vel, nnodes, elem1d->nodes);
    }
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    sarray_init_dbl(elem_rhs, nnodes);
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   VELOCITY CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the convection boundary flux addition to the SW 2D elemental residual. \n
     * \note  If the hydro flux is into the model (-), then the constituent flux must be assigned - it is set to 0 if not.
     * \note  If the hydro flux is out of the model, then the assigned constituent flux is ignored and calculated implicitly.
     *********************************************************************************************/
    if(hydro_bc_flag == BCT_VEL_NEU) {
        
        // get user input water flux
        int isers = str_value_hydro->isigma;
        double hydro_flux_normal = -1.0 * sseries_get_value(isers, mod->series_head, 0);

        if(hydro_flux_normal < -SMALL) {
            // get user input consitutent flux if given, otherwise set to 0
            if(con_bc_flag > NORMAL) {
                int isers = str_value_con->isigma;
                double user_c = sseries_get_value(isers, mod->series_head, 0) * c_inv;
                integrate_line_phi(djac, user_c * hydro_flux_normal * dt, elem_rhs);
            }
        } else {
            integrate_line_phi_f(djac, hydro_flux_normal * dt, elem_c, elem_rhs);
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("2D TRANSPORT BCT_VEL_NEU: ",NDONSEG, ie, elem1d->nodes, elem_rhs);
            Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
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
        integrate_line_phi_f_g_h(djac, dt, elem_vel_sed_nrml[0], elem_vel_sed_nrml[1], elem_head[0], elem_head[1], elem_c[0], elem_c[1], elem_rhs);
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("2D TRANSPORT BCT_OUTFLOW: ",NDONSEG, ie, elem1d->nodes, elem_rhs);
            Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
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
    else if (hydro_bc_flag == BCT_PRS_NEU || hydro_bc_flag == BCT_DIS_NEU) {
        double hydro_flux = one_2 * (elem_vel_sed_nrml[0] * elem_head[0] + elem_vel_sed_nrml[1] * elem_head[1]);
        if(hydro_flux >= 0.) {
            integrate_line_phi_f_g_h(djac, dt, elem_vel_sed_nrml[0], elem_vel_sed_nrml[1], elem_head[0], elem_head[1], elem_c[0], elem_c[1], elem_rhs);
        } else {
            // get user input consitutent flux if given, otherwise set to 0
            if(con_bc_flag == BCT_NEU) {
                int isers = str_value_con->isigma;
                double user_c = sseries_get_value(isers, mod->series_head, 0) * c_inv;
                integrate_line_phi_f_g(djac, user_c * dt, elem_head, elem_vel_sed_nrml, elem_rhs);
            }
        }
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("2D TRANSPORT BCT_PRS_NEU ",NDONSEG, ie, elem1d->nodes, elem_rhs);
            Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
        }
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                  DRY MASS ALONG A LINE
     *--------------------------------------------------------------------------------------------
     * Not sure what this is for. \n
     * \note
     *********************************************************************************************/
    else if(hydro_bc_flag == BCT_NEU_LOAD) {
        int isers = str_value_con->isigma;
        if(isers != UNSET_INT) {
            double user_c = -sseries_get_value(isers, mod->series_head, 0) * c_inv;
            integrate_line_phi(djac, user_c * dt, elem_rhs);
#ifdef _DEBUG
            if (DEBUG_LOCAL == ON) {
                rhs_1dof("2D TRANSPORT BCT_NEU_LOAD ",NDONSEG, ie, elem1d->nodes, elem_rhs);
                Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
            }
#endif
        }

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
        if (elem_head[i] < 0.0) elem_rhs[i] = 0.0;
        if (str_value_hydro->bc_flag == BCT_DIS_NEU) {
            if (elem_old_head[i] < 0.0) elem_rhs[i] = 0.0;
        }
        // note :: get node not element string here
        int node_string = mod->grid->node[ elem1d->nodes[i] ].string;
        if (node_string > NORMAL) {
            if (str_value_con->bc_flag == BCT_DIR || str_value_con->bc_flag == BCT_CEQ) {
                elem_rhs[i] = 0;
            }
        }
    }
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        rhs_1dof("2D TRANSPORT DIRICHLET BC ",NDONSEG, ie, elem1d->nodes, elem_rhs);
        Is_DoubleArray_Inf_or_NaN(elem_rhs ,nnodes ,__FILE__ ,__LINE__);
    }
#endif
    
#ifdef _DEBUG
    time_t time2;  time(&time2);
    TIME_IN_2D_TRANSPORT_BOUNDARY_RESID += difftime(time2,time1);
#endif
    

}
