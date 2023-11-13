/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes the 1D element residual for the diffusive wave model.
 *  \author    Charlie Berger, Ph.D.
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
 * @param[in]  perturbation   the Newton perturbation
 * @param[in]  perturb_node  the node to be perturbed
 * @param[in]  perturb_var   the variable to be perturbed
 * @param[in]  perturb_sign  the direction of Newton perturbation
 * @param[in]  PRESSURE      a flag to include pressure terms (TRUE/FALSE)
 * @param[in]  DEBUG         a debug flag
 *
 *  \details Integrates the weak, discrete outflow boundary terms: \n
 *  \f{eqnarray*}{
 *   R_{i,boundary}^e &=& dt * \bigg( \bcConv{}{e}{\phidd{i}}{(\velb{h} \,\depth{h})} \bigg)
 *  \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_diffusive_boundary_resid(SMODEL *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    
    int i, DEBUG_LOCAL = OFF;
    
    /*
     if (mod->grid->elem1d[ie].nodes[0]==0 ||
     mod->grid->elem1d[ie].nodes[1]==0 ||
     mod->grid->elem1d[ie].nodes[2]==0) {
     DEBUG_LOCAL = ON;
     } else {
     DEBUG_LOCAL = OFF;
     }
     */
    
    // aliases
    SSW_2D *sw2 = mod->sw->d2;
    SGRID *grid = mod->grid;
    SELEM_1D *elem1d = &(grid->elem1d[ie]);// be careful not to shallow copy here
    int nnodes = elem1d->nnodes;
    double dt = mod->dt;
    double djac = elem1d->djac;
    double gravity = mod->gravity;
    double rhs[nnodes];
    int string = elem1d->string;
    STR_VALUE *str_values = mod->str_values;

    /* Gajanan gkc for correctly calculate 1D elem velocities using gradients from attached 2D elements. */
    int ie2d = grid->elem1d[ie].elem2d;
    SELEM_2D *elem2d = &(grid->elem2d[ie2d]);
    
#ifdef _DEBUG
    assert(nnodes == NDONSEG); // only 1D element is a line segment
    assert(djac > 0); // only 1D element is a line segment
    assert(perturb_node == 0 || perturb_node == 1 || perturb_node == UNSET_INT);
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    double elem_head[nnodes];
    global_to_local_dbl(sw2->head, elem_head, elem1d->nodes, nnodes);
    if (perturb_var == PERTURB_H) elem_head[perturb_node] += perturb_sign * perturbation;
    
    double elem2d_head[elem2d->nnodes];
    global_to_local_dbl(sw2->head, elem2d_head, elem2d->nodes, elem2d->nnodes);
    if (perturb_var == PERTURB_H) {
        for (i=0; i<elem2d->nnodes; i++) {
            if (elem2d->nodes[i]==elem1d->nodes[perturb_node]){
                elem2d_head[i] += perturb_sign * perturbation;
                assert(elem2d_head[i]==elem_head[perturb_node]); /* Yes. Floating point equality check. Intentional.*/
                break;
            }
        }
    }
    double node2d_z[elem2d->nnodes];
    SNODE nodes2d[elem2d->nnodes];
    for (i=0; i<elem2d->nnodes; i++) {
        snode_copy(&(nodes2d[i]), mod->grid->node[elem2d->nodes[i]]);
        node2d_z[i] = nodes2d[i].z;
    }
    SVECT elem2d_nds[elem2d->nnodes];
    snode2svect(nodes2d, elem2d_nds, elem2d->nnodes);


    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    double slope = 0.;
    SVECT2D elem_vel[nnodes]; svect2d_init_array(elem_vel, nnodes);
    for (i=0; i<nnodes; i++) {
        if (elem_head[i] > 0){
            //Gajanan gkc warning: I do not think this is the correct way to calculate the slope here.
            //I believe the 2D element that the 1D element is attached to must be used to calculate
            //the slope. This is because dh/dx and dh/dz are affected even by the 3rd node of the 
            //attached 2D element.
            slope = fabs((elem_head[i] + grid->node[elem1d->nodes[i]].z) * elem1d->grad_shp[i]);
            slope *= slope;
            slope = pow(slope,0.5);
            slope = pow(slope,0.5);
            
#ifdef _DEBUG
            assert(fabs(mod->str_values[string].fterms.manningsn) > SMALL);
#endif
            //Gajanan gkc warning: Surely this must be wrong? u[n1]=v[n1] and u[n2]=v[n2], always according to the below statement.
            //elem_vel[i].x = mod->manning_units_constant * slope *pow(elem_head[i], 2./3.)/mod->str_values[string].fterms.manningsn;
            //elem_vel[i].y = mod->manning_units_constant * slope *pow(elem_head[i], 2./3.)/mod->str_values[string].fterms.manningsn;
        }
    }
    // Gajanan gkc warning: Temporary work-around for the above problem.
    // Adding the line below. Note that although velocity in the residual
    // calculation is likely better with the line below, the Jacobian
    // contribution when perturbing depth is not reflected in this fix.
    // So velocity ends up acting like a constant for sake of Jacobian
    // calculation which should not be the case since vel depends on depth.
    global_to_local_svect2d(sw2->vel, elem_vel, elem1d->nodes, nnodes);
    
    double elem2d_avg_depth = 0., area = 0.;
    int isTriangle = NO;  if (elem2d->nnodes == NDONTRI) isTriangle = YES;
    if (isTriangle == TRUE) {
        area = elem2d->djac;
        elem2d_avg_depth = integrate_triangle_f(1.,1.,elem2d_head); // cjt :: djac = area cancel here
    } else {
        area = integrate_quadrilateral_area(elem2d_nds, 1.);
        elem2d_avg_depth = integrate_quadrilateral_f(elem2d_nds,1./area,elem2d_head);
    }
    SVECT2D elem2d_vel[elem2d->nnodes], vv;
    vv = getDiffusiveWaveVelocities(elem2d->nnodes, elem2d_head, elem2d_avg_depth, node2d_z, elem2d->grad_shp, mod->str_values[string].fterms.manningsn, mod->manning_units_constant, elem2d_vel);
    //printf("\n1D: mannings constant: %20.10e \t mannings n: %20.10e",
    //           mod->manning_units_constant,mod->str_values[string].fterms.manningsn);
#ifdef _DEBUG
    //printf("\nElem1d[%i], part of Elem2d[%i]: vavg.x=%f, vavg.y=%f", ie, ie2d, vv.x, vv.y);
#endif

    elem_vel[0].x = vv.x;    elem_vel[1].x = vv.x;
    elem_vel[0].y = vv.y;    elem_vel[1].y = vv.y;
//    for (i=0;i<elem2d->nnodes;i++){
//        if (elem2d->nodes[i]==elem1d->nodes[0]){
//            //if(elem_vel[0].x != elem2d_vel[i].x || elem_vel[0].y != elem2d_vel[i].y){
//            //    printf("\nE1d[%4i] : Before[0] : vx = % .9e, vy = % .9e, After[0] : vx = % .9e, vy = % .9e", ie,
//            //            elem_vel[0].x, elem_vel[0].y, elem2d_vel[i].x, elem2d_vel[i].y);
//            //}
//            elem_vel[0].x = elem2d_vel[i].x; // vv.x; //
//            elem_vel[0].y = elem2d_vel[i].y; // vv.y; //
//        }
//        else if (elem2d->nodes[i]==elem1d->nodes[1]){
//            //if(elem_vel[1].x != elem2d_vel[i].x || elem_vel[1].y != elem2d_vel[i].y){
//            //    printf("\nE1d[%4i] : Before[1] : vx = % .9e, vy = % .9e, After[1] : vx = % .9e, vy = % .9e", ie,
//            //            elem_vel[1].x, elem_vel[1].y, elem2d_vel[i].x, elem2d_vel[i].y);
//            //}
//            elem_vel[1].x = elem2d_vel[i].x; // vv.x; //
//            elem_vel[1].y = elem2d_vel[i].y; // vv.y; //
//        }
//#ifdef _DEBUG
//        //printf("\n    vel[%i].x=%f, vel[%i].y=%f", i, elem2d_vel[i].x, i, elem2d_vel[i].y);
//#endif
//    }

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        printf("DIFFUSIVE WAVE 1D ELEM RESID :: ie: %d \t dt: %20.10f \t djac: %20.10f \n",ie,dt,djac);
        if (perturb_var == PERTURB_H) {
            printf("\t perturbing H || node[%d]: %d || perturbation: %20.10e\n",perturb_node,elem1d->nodes[perturb_node],perturb_sign*perturbation);
        }
        printf("slope: %20.10e \t mannings constant: %20.10e \t mannings n: %20.10e\n",
               slope,mod->manning_units_constant,mod->str_values[string].fterms.manningsn);
        printf("\n--------------------------------------------------------- \n");
        printScreen_debug2_dbl("elem_head", elem_head, nnodes, elem1d->nodes);
        printScreen_debug_svec2d("elem_vel", elem_vel, nnodes, elem1d->nodes);
    }
#endif
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    sarray_init_dbl(elem_rhs, nnodes);
    
    if (string > NORMAL) {
        
        /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         *                                    OUTFLOW CONTRIBUTION
         *--------------------------------------------------------------------------------------------
         * Calculates the outflow addition to the elemental residual. \n
         *
         * \note CJT \:: calculates implicit flow using implicitly calculated water elevation
         *********************************************************************************************/
        if (str_values[string].ol_flow.bc_flag == BCT_OUTFLOW) {
            if (elem_head[0] > 0 && elem_head[1] > 0) {
                integrate_line_phi_h_v_dot_n(djac, dt, elem_head, elem_vel, elem1d->nrml, elem_rhs);
#ifdef DEBUG
                if (DEBUG == ON || DEBUG_LOCAL == ON) {
                    rhs_1dof("1D DIFFUSIVE WAVE || OUTFLOW: ", nnodes, ie, elem1d->nodes, elem_rhs);
                }
#endif
            }
        }
        
        /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         *                                   VELOCITY CONTRIBUTION
         *--------------------------------------------------------------------------------------------
         * Calculates the velocity/discharge addition to the elemental residual. \n
         * Here user prescribes: \f$ q_n \f$ \n
         *
         * \note CJT \:: calculates explicitly perscribed normal velocity/discharge
         *********************************************************************************************/
        else if (str_values[string].ol_flow.bc_flag == BCT_VEL_NEU) {
            int isers = str_values[string].ol_flow.isigma;
            double flux = -1.0 * sseries_get_value(isers, mod->series_head, 0);
            
            if (elem_head[0] > 0 && elem_head[1] > 0) {
                integrate_line_phi(djac, dt * flux, elem_rhs);
#ifdef DEBUG
                if (DEBUG == ON || DEBUG_LOCAL == ON) {
                    rhs_1dof("1D DIFFUSIVE WAVE || NB VEL: ", nnodes, ie, elem1d->nodes, elem_rhs);
                }
#endif
            }
        }
        
        /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         *                             TAILWATER/PRESSURE CONTRIBUTION
         *--------------------------------------------------------------------------------------------
         * Calculates the tailwater or pressure addition to the elemental residual. \n
         *
         * \note CJT \:: calculates implicit flow using user provides water elevation
         *********************************************************************************************/
        else if (str_values[string].ol_flow.bc_flag == BCT_PRS_NEU) {
            double elem_depth[nnodes];
            int isers = str_values[string].ol_flow.isigma;
            double tail_elev = sseries_get_value(isers, mod->series_head,0);
            
            for (i=0; i<nnodes; i++) {
                elem_depth[i] = tail_elev - grid->node[elem1d->nodes[i]].z;
                elem_vel[i].x = mod->manning_units_constant * slope *pow(elem_depth[i], 2./3.)/mod->str_values[string].fterms.manningsn;
                elem_vel[i].y = mod->manning_units_constant * slope *pow(elem_depth[i], 2./3.)/mod->str_values[string].fterms.manningsn;
            }
            if (elem_depth[0] <= 0. || elem_depth[1] <= 0. || elem_head[0] <= 0. || elem_head[1] <= 0.) return;
            
            integrate_line_phi_h_v_dot_n(djac, dt, elem_depth, elem_vel, elem1d->nrml, elem_rhs);
#ifdef DEBUG
            if (DEBUG == ON || DEBUG_LOCAL == ON) {
                rhs_1dof("1D DIFFUSIVE WAVE || TAILWATER: ", nnodes, ie, elem1d->nodes, elem_rhs);
            }
#endif
            
        }
    }
}




