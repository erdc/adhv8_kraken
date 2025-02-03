/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  sw2_body_resid.c This file collections functions responsible for
 *          the 2D shallow water body contributions to the elemental residual.              */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//#include "global_header.h"
#include "adh.h"
// prototypes
void fe_sw2_temporal(int ie, SELEM_2D *elem2d, int nnodes, SVECT *elem_nds, double djac, double drying_lower_limit, double *elem_head, SVECT2D *elem_vel, double wd_factor, double dt_factor, double *elem_rhs, char *string, int DEBUG, int DEBUG_LOCAL, int wd_flag);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 2D shallow water elemental residual.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 * @param[in] mod (SMODEL_SUPER *) - the super model, contains pointer to grid, dependent vars, maps, and such
 * @param[in,out] elem_rhs (double *) - the 2D elemental residual array
 * @param[in]  ie  (int) - the elemental id
 * @param[in]  pertubation (double) - the F-D approximation size aka Newton pertubation
 * @param[in]  perturb_node (int) -  the index of node to be perturbed (local to element)
 * @param[in]  perturb_var (int)  - the index of the variable to be perturbed
 * @param[in]  perturb_sign (int)  - the direction of Newton perturbation (either -1 or 1)
 * @param[in]  DEBUG    (int) - a debug flag
 * \returns a possible error code
 *
 *  \details Solves the body integals of the following weak, discrete body terms of the 2D shallow water equation: \n
 *  \f{eqnarray*}{
 *  \weakSWDAcont{e}{i}{h} \\
 *  \weakSWMxDD{e}{i}{h} \\
 *  \weakSWMyDD{e}{i}{h}
 *  \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int fe_sw2_body_resid(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    //local debug flags
    int DEBUG_LOCAL = OFF;
    int debug_node = UNSET_INT;
    int debug_pe = 0;
    int DEBUG_PICKETS = OFF;
   
    //convenitent if looking at only a specific node
    if (mod->grid->elem2d[ie].nodes[0] == debug_node-1 ||
        mod->grid->elem2d[ie].nodes[1] == debug_node-1 ||
        mod->grid->elem2d[ie].nodes[2] == debug_node-1) {
#ifdef _MESSG
        int ierr_code, myid=0;
        ierr_code = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        if (myid == debug_pe) {
#endif
      DEBUG = ON;
#ifdef _MESSG
        }
#endif
    }
#ifdef _DEBUG
    if (DEBUG == ON) DEBUG_LOCAL = ON;
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    time_t time1;  time(&time1);
#endif
    //local variables
    int i;
    int PRESSURE_FLAG = ON;
    double t1, t2, t3, t4, t5;   
    // aliases for convenience
    SSW *sw2 = mod->sw;
    SGRID *grid = mod->grid;
    SELEM_2D *elem2d = &(grid->elem2d[ie]); // be careful not to shallow copy here
    SVECT2D *grad_shp = elem2d->grad_shp;
    int nnodes = elem2d->nnodes;
    int ndof = nnodes*3; //this is a 3 dof problem
    double alpha = mod->tau_temporal; //where to put this? this is specific to choice of time stepping scheme
    double dt = (*mod->dt)/(mod->nsubsteps); //maybe save this instead of computing every time but it can change
    double djac = elem2d->djac;
    double g = mod->gravity; //where to put this
    double drying_lower_limit = sw2->drying_lower_limit;
    int imat = elem2d->mat; // this will need to be mod->coverage->coverage_2d[SW2_INDEX][ie]
    STR_VALUE str_values = mod->str_values[elem2d->string]; //still in use?
    int wd_flag = sw2->dvar.elem_flags[sw2->WD_FLAG][ie]; //WARNING, ie is not always correct, could be sw2->dvar.dvar_elem_map[ie] need to use Coreys map eventually? how to avoid conditional. a function pointer?
    //only make for quadrilateral elements?
    SQUAD *quad = grid->quad_rect; // for quadrilateral quadrature calculations
    int isTriangle = NO;  if (nnodes == NDONTRI) isTriangle = YES;
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        if (isTriangle == TRUE) assert(djac>SMALL);
        assert(alpha > -1E-6 && alpha < 1.0000001 );
        assert(imat > -1);
        assert(perturb_var == PERTURB_NONE || perturb_var == PERTURB_H || perturb_var == PERTURB_U || perturb_var == PERTURB_V);
    }
#endif    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GRID VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/    
    double node_z[nnodes];
    SNODE nodes[nnodes];
    for (i=0; i<nnodes; i++) {
        snode_copy(&(nodes[i]), grid->node[elem2d->nodes[i]]);
        node_z[i] = nodes[i].z;
    }
    SVECT elem_nds[nnodes];
    snode2svect(nodes, elem_nds, nnodes);
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 
    double elem_head[nnodes]; SVECT2D elem_vel[nnodes];
    //for now use cg map but use Corey's index maps later on
    global_to_local_dbl_cg(elem_head, mod->sol, elem2d->nodes, nnodes, PERTURB_H, mod->dof_map_local, mod->node_physics_mat);
    global_to_local_SVECT2D_cg(elem_vel, mod->sol, elem2d->nodes, nnodes, PERTURB_U, PERTURB_V, mod->dof_map_local, mod->node_physics_mat);
    if (perturb_var == PERTURB_H) {
        elem_head[perturb_node] += perturb_sign * perturbation;
    } else if (perturb_var == PERTURB_U) {
        PRESSURE_FLAG = OFF;
        elem_vel[perturb_node].x += perturb_sign * perturbation;
    } else if (perturb_var == PERTURB_V) {
        PRESSURE_FLAG = OFF;
        elem_vel[perturb_node].y += perturb_sign * perturbation;
    }
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    double elem_old_head[nnodes], elem_older_head[nnodes];
    global_to_local_dbl_cg(elem_old_head, mod->sol_old, elem2d->nodes, nnodes, PERTURB_H, mod->dof_map_local, mod->node_physics_mat);
    global_to_local_dbl_cg(elem_older_head, mod->sol_older, elem2d->nodes, nnodes, PERTURB_H, mod->dof_map_local, mod->node_physics_mat);  
    SVECT2D elem_old_vel[nnodes], elem_older_vel[nnodes];
    global_to_local_SVECT2D_cg(elem_old_vel, mod->sol_old, elem2d->nodes, nnodes, PERTURB_U, PERTURB_V, mod->dof_map_local, mod->node_physics_mat);
    global_to_local_SVECT2D_cg(elem_older_vel, mod->sol_older, elem2d->nodes, nnodes, PERTURB_U, PERTURB_V, mod->dof_map_local, mod->node_physics_mat);  
    double elem_density[nnodes]; sarray_init_value_dbl(elem_density,nnodes,0.);
    if (FALSE){//(mod->flag.BAROCLINIC == 1) { //where to put this?
        //global_to_local_dbl(sw2->nd, elem_density, elem2d->nodes, nnodes);
        //use array mapping
        global_to_local_dbl_cg_arr_map(elem_density, elem2d->nodes, nnodes, sw2->dvar.dvar_node_map, sw2->dvar.nodal_dvar[sw2->DENSITY]);
        for (i=0; i<nnodes; i++) {
            if (elem_head[i] < 0) elem_density[i] = 0;
           }
    }  
    // local external pressure (atmospheric, lids, etc) :: CJT :: should we use local density here?
    double elem_pressure[nnodes]; sarray_init_value_dbl(elem_pressure,nnodes,0.); // cjt :: used in SUPG only right now
    if (PRESSURE_FLAG == ON) {
        int id = UNSET_INT;
        for (i = 0; i < nnodes; i++) {
            if (nodes[i].string > NORMAL) { //bug?? switched to nodes alias not gri
                if (mod->str_values[nodes[i].string].ol_flow.bc_flag == BCT_LID_DFT) {
                    id = mod->str_values[nodes[i].string].ol_flow.iu_0;
                    t1 = sseries_get_value(id, mod->series_head,0);
                    elem_pressure[i] = mod->density * g * t1;
                } else if (mod->str_values[nodes[i].string].ol_flow.bc_flag == BCT_LID_ELV) {
                    id = mod->str_values[nodes[i].string].ol_flow.iu_0;
                    t1 = sseries_get_value(id, mod->series_head,0);
                    t1 -= elem_nds[i].z;
                    if (elem_head[i] > t1) elem_pressure[i] = 2. * mod->density * g * (elem_head[i] - t1);
                } else if (mod->str_values[nodes[i].string].ol_flow.bc_flag == BCT_LID_DEP) {
                    id = mod->str_values[nodes[i].string].ol_flow.iu_0;
                    t1 = sseries_get_value(id, mod->series_head,0);
                if (elem_head[i] > t1) elem_pressure[i] = 2. * mod->density * g * (elem_head[i] - t1);
                }
            }
        }
    }
    // local water elevation
    double elem_wse[nnodes]; sarray_add_dbl(elem_wse, elem_head, node_z, nnodes);
    // velocity components for convienence
    double elem_u[nnodes], elem_v[nnodes];
    dumpVector2D(elem_vel, nnodes, elem_u, elem_v);
    double elem_old_u[nnodes], elem_old_v[nnodes];
    dumpVector2D(elem_old_vel, nnodes, elem_old_u, elem_old_v);   
    // non-wet/dry averages, area and function/vector gradients
    double area = 0., elem_h_avg = 0., elem_h_old_avg = 0., elem_density_avg = 0.;
    double elem_u_avg = 0., elem_v_avg = 0., elem_u_old_avg = 0., elem_v_old_avg = 0.;
    SVECT2D elem_grad_z, elem_grad_head, elem_grad_u, elem_grad_v, elem_grad_wse, elem_grad_prs, elem_grad_density;
    svect2d_init(&elem_grad_z); svect2d_init(&elem_grad_head); svect2d_init(&elem_grad_u); svect2d_init(&elem_grad_v);
    svect2d_init(&elem_grad_wse); svect2d_init(&elem_grad_prs); svect2d_init(&elem_grad_density);
    if (isTriangle == TRUE) {
        area = djac;
        // elemental averages - djacs cancel here
        elem_h_avg =       integrate_triangle_f(1.,1.,elem_head);
        elem_h_old_avg =   integrate_triangle_f(1.,1.,elem_old_head);
        elem_u_avg =       integrate_triangle_f(1.,1.,elem_u);
        elem_v_avg =       integrate_triangle_f(1.,1.,elem_v);
        elem_u_old_avg =   integrate_triangle_f(1.,1.,elem_old_u);
        elem_v_old_avg =   integrate_triangle_f(1.,1.,elem_old_v);
        elem_density_avg = integrate_triangle_f(1.,1.,elem_density);
        // function gradients
        grad2d_phi_f(grad_shp, node_z, &elem_grad_z, nnodes);
        grad2d_phi_f(grad_shp, elem_head, &elem_grad_head, nnodes);
        grad2d_phi_f(grad_shp, elem_pressure, &elem_grad_prs, nnodes);
        grad2d_phi_f(grad_shp, elem_density, &elem_grad_density, nnodes);
        grad2d_phi_dot_v(grad_shp, elem_vel, &elem_grad_u, &elem_grad_v, nnodes);
        elem_grad_wse = svect2d_add(elem_grad_z, elem_grad_head);       
    } else {
        area = integrate_quadrilateral_area(elem_nds, 1.);        
        // elemental averages
        elem_h_avg =       integrate_quadrilateral_f(elem_nds,1./area,elem_head);
        elem_h_old_avg =   integrate_quadrilateral_f(elem_nds,1./area,elem_old_head);
        elem_u_avg =       integrate_quadrilateral_f(elem_nds,1./area,elem_u);
        elem_v_avg =       integrate_quadrilateral_f(elem_nds,1./area,elem_v);
        elem_u_old_avg =   integrate_quadrilateral_f(elem_nds,1./area,elem_old_u);
        elem_v_old_avg =   integrate_quadrilateral_f(elem_nds,1./area,elem_old_v);
        elem_density_avg = integrate_quadrilateral_f(elem_nds,1./area,elem_density);
        elem_grad_u =  integrate_quadrilateral_gradF(elem_nds, 1./area, elem_u); // assume grad_u = grad_u_avg here
        elem_grad_v =  integrate_quadrilateral_gradF(elem_nds, 1./area, elem_v); // assume grad_v = grad_v_avg here
    }
    SVECT2D elem_vel_avg, elem_old_vel_avg;
    elem_vel_avg.x = elem_u_avg; elem_old_vel_avg.x = elem_u_old_avg;
    elem_vel_avg.y = elem_v_avg; elem_old_vel_avg.y = elem_v_old_avg;
    double elem_vel_avg_mag = svect2d_mag(elem_vel_avg);
    double elem_old_vel_avg_mag = svect2d_mag(elem_old_vel_avg);    
    // non-wetting and drying unit normals
    double vst_i = elem_u_avg / MAX(elem_vel_avg_mag, 1.E-8);
    double vst_j = elem_v_avg / MAX(elem_vel_avg_mag, 1.E-8);    
    // wet/dry factors and elementally averages
    double factor = 1., factor_old = 1., factor_older = 1.;
    double elem_h_wd_avg = elem_h_avg;
    double elem_h_old_wd_avg = elem_h_old_avg;
    SVECT2D elem_vel_wd_avg, elem_old_vel_wd_avg;
    elem_vel_wd_avg.x = elem_u_avg, elem_old_vel_wd_avg.x = elem_u_old_avg;
    elem_vel_wd_avg.y = elem_v_avg, elem_old_vel_wd_avg.y = elem_v_old_avg;
    if (isTriangle == YES) {      
        fe_sw2_wd_average(elem_nds, elem_head, elem_vel, elem_head, djac, &elem_h_wd_avg, &elem_vel_wd_avg);
        fe_sw2_wd_average(elem_nds, elem_old_head, elem_old_vel, elem_old_head, djac, &elem_h_old_wd_avg, &elem_old_vel_wd_avg); 
        factor =       fe_sw2_wet_dry_factor(elem_nds, elem_head,       djac);
        factor_old =   fe_sw2_wet_dry_factor(elem_nds, elem_old_head,   djac);
        factor_older = fe_sw2_wet_dry_factor(elem_nds, elem_older_head, djac);
    }
    double elem_vel_wd_avg_mag     = svect2d_mag(elem_vel_wd_avg);
    double elem_old_vel_wd_avg_mag = svect2d_mag(elem_old_vel_wd_avg);  
    // friction variables
    double h_fric = MAX(elem_h_wd_avg, SMALL);
    double roughness = 0.03;//fe_sw2_get_roughness(mod, h_fric, elem_vel_wd_avg_mag, elem2d->string, UNUSED);
    // loads the diffusion eddy tensor
    double ev_st = 0., ev_tr = 0.;
    STENSOR2DAI ev; stensor2dai_init(&ev);
    //TO DO:
    //hard set to 0 for now
    ev.xx = 0.0;
    ev.yy = 0.0;
    ev.xy = 0.0;
    ev.yx = 0.0;
    ev_st = fabs(ev.xx - ev.yy);
    ev_tr = MIN(ev.xx, ev.yy);
//    if (mod->mat[imat].sw->EVSF == YES) {  // user supplied eddy viscosity, should be in coverages now?
//        ev.xx = mod->mat[imat].sw->ev.xx;
//        ev.yy = mod->mat[imat].sw->ev.yy;
//        ev.xy = mod->mat[imat].sw->ev.xy;
//        ev.yx = mod->mat[imat].sw->ev.xy;
//        ev.xx += mod->viscosity;
//        ev.yy += mod->viscosity;
//        ev.xy += mod->viscosity;
//        ev.yx += mod->viscosity;    
//        ev_st = fabs(ev.xx - ev.yy);
//        ev_tr = MIN(ev.xx, ev.yy);
//    } else if (mod->mat[imat].sw->EEVF == YES) {  // use estimated eddy viscosity
//        if (isTriangle == TRUE) {
//            fe_sw2_get_EEVF(mod->mat[imat].sw->eev_mode, mod->mat[imat].sw->eev_coef, g, drying_lower_limit, h_fric, roughness, area, elem_grad_u, elem_grad_v, elem_vel_avg_mag, wd_flag, &ev_st, &ev_tr);
//        } else {
//            ///Users/rditlcjt/Dropbox/Science/adh/repo_gitlab/src/adh/fe/sw2/fe_sw2_body_resid.c use elementally averaged gradients here
//            SVECT2D elem_grad_u_avg = integrate_quadrilateral_gradF(elem_nds, 1./area, elem_u);
//            SVECT2D elem_grad_v_avg = integrate_quadrilateral_gradF(elem_nds, 1./area, elem_v);      
//            fe_sw2_get_EEVF(mod->mat[imat].sw->eev_mode, mod->mat[imat].sw->eev_coef, g, drying_lower_limit, h_fric, roughness, area, elem_grad_u_avg, elem_grad_v_avg, elem_vel_avg_mag, wd_flag, &ev_st, &ev_tr);
//        }
//        ev.xx = ev_st;
//        ev.xy = ev_tr;
//        ev.yx = ev_tr;
//        ev.yy = ev_st;
//    }    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        printf("\nSHALLOW WATER 2D ELEM RESID :: ie: %d \t dt: %20.10f \t area: %20.10f \t wd_flag: %d \t djac: %30.20e : %30.20e %30.20e %30.20e",ie,dt,area,wd_flag,djac,mod->grid->node[elem2d->nodes[2]].x,mod->grid->node[elem2d->nodes[2]].y,mod->grid->node[elem2d->nodes[2]].z);
        if (perturb_var == PERTURB_H) {
            printf("\t perturbing H  || node: %d || perturbation: %30.20e\n",nodes[perturb_node].id,perturb_sign*perturbation);
        } else if (perturb_var == PERTURB_U) {
            printf("\t perturbing U  || node: %d || perturbation: %30.20e\n",nodes[perturb_node].id,perturb_sign*perturbation);
        } else if (perturb_var == PERTURB_V) {
            printf("\t perturbing V  || node: %d || perturbation: %30.20e\n",nodes[perturb_node].id,perturb_sign*perturbation);
        }
        //selem2d_printScreen(elem2d);        
        printf("\n------------------ ELEMENTAL AVERAGES --------------------\n");
        printf("elem_h_avg: %30.20e elem_h_wd_avg: %30.20e \n",elem_h_avg,elem_h_wd_avg);
        printf("elem_h_old_avg: %20.10e  elem_h_wd_avg: %20.10e \n",elem_h_old_avg ,elem_h_old_wd_avg);
        printf("elem_u_avg: %20.10e  elem_u_wd_avg: %20.10e \n",elem_vel_avg.x,elem_vel_wd_avg.x);
        printf("elem_v_avg: %20.10e  elem_v_wd_avg: %20.10e \n",elem_vel_avg.y,elem_vel_wd_avg.y);
        printf("elem_u_old_avg: %20.10e  elem_u_old_wd_avg: %20.10e \n",elem_old_vel_avg.x,elem_old_vel_wd_avg.x);
        printf("elem_v_old_avg: %20.10e  elem_v_old_wd_avg: %20.10e \n",elem_old_vel_avg.y,elem_old_vel_wd_avg.y);
        printf("elem_density_avg: %20.10e\n",elem_density_avg);
        printf("elem_vel_avg_mag: %20.10e \t elem_vel_wd_avg_mag: %20.10e \n",elem_vel_avg_mag,elem_vel_wd_avg_mag);
        printf("elem_old_vel_avg_mag: %20.10e \t elem_old_vel_wd_avg_mag: %20.10e \n",elem_old_vel_avg_mag,elem_old_vel_wd_avg_mag);        
        printf("\n------------------ FUNCTION GRADIENTS --------------------\n");
        printf("elem_grad_z: {%20.10e, %20.10e}\n",elem_grad_z.x,elem_grad_z.y);
        printf("elem_grad_head: {%20.10e, %20.10e}\n",elem_grad_head.x,elem_grad_head.y);
        printf("elem_grad_prs: {%20.10e, %20.10e}\n",elem_grad_prs.x,elem_grad_prs.y);
        printf("elem_grad_density: {%20.10e, %20.10e}\n",elem_grad_density.x,elem_grad_density.y);
        printf("elem_grad_vel.x: {%20.10e, %20.10e}\n",elem_grad_u.x,elem_grad_u.y);
        printf("elem_grad_vel.y: {%20.10e, %20.10e}\n",elem_grad_v.x,elem_grad_v.y);
        printf("elem_grad_wse: {%20.10e, %20.10e}\n",elem_grad_wse.x,elem_grad_wse.y);        
        printf("\n------------------- LINEAR VARIABLES ---------------------\n");
        printScreen_debug2_dbl("elem_head", elem_head, nnodes, elem2d->nodes);
        printScreen_debug2_dbl("elem_old_head", elem_old_head, nnodes, elem2d->nodes);
        printScreen_debug2_dbl("elem_older_head", elem_older_head, nnodes, elem2d->nodes);
        printScreen_debug_svec2d("elem_vel", elem_vel, nnodes, elem2d->nodes);
        printScreen_debug_svec2d("elem_old_vel", elem_old_vel, nnodes, elem2d->nodes);
        printScreen_debug_svec2d("elem_older_vel", elem_older_vel, nnodes, elem2d->nodes);       
        printf("\n--------------------- OTHER STUFF ------------------------\n");
        printf("roughness: %20.10e\n",roughness);
        printf("wet-dry factors: current: %30.20e \t old: %30.20e \t older: %30.20e\n",factor,factor_old,factor_older);
        printf("eddy viscosity tensor: "); stensor2dai_printScreen(ev);
    }
#endif  
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                  FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    //dof3_init_array(elem_rhs,nnodes); //now a double array and we do not need to 0!!!
    //keep temp rhs for convenience and debugging
    double rhs[ndof];
    sarray_init_dbl(rhs,ndof);
    //adopting the convention:
    //continuity
    //momentum x
    //momentum y
    double rhs_c_eq[nnodes], rhs_x_eq[nnodes], rhs_y_eq[nnodes];
    double vars[2]; // for passing doubles through wet-dry routine  
    // this is used to store later for transport
    //NEED TO REINCORPORATE
    for (i=0; i<nnodes; i++) {
        //should be
        sw2->elem_rhs_dacont_extra_terms[ie][i] = 0.;
    }    
   /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                 SHOCK CAPTURING CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates residual contributions from momentum and mass shock capturing. \n
     * \note CJT \:: wd_flag = 0 (element fully wet), 1 (element partially dry), 2 (element full dry)
     * \note CJT \:: shock capturing terms only contribute for partially or fully dry elements
     * \note CJT \:: formulation?
     *********************************************************************************************/    
    if (isTriangle == TRUE) {
        if (wd_flag == 1 || wd_flag == 2) {   
                    
            double elem_wse_grad_mag = svect2d_mag_safe(elem_grad_wse); //warning, this adds an epsilon to avoid 0
            double elem_u_grad_mag   = svect2d_mag_safe(elem_grad_u);
            double elem_v_grad_mag   = svect2d_mag_safe(elem_grad_v);            
            double sh_cap_coef_u = one_2 * sw2->drying_upper_limit * djac * elem_u_grad_mag;
            double sh_cap_coef_v = one_2 * sw2->drying_upper_limit * djac * elem_v_grad_mag;            
            double shock_avg_depth = sw2->drying_upper_limit;
            if (shock_avg_depth < sw2->drying_upper_limit) shock_avg_depth = sw2->drying_upper_limit;
            double temp = (elem_vel_avg_mag / sqrt(g * sw2->drying_upper_limit)); 
            if (temp > 1) temp = 1.;            
            double sh_cap_coef = one_2 * elem_vel_avg_mag * djac * sw2->drying_upper_limit * elem_wse_grad_mag / (shock_avg_depth * shock_avg_depth + SMALL);            
            SVECT2D grad_unorm, grad_vnorm, norm;
            grad_unorm.x = elem_grad_u.x / elem_u_grad_mag;  grad_unorm.y = elem_grad_u.y / elem_u_grad_mag;
            grad_vnorm.x = elem_grad_v.x / elem_v_grad_mag;  grad_vnorm.y = elem_grad_v.y / elem_v_grad_mag;            
            for (i=0; i<nnodes; i++) {
                // mass contributions, every 3rd equation is continuity
                rhs[i*3] = dt * one_2 * g * djac * ((sw2->drying_lower_limit * elem_grad_wse.x) * grad_shp[i].x + (sw2->drying_lower_limit * elem_grad_wse.y) * grad_shp[i].y);        
                // momentum contribution
                norm.x = svect2d_dotp(grad_unorm, elem_grad_u);
                norm.y = svect2d_dotp(grad_vnorm, elem_grad_v);
                //x momentum
                rhs[i*3+1] = sh_cap_coef_u * dt * djac * norm.x * svect2d_dotp(grad_unorm, grad_shp[i]);
                //y momentum
                rhs[i*3+2] = sh_cap_coef_v * dt * djac * norm.y * svect2d_dotp(grad_vnorm, grad_shp[i]);
            }            
            for (i=0; i<nnodes; i++) {
                //NEED TO FIX!!!
                sw2->elem_rhs_dacont_extra_terms[ie][i] += rhs[i*3]; // cjt :: store
                elem_rhs[i*3] += rhs[i*3];
                elem_rhs[i*3+1] += rhs[i*3+1];
                elem_rhs[i*3+2] += rhs[i*3+2];
            }            
#ifdef _DEBUG
            if (DEBUG == ON || DEBUG_LOCAL == ON) {
                printScreen_rhs_3dof("2D SW || SHOCK CAPTURING",nnodes, ie, elem2d->nodes, elem_rhs);
                Is_DoubleArray_Inf_or_NaN(elem_rhs, nnodes*3 ,__FILE__ ,__LINE__);
                printf("grad_wse: %20.10e \t %20.10e\n",elem_grad_wse.x,elem_grad_wse.y);
                printf("grad_u: %20.10e \t %20.10e\n",elem_grad_u.x,elem_grad_u.y);
                printf("grad_v: %20.10e \t %20.10e\n",elem_grad_v.x,elem_grad_v.y);
                printf("norm.x: %20.10e \t norm.y: %20.10e \t sh_cap_coef_u: %20.10e \t sh_cap_coef_v: %20.10e \n",norm.x,norm.y,sh_cap_coef_u,sh_cap_coef_v);
            }
#endif
        }

    }    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    ADVECTION CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the advection addition to the SW body terms. \n
     *  \f{eqnarray*}{
     *   \residDA{i}{}{c}   &=& -  dt * \bodyConv{\,2d}{e}{\phidd{i}}{(\velb{h} \, \depth{h})} \\
     *   \residDA{i}{}{mx}  &=& -  dt * \bodyConv{\,2d}{e}{\phidd{i}}{(\velb{h} \, \ub{h}\depth{h})} \\
     *   \residDA{i}{}{my}  &=& -  dt * \bodyConv{\,2d}{e}{\phidd{i}}{(\velb{h} \, \vb{h}\depth{h})}
     *  \f}
     * \note GS  \:: Wrap in wet/dry integration
     * \note CJT \:: Until we include quadrilaterals in wetting and drying, this is only for triangles
     *********************************************************************************************/    
    sarray_init_dbl(rhs,ndof);
    if (isTriangle == YES) {
        SVECT2D elem_vel_wd[nnodes]; svect2d_copy_array(elem_vel_wd, elem_vel, nnodes);
        vars[0] = dt;
        vars[1] = -1;
        t1 = fe_sw2_wet_dry_wrapper(rhs, elem_nds, elem_head, grad_shp, elem_vel, elem_vel_wd, NULL, NULL, djac, ON, DEBUG, vars, fe_sw2_wd_convection_triangle);
    } else {
        sarray_init_dbl(rhs_c_eq,nnodes);
        sarray_init_dbl(rhs_x_eq,nnodes);
        sarray_init_dbl(rhs_y_eq,nnodes);
        integrate_quadrilateral_gradPhi_dot_f_v(elem_nds, -dt, elem_head, elem_vel, rhs_c_eq);
        integrate_quadrilateral_gradPhi_dot_f_g_v(elem_nds, -dt, elem_head, elem_u, elem_vel, rhs_x_eq);
        integrate_quadrilateral_gradPhi_dot_f_g_v(elem_nds, -dt, elem_head, elem_v, elem_vel, rhs_y_eq);
        for (i=0; i<nnodes; i++) {
            rhs[i*3] = rhs_c_eq[i];
            rhs[i*3+1] = rhs_x_eq[i];
            rhs[i*3+2] = rhs_y_eq[i];
        }
    }
    for (i=0; i<ndof; i++) {
        elem_rhs[i] += rhs[i];
    }
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        printScreen_rhs_3dof("2D SW || ADVECTION",nnodes, ie, elem2d->nodes, rhs);
    }
#endif 
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    TEMPORAL CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates the temporal addition to the elemental residual. \n
     *  \f{eqnarray*}{
     *   \residDA{i}{}{c}   &=& dt * \bodyTime{\,2d}{e}{\phidd{i}}{\depth{h}} \\
     *   \residDA{i}{}{mx}  &=& dt * \bodyTime{\,2d}{e}{\phidd{i}}{\,\ub{h}\depth{h}} \\
     *   \residDA{i}{}{my}  &=& dt * \bodyTime{\,2d}{e}{\phidd{i}}{\,\vb{h}\depth{h}}
     *  \f}
     * \note GLB \:: consistent mass integration of continuity so only wet portion is itegrated, momentum temporal is lumped
     * \note CJT \:: consistent mass temporal integral:
     * \f$ \bodyTime{\,2d}{e}{\phidd{i}}{\depth{h}_1 \, \phidd{1} + \depth{h}_2 \, \phidd{2} + \depth{h}_3 \, \phidd{3}} \f$
     * \note CJT \:: lumped & grouped momentum temporal integral:
     * \f$ \bodyTime{\,2d}{e}{\phidd{i}}{\depth{h}_i \, \ub{h}_i} \f$
     *********************************************************************************************/    
    // ++ t(n+1) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_sw2_temporal(ie, elem2d, nnodes, elem_nds, djac, drying_lower_limit, elem_head, elem_vel, factor,
                    (1. + alpha / 2.), elem_rhs, "2D SW || TEMPORAL T(N+1)", DEBUG, DEBUG_LOCAL, UNSET_INT);    
    // ++ t(n) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_sw2_temporal(ie, elem2d, nnodes, elem_nds, djac, drying_lower_limit, elem_old_head, elem_old_vel, factor_old,
                    (-1.) * (1. + alpha), elem_rhs, "2D SW || TEMPORAL T(N)", DEBUG, DEBUG_LOCAL, wd_flag);    
    // ++ t(n-1) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_sw2_temporal(ie, elem2d, nnodes, elem_nds, djac, drying_lower_limit, elem_older_head, elem_older_vel, factor_older,
                    (alpha / 2.), elem_rhs, "2D SW || TEMPORAL T(N-1)", DEBUG, DEBUG_LOCAL, wd_flag);   
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    DIFFUSION CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates the diffusion addition to the elemental residual. \n
     *  \f{eqnarray*}{
     *   \residDA{i}{}{mx}  &=& dt * \bodyDiffusion{\,2d}{e}{\phidd{i}}{\depth{h}\,\diffTensorDA{x}{h}{\diffTensorSW}{}}  \\
     *   \residDA{i}{}{my}  &=& dt * \bodyDiffusion{\,2d}{e}{\phidd{i}}{\depth{h}\,\diffTensorDA{y}{h}{\diffTensorSW}{}}
     *  \f}
     * \note CJT \:: currently uses elementally averaged depth
     *
     *********************************************************************************************/    
    // get elementally averaged depth
    double temp = elem_h_wd_avg;
    if (elem_h_wd_avg < sw2->drying_upper_limit) {
        temp = sw2->drying_upper_limit;
    }    
    sarray_init_dbl(rhs_x_eq,nnodes);
    sarray_init_dbl(rhs_y_eq,nnodes);    
    if (isTriangle == TRUE) {        
        // sigma_xx = 2 ev.xx du_dx :: sigma_xy = ev.xy (du_dy + dv_dx) :: sigma_yx = sigma_xy :: sigma_yy = 2 ev.yy dv_dy
        SVECT2D diffusive_flux_x, diffusive_flux_y;
        diffusive_flux_x.x = 2. * ev_st * vst_i * (vst_i * elem_grad_u.x + vst_j * elem_grad_u.y) + 2. * ev_tr * elem_grad_u.x;
        diffusive_flux_x.y = ev_st * vst_j * (vst_j * elem_grad_u.y + vst_i * elem_grad_u.x + vst_i * elem_grad_v.x + vst_j * elem_grad_v.y) + ev_tr * (elem_grad_u.y + elem_grad_v.x);
        diffusive_flux_y.x = ev_st * vst_i * (vst_j * elem_grad_u.y + vst_i * elem_grad_u.x + vst_i * elem_grad_v.x + vst_j * elem_grad_v.y) + ev_tr * (elem_grad_u.y + elem_grad_v.x);
        diffusive_flux_y.y = 2. * ev_st * vst_j * (vst_j * elem_grad_v.y + vst_i * elem_grad_v.x) + 2. * ev_tr * elem_grad_v.y;        
        integrate_triangle_gradPhi_dot_vbar(grad_shp, djac, temp, diffusive_flux_x, rhs_x_eq);
        integrate_triangle_gradPhi_dot_vbar(grad_shp, djac, temp, diffusive_flux_y, rhs_y_eq);        
    } else {        
        integrate_quadrilateral_gradPhi_dot_Dv(elem_nds, quad, elem_h_wd_avg, ev, elem_vel, rhs_x_eq, rhs_y_eq);
    }    
    t1 = dt * factor;
    for (i=0; i<nnodes; i++) {
        elem_rhs[i*3+1]+= t1 * rhs_x_eq[i];
        elem_rhs[i*3+2]+= t1 * rhs_y_eq[i];
    }
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        for (i=0; i<nnodes; i++) {
            rhs[i*3] = 0.;
            rhs[i*3+1]= t1 * rhs_x_eq[i];
            rhs[i*3+2]= t1 * rhs_y_eq[i];
        }
        printScreen_rhs_3dof("2D SW || DIFFUSION",nnodes, ie, elem2d->nodes, rhs);
    }
#endif    
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    FRICTION CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates the friction addition to the elemental residual. \n
     *  \f{eqnarray*}{
     *   \residDA{i}{}{mx}  &=& dt * \bodyFriction{\,2d}{e}{\phidd{i}}{\ub{h}}{\velb{h}}   \\
     *   \residDA{i}{}{my}  &=& dt * \bodyFriction{\,2d}{e}{\phidd{i}}{\vb{h}}{\velb{h}}
     *  \f}
     *
     *********************************************************************************************/    
    double resistance = roughness / (2.0 * g);
    double elem_friction_x[nnodes], elem_friction_y[nnodes];
    for (i=0; i<nnodes; i++) {
        t1 = svect2d_mag(elem_vel[i]);
        elem_friction_x[i] = elem_vel[i].x * t1;
        elem_friction_y[i] = elem_vel[i].y * t1;
    }
    sarray_init_dbl(rhs_x_eq, nnodes); sarray_init_dbl(rhs_y_eq, nnodes);
    t1 = dt * g * resistance * factor;
    if (isTriangle == TRUE) {
        integrate_triangle_phi_f(djac, t1, elem_friction_x, rhs_x_eq); //NOTE: NOT CONSISTENT INTEGRATION
        integrate_triangle_phi_f(djac, t1, elem_friction_y, rhs_y_eq);
    } else {
        integrate_quadrilateral_phi_f(elem_nds, t1, elem_friction_x, rhs_x_eq);
        integrate_quadrilateral_phi_f(elem_nds, t1, elem_friction_y, rhs_y_eq);
    }    
    for (i=0; i<nnodes; i++) {
        elem_rhs[i*3+1] += rhs_x_eq[i];
        elem_rhs[i*3+2] += rhs_y_eq[i];
    }
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        for (i=0; i<nnodes; i++) {
            rhs[i*3] = 0.;
            rhs[i*3+1] = rhs_x_eq[i];
            rhs[i*3+2] = rhs_y_eq[i];
        }
        printScreen_rhs_3dof("2D SW || FRICTION",nnodes, ie, elem2d->nodes, rhs);
        printf("roughness: %20.10f \t resistance: %20.10f t1: %20.10f \n",roughness, resistance,t1);
    }
#endif  
    // calculate bed stress for PG
    SVECT2D elem_bed_stress;
    elem_bed_stress.x = resistance * elem_vel_wd_avg.x * elem_vel_wd_avg_mag;
    elem_bed_stress.y = resistance * elem_vel_wd_avg.y * elem_vel_wd_avg_mag;
        
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    VORTICITY CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates the vorticity addition to the elemental residual. \n
     *  \f{eqnarray*}{
     *   \residDA{i}{}{mx}  &=& dt * \bodySource{\,2d}{e}{\phidd{i}}{stream_x}  \\
     *   \residDA{i}{}{my}  &=& dt * \bodySource{\,2d}{e}{\phidd{i}}{stream_y}
     *  \f}
     *
     *  \note CJT \:: currently only for triangles, but can use elementally averaged gradients and include quads
     *********************************************************************************************/    
    //Ask corey about vorticity_id
//    if (mod->vorticity_id > UNSET_INT && nnodes == NDONTRI) {
//        SVECT2D stream_accel_vect; svect2d_init(&stream_accel_vect);
//        double elem_c[nnodes];
//        //need to change this to get the correct constituent from transport
//        for (i=0; i<nnodes; i++) elem_c[i] =  mod->con[mod->vorticity_id].concentration[elem2d->nodes[i]];        
//        // 0 is the routine flag, 0 indicates that hydro is calling
//        stream_accel_vect = tl_bendway_correction(grad_shp, elem_old_vel, elem_old_head, elem_c, mod->con[mod->vorticity_id].property[1],
//                                                  mod->con[mod->vorticity_id].property[2], sw2->drying_lower_limit, roughness, mod->density, 0, ie);      
//        for (i=0; i<nnodes; i++) {
//            elem_rhs[i*3+1] += one_3 * dt * stream_accel_vect.x * djac;
//            elem_rhs[i*3+2] += one_3 * dt * stream_accel_vect.y * djac;
//        }
//#ifdef _DEBUG
//        if (DEBUG == ON || DEBUG_LOCAL == ON) {
//            printScreen_rhs_3dof("2D SW || VORTICITY",nnodes, ie, elem2d->nodes, elem_rhs);
//        }
//#endif
//    }    
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    CORIOLIS CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates the coriolis addition to the elemental residual. \n
     *  \f{eqnarray*}{
     *   \residDA{i}{}{mx}  &=& - dt * \bodyCoriolis{\,2d}{e}{\phidd{i}}{\depth{h}\,\vb{h}}   \\
     *   \residDA{i}{}{my}  &=& + dt * \bodyCoriolis{\,2d}{e}{\phidd{i}}{\depth{h}\,\ub{h}}
     *  \f}
     *
     *********************************************************************************************/    
    sarray_init_dbl(rhs_x_eq, nnodes); sarray_init_dbl(rhs_y_eq, nnodes);
    //double coriolis_speed = get_coriolis_angular_speed(mod->mat[imat].sw->coriolis); //need to fix
    double coriolis_speed = 0.0;
    t1 = dt * factor * coriolis_speed;
    if (isTriangle == TRUE) {
        integrate_triangle_phi_f_g(djac, -t1, elem_head, elem_v, rhs_x_eq); //no wd wrapper?
        integrate_triangle_phi_f_g(djac, +t1, elem_head, elem_u, rhs_y_eq);
    } else {
        integrate_quadrilateral_phi_f_g(elem_nds, -t1, elem_head, elem_v, rhs_x_eq);
        integrate_quadrilateral_phi_f_g(elem_nds, +t1, elem_head, elem_u, rhs_y_eq);
    }
    for (i=0; i<nnodes; i++) {
        elem_rhs[i*3+1] += rhs_x_eq[i];
        elem_rhs[i*3+2] += rhs_y_eq[i];
    }
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        for (i=0; i<nnodes; i++) {
            rhs[i*3] = 0.;
            rhs[i*3+1] = rhs_x_eq[i];
            rhs[i*3+2] = rhs_y_eq[i];
        }
        printScreen_rhs_3dof("2D SW || CORIOLIS",nnodes, ie, elem2d->nodes, elem_rhs);
    }
#endif    
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                 WIND AND WAVE CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates the wind and wave addition to the elemental residual. \n
     *  \f{eqnarray*}{
     *   \residDA{i}{}{mx}  &=& -dt * \bodyStress{\,2d}{e}{\phidd{i}}{\lrpb{\tau_{winds,x}^{\,h} + \tau_{waves,x}^{\,h}}}   \\
     *   \residDA{i}{}{my}  &=& -dt * \bodyStress{\,2d}{e}{\phidd{i}}{\lrpb{\tau_{winds,y}^{\,h} + \tau_{waves,y}^{\,h}}}
     *  \f}
     *
     *********************************************************************************************/
    //Mark, work needs to be done here
    SVECT2D total_stress; svect2d_init(&total_stress);
    //if (wd_flag == 0) { // only apply on fully wet element        
//        SVECT2D wind_stress;  svect2d_init(&wind_stress);
//        SVECT2D wave_stress;  svect2d_init(&wave_stress);        
//        if (mod->flag.WIND && wd_flag == 0) {
//            int nws = 0;
//#ifdef _WINDLIB
//            nws = mod->windlib->nws;
//#endif
//            wind_stress = swind_elem2d_get_local_stress(2, sw2->winds, elem2d, elem2d->nodes, elem_h_wd_avg, g,
//                mod->density, mod->mat[imat].sw->windatt, mod->mat[imat].sw->wind_flag, mod->flag.WIND_LIBRARY, nws);
//        }        
//        if (mod->flag.WAVE && wd_flag == 0) {
//            wave_stress = swave_elem2d_get_local_stress(2, mod->flag.CSTORM_WSID, mod->sw->d2->waves, *elem2d, elem2d->nodes);
//        }        
//        total_stress.x = wind_stress.x + wave_stress.x;
//        total_stress.y = wind_stress.y + wave_stress.y;        
//        sarray_init_dbl(rhs_x_eq, nnodes); sarray_init_dbl(rhs_y_eq, nnodes);
//        if (isTriangle == TRUE) {
//            integrate_triangle_phi(djac, total_stress.x, rhs_x_eq);
//            integrate_triangle_phi(djac, total_stress.y, rhs_y_eq);
//        } else {
//            integrate_quadrilateral_phi(elem_nds, total_stress.x, rhs_x_eq);
//            integrate_quadrilateral_phi(elem_nds, total_stress.y, rhs_y_eq);
//        }        
//        for (i=0; i<nnodes; i++) {
//            elem_rhs[i*3+1] -= dt * rhs_x_eq[i];
//            elem_rhs[i*3+2] -= dt * rhs_y_eq[i];
//        }
//#ifdef _DEBUG
//        if (DEBUG == ON || DEBUG_LOCAL == ON) {
//            for (i=0; i<nnodes; i++) {
//                rhs[i*3] = 0.;
//                rhs[i*3+1] = - dt * rhs_x_eq[i];
//                rhs[i*3+2] = - dt * rhs_y_eq[i];
//            }
//            printScreen_rhs_3dof("2D SW || SURFACE_STRESS",nnodes, ie, elem2d->nodes, rhs);
//        }
//#endif
//    }    
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                 RAINFALL CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates the rainfall addition to the elemental residual. \n
     * \f$ \residDA{i}{}{c} = dt * \bodySource{\,2d}{e}{\phidd{i}}{R} \f$
     *
     *********************************************************************************************/
    double elem_water_source = 0.;
    //need to figure this out too
//    if (mod->str_values[elem2d->string].ol_flow.bc_flag == BCT_WATER_SOURCE) {
//        int isers = mod->str_values[elem2d->string].ol_flow.isigma;
//        elem_water_source = sseries_get_value(isers, mod->series_head, 0);        
//        sarray_init_dbl(rhs_c_eq, nnodes);
//        if (isTriangle == TRUE) {
//            integrate_triangle_phi(djac, elem_water_source, rhs_c_eq);
//        } else {
//            integrate_quadrilateral_phi(elem_nds, elem_water_source, rhs_c_eq);
//        }        
//        t1 = dt * factor;
//        for (i=0; i<nnodes; i++) {
//            elem_rhs[i*3] -= t1 * rhs_c_eq[i];
//        }
//#ifdef _DEBUG
//        if (DEBUG == ON || DEBUG_LOCAL == ON) {
//            for (i=0; i<nnodes; i++) {
//                rhs[i*3] = -t1 * rhs_c_eq[i];
//                rhs[i*3+1] = 0;
//                rhs[i*3+2] = 0;
//            }
//            printScreen_rhs_3dof("2D SW || SOURCE",nnodes, ie, elem2d->nodes, rhs);
//        }
//#endif
//    }    
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                 PRESSURE CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates the pressure additions to the elemental residual. \n
     *********************************************************************************************/
    if (PRESSURE_FLAG == ON) {        
        /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n
         * Calculates the body pressure on element. \n
         *  \f{eqnarray*}{
         *   \residDA{i}{}{mx}  &=& -dt * \bodyPressureDD{\,2d}{e}{\phidd{i}}{x}{(\depth{h})}    \\
         *   \residDA{i}{}{my}  &=& -dt * \bodyPressureDD{\,2d}{e}{\phidd{i}}{y}{(\depth{h})}
         *  \f}
         *********************************************************************************************/
        sarray_init_dbl(rhs, ndof);
        t1 = -one_2 * g * dt;
        if (isTriangle == TRUE) {
            vars[0] = 1;
            vars[1] = t1;
            t2 = fe_sw2_wet_dry_wrapper(rhs, elem_nds, elem_head, grad_shp, NULL, NULL, NULL, NULL, djac, ON, DEBUG, vars, fe_sw2_wd_pressure_triangle);
        } else {
            sarray_init_dbl(rhs_x_eq,nnodes);
            sarray_init_dbl(rhs_y_eq,nnodes);
            integrate_quadrilateral_dphi_f_f(elem_nds, t1, elem_head, rhs_x_eq, rhs_y_eq);
            for (i=0; i<nnodes; i++) {
                rhs[i*3] = 0.;
                rhs[i*3+1] = rhs_x_eq[i];
                rhs[i*3+2] = rhs_y_eq[i];
            }
        }
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON ){//|| debug.pressure || debug.rhs) {
            printScreen_rhs_3dof("2D SW || PRESSURE",nnodes, ie, elem2d->nodes, rhs);
        }
#endif
        sarray_add_replace_dbl(elem_rhs, rhs, ndof);        
        /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n
         * Calculates the bathymetry gradient integral on element. \n
         *  \f{eqnarray*}{
         *   \residDA{i}{}{mx}  &=& -dt * \bodySource{\,2d}{e}{\phidd{i}}{ g \, \depth{h} \frac{ \partial z_\bed{}}{\partial x} }    \\
         *   \residDA{i}{}{my}  &=& -dt * \bodySource{\,2d}{e}{\phidd{i}}{ g \, \depth{h} \frac{ \partial z_\bed{}}{\partial y} }
         *  \f}
         *********************************************************************************************/
        sarray_init_dbl(rhs, ndof);
        t1 = g * dt;
        if (isTriangle == TRUE) {
            vars[0] = 1;
            vars[1] = t1;
            t2 = fe_sw2_wet_dry_wrapper(rhs, elem_nds, elem_head, NULL, &elem_grad_z, NULL, NULL, NULL, djac, ON, DEBUG, vars, fe_sw2_wd_bodyForce_triangle);
        } else {
            sarray_init_dbl(rhs_x_eq,nnodes);
            sarray_init_dbl(rhs_y_eq,nnodes);
            integrate_quadrilateral_phi_h_df(elem_nds, t1, elem_head, node_z, rhs_x_eq, rhs_y_eq);
            for (i=0; i<nnodes; i++) {
                rhs[i*3] = 0.;
                rhs[i*3+1] = rhs_x_eq[i];
                rhs[i*3+2] = rhs_y_eq[i];
            }
        }
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON){// || debug.pressure || debug.rhs) {
            printScreen_rhs_3dof("2D SW || BODY FORCE",nnodes, ie, elem2d->nodes, rhs);
        }
#endif
        sarray_add_replace_dbl(elem_rhs, rhs, ndof);       
        /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n
         * Calculates the boundary integral on element edges. This is to auto-find no-flow boundaries.  \n
         *  \f{eqnarray*}{
         *   \residDA{i}{}{mx}  &=& dt * \bcPressureDD{}{e}{\phidd{i}}{(\depth{h})}{x}     \\
         *   \residDA{i}{}{my}  &=& dt * \bcPressureDD{}{e}{\phidd{i}}{(\depth{h})}{y}
         *  \f}
         *********************************************************************************************/
        sarray_init_dbl(rhs, ndof);
        t1 = one_2 * dt * g;
        if (isTriangle == TRUE) {
            vars[0] = 1;
            vars[1] = t1;
            t2 = fe_sw2_wet_dry_wrapper(rhs, elem_nds, elem_head, NULL, NULL, NULL, NULL, NULL, djac, ON, DEBUG, vars, fe_sw2_wd_boundaryPressure_triangle);
        } else {
            int nd1, nd2, iedge;
            double delx, dely;		/* the signed edge x and y lengths */            
            // loop over edges to evaluate edge integrals
            for (iedge=0; iedge<nnodes; iedge++) {
                nd1 = grid->nd_on_QuadEdge[iedge][0]; //doesnt seem to exist anymore, can add back small array
                nd2 = grid->nd_on_QuadEdge[iedge][1];                
                delx = elem_nds[nd1].x - elem_nds[nd2].x;
                dely = elem_nds[nd2].y - elem_nds[nd1].y;                
                sarray_init_dbl(rhs_c_eq,nnodes);
                integrate_line_phi_h_h(1., t1, elem_head[nd1], elem_head[nd2], rhs_c_eq); // set djac to 1 to multiply below
                rhs[nd1*3+1] += dely * rhs_c_eq[0];
                rhs[nd2*3+1] += dely * rhs_c_eq[1];
                rhs[nd1*3+2] += delx * rhs_c_eq[0];
                rhs[nd2*3+2] += delx * rhs_c_eq[1];
            }
        }
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON ){// debug.pressure || debug.rhs) {
            printScreen_rhs_3dof("2D SW || BOUNDARY PRESSURE",nnodes, ie, elem2d->nodes, rhs);
        }
#endif
        sarray_add_replace_dbl(elem_rhs, rhs, ndof);                 
        /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         *                                BAROCLINIC PRESSURE CONTRIBUTION
         *-------------------------------------------------------------------------------------------
         * Calculates the density pressure additions to the elemental residual. \n
         *********************************************************************************************/
        //NEED TO FIGURE OUT FLAGS 
        if (FALSE){//(mod->flag.BAROCLINIC != OFF) {            
            /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n
             * Calculates the density driven body pressure on element. \n
             *  \f{eqnarray*}{
             *   \residDA{i}{}{mx}  &=& -dt * \bodyDensityPressureDD{\,2d}{e}{\phidd{i}}{x}{(\depth{h})}    \\
             *   \residDA{i}{}{my}  &=& -dt * \bodyDensityPressureDD{\,2d}{e}{\phidd{i}}{y}{(\depth{h})}
             *  \f}
             *********************************************************************************************/
            sarray_init_dbl(rhs, ndof);
            t1 = -one_2 * g * dt;
            if (isTriangle == TRUE) {
                vars[0] = 1;
                vars[1] = t1;
                t2 = fe_sw2_wet_dry_wrapper(rhs, elem_nds, elem_head, grad_shp, NULL, NULL, elem_density, NULL, djac, ON, DEBUG, vars, fe_sw2_wd_densityPressure_triangle);
            } else {
                sarray_init_dbl(rhs_x_eq,nnodes);
                sarray_init_dbl(rhs_y_eq,nnodes);
                integrate_quadrilateral_dphi_f_h_h(elem_nds, t1, elem_density, elem_head, rhs_x_eq, rhs_y_eq);
                for (i=0; i<nnodes; i++) {
                    rhs[i*3] = 0.;
                    rhs[i*3+1] = rhs_x_eq[i];
                    rhs[i*3+2] = rhs_y_eq[i];
                }
            }
#ifdef _DEBUG
            if (DEBUG == ON || DEBUG_LOCAL == ON){ //|| debug.pressure || debug.rhs) {
                printScreen_rhs_3dof("2D SW || DENSITY PRESSURE",nnodes, ie, elem2d->nodes, rhs);
            }
#endif
            sarray_add_replace_dbl(elem_rhs, rhs, ndof);     
            /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n
             * Calculates the density bathymetry gradient integral on element. \n
             *  \f{eqnarray*}{
             *   \residDA{i}{}{mx}  &=& -dt * \bodySource{\,2d}{e}{\phidd{i}}{ g \, \depth{h} \frac{ \partial z_\bed{}}{\partial x} \, \frac{\bigtriangleup \rho}{\rho_o}}    \\
             *   \residDA{i}{}{my}  &=& -dt * \bodySource{\,2d}{e}{\phidd{i}}{ g \, \depth{h} \frac{ \partial z_\bed{}}{\partial y} \, \frac{\bigtriangleup \rho}{\rho_o}}
             *  \f}
             *********************************************************************************************/
            sarray_init_dbl(rhs, ndof);
            t1 = g * dt;
            if (isTriangle == TRUE) {
                vars[0] = 1;
                vars[1] = t1;
                t2 = fe_sw2_wet_dry_wrapper(rhs, elem_nds, elem_head, NULL, &elem_grad_z, NULL, elem_density, NULL, djac, ON, DEBUG, vars, fe_sw2_wd_densityBodyForce_triangle);
            } else {
                sarray_init_dbl(rhs_x_eq,nnodes);
                sarray_init_dbl(rhs_y_eq,nnodes);
                integrate_quadrilateral_phi_h_g_df(elem_nds, t1, elem_density, elem_head, node_z, rhs_x_eq, rhs_y_eq);
                for (i=0; i<nnodes; i++) {
                    rhs[i*3] = 0.;
                    rhs[i*3+1] = rhs_x_eq[i];
                    rhs[i*3+2] = rhs_y_eq[i];
                }
            }
#ifdef _DEBUG
            if (DEBUG == ON || DEBUG_LOCAL == ON){// || debug.pressure || debug.rhs) {
                printScreen_rhs_3dof("2D SW || DENSITY BODY FORCE PRESSURE",nnodes, ie, elem2d->nodes, rhs);
            }
#endif
            sarray_add_replace_dbl(elem_rhs, rhs, ndof);             
            /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n
             * Calculates the density driven pressure boundary integral on element edges. This is to auto-find no-flow boundaries.  \n
             *  \f{eqnarray*}{
             *   \residDA{i}{}{mx}  &=& dt * \bcDensityPressureDD{}{e}{\phidd{i}}{(\depth{h})}{x}     \\
             *   \residDA{i}{}{my}  &=& dt * \bcDensityPressureDD{}{e}{\phidd{i}}{(\depth{h})}{y}
             *  \f}
             *********************************************************************************************/
            sarray_init_dbl(rhs, ndof);
            t1 = one_2 * g * dt;
            if (isTriangle == TRUE) {
                vars[0] = 1;
                vars[1] = t1;
                t2 = fe_sw2_wet_dry_wrapper(rhs, elem_nds, elem_head, NULL, NULL, NULL, elem_density, NULL, djac, ON, DEBUG, vars, fe_sw2_wd_densityBoundaryPressure_triangle);
            } else {
                int nd1, nd2, iedge;
                double delx, dely;		/* the signed edge x and y lengths */                
                // loop over edges to evaluate edge integrals
                for (iedge=0; iedge<nnodes; iedge++) {
                    nd1 = grid->nd_on_QuadEdge[iedge][0];
                    nd2 = grid->nd_on_QuadEdge[iedge][1];                    
                    delx = elem_nds[nd1].x - elem_nds[nd2].x;
                    dely = elem_nds[nd2].y - elem_nds[nd1].y;                   
                    sarray_init_dbl(rhs_c_eq,nnodes);
                    integrate_line_phi_h_h_g(1., t1, elem_head[nd1], elem_head[nd2], elem_density[nd1], elem_density[nd2], rhs_c_eq); // set djac to 1 to multiply below
                    rhs[nd1*3+1] += dely * rhs_c_eq[0];
                    rhs[nd2*3+1] += dely * rhs_c_eq[1];
                    rhs[nd1*3+2] += delx * rhs_c_eq[0];
                    rhs[nd2*3+2] += delx * rhs_c_eq[1];
                }
            }
#ifdef _DEBUG
            if (DEBUG == ON || DEBUG_LOCAL == ON){// || debug.pressure || debug.rhs) {
                printScreen_rhs_3dof("2D SW || DENSITY BOUNDARY PRESSURE",nnodes, ie, elem2d->nodes, rhs);
            }
#endif
            sarray_add_replace_dbl(elem_rhs, rhs, ndof);  
        }//end of BAROCLINIC
    } //end of PRESSURE     
    /*!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                      SUPG CONTRIBUTION
     *-------------------------------------------------------------------------------------------
     * Calculates the SUPG  addition to the elemental residual. \n
     *  \f{eqnarray*}{
     *  \weakSwCsupgDD{2d}{e}{i} {\ub{h}} {\depth{h}} {\vb{h}} {\depth{h}} {c} \\
     *  \weakSwMXsupgDD{2d}{e}{i}{g\,\depth{h}}{\ub{h} \depth{h}}{}{\vb{h} \depth{h}}{mx} \\
     *  \weakSwMYsupgDD{2d}{e}{i}{g\,\depth{h}}{\vb{h} \depth{h}}{}{\ub{h} \depth{h}}{my}
     *  \f}
     *
     * \note CJT \:: For tau, use method flag == 2, the elemental area to get elemental length estimate
     * \note CJT \:: Not sure what g_factor is in the tau calculation and why area is X factor_old
     *********************************************************************************************/
    if (is_double_small(elem_h_wd_avg) == NO) {        
        sarray_init_dbl(rhs_c_eq, nnodes);
        sarray_init_dbl(rhs_x_eq, nnodes);
        sarray_init_dbl(rhs_y_eq, nnodes);        
        // get SUPG tau :: tau = alpha * le / sqrt(u^2 + v^2 + g*H)
        // CJT \::
        double c_eq = 0., x_eq = 0., y_eq = 0.;
        double g_factor = (1+elem_density_avg) * g;
        double grad_shp_x[nnodes]; for (i=0; i<nnodes; i++) {grad_shp_x[i] = grad_shp[i].x;}
        double grad_shp_y[nnodes]; for (i=0; i<nnodes; i++) {grad_shp_y[i] = grad_shp[i].y;}
        double elem_h_wd_avg_denom = MAX(elem_h_wd_avg, NOT_QUITE_SMALL);
        double tau = fe_get_supg_tau_sw(nnodes, elem_nds, g_factor, elem_h_wd_avg_denom, elem_vel_wd_avg.x, elem_vel_wd_avg.y,0.,grad_shp_x,grad_shp_y,NULL,area * factor_old, sw2->tau_pg, 2, 2, 2);//need to figure out         
        // calculate constant strong equation set. This may mean elementally averaging each term.
        if (isTriangle == TRUE) {
            //replaced            
            double dfdt[ndof];
            double advection[ndof]; 
            int redistribute_flag = OFF; // do not redistribute in SUPG            
            // elementally averaged wet/dry temporal derivatives
            double dhdt = 0., dudt = 0., dvdt = 0.;
            sarray_init_dbl(dfdt, ndof); // cjt :: note: older version of AdH all use current velocity here ... looks wrong
            sarray_init_dbl(advection, ndof);
            vars[0] = 1;
            vars[1] = 1;
            factor_old = fe_sw2_wet_dry_wrapper(dfdt, elem_nds, elem_old_head, NULL, NULL, elem_old_vel, NULL, elem_old_head, djac, redistribute_flag, DEBUG, vars, fe_sw2_wd_integrate_triangle_f); //fe_sw2_wd_average_tri);
            //dudt = -dfdt[0].x_eq;  dvdt = -dfdt[0].y_eq; dhdt = -dfdt[0].c_eq;
            dudt = -dfdt[1];  dvdt = -dfdt[2]; dhdt = -dfdt[0];            
            sarray_init_dbl(dfdt, ndof);
            vars[0] = 1;
            vars[1] = 1;
            factor = fe_sw2_wet_dry_wrapper(dfdt, elem_nds, elem_head, NULL, NULL, elem_vel, NULL, elem_head, djac, redistribute_flag, DEBUG, vars, fe_sw2_wd_integrate_triangle_f); //fe_sw2_wd_average_tri);
            dudt += dfdt[1];  dvdt += dfdt[2];  dhdt += dfdt[0];           
            // CJT :: tried this, doesnt work for angle dam for some reason ...
            //dhdt = elem_h_wd_avg - elem_h_old_wd_avg;
            //dudt = elem_vel_wd_avg.x - elem_old_vel_wd_avg.x;
            //dvdt = elem_vel_wd_avg.y - elem_old_vel_wd_avg.y;            
            // wet/dry convection terms :: these are already constant since they a gradients
            vars[0] = dt;
            vars[1] = 1;
            factor = fe_sw2_wet_dry_wrapper(advection, elem_nds, elem_head, grad_shp, NULL, elem_vel, NULL, elem_head, djac, redistribute_flag, DEBUG, vars, fe_sw2_wd_gls_convection_triangle);         
            // pressure terms ::  these are already constant since they a gradients
            //  g * dh/dx * dt + g * dz_dx * dt + dpress/dx * dt
            double x_pressure = g_factor * (elem_grad_head.x + elem_grad_z.x) * dt * factor + elem_grad_prs.x * dt * factor;
            double y_pressure = g_factor * (elem_grad_head.y + elem_grad_z.y) * dt * factor + elem_grad_prs.y * dt * factor;            
            // friction terms :: cjt :: has factor but not wet-dry wrapped??
            double x_friction = g * elem_bed_stress.x * dt * factor;
            double y_friction = g * elem_bed_stress.y * dt * factor;           
            // coriolis terms 
            //NEED TO ADD BACK
            double x_coriolis =0.0;// -elem_vel_wd_avg.y * coriolis_speed;
            double y_coriolis = 0.0;// elem_vel_wd_avg.x * coriolis_speed;

            // constant, strong equation set :: cjt :: trunk has factor on total stress and coriolis, but there are not wet-dry wrapped
            c_eq = dhdt + advection[0] + factor * elem_water_source * djac * dt;
            x_eq = dudt + advection[1] + (x_pressure + x_friction/elem_h_wd_avg_denom) * djac + (x_coriolis - total_stress.x / elem_h_wd_avg_denom) * djac * factor * dt;
            y_eq = dvdt + advection[2] + (y_pressure + y_friction/elem_h_wd_avg_denom) * djac + (y_coriolis - total_stress.y / elem_h_wd_avg_denom) * djac * factor * dt;           
            //old adh way :: this does make a subtle difference in codes when comparing *************
            //                        double tx_v[3] = { 0., 0., 0.};
            //                        double ty_v[3] = { 0., 0., 0.};
            //                        double txy[2][9] = { {0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 0.} };
            //                        double celerity_squared = g_factor * elem_h_wd_avg_denom;
            //                        if (elem_h_wd_avg_denom <= 0.0) celerity_squared = 0.;
            //                        double value = sqrt(pow(elem_vel_wd_avg.x,2) + pow(elem_vel_wd_avg.y,2) + celerity_squared + SMALL);
            //                        txy[0][0] = elem_vel_wd_avg.x / value;
            //                        txy[0][1] = elem_h_wd_avg / value;
            //                        txy[0][2] = 0.;
            //                        txy[0][3] = g_factor / value;
            //                        txy[0][4] = elem_vel_wd_avg.x / value;
            //                        txy[0][5] = 0.;
            //                        txy[0][6] = 0.;
            //                        txy[0][7] = 0.;
            //                        txy[0][8] = elem_vel_wd_avg.x / value;
            //            
            //                        txy[1][0] = elem_vel_wd_avg.y / value;
            //                        txy[1][1] = 0.;
            //                        txy[1][2] = elem_h_wd_avg / value;
            //                        txy[1][3] = 0.;
            //                        txy[1][4] = elem_vel_wd_avg.y / value;
            //                        txy[1][5] = 0.;
            //                        txy[1][6] = g_factor / value;
            //                        txy[1][7] = 0.;
            //                        txy[1][8] = elem_vel_wd_avg.y / value;
            //            
            //            
            //                        double delx = mod->tau_pg * sqrt(djac * factor_old);
            //                        double dely = delx;
            //                        tx_v[0] = delx * (txy[0][0] * c_eq + txy[0][1] * x_eq + txy[0][2] * y_eq);
            //                        tx_v[1] = delx * (txy[0][3] * c_eq * elem_h_wd_avg + txy[0][4] * x_eq * elem_h_wd_avg + txy[0][5] * y_eq * elem_h_wd_avg);  /* COME BACK */
            //                        tx_v[2] = delx * (txy[0][6] * c_eq * elem_h_wd_avg + txy[0][7] * x_eq * elem_h_wd_avg + txy[0][8] * y_eq * elem_h_wd_avg);  /* COME BACK */
            //                        ty_v[0] = dely * (txy[1][0] * c_eq + txy[1][1] * x_eq + txy[1][2] * y_eq);
            //                        ty_v[1] = dely * (txy[1][3] * c_eq * elem_h_wd_avg + txy[1][4] * x_eq * elem_h_wd_avg + txy[1][5] * y_eq * elem_h_wd_avg);  /* COME BACK */
            //                        ty_v[2] = dely * (txy[1][6] * c_eq * elem_h_wd_avg + txy[1][7] * x_eq * elem_h_wd_avg + txy[1][8] * y_eq * elem_h_wd_avg);  /* COME BACK */
            //            
            //                        double const_contrib = 0.;
            //            
            //                        const_contrib = tx_v[0] * 1.0;
            //                        rhs_c_eq[0] = const_contrib * grad_shp[0].x;
            //                        rhs_c_eq[1] = const_contrib * grad_shp[1].x;
            //                        rhs_c_eq[2] = const_contrib * grad_shp[2].x;
            //            
            //                        const_contrib = tx_v[1] * 1.0;
            //                        rhs_x_eq[0] = const_contrib * grad_shp[0].x;
            //                        rhs_x_eq[1] = const_contrib * grad_shp[1].x;
            //                        rhs_x_eq[2] = const_contrib * grad_shp[2].x;
            //            
            //                        const_contrib = tx_v[2] * 1.0;
            //                        rhs_y_eq[0] = const_contrib * grad_shp[0].x;
            //                        rhs_y_eq[1] = const_contrib * grad_shp[1].x;
            //                        rhs_y_eq[2] = const_contrib * grad_shp[2].x;
            //            
            //                        const_contrib = ty_v[0] * 1.0;
            //                        rhs_c_eq[0] += const_contrib * grad_shp[0].y;
            //                        rhs_c_eq[1] += const_contrib * grad_shp[1].y;
            //                        rhs_c_eq[2] += const_contrib * grad_shp[2].y;
            //            
            //                        const_contrib = ty_v[1] * 1.0;
            //                        rhs_x_eq[0] += const_contrib * grad_shp[0].y;
            //                        rhs_x_eq[1] += const_contrib * grad_shp[1].y;
            //                        rhs_x_eq[2] += const_contrib * grad_shp[2].y;
            //            
            //                        const_contrib = ty_v[2] * 1.0;
            //                        rhs_y_eq[0] += const_contrib * grad_shp[0].y;
            //                        rhs_y_eq[1] += const_contrib * grad_shp[1].y;
            //                        rhs_y_eq[2] += const_contrib * grad_shp[2].y;
            //            
            //#ifdef _DEBUG
            //                        if (DEBUG == ON || DEBUG_LOCAL == ON) {
            //                            printf("delx: %30.20e \t value: %30.20e \t h_avg: %30.20e \t u_avg: %30.20e \t v_avg: %30.20e\n",delx,value,elem_h_wd_avg,elem_vel_wd_avg.x,elem_vel_wd_avg.y);
            //                            printf("tau: %30.20e \t elem_h_wd_avg: %30.20e \t elem_vel_wd_avg.x: %30.20e \t elem_vel_wd_avg.y: %30.20e \t g_factor: %30.20e\n",tau,elem_h_wd_avg,elem_vel_wd_avg.x,elem_vel_wd_avg.y,g_factor);
            //                            printf("dhdt: %30.20e \t dudt: %30.20e \t dvdt: %30.20e \n",dhdt,dudt,dvdt);
            //                            printf("advection:  c: %30.20e \t x: %30.20e \t y: %30.20e\n",advection[0].c_eq,advection[0].x_eq,advection[0].y_eq);
            //                            printf("coriolis:  x: %30.20e \t y: %30.20e\n",x_coriolis,y_coriolis);
            //                            printf("pressure:  x: %30.20e \t y: %30.20e\n",x_pressure,y_pressure);
            //                            printf("friction:  x: %30.20e \t y: %30.20e\n",x_friction,y_friction);
            //                            printf("grad_x :: tau_c_c: %30.20e \t tau_c_x: %30.20e \t tau_c_y: %30.20e\n",tau * elem_vel_wd_avg.x,tau * elem_h_wd_avg,tau * 0);
            //                            printf("grad_x :: tau_x_c: %30.20e \t tau_x_x: %30.20e \t tau_x_y: %30.20e\n",tau * g_factor * elem_h_wd_avg ,tau * elem_vel_wd_avg.x * elem_h_wd_avg,tau * 0);
            //                            printf("grad_x :: tau_y_c: %30.20e \t tau_y_x: %30.20e \t tau_y_y: %30.20e\n",tau * 0 ,tau * 0,tau * elem_vel_wd_avg.x * elem_h_wd_avg);
            //            
            //                            printf("grad_y :: tau_c_c: %30.20e \t tau_c_x: %30.20e \t tau_c_y: %30.20e\n",tau * elem_vel_wd_avg.y, tau * 0, tau * elem_h_wd_avg);
            //                            printf("grad_y :: tau_x_c: %30.20e \t tau_x_x: %30.20e \t tau_x_y: %30.20e\n",tau * 0 ,tau * elem_vel_wd_avg.y * elem_h_wd_avg,tau * 0);
            //                            printf("grad_y :: tau_y_c: %30.20e \t tau_y_x: %30.20e \t tau_y_y: %30.20e\n",tau * g_factor * elem_h_wd_avg ,tau * 0,tau * elem_vel_wd_avg.y * elem_h_wd_avg);
            //                            printf("tx_v: %30.20e \t %30.20e \t %30.20e\n",tx_v[0],tx_v[1],tx_v[2]);
            //                            printf("ty_v: %30.20e \t %30.20e \t %30.20e\n",ty_v[0],ty_v[1],ty_v[2]);
            //                            printf("c_eq: %30.20e \t x_eq: %30.20e \t y_eq: %30.20e\n",c_eq,x_eq,y_eq);
            //                        }
            //#endif            
            for (i=0; i<nnodes; i++) {
                rhs_c_eq[i] = tau * (grad_shp[i].x * (elem_vel_wd_avg.x * c_eq + elem_h_wd_avg * x_eq) +
                                     grad_shp[i].y * (elem_vel_wd_avg.y * c_eq + elem_h_wd_avg * y_eq));
                rhs_x_eq[i] = tau * (grad_shp[i].x * (g_factor * elem_h_wd_avg * c_eq + elem_vel_wd_avg.x * elem_h_wd_avg * x_eq) +
                                     grad_shp[i].y * (elem_vel_wd_avg.y * elem_h_wd_avg * x_eq));
                rhs_y_eq[i] = tau * (grad_shp[i].y * (g_factor * elem_h_wd_avg * c_eq + elem_vel_wd_avg.y * elem_h_wd_avg * y_eq) +
                                     grad_shp[i].x * (elem_vel_wd_avg.x * elem_h_wd_avg * y_eq));                
                rhs_c_eq[i] = (grad_shp[i].x * (tau * elem_vel_wd_avg.x * c_eq + tau * elem_h_wd_avg * x_eq) +
                               grad_shp[i].y * (tau * elem_vel_wd_avg.y * c_eq + tau * elem_h_wd_avg * y_eq));
                rhs_x_eq[i] = (grad_shp[i].x * (tau * g_factor * elem_h_wd_avg * c_eq + tau * elem_vel_wd_avg.x * elem_h_wd_avg * x_eq) +
                               grad_shp[i].y * (tau * elem_vel_wd_avg.y * elem_h_wd_avg * x_eq));
                rhs_y_eq[i] = (grad_shp[i].y * (tau * g_factor * elem_h_wd_avg * c_eq + tau * elem_vel_wd_avg.y * elem_h_wd_avg * y_eq) +
                               grad_shp[i].x * (tau * elem_vel_wd_avg.x * elem_h_wd_avg * y_eq));
            }   

        } else {            
            // elementally averaged gradients
            SVECT2D elem_avg_grad_head = integrate_quadrilateral_gradF(elem_nds, 1./area, elem_head);
            SVECT2D elem_avg_grad_z =    integrate_quadrilateral_gradF(elem_nds, 1./area, node_z);
            SVECT2D elem_avg_grad_p =    integrate_quadrilateral_gradF(elem_nds, 1./area, elem_pressure);
            SVECT2D elem_avg_grad_u =    integrate_quadrilateral_gradF(elem_nds, 1./area, elem_u);
            SVECT2D elem_avg_grad_v =    integrate_quadrilateral_gradF(elem_nds, 1./area, elem_v);            
            // elementally average temporal terms
            double dhdt = elem_h_avg - elem_h_old_avg;
            double dudt = elem_u_avg - elem_u_old_avg;
            double dvdt = elem_v_avg - elem_v_old_avg;            
            // friction terms
            double x_friction = g * elem_bed_stress.x * dt / elem_h_wd_avg_denom;
            double y_friction = g * elem_bed_stress.y * dt / elem_h_wd_avg_denom;            
            // elementally averaged coriolis terms
            double x_coriolis = -elem_v_avg * coriolis_speed * dt;
            double y_coriolis =  elem_u_avg * coriolis_speed * dt;            
            // elementally averaged pressure terms :: g * dh/dx * dt + g * dz_dx * dt + dpress/dx * dt
            double x_pressure = g_factor * (elem_avg_grad_head.x + elem_avg_grad_z.x) * dt + elem_avg_grad_p.x * dt;
            double y_pressure = g_factor * (elem_avg_grad_head.y + elem_avg_grad_z.y) * dt + elem_avg_grad_p.y * dt;            
            // elementally averaged advection terms (cjt :: could average the whole term together as well)
            double c_advection = (elem_u_avg * elem_avg_grad_head.x + elem_h_avg * elem_avg_grad_u.x +
                                  elem_v_avg * elem_avg_grad_head.y + elem_h_avg * elem_avg_grad_v.y) * dt;
            double x_advection = (elem_u_avg * elem_avg_grad_u.x + elem_v_avg * elem_avg_grad_u.y) * dt;
            double y_advection = (elem_u_avg * elem_avg_grad_v.x + elem_v_avg * elem_avg_grad_v.y) * dt;            
            // elementally constant equations
            c_eq = dhdt + c_advection - elem_water_source * dt;
            x_eq = dudt + x_advection + x_pressure + x_friction + x_coriolis - total_stress.x / elem_h_wd_avg_denom;
            y_eq = dvdt + y_advection + y_pressure + y_friction + y_coriolis - total_stress.y / elem_h_wd_avg_denom;            
            // integrate[ tau * (grad(phi) dot elem_vel_wd_avg.y) * constant_strong_eq ]
            SVECT2D term;
            term.x = elem_vel_wd_avg.x * c_eq + elem_h_wd_avg * x_eq;
            term.y = elem_vel_wd_avg.y * c_eq + elem_h_wd_avg * y_eq;
            integrate_quadrilateral_gradPhi_dot_vbar(elem_nds, dt * tau, term, rhs_c_eq);            
            term.x = g_factor * g * elem_h_wd_avg * c_eq + elem_vel_wd_avg.x * elem_h_wd_avg * x_eq;
            term.y = elem_vel_wd_avg.y * elem_h_wd_avg * x_eq;
            integrate_quadrilateral_gradPhi_dot_vbar(elem_nds, dt * tau, term, rhs_x_eq);            
            term.x = g_factor * g * elem_h_wd_avg * c_eq + elem_vel_wd_avg.y * elem_h_wd_avg * y_eq;
            term.y = elem_vel_wd_avg.x * elem_h_wd_avg * y_eq;
            integrate_quadrilateral_gradPhi_dot_vbar(elem_nds, dt * tau, term, rhs_y_eq);            
        }
        for (i=0; i<nnodes; i++) {
            //NEED TO ADD BACK IN
            sw2->elem_rhs_dacont_extra_terms[ie][i] += rhs_c_eq[i];
            elem_rhs[i*3] += rhs_c_eq[i];
            elem_rhs[i*3+1] += rhs_x_eq[i];
            elem_rhs[i*3+2] += rhs_y_eq[i];
        }
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) {
            for (i=0; i<nnodes; i++) {
                rhs[i*3] = rhs_c_eq[i];
                rhs[i*3+1] = rhs_x_eq[i];
                rhs[i*3+2] = rhs_y_eq[i];
            }
            printScreen_rhs_3dof("2D SW || SUPG",nnodes, ie, elem2d->nodes, rhs);
            Is_DoubleArray_Inf_or_NaN(rhs_c_eq, nnodes ,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(rhs_x_eq, nnodes ,__FILE__ ,__LINE__);
            Is_DoubleArray_Inf_or_NaN(rhs_y_eq, nnodes ,__FILE__ ,__LINE__);
        }
#endif

    }    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    SET DIRICHLET BCS
     *--------------------------------------------------------------------------------------------
     * Sets the Dirichlet boundary conditions. \n
     *********************************************************************************************/
    //Don't think we need this anymore??
//    int istring = UNSET_INT;
//    for (i=0; i<nnodes; i++) {
//        istring = grid->node[elem2d->nodes[i]].string;
//        if (istring > NORMAL) {
//            if (mod->str_values[istring].ol_flow.bc_flag == BCT_PRS_DIR) {
//                elem_rhs[i].c_eq = 0.0;
//            }
//            else if(mod->str_values[istring].ol_flow.bc_flag == BCT_VEL_DIR) {
//                elem_rhs[i].x_eq = 0.0;
//                elem_rhs[i].y_eq = 0.0;
//            }
//            else if(mod->str_values[istring].ol_flow.bc_flag == BCT_VEL_PRS_DIR) {
//                elem_rhs[i].c_eq = 0.0;
//                elem_rhs[i].x_eq = 0.0;
//                elem_rhs[i].y_eq = 0.0;
//            }
//        }
//    }
//    
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        printScreen_rhs_3dof("2D SW ||TOTAL RHS",nnodes, ie, elem2d->nodes, elem_rhs);
        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    }    
    time_t time2;  time(&time2);
    //TIME_IN_2D_SW_BODY_RESID += difftime(time2,time1);
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 2D body temporal contributions to the shallow water equations.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  \note CJT\:: this calculation groups hu,v in the momentum temporal terms
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_temporal(int ie, SELEM_2D *elem2d, int nnodes, SVECT *elem_nds, double djac, double drying_lower_limit, double *elem_head, SVECT2D *elem_vel, double wd_factor, double dt_factor, double *elem_rhs, char *string, int DEBUG, int DEBUG_LOCAL, int wd_flag) {    
    int i;
    double con = dt_factor * wd_factor;
    double rhs[nnodes*3]; sarray_init_dbl(rhs,nnodes*3);
    double uh[nnodes], vh[nnodes];
    double rhs_c_eq[nnodes], rhs_x_eq[nnodes], rhs_y_eq[nnodes];
    double vars[2];    
    for (i=0; i<nnodes; i++) {
        uh[i] = elem_vel[i].x * MAX(elem_head[i], drying_lower_limit);; // CJT :: group momentum terms
        vh[i] = elem_vel[i].y * MAX(elem_head[i], drying_lower_limit);; // CJT :: group momentum terms
        rhs_c_eq[i] = 0.; rhs_x_eq[i] = 0.; rhs_y_eq[i] = 0.;
    }   
    if (nnodes == NDONTRI) {
        // consistent mass matrix for continuity :: cjt :: gives slightly different results than trunk, but ok, its just the simplified integration routine
        vars[0] = 1.0;
        vars[1] = dt_factor;        
        double t1 = fe_sw2_wet_dry_wrapper(rhs, elem_nds, elem_head, NULL, NULL, NULL, NULL, NULL, djac, ON, DEBUG, vars, fe_sw2_wd_continuity_temporal_triangle);        
        // lumped AND grouped momentum mass matrix
        integrate_triangle_phi_f_lump(djac, con, uh, rhs_x_eq);
        integrate_triangle_phi_f_lump(djac, con, vh, rhs_y_eq);        
        for (i=0; i<nnodes; i++) {
            if (wd_flag != 2) {
                rhs_c_eq[i] = rhs[i*3];
                elem_rhs[i*3] += rhs[i*3];
            }
            elem_rhs[i*3+1] += rhs_x_eq[i];
            elem_rhs[i*3+2] += rhs_y_eq[i];
        }        
    } else {
        // consistent mass matrix for continuity
        integrate_quadrilateral_phi_f(elem_nds, 1., elem_head, rhs_c_eq);        
        // lumped AND grouped momentum mass matrix
        integrate_quadrilateral_phi_f_lump(elem_nds, 1., uh, rhs_x_eq);
        integrate_quadrilateral_phi_f_lump(elem_nds, 1., vh, rhs_y_eq);       
        for (i=0; i<nnodes; i++) {
            rhs_c_eq[i] = dt_factor * rhs_c_eq[i];
            rhs_x_eq[i] = con * rhs_x_eq[i];
            rhs_y_eq[i] = con * rhs_y_eq[i];
            elem_rhs[i*3] += rhs_c_eq[i];
            elem_rhs[i*3+1] += rhs_x_eq[i];
            elem_rhs[i*3+2]  += rhs_y_eq[i];
        }        
    }
#ifdef _DEBUG
    for (i=0; i<nnodes; i++) {
        rhs[i*3] = rhs_c_eq[i];
        rhs[i*3+1] = rhs_x_eq[i];
        rhs[i*3+2] = rhs_y_eq[i];
    }
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        printScreen_rhs_3dof(string,nnodes, ie, elem2d->nodes, rhs);
        Is_DoubleArray_Inf_or_NaN(rhs_c_eq, nnodes ,__FILE__ ,__LINE__);
        Is_DoubleArray_Inf_or_NaN(rhs_x_eq, nnodes ,__FILE__ ,__LINE__);
        Is_DoubleArray_Inf_or_NaN(rhs_y_eq, nnodes ,__FILE__ ,__LINE__);
    }
#endif
}

