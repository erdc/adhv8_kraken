/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_3d_transport_body_resid.c This file collections functions responsible for
 *          the 3D shallow water transport body contributions to the elemental residual.    */
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

static SGRID *grid = NULL;
static SQUAD *quad = NULL;
static int icol = 0, nnodes = 0, nnodes_quad = 0, isElementTetrahedron = TRUE, imat = UNSET_INT;
static SELEM_3D *elem3d = NULL;
static SELEM_2D *elem2d_sur = NULL;
static SELEM_2D *elem2d_bed = NULL;
static SNODE elem_nodes[NDONPRISM];
static double g = 0., dt = 0., alpha = 0.;
static double elem_c[NDONPRISM], elem_c_old[NDONPRISM], elem_c_older[NDONPRISM], elem_density[NDONPRISM];
static double elem_dpl[NDONPRISM], elem_dpl_old[NDONPRISM], elem_dpl_older[NDONPRISM];
static double elem_bed_dpl[NDONPRISM], elem_bed_dpl_old[NDONPRISM], elem_bed_dpl_older[NDONPRISM];
static SVECT elem_nds[NDONPRISMQUAD], elem_nds_old[NDONPRISMQUAD], elem_nds_older[NDONPRISMQUAD];
static SVECT elem_vel[NDONPRISM], elem_rel_vel[NDONPRISM], elem_grid_vel[NDONPRISM];
static double elem_u[NDONPRISM], elem_v[NDONPRISM], elem_w[NDONPRISM];
static double elem_ur[NDONPRISM], elem_vr[NDONPRISM], elem_wr[NDONPRISM];
static SVECT grad_phi[NDONPRISM], grad_phi_old[NDONPRISM], grad_phi_older[NDONPRISM];
static double elem_volume = 0., elem_volume_old = 0., elem_volume_older = 0., elem_djac = 0., elem_djac_old = 0., elem_djac_older = 0.;
static SVECT elem_grad_u, elem_grad_v, elem_grad_w, elem_grad_c, elem_grad_density;
static double elem_avg_c = 0., elem_avg_c_old = 0., elem_avg_c_older = 0., elem_avg_density = 0., elem_avg_u = 0., elem_avg_v = 0., elem_avg_w = 0.;
static double elem_avg_ur = 0., elem_avg_vr = 0., elem_avg_wr = 0.;
static SVECT  elem_avg_grad_ur,  elem_avg_grad_vr,  elem_avg_grad_wr, elem_avg_grad_density, elem_avg_grad_c;
static double elem_avg_strong_resid = 0., elem_avg_depth = 0., elem_avg_vmag_safe = 0.;
static double molecular_diffusion = 0., turbulent_diffusion = 0., vertical_diffusion = 0., tur_viscosity = 0., total_diffusion = 0.;
static double dist_above_bed = 0., average_bed_z = 0., average_elem_z = 0.;


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// file prototypes
int fe_3d_transport_body_test_for_NAN_or_INF();
void fe_3d_transport_body_debug(double perturbation, int perturb_node, int perturb_var, int perturb_sign);
void fe_3d_trns_temporal(SELEM_3D *elem3d, SVECT *elem_nds, double djac, double *elem_c, double dt_factor, double *elem_rhs, char *string,
                         double perturbation, int perturb_node, int perturb_var, int perturb_sign);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 3D general transport body contributions to the elemental residual.
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
 * \details Solves the following weak, discrete body terms of the diffusive wave equation: \n
 * \f[  R_{i,body}^{\,e} = \bodyTime{\,3d}{e}{\phiddd{i}}{\, \ctrns{h}{j}{trns}} \,-
 *                     \bodyConv{\,3d}{e}{\phiddd{i}}{ (\velc{h}{trns} \, \ctrns{h}{j}{trns}) } \,+
 *                     \bodyDiffusion{\,3d}{e}{\phiddd{i}}{\diffTensor{trns}{h}{\diffTensorTRN}{j,}} \,-
 *                     \bodySource{\,3d}{e}{\phiddd{i}}{\srcGen{trns}{h}{j,}} \,+
 *                     \bodySUPG{\supg{i}{trns}{e}}
 * \f]
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_3d_transport_body_resid(SMODEL *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    
    int i;
    
    // DO NOT TAMPER WITH +++++++++++++++++++++++++++
    // for time-dependent debugging
    if (mod->t_prev > DEBUG_TIME) DEBUG_LOCAL = ON;
    
    // for node-dependent debugging
    for (i=0; i<mod->grid->elem3d[ie].nnodes; i++) {
        if (mod->grid->elem3d[ie].nodes[i] == DEBUG_NODE_ID) {
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
    g = mod->gravity;
    dt = mod->dt;
    alpha = mod->tau_temporal;
    
    // decide which hydro to use
    SSW_3D *sw = NULL;
    SNS_3D *ns = NULL;
    SMAT_SW mat_SW;
    SMAT_NS mat_NS;
    
    imat = mod->grid->elem3d[ie].mat;
    if (mod->flag.SW_FLOW == ON) {
        sw = mod->sw->d3;
        mat_SW = *(mod->mat[imat].sw);
    } else if (mod->flag.NS_FLOW == ON) {
        ns = mod->ns->d3;
        mat_NS = *(mod->mat[imat].ns);
    }
    
    // decide with transport constituent to use
    int index = UNSET_INT;
    SMAT_TRN mat_trans;
    SCON *con = NULL;
#ifdef _SEDIMENT
    SSUSLOAD *susload = NULL;
#endif
    double *c = NULL, *c_old = NULL, *c_older = NULL, *source = NULL, *sink = NULL, *mfcf = NULL;
    SVECT2D *vcf = NULL;
    if (mod->is_sediment_running) {
#ifdef SEDIMENT
        index = mod->ised;
        susload = &(mod->sed->susload[index]);
        mat_trans = mod->mat[imat].sed[index];
        c = susload->c;
        c_old = susload->c_old;
        c_older = susload->c_older;
        source = susload->source;
        sink = susload->sink;
        mfcf = susload->mfcf;
        vcf = susload->vcf;
#endif
    } else {
        index = mod->itrns;
        con = &(mod->con[index]);
        mat_trans = mod->mat[imat].trn[index];
        c = con->concentration;
        c_old = con->old_concentration;
        c_older = con->older_concentration;
        source = con->source;
        sink = con->sink;
        mfcf = con->mfcf;
        vcf = con->vcf;
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // GRID VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    grid = mod->grid;
    elem3d = &(grid->elem3d[ie]);
    
    icol = UNSET_INT;
    elem2d_sur = NULL;
    elem2d_bed = NULL;
    if (grid->type == COLUMNAR) {
        icol = grid->elem3d[ie].icol;
        elem2d_sur = &(grid->elem2d[grid->elem2d_sur[icol]]);
        elem2d_bed = &(grid->elem2d[grid->elem2d_bed[icol]]);
    }
    
    nnodes = elem3d->nnodes;
    nnodes_quad = elem3d->nnodes_quad;
    isElementTetrahedron = TRUE;
    if (nnodes == NDONPRISM) isElementTetrahedron = FALSE;
    for (i=0; i<nnodes; i++) {snode_copy(&(elem_nodes[i]), grid->node[elem3d->nodes[i]]);}
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    sarray_init_dbl(elem_c,NDONPRISM);
    global_to_local_dbl(c, elem_c, elem3d->nodes, nnodes);
    if (perturb_var == PERTURB_C) elem_c[perturb_node] += perturb_sign * perturbation;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // constituent concentrations
    sarray_init_dbl(elem_c_old,NDONPRISM);
    sarray_init_dbl(elem_c_older,NDONPRISM);
    global_to_local_dbl(c_old, elem_c_old, elem3d->nodes, nnodes);
    global_to_local_dbl(c_older, elem_c_older, elem3d->nodes, nnodes);
    
    // displacements
    double *dpl, *old_dpl, *older_dpl, *bed_dpl, *old_bed_dpl, *older_bed_dpl;
    if (mod->flag.SW_FLOW == ON) {
        dpl = sw->displacement; old_dpl = sw->old_displacement; older_dpl = sw->older_displacement;
    } else if (mod->flag.NS_FLOW == ON) {
        dpl = ns->displacement; old_dpl = ns->old_displacement; older_dpl = ns->older_displacement;
    }
    sarray_init_dbl(elem_dpl,NDONPRISM);        sarray_init_dbl(elem_bed_dpl,NDONPRISM);
    sarray_init_dbl(elem_dpl_old,NDONPRISM);    sarray_init_dbl(elem_bed_dpl_old,NDONPRISM);
    sarray_init_dbl(elem_dpl_older,NDONPRISM);  sarray_init_dbl(elem_bed_dpl_older,NDONPRISM);
    global_to_local_dbl(dpl,       elem_dpl,       elem3d->nodes, nnodes);
    global_to_local_dbl(old_dpl,   elem_dpl_old,   elem3d->nodes, nnodes);
    global_to_local_dbl(older_dpl, elem_dpl_older, elem3d->nodes, nnodes);
#ifdef _SEDIMENT
    // cjt :: add bed displacement to total displacement
    if (mod->flag.SW_FLOW == ON) {
        bed_dpl = sw->bed_displacement; old_bed_dpl = sw->old_bed_displacement; older_bed_dpl = sw->older_bed_displacement;
    } else (mod->flag.NS_FLOW == ON) {
        bed_dpl = ns->bed_displacement; old_bed_dpl = ns->old_bed_displacement; older_bed_dpl = ns->older_bed_displacement;
    }
    global_to_local_dbl(bed_dpl,       elem_bed_dpl,       elem3d->nodes, nnodes);
    global_to_local_dbl(old_bed_dpl,   elem_bed_dpl_old,   elem3d->nodes, nnodes);
    global_to_local_dbl(older_bed_dpl, elem_bed_dpl_older, elem3d->nodes, nnodes);
    sarray_add_replace_dbl(elem_dpl,       elem_bed_dpl,       nnodes);
    sarray_add_replace_dbl(elem_dpl_old,   elem_bed_dpl_old,   nnodes);
    sarray_add_replace_dbl(elem_dpl_older, elem_bed_dpl_older, nnodes);
#endif
    
    // nodal vectors with displacements
    svect_init_array(elem_nds,NDONPRISM);
    svect_init_array(elem_nds_old,NDONPRISM);
    svect_init_array(elem_nds_older,NDONPRISM);
    elem_get_midpt_locations2(elem_nodes, elem_dpl,       elem_nds,       nnodes, nnodes_quad, elem3d->edges);
    elem_get_midpt_locations2(elem_nodes, elem_dpl_old,   elem_nds_old,   nnodes, nnodes_quad, elem3d->edges);
    elem_get_midpt_locations2(elem_nodes, elem_dpl_older, elem_nds_older, nnodes, nnodes_quad, elem3d->edges);
    
    // velocities
    svect_init_array(elem_vel,NDONPRISM);
    if (mod->flag.SW_FLOW == ON) {
        global_to_local_svect(sw->vel, elem_vel, elem3d->nodes, nnodes);
    } else if (mod->flag.NS_FLOW == ON) {
        global_to_local_svect(ns->vel, elem_vel, elem3d->nodes, nnodes);
    }
    dumpVector(elem_vel, nnodes, elem_u, elem_v, elem_w);
#ifdef _SEDIMENT
    // Since sediment only exchanges with the bed at the bed nodes, we have to force it to fall above (cjt)
    if (mod->is_sediment_running) {
        if (mod->sed->grain[index].type == CLA) {
            for (inode=0; inode<nnodes; inode++) {
                elem_vel[inode].z -= mod->sed->grain[index].clay.settling_velocity;
                elem_w[inode] = elem_vel[inode].z;
            }
        }
    }
#endif
    
    // grid and relative velocities
    svect_init_array(elem_grid_vel, nnodes);
    svect_init_array(elem_rel_vel,NDONPRISM);
    for (i=0; i<nnodes; i++) {elem_grid_vel[i].z = (elem_dpl[i] - elem_dpl_old[i])/dt;}
    svect_subtract_array2(elem_rel_vel, elem_vel, elem_grid_vel, nnodes);
    dumpVector(elem_rel_vel, nnodes, elem_ur, elem_vr, elem_wr);
    
    // get baroclinic density
    sarray_init_value_dbl(elem_density, nnodes, 1.);
    if (mod->flag.BAROCLINIC == 1) {
        global_to_local_dbl(sw->density, elem_density, elem3d->nodes, nnodes);
    }
    // normalize density (only used in turbulence)
    if (mod->flag.BAROCLINIC > 0) {
        ELEM3D_SCALE(elem_density, 1. / mod->density);
    }
    
    // calculate tetrehedral jacobians and elemental averages
    elem_volume = 0.; elem_volume_old = 0.; elem_volume_older = 0.; elem_djac = 0.; elem_djac_old = 0.; elem_djac_older = 0.;
    elem_avg_c = 0.; elem_avg_c_old = 0.; elem_avg_c_older = 0.; elem_avg_density = 0.; elem_avg_u = 0.; elem_avg_v = 0.; elem_avg_w = 0.;
    elem_avg_ur = 0.; elem_avg_vr = 0.; elem_avg_wr = 0.;
    elem_avg_strong_resid = 0.; elem_avg_depth = 0.; elem_avg_vmag_safe = 0.;
    molecular_diffusion = 0.; turbulent_diffusion = 0.; vertical_diffusion = 0.; tur_viscosity = 0.; total_diffusion = 0.;
    dist_above_bed = 0.; average_bed_z = 0.; average_elem_z = 0.;
    svect_init_array(grad_phi,NDONPRISM); svect_init_array(grad_phi_old,NDONPRISM); svect_init_array(grad_phi_older,NDONPRISM);
    svect_init(&elem_grad_u); svect_init(&elem_grad_v); svect_init(&elem_grad_w); svect_init(&elem_grad_density);
    svect_init(&elem_avg_grad_ur); svect_init(&elem_avg_grad_vr); svect_init(&elem_avg_grad_wr);
    svect_init(&elem_avg_grad_density); svect_init(&elem_avg_grad_c);
    
    if (isElementTetrahedron == TRUE) { // djacs cancel here (djac = volume);
        quad = grid->quad_tet;
        
        elem_djac = get_tet_linear_djac_gradPhi2(NULL, elem_nds, grad_phi);
        elem_djac_old = get_tet_linear_djac_gradPhi2(NULL, elem_nds_old, grad_phi_old);
        elem_djac_older = get_tet_linear_djac_gradPhi2(NULL, elem_nds_older, grad_phi_older);
        
        elem_volume = elem_djac;
        
        grad_phi_dot_v(grad_phi, elem_vel, &elem_grad_u, &elem_grad_v, &elem_grad_w, nnodes);
        grad_phi_f(grad_phi, elem_c, &elem_grad_c, nnodes);
        grad_phi_f(grad_phi, elem_density, &elem_grad_density, nnodes);
        
        elem_avg_c = integrate_tetrahedron_f(1.,1.,elem_c);
        elem_avg_c_old = integrate_tetrahedron_f(1.,1.,elem_c_old);
        elem_avg_c_older = integrate_tetrahedron_f(1.,1.,elem_c_older);
        elem_avg_density = integrate_tetrahedron_f(1.,1.,elem_density);
        elem_avg_u = integrate_tetrahedron_f(1.,1.,elem_u);
        elem_avg_v = integrate_tetrahedron_f(1.,1.,elem_v);
        elem_avg_w = integrate_tetrahedron_f(1.,1.,elem_w);
        elem_avg_ur = integrate_tetrahedron_f(1.,1.,elem_ur);
        elem_avg_vr = integrate_tetrahedron_f(1.,1.,elem_vr);
        elem_avg_wr = integrate_tetrahedron_f(1.,1.,elem_wr);
        elem_avg_grad_ur = elem_grad_u;
        elem_avg_grad_vr = elem_grad_v;
        elem_avg_grad_wr = elem_grad_w;
        elem_avg_grad_density = elem_grad_density;
        elem_avg_grad_c = elem_grad_c;
    } else {
        quad = grid->quad_prism;
        elem_volume = get_triprism_volume(elem_nds);
        elem_volume_old = get_triprism_volume(elem_nds_old);
        elem_volume_older = get_triprism_volume(elem_nds_older);
        
        elem_avg_c =             integrate_triPrism_f(elem_nds, 1./elem_volume, elem_c);
        elem_avg_c_old =         integrate_triPrism_f(elem_nds_old, 1./elem_volume_old, elem_c_old);
        elem_avg_c_older =       integrate_triPrism_f(elem_nds_older, 1./elem_volume_older, elem_c_older);
        elem_avg_density =       integrate_triPrism_f(elem_nds, 1./elem_volume, elem_density);
        elem_avg_u =             integrate_triPrism_f(elem_nds, 1./elem_volume, elem_u);
        elem_avg_v =             integrate_triPrism_f(elem_nds, 1./elem_volume, elem_v);
        elem_avg_w =             integrate_triPrism_f(elem_nds, 1./elem_volume, elem_w);
        elem_avg_ur =            integrate_triPrism_f(elem_nds, 1./elem_volume, elem_ur);
        elem_avg_vr =            integrate_triPrism_f(elem_nds, 1./elem_volume, elem_vr);
        elem_avg_wr =            integrate_triPrism_f(elem_nds, 1./elem_volume, elem_wr);
        elem_avg_grad_ur =       integrate_triPrism_df_full(elem_nds, 1./elem_volume, elem_ur);
        elem_avg_grad_vr =       integrate_triPrism_df_full(elem_nds, 1./elem_volume, elem_vr);
        elem_avg_grad_wr =       integrate_triPrism_df_full(elem_nds, 1./elem_volume, elem_wr);
        elem_avg_grad_density =  integrate_triPrism_df_full(elem_nds, 1./elem_volume, elem_density);
        elem_avg_grad_c =        integrate_triPrism_df_full(elem_nds, 1./elem_volume, elem_c);
    }
    SVECT elem_avg_vel_rel;
    elem_avg_vel_rel.x = elem_avg_ur;
    elem_avg_vel_rel.y = elem_avg_vr;
    elem_avg_vel_rel.z = elem_avg_wr;
    double elem_avg_c_32dt = get_second_order(elem_avg_c, elem_avg_c_old);
    double elem_avg_c_12dt = get_second_order(elem_avg_c_old, elem_avg_c_older);
    
    elem_avg_strong_resid = (elem_avg_c_32dt - elem_avg_c_12dt)/dt + elem_avg_grad_c.x * elem_avg_ur + elem_avg_grad_c.y * elem_avg_vr + elem_avg_grad_c.z * elem_avg_wr;
    //double elem_avg_strong_resid = (elem_avg_c_32dt - elem_avg_c_12dt)/dt + (elem_avg_c * ( elem_avg_grad_ur.x +  elem_avg_grad_vr.y +  elem_avg_grad_wr.z) + svect_dotp(elem_avg_vel_rel,elem_avg_grad_c));
    
    // calculate depths for turbulence calculations
    if (mod->flag.SW_FLOW == ON) {
        elem_avg_depth = tl_find_avg_column_depth(grid, elem2d_sur->nodes, sw->displacement);
    } else {
        // CJT :: NEED A WAY TO GET AT ELEMENTAL DEPTH FOR UNSTRUCTURED FLOWS
    }
    elem_avg_vmag_safe = sqrt(pow(elem_avg_u,2) + pow(elem_avg_v,2) + pow(elem_avg_w,2) + NOT_QUITE_SMALL);
    
    // calculate the distance above the bed of the element centroid (cjt :: do it this way so we don't have to calculate pressure)
    if (mod->flag.SW_FLOW == ON) {
#ifdef _DEBUG
        assert(elem2d_bed->nnodes == NDONTRI);
#endif
        average_bed_z = one_3 * (grid->node[elem2d_bed->nodes[0]].z + grid->node[elem2d_bed->nodes[1]].z+ grid->node[elem2d_bed->nodes[2]].z);
        average_elem_z = 0;
        for (i=0; i<nnodes; i++) {average_elem_z += (elem_nds[i].z + elem_dpl[i]);}
        average_elem_z /= (double) nnodes;
        dist_above_bed = average_elem_z - average_bed_z;
    }
    
    // calculate the diffusion tensor
    molecular_diffusion = mod->viscosity;
    turbulent_diffusion = mat_trans.d_m;
    vertical_diffusion = mat_trans.d_m;
    tur_viscosity = 0;
    
    if (mod->flag.SW_FLOW == ON) {
        if (mat_SW.turbulence_model_xy == 0) {
            tur_viscosity = tur_smag( elem_avg_grad_ur.x,  elem_avg_grad_ur.y,  elem_avg_grad_vr.x,  elem_avg_grad_vr.y, elem_volume, mat_SW.smag_coeff);
        }
        if (mat_SW.turbulence_model_xy == 1) {
            tur_viscosity = tur_ws(elem_avg_vmag_safe, elem_avg_depth, mat_SW.smag_coeff);
        }
        if (mat_SW.turbulence_model_z == 1) {
            vertical_diffusion += tur_MY_2(g, elem_avg_depth, dist_above_bed,  elem_avg_grad_ur.z,  elem_avg_grad_vr.z, elem_avg_grad_density.z, elem_avg_density, mat_SW.supression_func, 2);
        }
    }
    turbulent_diffusion += tur_viscosity;
    total_diffusion = molecular_diffusion + turbulent_diffusion; //horizontal
    
    // Store elemental eddy viscosity
    grid->trn_diff[ie] = total_diffusion;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    assert(ie == elem3d->id);
    assert(alpha > -1E-6 && alpha < 1.1);
    assert(imat >= 0);
    assert(perturb_var == PERTURB_NONE || perturb_var == PERTURB_C);
    if (isElementTetrahedron == TRUE) {
        assert(elem_djac > 0);
        assert(elem_djac_old > 0);
        assert(elem_djac_older > 0);
    } else {
        assert(elem_volume > 0);
        assert(elem_volume_old > 0);
        assert(elem_volume_older > 0);
    }
    
    // if there is a NAN or INF, print arrays and possibly exit
    int print_debug_info = fe_3d_transport_body_test_for_NAN_or_INF();
    if (print_debug_info != NO) {
        tl_check_all_pickets(__FILE__,__LINE__);
        fe_3d_transport_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
    
    if (DEBUG_LOCAL == ON) fe_3d_transport_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    
#endif
    
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                  FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    double rhs[nnodes];
    sarray_init_dbl(elem_rhs, nnodes);
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                   TEMPORAL CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the temporal addition to the 3D transport elemental residual. \n
     * \note
     *********************************************************************************************/
    
    // ++ t(n+1) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_3d_trns_temporal(elem3d, elem_nds, elem_djac, elem_c, (1. + alpha / 2.), elem_rhs, "3D TRANSPORT || TEMPORAL T(N+1): ",
                        perturbation, perturb_node, perturb_var, perturb_sign);
    
    // ++ t(n) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_3d_trns_temporal(elem3d, elem_nds_old, elem_djac_old, elem_c_old, (-1.) * (1. + alpha), elem_rhs, "3D TRANSPORT || TEMPORAL T(N): ",
                        perturbation, perturb_node, perturb_var, perturb_sign);
    
    // ++ t(n-1) terms +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fe_3d_trns_temporal(elem3d, elem_nds_older, elem_djac_older, elem_c_older, (alpha / 2.), elem_rhs, "3D TRANSPORT || TEMPORAL T(N-1): ",
                        perturbation, perturb_node, perturb_var, perturb_sign);
    
    // old AdH 3D way
    //    int inode;
    //    double elem_rhs1[nnodes]; sarray_init_dbl(elem_rhs1, nnodes);
    //    double elem_rhs2[nnodes]; sarray_init_dbl(elem_rhs2, nnodes);
    //    double elem_rhs3[nnodes]; sarray_init_dbl(elem_rhs3, nnodes);
    //    fe_3d_trns_temporal(elem3d, elem_nds, elem_djac, elem_c, 1, elem_rhs1, "3D TRANSPORT || TEMPORAL T(N+1): ", DEBUG_LOCAL);
    //    fe_3d_trns_temporal(elem3d, elem_nds_old, elem_djac_old, elem_c_old, 1, elem_rhs2, "3D TRANSPORT || TEMPORAL T(N): ", DEBUG_LOCAL);
    //    fe_3d_trns_temporal(elem3d, elem_nds_older, elem_djac_older, elem_c_older, 1, elem_rhs3, "3D TRANSPORT || TEMPORAL T(N-1): ", DEBUG_LOCAL);
    //    double new_vol_c[nnodes]; sarray_init_dbl(new_vol_c, nnodes);
    //    double old_vol_c[nnodes]; sarray_init_dbl(old_vol_c, nnodes);
    //    ELEM3D_GET_SECONDORDER(new_vol_c, elem_rhs1, elem_rhs2);
    //    ELEM3D_GET_SECONDORDER(old_vol_c, elem_rhs2, elem_rhs3);
    //    for (inode=0; inode<nnodes; inode++) {
    //        elem_rhs[inode] += new_vol_c[inode] - old_vol_c[inode];
    //    }
    
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                     ADVECTION CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the advection addition to the 3D transport elemental residual. \n
     * \note
     *********************************************************************************************/
    
    // Galerkin contribution
    sarray_init_dbl(rhs, nnodes);
    if (isElementTetrahedron == TRUE) {
        integrate_tetrahedron_gradPhi_dot_f_v(grad_phi, elem_djac, -dt, elem_c, elem_rel_vel, rhs);
    } else {
        integrate_triPrism_gradPhi_dot_f_v(elem_nds, -dt, elem_c, elem_rel_vel, rhs);
        
        //                int iqp, quad_order = 3; // quadrature order
        //                SQUAD_PT *qp = NULL;
        //                for (iqp=0; iqp<quad[quad_order].n; iqp++) {
        //                    qp = &(quad[quad_order].pt[iqp]);
        //
        //                    // evalute triangular prism djac and cartesian shape function gradients at quadrature point
        //                    qp->djac = get_triprism_linear_djac_gradPhi(qp->xhat, qp->yhat, qp->zhat, elem_nds, qp->grad_shp);
        //
        //                    // evaluate transport equation at quadrature points
        //                    double qp_c = SQUAD_get_function(qp, elem_c, nnodes);
        //                    double qp_u = SQUAD_get_function(qp, elem_u, nnodes);
        //                    double qp_v = SQUAD_get_function(qp, elem_v, nnodes);
        //                    double qp_w = SQUAD_get_function(qp, elem_w, nnodes);
        //                    double t1 = -dt * qp->djac * qp->w * qp_c;
        //                    for (i=0; i<nnodes; i++) {
        //                        rhs[i] += t1 * (qp->grad_shp[i].x * qp_u + qp->grad_shp[i].y * qp_v + qp->grad_shp[i].z * qp_w);
        //                    }
        //                }
    }
    for (i=0; i<nnodes; i++) {elem_rhs[i] += rhs[i];}
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        rhs_1dof("3D TRANSPORT || ADVECTION: ",nnodes, ie, elem3d->nodes, rhs);
    }
    
    if (Is_DoubleArray_Inf_or_NaN_noExit(rhs,nnodes,__FILE__ ,__LINE__) != NO) {
        tl_check_all_pickets(__FILE__,__LINE__);
        rhs_1dof("3D TRANSPORT || ADVECTION: ", nnodes, elem3d->id, elem3d->nodes, rhs);
        printf("\n");
        fe_3d_transport_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                 DIFFUSION CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the diffusion addition to the 3D transport elemental residual. \n
     * \note SUPG does not contribute here, since first order elements are used
     *********************************************************************************************/
    STENSOR d_total; stensor_init(&d_total);
    d_total.xx = total_diffusion;
    d_total.yy = total_diffusion;
    d_total.zz = vertical_diffusion;
    d_total.xy = 0.0;
    d_total.xz = 0.0;
    d_total.yz = 0.0;
    
    
    sarray_init_dbl(rhs, nnodes);
    if (isElementTetrahedron == TRUE) {
        SVECT diffusive_flux; VT_3D_TENS_VECT_PROD(diffusive_flux, d_total, elem_grad_c);
        integrate_tetrahedron_gradPhi_dot_vcon(grad_phi, elem_djac, dt, diffusive_flux, rhs);
    } else {
        integrate_triPrism_gradPhi_dot_DgradF(elem_nds, quad, dt, d_total, elem_c, rhs);
        
        //        double rhs_x[NDONPRISM], rhs_y[NDONPRISM], rhs_z[NDONPRISM];
        //        SVECT diffusive_flux; VT_3D_TENS_VECT_PROD(diffusive_flux, d_total, elem_avg_grad_c);
        //        integrate_triPrism_dphi(elem_nds, dt * diffusive_flux.x, rhs_x, 1);
        //        integrate_triPrism_dphi(elem_nds, dt * diffusive_flux.y, rhs_y, 2);
        //        integrate_triPrism_dphi(elem_nds, dt * diffusive_flux.z, rhs_z, 3);
        //        for (i=0; i<nnodes; i++) rhs[i] = (rhs_x[i] + rhs_y[i] + rhs_z[i]);
    }
    for (i=0; i<nnodes; i++) {elem_rhs[i] += rhs[i];}
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        printf("diffusion tensor: xx: %20.10f yy: %20.10f zz: %20.10f \n",d_total.xx, d_total.yy, d_total.zz);
        rhs_1dof("3D TRANSPORT || DIFFUSION - GALERKIN: ",nnodes, ie, elem3d->nodes, rhs);
    }
    
    if (Is_DoubleArray_Inf_or_NaN_noExit(rhs,nnodes,__FILE__ ,__LINE__) != NO) {
        tl_check_all_pickets(__FILE__,__LINE__);
        printf("diffusion tensor: xx: %20.10f yy: %20.10f zz: %20.10f \n",d_total.xx, d_total.yy, d_total.zz);
        rhs_1dof("3D TRANSPORT || DIFFUSION - GALERKIN: ", nnodes, elem3d->id, elem3d->nodes, rhs);
        printf("\n");
        fe_3d_transport_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    SUPG CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the streamline upwind Petrov-Galerkin addition to the 3D transport elemental residual. \n
     * \note
     *********************************************************************************************/
    
    sarray_init_dbl(rhs, nnodes);
    double tau = 0.;
    if (isElementTetrahedron == TRUE) {
        
        double grad_phi_x[nnodes], grad_phi_y[nnodes], grad_phi_z[nnodes];
        for (i=0; i<nnodes; i++) {
            grad_phi_x[i] = grad_phi[i].x;
            grad_phi_y[i] = grad_phi[i].y;
            grad_phi_z[i] = grad_phi[i].z;
        }
        int tau_method_flag = 1; // flag for calculating tau :: Tezduyar and Park
        int le_method_flag = UNSET_INT; // flag for calculating elemental length :: not used yet
        tau = fe_get_supg_tau(nnodes, elem_nds, total_diffusion, elem_avg_ur, elem_avg_vr, elem_avg_wr,
                              grad_phi_x, grad_phi_y, grad_phi_z, elem_volume, mod->tau_pg, 3, tau_method_flag, le_method_flag);
        
        
        // direction integration :: assumes an average concentration and velocity in the strong residual
        integrate_tetrahedron_gradPhi_dot_v(grad_phi, elem_djac, dt * tau * elem_avg_strong_resid, elem_rel_vel, rhs);
        
        // use quadrature :: does not assume an average concentration or velocity
        //        int iqp, quad_order = 1; // quadrature order
        //        SQUAD_PT *qp = NULL;
        //        for (iqp=0; iqp<quad[quad_order].n; iqp++) {
        //            qp = &(quad[quad_order].pt[iqp]);
        //
        //            // evalute triangular prism djac and cartesian shape function gradients at quadrature point
        //            qp->djac = elem_volume;
        //
        //            // evaluate transport equation at quadrature points
        //            double qp_ur = 0., qp_vr = 0., qp_wr = 0., qp_c = 0., qp_c_old = 0., qp_c_older = 0.;
        //            qp_ur =      SQUAD_get_function(qp, elem_ur, nnodes);
        //            qp_vr =      SQUAD_get_function(qp, elem_vr, nnodes);
        //            qp_wr =      SQUAD_get_function(qp, elem_wr, nnodes);
        //            qp_c =       SQUAD_get_function(qp, elem_c, nnodes);
        //            qp_c_old =   SQUAD_get_function(qp, elem_c_old, nnodes);
        //            qp_c_older = SQUAD_get_function(qp, elem_c_older, nnodes);
        //            double qp_c32dt = get_second_order(qp_c, qp_c_old);
        //            double qp_c12dt = get_second_order(qp_c_old, qp_c_older);
        //            double qp_Rcont = (qp_c32dt - qp_c12dt)/dt + (elem_grad_c.x * qp_ur + elem_grad_c.y * qp_vr + elem_grad_c.y * qp_wr);;
        //            double t1 = tau * dt * qp->djac * qp->w * qp_Rcont;
        //            for (i=0; i<nnodes; i++) {
        //                rhs[i] += t1 * svect_dotp(elem_avg_vel_rel, qp->grad_shp[i]);
        //            }
        //        }
    } else {
        
        int tau_method_flag = 1; // flag for calculating tau :: Tezduyar and Park
        int le_method_flag = UNSET_INT; // flag for calculating elemental length :: not used yet
        tau = fe_get_supg_tau(nnodes, elem_nds, total_diffusion, elem_avg_ur, elem_avg_vr, elem_avg_wr,
                              NULL, NULL, NULL, elem_volume, mod->tau_pg, 3, tau_method_flag, le_method_flag);
#ifdef _DEBUG
        assert(tau > -1e-6);
#endif
        //tau *= 0.0001;
        
        // evalute triangular prism djac and cartesian shape function gradients at quadrature point
        integrate_triPrism_gradPhi_dot_v(elem_nds, dt * tau * elem_avg_strong_resid, elem_rel_vel, rhs);
        
        // use quadrature
        int iqp, quad_order = 3; // quadrature order
        SQUAD_PT *qp = NULL;
        
        for (iqp=0; iqp<quad[quad_order].n; iqp++) {
            qp = &(quad[quad_order].pt[iqp]);
            
            // evalute triangular prism djac and cartesian shape function gradients at quadrature point
            qp->djac = get_triprism_linear_djac_gradPhi(qp->xhat, qp->yhat, qp->zhat, elem_nds, qp->grad_shp);
            
            // evaluate continuity equation at quadrature points
            double qp_grad_vr = 0., qp_c = 0., qp_c_old = 0., qp_c_older = 0.;
            SVECT qp_grad_c, qp_vr;
            svect_init(&qp_grad_c); svect_init(&qp_vr);
            
            for (i=0; i<nnodes; i++) {
                qp_grad_vr += elem_rel_vel[i].x * qp->grad_shp[i].x;
                qp_grad_vr += elem_rel_vel[i].y * qp->grad_shp[i].y;
                qp_grad_vr += elem_rel_vel[i].z * qp->grad_shp[i].z;
                
                qp_grad_c.x += elem_c[i] * qp->grad_shp[i].x;
                qp_grad_c.y += elem_c[i] * qp->grad_shp[i].y;
                qp_grad_c.z += elem_c[i] * qp->grad_shp[i].z;
            }
            
            qp_vr.x =    SQUAD_get_function(qp, elem_ur, nnodes);
            qp_vr.y =    SQUAD_get_function(qp, elem_vr, nnodes);
            qp_vr.z =    SQUAD_get_function(qp, elem_wr, nnodes);
            qp_c =       SQUAD_get_function(qp, elem_c, nnodes);
            qp_c_old =   SQUAD_get_function(qp, elem_c_old, nnodes);
            qp_c_older = SQUAD_get_function(qp, elem_c_older, nnodes);
            double qp_c32dt = get_second_order(qp_c, qp_c_old);
            double qp_c12dt = get_second_order(qp_c_old, qp_c_older);
            double qp_Rcont = (qp_c32dt - qp_c12dt)/dt + (qp_c * qp_grad_vr + svect_dotp(qp_vr,qp_grad_c));
            double t1 = tau * dt * qp->djac * qp->w * qp_Rcont;
            for (i=0; i<nnodes; i++) {
                rhs[i] += t1 * svect_dotp(elem_avg_vel_rel, qp->grad_shp[i]);
                //rhs[i] += t1 * (qp_vr.x * qp->grad_shp[i].x + qp_vr.y * qp->grad_shp[i].y + qp_vr.z * qp->grad_shp[i].z);
            }
        }
    }
    for (i=0; i<nnodes; i++) {elem_rhs[i] += rhs[i];}
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        printf("tau: %20.10e\n",tau);
        rhs_1dof("3D TRANSPORT || SUPG: ",nnodes, ie, elem3d->nodes, rhs);
    }
    
    if (Is_DoubleArray_Inf_or_NaN_noExit(rhs,nnodes,__FILE__ ,__LINE__) != NO) {
        tl_check_all_pickets(__FILE__,__LINE__);
        printf("tau: %20.10e\n",tau);
        rhs_1dof("3D TRANSPORT || SUPG: ", nnodes, elem3d->id, elem3d->nodes, rhs);
        printf("\n");
        fe_3d_transport_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
#endif
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                             DISCONTINUITY CAPTURING CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the streamline upwind Petrov-Galerkin addition to the 3D transport elemental residual. \n
     * \note
     *********************************************************************************************/
    
    // DC contributions (CJT :: leave this out for now (we don't use)... come back to it later)
    
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                               EXTRA SW CONTINUITY TERM CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Adds any extra 3D SW depth-averged continuity and 3D continuity terms to the transport elemental residual. \n
     * note: CJT \:: already multiplied by dt
     *********************************************************************************************/
    if (mod->flag.SW_FLOW == ON && debug.no_hydro == OFF) {
        for (i=0; i<nnodes; i++) {
            rhs[i] = elem_avg_c * (mod->sw->d3->elem_rhs_supg_dacont[i][ie] + mod->sw->d3->elem_rhs_supg_cont[i][ie]);
            elem_rhs[i] += rhs[i];
        }
        
#ifdef _DEBUG
        if (DEBUG_LOCAL == ON) {
            rhs_1dof("3D TRANSPORT || SW CONTINUITY: ",nnodes, ie, elem3d->nodes, rhs);
            Is_DoubleArray_Inf_or_NaN(rhs ,nnodes ,__FILE__ ,__LINE__);
        }
        
        if (Is_DoubleArray_Inf_or_NaN_noExit(rhs,nnodes,__FILE__ ,__LINE__) != NO) {
            tl_check_all_pickets(__FILE__,__LINE__);
            rhs_1dof("3D TRANSPORT || SW CONTINUITY: ", nnodes, elem3d->id, elem3d->nodes, rhs);
            printf("\n");
            fe_3d_transport_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
            if (DEBUG_EXIT_ON_NAN) exit(-1);
        }
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                    SET DIRICHLET BCS
     *--------------------------------------------------------------------------------------------
     * Sets the Dirichlet boundary conditions. \n
     *********************************************************************************************/
    int istring = UNSET_INT;
    INPUT_DATA *str_value;
    for (i=0; i<nnodes; i++) {
        istring = grid->node[elem3d->nodes[i]].string;
        if (mod->is_sediment_running) {
            str_value = &(mod->str_values[istring].sed[mod->ised]);
        } else {
            str_value = &(mod->str_values[istring].trans[mod->itrns]);
        }
        if (istring > NORMAL) {
            if (str_value->bc_flag == BCT_DIR || str_value->bc_flag == BCT_CEQ) {
                elem_rhs[i] = 0.0;
            }
        }
    }
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        rhs_1dof("3D TRANSPORT || TOTAL RHS AFER DB: ",nnodes, ie, elem3d->nodes, elem_rhs);
        if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    }
    
    if (Is_DoubleArray_Inf_or_NaN_noExit(rhs,nnodes,__FILE__ ,__LINE__) != NO) {
        tl_check_all_pickets(__FILE__,__LINE__);
        rhs_1dof("3D TRANSPORT || TOTAL RHS AFER DB: ", nnodes, elem3d->id, elem3d->nodes, rhs);
        printf("\n");
        fe_3d_transport_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
    
    time_t time2;  time(&time2);
    TIME_IN_3D_TRANSPORT_BODY_RESID += difftime(time2,time1);
#endif
    
    return;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 3D transport body temporal contributions to the shallow water equations.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_3d_trns_temporal(SELEM_3D *elem3d, SVECT *elem_nds, double djac, double *elem_c, double dt_factor, double *elem_rhs, char *string,
                         double perturbation, int perturb_node, int perturb_var, int perturb_sign) {
    int i;
    int nnodes = elem3d->nnodes;
    double rhs[nnodes]; sarray_init_dbl(rhs, nnodes);
    if (nnodes == NDONTET) {
        integrate_tetrahedron_phi_f(djac, dt_factor, elem_c, rhs);
    } else {
        integrate_triPrism_phi_f(elem_nds, dt_factor, elem_c, rhs);
    }
    for (i=0; i<nnodes; i++) {elem_rhs[i] += rhs[i];}
    
#ifdef _DEBUG
    if (DEBUG_LOCAL == ON) {
        rhs_1dof(string, nnodes, elem3d->id, elem3d->nodes, rhs);
    }
    
    if (Is_DoubleArray_Inf_or_NaN_noExit(rhs,nnodes,__FILE__ ,__LINE__) != NO) {
        tl_check_all_pickets(__FILE__,__LINE__); printf("\n");
        rhs_1dof(string, nnodes, elem3d->id, elem3d->nodes, rhs);
        fe_3d_transport_body_debug(perturbation, perturb_node, perturb_var, perturb_sign);
        if (DEBUG_EXIT_ON_NAN) exit(-1);
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Checks 3D transport body variables for NANS or INFS.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int fe_3d_transport_body_test_for_NAN_or_INF() {
    
    int print_debug_info = NO;

    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_c, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_c has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_c_old, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_c_old has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_c_older, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_c_older has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_dpl, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_dpl has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_dpl_old, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_dpl_old has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_DoubleArray_Inf_or_NaN_noExit(elem_dpl_older, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_dpl_older has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_vel, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_ve; has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_grid_vel, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_grid_vel has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_rel_vel, nnodes ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_rel_vel has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_nds, nnodes_quad ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_nds has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_nds_old, nnodes_quad ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_nds_old has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_vectorArray_Inf_or_NaN_noExit(elem_nds_older, nnodes_quad ,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_nds_older has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_c,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_c has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_c_old,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_c_old has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_c_older,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_c_older has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_density,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_density has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_u,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_u has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_v,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_v has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_w,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_w has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_ur,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_ur has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_vr,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_vr has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_wr,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_wr has either a NAN or INF, printing debug info \n\n");
    
    print_debug_info = Is_Vector_Inf_or_NaN_noExit(elem_avg_grad_density,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_grad_density has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Vector_Inf_or_NaN_noExit(elem_avg_grad_c,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_grad_c has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(average_elem_z,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("average_elem_z has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(average_bed_z,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("average_bed_z has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(dist_above_bed,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("dist_above_bed has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_depth,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_depth has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_vmag_safe,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_vmag_safe has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(molecular_diffusion,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("molecular_diffusion has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(turbulent_diffusion,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("turbulent_diffusion has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(vertical_diffusion,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("vertical_diffusion has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(total_diffusion,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("total_diffusion has either a NAN or INF, printing debug info \n\n");
    print_debug_info = Is_Double_Inf_or_NaN_noExit(elem_avg_strong_resid,__FILE__ ,__LINE__);
    if (print_debug_info != NO) printf("elem_avg_strong_resid has either a NAN or INF, printing debug info \n\n");
    
    return print_debug_info;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Prints to screen DEBUG info for the 3D transport body residual.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_3d_transport_body_debug(double perturbation, int perturb_node, int perturb_var, int perturb_sign) {
    
    printf("SW-3D TRANSPORT BODY ELEM RESID :: ie: %d \t dt: %20.10f \t volume: %20.10f",elem3d->id,dt,elem_volume);
    if (perturb_var == PERTURB_C) {
        printf("\t perturbing C || node: %d || perturbation: %20.10e\n",elem3d->nodes[perturb_node],perturb_sign*perturbation);
    }
    if (isElementTetrahedron == TRUE) {
        printf("\ntetrahedral djacs :: djac: %30.20f \t old_djac: %30.20f \t older_djac: %30.20f \n",elem_djac,elem_djac_old,elem_djac_older);
    }
    //selem3d_printScreen(elem3d);
    printScreen_debug_vec("node locations: ",elem_nds, nnodes);
    printScreen_debug_vec("node locations old: ",elem_nds_old, nnodes);
    printScreen_debug_vec("node locations older: ",elem_nds_older, nnodes);
    printScreen_debug_vec("grad_phi: ",grad_phi, nnodes);
    printScreen_debug_vec("grad_phi_old: ",grad_phi_old, nnodes);
    printScreen_debug_vec("grad_phi_older: ",grad_phi_older, nnodes);
    printScreen_debug2_dbl("elem_c", elem_c, nnodes, elem3d->nodes);
    printScreen_debug2_dbl("elem_c_old", elem_c_old, nnodes, elem3d->nodes);
    printScreen_debug2_dbl("elem_c_older", elem_c_older, nnodes, elem3d->nodes);
    printScreen_debug2_dbl("elem_displacement", elem_dpl, nnodes, elem3d->nodes);
    printScreen_debug2_dbl("elem_displacement_old", elem_dpl_old, nnodes, elem3d->nodes);
    printScreen_debug2_dbl("elem_displacement_older", elem_dpl_older, nnodes, elem3d->nodes);
    printScreen_debug_svect("elem_vel", elem_vel, nnodes, elem3d->nodes);
    printScreen_debug_svect("elem_grid_vel", elem_grid_vel, nnodes, elem3d->nodes);
    printScreen_debug_svect("elem_rel_vel", elem_rel_vel, nnodes, elem3d->nodes);
    printScreen_debug2_dbl("elem_density", elem_density, nnodes, elem3d->nodes);
    printf("elem_avg_c: %20.10f \t elem_avg_c_old: %20.10f \t elem_avg_c_older: %20.10f\n",elem_avg_c,elem_avg_c_old,elem_avg_c_older);
    printf("elem_avg_density: %20.10f\n",elem_avg_density);
    printf("elem_avg_u: %20.10f \t elem_avg_v: %20.10f \t elem_avg_w: %20.10f\n",elem_avg_u, elem_avg_v, elem_avg_w);
    printf("elem_avg_grad_ur: ddx = %20.10f \t ddy = %20.10f \t ddz = %20.10f\n", elem_avg_grad_ur.x, elem_avg_grad_ur.y, elem_avg_grad_ur.z);
    printf("elem_avg_grad_vr: ddx = %20.10f \t ddy = %20.10f \t ddz = %20.10f\n", elem_avg_grad_vr.x, elem_avg_grad_vr.y, elem_avg_grad_vr.z);
    printf("elem_avg_grad_wr: ddx = %20.10f \t ddy = %20.10f \t ddz = %20.10f\n", elem_avg_grad_wr.x, elem_avg_grad_wr.y, elem_avg_grad_wr.z);
    printf("elem_avg_grad_density: ddx = %20.10f \t ddy = %20.10f \t ddz = %20.10f\n",elem_avg_grad_density.x,elem_avg_grad_density.y,elem_avg_grad_density.z);
    
    printf("average_elem_z: %20.10f \t average_bed_z: %20.10f \t dist_above_bed: %20.10f\n",average_elem_z,average_bed_z,dist_above_bed);
    printf("elem_avg_depth: %20.10f\n",elem_avg_depth);
    printf("elem_avg_vmag_safe: %20.10f\n",elem_avg_vmag_safe);
    printf("molecular_diffusion: %20.10f \t turbulent_diffusion: %20.10f \t vertical_diffusion: %20.10f \t total_diffusion: %20.10f \n",molecular_diffusion,turbulent_diffusion,vertical_diffusion,total_diffusion);
    printf("elem_avg_strong_resid: %20.10f\n",elem_avg_strong_resid);
}
