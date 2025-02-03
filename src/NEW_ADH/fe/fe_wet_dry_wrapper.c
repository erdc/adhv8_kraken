/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_wet_dry_wrapper.c This file collects the 2D shallow water wet/dry wrapper.
 *   All function arguements assumed to be consistent with function prototype.              */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the wet-dry integration factor :: essentially a stripped down version of fe_sw2_wet_dry_wrapper
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gary Brown
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \date      October, 2017
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in]  (required) x the 2d element nodal coordinates
 * @param[in]  (required) h the 2d element nodal depths - this is wet-dry averaged
 * @param[in]  (required) djac the 2d element area
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double fe_sw2_wet_dry_factor(SVECT *x, double *h, double djac) {
    
    int i, is_dry, total_vtx, passcode, only_vtx;
    double elem_djac, factor;
    
    passcode = 0; total_vtx = 0;
    
    for(i=0; i<3; i++) {
        is_dry = (h[i] <= 0.0) ? 1 : 0;
        total_vtx += is_dry;
        passcode += (is_dry << i);
    }
    
    if(total_vtx == 3) return 0.0;
    else if(total_vtx == 0) return 1.0;
    
    if(total_vtx != 1) passcode = (7 - passcode);
    
    only_vtx = passcode / 2;
    
    double xsi, new_h[3];
    SVECT2D new_x[3];
    for(i=0; i<3; i++) {
        if (i != only_vtx) {
            xsi = h[i] / (h[i] - h[only_vtx]);
            FIXED_POS(x[i].x, x[only_vtx].x, xsi, new_x[i].x);
            FIXED_POS(x[i].y, x[only_vtx].y, xsi, new_x[i].y);
            new_h[i] = 0;
        } else {
            new_x[i].x = x[i].x;
            new_x[i].y = x[i].y;
            new_h[i] = h[i];
        }
    }
    
    elem_djac = djac;
    TRI_AREA(new_x[0], new_x[1], new_x[2], djac);
    
    if(djac / elem_djac < SMALL) return (double)(total_vtx % 2);
    factor = (total_vtx >> 1) - (total_vtx % 2);
    factor *= djac / elem_djac;
    return factor + ((factor < 0) ? 1.0 : 0.0);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Wet/Dry wrapper for the 2D depth-averaged continuity equation temporal integral
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gary Brown
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \date      October, 2017
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elem_rhs the 2D elemental residual array contribution with wet/dry inclusion
 * @param[in]  (required) x the 2d element nodal coordinates
 * @param[in]  (required) h the 2d element nodal depths - this is wet-dry averaged
 * @param[in]  (optional) grad_phi the elemental shape function gradients in configuration space, NULL if not included
 * @param[in]  (optional) v the element nodal velocities, NULL if not included
 * @param[in]  (optional) grad_z the elemental bathymetry gradients, NULL if not included
 * @param[in]  (optional) f generic scalar function which is integrated over the whole element
 * @param[in]  (optional) f_wet_dry a generic scalar function which is wet-dry integrated
 * @param[in]  (required) djac the 2d element area
 * @param[in]  (required) dt timestep
 * @param[in]  (required) c an integration constant
 * @param[in]  (required) void(*fe_sw_func) a function that performs the consisitent mass integration
 *
 * \note CJT \:: only handles triangles right now
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double fe_sw2_wet_dry_wrapper(double *elem_rhs, SVECT *x, double *h, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, int redistribute_flag, int DEBUG, double *vars, void (*fe_sw_func) (SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *elem_rhs)) {
    
#ifdef _DEBUG
    if (DEBUG) {
        assert(elem_rhs != NULL);
        assert(x != NULL);
        assert(h != NULL);
    }
#endif
    
    int i, is_dry, total_vtx, passcode, only_vtx;
    
    // here's where the magic happens - passcode is a number that will be used
    // in both its binary and decimal representations to flag which vertices are dry
    // total_vtx is the number of dry vertices (either 0, 1, 2, or 3)
    passcode = 0;
    total_vtx = 0;
    
    // checking to see if the vertices are dry or not
    // if a vertex is dry, increment total_vtx and set the corresponding
    // bit of the binary representation of passcode to '1' by using the binary
    // shift operator
    for (i = 0; i < NDONTRI; i++) {
        is_dry = (h[i] <= 0.0) ? 1 : 0;
        total_vtx += is_dry;
        passcode += (is_dry << i);
    }
    
    // if total_vtx = 3, we're all dry and no integration needs to take place
    // if total vtx = 0, we're all wet and the entire element may be integrated as normal
    // otherwise (total_vtx=1 or 2), we can prepare to create a new element by performing
    // mathematical magic using the decimal representation of passcode to find the 'unique' element
    // i.e., that which is the only one of its kind (dry or wet) in the original element
    
    if (total_vtx == NDONTRI) {
#ifdef _DEBUG
        if (DEBUG) printf("fe_sw2_wet_dry_wrapper :: element fully dry\n");
#endif
        return 0;
    } else if (total_vtx == 0) {
        fe_sw_func(x, h, grad_phi, v, v_wet_dry, f, f_wet_dry, djac, vars, elem_rhs);
#ifdef _DEBUG
        if (DEBUG) printf("fe_sw2_wet_dry_wrapper :: element fully wet\n");
#endif
        return 1.0;
    }
    
    
    // if the new element is completely dry (total_vtx=1), we integrate the entire original element
    // so that we may subtract the contribution of the new dry element later
    // CJT :: don't mean nearly dry, not completely here?
    // CJT :: should this not be checked for 1 *or* 2 ??
    if (total_vtx == 1) {
        fe_sw_func(x, h, grad_phi, v, v_wet_dry, f, f_wet_dry, djac, vars, elem_rhs);
    } else {
        passcode = (7 - passcode);
    }
    only_vtx = passcode / 2;
    
    
    // at this point integration has been performed such that elem_rhs = 0 (total_vtx = 3) or its value
    // the entire original element had been wet (total_vtx = 0 or 1)
    // our next task is to find the remaining 2 vertices of the new element
    // these will lie along the sides of the original element between only_vtx and each of the other two
    // vertices at the point where h=0
    double xsi, new_h[NDONTRI] = {0.0, 0.0, 0.0};
    SVECT new_x[NDONTRI];
    double *new_f_wet_dry  = NULL;
    if (f_wet_dry != NULL) {
        new_f_wet_dry = (double *) tl_alloc(sizeof(double), NDONTRI);
        for (i=0; i<NDONTRI; i++) new_f_wet_dry[i] = 0.0;
    }
    SVECT2D *new_v_wet_dry = NULL;
    if (v_wet_dry != NULL) {
        new_v_wet_dry = (SVECT2D *) tl_alloc(sizeof(SVECT2D), NDONTRI);
        for (i=0; i<NDONTRI; i++) {new_v_wet_dry[i].x = 0.0; new_v_wet_dry[i].y = 0.0;}
    }
    
    
    // loop through the vertices and find only_vtx then calculate position of the 2 new verticies
    for (i = 0; i < NDONTRI; i++) {
        if (i != only_vtx) {
            
            // interpolate depth
            xsi = h[i] / (h[i] - h[only_vtx]);
            
            // get wet-dry coordinates
            FIXED_POS(x[i].x, x[only_vtx].x, xsi, new_x[i].x);
            FIXED_POS(x[i].y, x[only_vtx].y, xsi, new_x[i].y);
            new_x[i].z = x[i].z;
            
            // get wet-dry vector
            if (v_wet_dry != NULL) {
                FIXED_POS(v_wet_dry[i].x, v_wet_dry[only_vtx].x, xsi, new_v_wet_dry[i].x);
                FIXED_POS(v_wet_dry[i].y, v_wet_dry[only_vtx].y, xsi, new_v_wet_dry[i].y);
            }
            
            // get wet-dry function
            if (f_wet_dry != NULL) FIXED_POS(f_wet_dry[i], f_wet_dry[only_vtx], xsi, new_f_wet_dry[i]);
            
        } else {
            new_x[i] = x[i];
            new_h[i] = h[i];
            if (f_wet_dry != NULL) new_f_wet_dry[i] = f_wet_dry[i];
            if (v_wet_dry != NULL) {
                new_v_wet_dry[i].x = v_wet_dry[i].x;
                new_v_wet_dry[i].y = v_wet_dry[i].y;
            }
        }
    }
#ifdef _DEBUG
    if (DEBUG) {
        printf("fe_sw2_wet_dry_wrapper :: new_x[0] :: %30.20e %30.20e %30.20e \n",new_x[0].x,new_x[0].y,new_x[0].z);
        printf("fe_sw2_wet_dry_wrapper :: new_x[1] :: %30.20e %30.20e %30.20e \n",new_x[1].x,new_x[1].y,new_x[1].z);
        printf("fe_sw2_wet_dry_wrapper :: new_x[2] :: %30.20e %30.20e %30.20e \n",new_x[2].x,new_x[2].y,new_x[2].z);
        printf("fe_sw2_wet_dry_wrapper :: new_h[0] :: %30.20e %30.20e %30.20e \n",new_h[0],new_h[0],new_h[0]);
        if (new_f_wet_dry != NULL) {
            printf("fe_sw2_wet_dry_wrapper :: new_f_wet_dry[0] :: %30.20e %30.20e %30.20e \n",new_f_wet_dry[0],new_f_wet_dry[0],new_f_wet_dry[0]);
        }
        if (new_v_wet_dry != NULL) {
            printf("fe_sw2_wet_dry_wrapper :: new_v_wet_dry[0] :: %30.20e %30.20e\n",new_v_wet_dry[0].x,new_v_wet_dry[0].y);
            printf("fe_sw2_wet_dry_wrapper :: new_v_wet_dry[1] :: %30.20e %30.20e\n",new_v_wet_dry[1].x,new_v_wet_dry[1].y);
            printf("fe_sw2_wet_dry_wrapper :: new_v_wet_dry[2] :: %30.20e %30.20e\n",new_v_wet_dry[2].x,new_v_wet_dry[2].y);
        }
    }
#endif
    
    // compute the area integration of the tri composed by new_x
    // - CJT :: it is possible to have only a very small % of the element wet, causing gradient errors for this routine.
    // - CJT :: in the case, we just assume the element is completely dry (doesn't this give mass errors?)
    double elem_djac = djac; // store this to get area ratio
    TRI_AREA(new_x[0], new_x[1], new_x[2], djac);
    if (djac / elem_djac < SMALL) {
        if (new_f_wet_dry != NULL) new_f_wet_dry = (double *)  tl_free(sizeof(double),  NDONTRI, new_f_wet_dry);
        if (new_v_wet_dry != NULL) new_v_wet_dry = (SVECT2D *) tl_free(sizeof(SVECT2D), NDONTRI, new_v_wet_dry);
        return (double) (total_vtx % 2);
    }
    
    // get new shape function gradients using new jacobian
    SVECT2D *new_grad_phi = NULL;
    if (grad_phi != NULL) {
        new_grad_phi = (SVECT2D *) tl_alloc(sizeof(SVECT2D), NDONTRI);
        GRAD_SHAPE(new_x[0], new_x[1], new_x[2], djac, new_grad_phi);
    }
    
    // Now integrate using wet/dry values
    double new_elem_rhs[NDONTRI*3] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    fe_sw_func(new_x, new_h, new_grad_phi, v, new_v_wet_dry, f, new_f_wet_dry, djac, vars, new_elem_rhs);
    
#ifdef _DEBUG
    if (DEBUG) {
        if (new_f_wet_dry != NULL) {
            printf("fe_sw2_wet_dry_wrapper :: new_f_wet_dry[0] :: %30.20e %30.20e %30.20e \n",new_f_wet_dry[0],new_f_wet_dry[0],new_f_wet_dry[0]);
        }
        if (new_grad_phi != NULL) {
            printf("fe_sw2_wet_dry_wrapper :: new_grad_phi[0] :: %30.20e %30.20e \n",new_grad_phi[0].x,new_grad_phi[0].y);
            printf("fe_sw2_wet_dry_wrapper :: new_grad_phi[1] :: %30.20e %30.20e \n",new_grad_phi[1].x,new_grad_phi[1].y);
            printf("fe_sw2_wet_dry_wrapper :: new_grad_phi[2] :: %30.20e %30.20e \n",new_grad_phi[2].x,new_grad_phi[2].y);
        }
        printf("fe_sw2_wet_dry_wrapper :: new_elem_rhs.c_eq before redistribution :: %30.20e %30.20e %30.20e \n",new_elem_rhs[0],new_elem_rhs[3],new_elem_rhs[6]);
        printf("fe_sw2_wet_dry_wrapper :: new_elem_rhs.x_eq before redistribution :: %30.20e %30.20e %30.20e \n",new_elem_rhs[1],new_elem_rhs[4],new_elem_rhs[7]);
        printf("fe_sw2_wet_dry_wrapper :: new_elem_rhs.y_eq before redistribution :: %30.20e %30.20e %30.20e \n",new_elem_rhs[2],new_elem_rhs[5],new_elem_rhs[8]);
    }
#endif
    
    double factor = (total_vtx >> 1) - (total_vtx % 2);
    if (redistribute_flag == ON) {    /* neither for avg term nor for PG */
        RE_DISTRIBUTE(x, new_x, new_elem_rhs, only_vtx);
    }
    
#ifdef _DEBUG
    if (DEBUG) {
        printf("fe_sw2_wet_dry_wrapper :: factor before redistribution: %30.20e \n",factor);
        printf("fe_sw2_wet_dry_wrapper :: new_elem_rhs.c_eq after redistribution :: %30.20e %30.20e %30.20e \n",new_elem_rhs[0],new_elem_rhs[3],new_elem_rhs[6]);
        printf("fe_sw2_wet_dry_wrapper :: new_elem_rhs.x_eq after redistribution :: %30.20e %30.20e %30.20e \n",new_elem_rhs[1],new_elem_rhs[4],new_elem_rhs[7]);
        printf("fe_sw2_wet_dry_wrapper :: new_elem_rhs.y_eq after redistribution :: %30.20e %30.20e %30.20e \n",new_elem_rhs[2],new_elem_rhs[5],new_elem_rhs[8]);
    }
#endif
    
    for (i = 0; i < NDONTRI; i++) {
        elem_rhs[i*3] += factor * new_elem_rhs[i*3];
        elem_rhs[i*3+1] += factor * new_elem_rhs[i*3+1];
        elem_rhs[i*3+2] += factor * new_elem_rhs[i*3+2];
    }
    factor *= djac / elem_djac;
    
#ifdef _DEBUG
    if (DEBUG) {
        printf("fe_sw2_wet_dry_wrapper :: elem_rhs.c_eq :: %30.20e %30.20e %30.20e \n",elem_rhs[0],elem_rhs[3],elem_rhs[6]);
        printf("fe_sw2_wet_dry_wrapper :: elem_rhs.x_eq :: %30.20e %30.20e %30.20e \n",elem_rhs[1],elem_rhs[4],elem_rhs[7]);
        printf("fe_sw2_wet_dry_wrapper :: elem_rhs.y_eq :: %30.20e %30.20e %30.20e \n",elem_rhs[2],elem_rhs[5],elem_rhs[8]);
        printf("fe_sw2_wet_dry_wrapper :: factor returned: %30.20e \n",factor + ((factor < 0) ? 1.0 : 0.0) );
    }
#endif
    
    if (new_f_wet_dry != NULL) new_f_wet_dry = (double *)  tl_free(sizeof(double),  NDONTRI, new_f_wet_dry);
    if (new_v_wet_dry != NULL) new_v_wet_dry = (SVECT2D *) tl_free(sizeof(SVECT2D), NDONTRI, new_v_wet_dry);
    if (new_grad_phi  != NULL) new_grad_phi =  (SVECT2D *) tl_free(sizeof(SVECT2D), NDONTRI, new_grad_phi);
    
    return factor + ((factor < 0) ? 1.0 : 0.0);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Wet/Dry wrapper for the 2D depth-averaged continuity equation temporal integral
 *  \author    Mark Loveland
 *
 * @param[inout] elem_rhs the 2D elemental residual array contribution with wet/dry inclusion
 * @param[in]  (required) x the 2d element nodal coordinates
 * @param[in]  (required) h the 2d element nodal depths - this is wet-dry averaged
 * @param[in]  (optional) grad_phi the elemental shape function gradients in configuration space, NULL if not included
 * @param[in]  (optional) v the element nodal velocities, NULL if not included
 * @param[in]  (optional) grad_z the elemental bathymetry gradients, NULL if not included
 * @param[in]  (optional) f generic scalar function which is integrated over the whole element
 * @param[in]  (optional) f_wet_dry a generic scalar function which is wet-dry integrated
 * @param[in]  (required) djac the 2d element area
 * @param[in]  (required) dt timestep
 * @param[in]  (required) c an integration constant
 * @param[in]  (required) void(*fe_sw_func) a function that performs the consisitent mass integration
 *
 * 
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double fe_sw2_wet_dry_wrapper_1d(double *elem_rhs, SVECT *x, double *h, SVECT2D *nrml, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, int redistribute_flag, int DEBUG, double *vars, void (*fe_sw_func) (SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *elem_rhs)) {
    
#ifdef _DEBUG
    if (DEBUG) {
        assert(elem_rhs != NULL);
        assert(x != NULL);
        assert(h != NULL);
    }
#endif
    
    int i, is_dry, total_vtx, passcode, only_vtx;
    //hardcode for now, change later
    int nnodes = 2;
    // here's where the magic happens - passcode is a number that will be used
    // in both its binary and decimal representations to flag which vertices are dry
    // total_vtx is the number of dry vertices (either 0, 1, 2, or 3)
    passcode = 0;
    total_vtx = 0;
    
    // checking to see if the vertices are dry or not
    // if a vertex is dry, increment total_vtx and set the corresponding
    // bit of the binary representation of passcode to '1' by using the binary
    // shift operator
    for (i = 0; i < nnodes; i++) {
        is_dry = (h[i] <= 0.0) ? 1 : 0;
        total_vtx += is_dry;
        passcode += (is_dry << i);
    }
    // passcode is 0 if both are wer
    // passcode is 1 if node 0 is dry and node 1 is wet
    // passcode is 2 if node 0 is wet and node 1 is dry
    // passcode is 3 if both nodes are dry 
    // if total_vtx = 2, we're all dry and no integration needs to take place
    // if total vtx = 0, we're all wet and the entire element may be integrated as normal
    // otherwise (total_vtx=1 ), we can prepare to create a new element by performing
    // mathematical magic using the decimal representation of passcode to find the 'unique' element
    // i.e., that which is the only one of its kind (dry or wet) in the original element
    
    if (total_vtx == nnodes) {
#ifdef _DEBUG
        if (DEBUG) printf("fe_sw2_wet_dry_wrapper :: element fully dry\n");
#endif
        return 0;
    } else if (total_vtx == 0) {
        fe_sw_func(x, h, nrml, v, v_wet_dry, f, f_wet_dry, djac, vars, elem_rhs);
#ifdef _DEBUG
        if (DEBUG) printf("fe_sw2_wet_dry_wrapper :: element fully wet\n");
#endif
        return 1.0;
    }
    
    
    // if the new element has one dry node (total_vtx=1), we integrate the entire original element
    // so that we may subtract the contribution of the new dry element later
    //for 1D i dont think we will need this
    //if (total_vtx == 1) {
    //    fe_sw_func(x, h, grad_phi, v, v_wet_dry, f, f_wet_dry, djac, vars, elem_rhs);
    //} else {
    //    passcode = (7 - passcode);
    //}

    //identifies wet node. either 0 or 1
    only_vtx = passcode % 2;
    
    
    // at this point integration has been performed such that elem_rhs = 0 (total_vtx = 2) or whole wet
    // our next task is to find the remaining 2 vertices of the new element
    // these will lie along the sides of the original element between only_vtx and each of the other two
    // vertices at the point where h=0
    double xsi;
    double new_h[nnodes];
    new_h[0]=0;
    new_h[1]=0;
    SVECT new_x[nnodes];
    double *new_f_wet_dry  = NULL;
    if (f_wet_dry != NULL) {
        new_f_wet_dry = (double *) tl_alloc(sizeof(double), nnodes);
        for (i=0; i<nnodes; i++) new_f_wet_dry[i] = 0.0;
    }
    SVECT2D *new_v_wet_dry = NULL;
    if (v_wet_dry != NULL) {
        new_v_wet_dry = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
        for (i=0; i<nnodes; i++) {new_v_wet_dry[i].x = 0.0; new_v_wet_dry[i].y = 0.0;}
    }
    
    
    // loop through the vertices and find only_vtx then calculate position of the 2 new verticies
    for (i = 0; i < nnodes; i++) {
        if (i != only_vtx) {
            
            // interpolate depth
            xsi = h[i] / (h[i] - h[only_vtx]);
            
            // get wet-dry coordinates
            FIXED_POS(x[i].x, x[only_vtx].x, xsi, new_x[i].x);
            FIXED_POS(x[i].y, x[only_vtx].y, xsi, new_x[i].y);
            new_x[i].z = x[i].z;
            
            // get wet-dry vector
            if (v_wet_dry != NULL) {
                FIXED_POS(v_wet_dry[i].x, v_wet_dry[only_vtx].x, xsi, new_v_wet_dry[i].x);
                FIXED_POS(v_wet_dry[i].y, v_wet_dry[only_vtx].y, xsi, new_v_wet_dry[i].y);
            }
            
            // get wet-dry function
            if (f_wet_dry != NULL) FIXED_POS(f_wet_dry[i], f_wet_dry[only_vtx], xsi, new_f_wet_dry[i]);
            
        } else {
            new_x[i] = x[i];
            new_h[i] = h[i];
            if (f_wet_dry != NULL) new_f_wet_dry[i] = f_wet_dry[i];
            if (v_wet_dry != NULL) {
                new_v_wet_dry[i].x = v_wet_dry[i].x;
                new_v_wet_dry[i].y = v_wet_dry[i].y;
            }
        }
    }
#ifdef _DEBUG
    if (DEBUG) {
        printf("fe_sw2_wet_dry_wrapper :: new_x[0] :: %30.20e %30.20e %30.20e \n",new_x[0].x,new_x[0].y,new_x[0].z);
        printf("fe_sw2_wet_dry_wrapper :: new_x[1] :: %30.20e %30.20e %30.20e \n",new_x[1].x,new_x[1].y,new_x[1].z);
        printf("fe_sw2_wet_dry_wrapper :: new_x[2] :: %30.20e %30.20e %30.20e \n",new_x[2].x,new_x[2].y,new_x[2].z);
        printf("fe_sw2_wet_dry_wrapper :: new_h[0] :: %30.20e %30.20e %30.20e \n",new_h[0],new_h[0],new_h[0]);
        if (new_f_wet_dry != NULL) {
            printf("fe_sw2_wet_dry_wrapper :: new_f_wet_dry[0] :: %30.20e %30.20e %30.20e \n",new_f_wet_dry[0],new_f_wet_dry[0],new_f_wet_dry[0]);
        }
        if (new_v_wet_dry != NULL) {
            printf("fe_sw2_wet_dry_wrapper :: new_v_wet_dry[0] :: %30.20e %30.20e\n",new_v_wet_dry[0].x,new_v_wet_dry[0].y);
            printf("fe_sw2_wet_dry_wrapper :: new_v_wet_dry[1] :: %30.20e %30.20e\n",new_v_wet_dry[1].x,new_v_wet_dry[1].y);
            printf("fe_sw2_wet_dry_wrapper :: new_v_wet_dry[2] :: %30.20e %30.20e\n",new_v_wet_dry[2].x,new_v_wet_dry[2].y);
        }
    }
#endif
    
    // compute the area integration of the tri composed by new_x
    // - CJT :: it is possible to have only a very small % of the element wet, causing gradient errors for this routine.
    // - CJT :: in the case, we just assume the element is completely dry (doesn't this give mass errors?)
    double elem_djac = djac; // store this to get area ratio

    //This is the triangle version
    //TRI_AREA(new_x[0], new_x[1], new_x[2], djac);
    //Mark, is this just area or is it inverse of area?
    djac=DIST_2D(new_x[0], new_x[1]);
    
    if (djac / elem_djac < SMALL) {
        if (new_f_wet_dry != NULL) new_f_wet_dry = (double *)  tl_free(sizeof(double),  NDONTRI, new_f_wet_dry);
        if (new_v_wet_dry != NULL) new_v_wet_dry = (SVECT2D *) tl_free(sizeof(SVECT2D), NDONTRI, new_v_wet_dry);
        return (double) (total_vtx % 2);
    }
    

    //Mark, do we need to do this?? Dont think so
    // get new shape function gradients using new jacobian
    //SVECT2D *new_grad_phi = NULL;
    //    if (grad_phi != NULL) {
    //        new_grad_phi = (SVECT2D *) tl_alloc(sizeof(SVECT2D), NDONTRI);
    //        GRAD_SHAPE(new_x[0], new_x[1], new_x[2], djac, new_grad_phi);
    //    }
    
    // Now integrate using wet/dry values
    double new_elem_rhs[nnodes*3];
    sarray_init_dbl(new_elem_rhs, nnodes*3);
    fe_sw_func(new_x, new_h, nrml, v, new_v_wet_dry, f, new_f_wet_dry, djac, vars, new_elem_rhs);
    
#ifdef _DEBUG
    if (DEBUG) {
        if (new_f_wet_dry != NULL) {
            printf("fe_sw2_wet_dry_wrapper :: new_f_wet_dry[0] :: %30.20e %30.20e %30.20e \n",new_f_wet_dry[0],new_f_wet_dry[0],new_f_wet_dry[0]);
        }
        printf("fe_sw2_wet_dry_wrapper :: new_elem_rhs.c_eq before redistribution :: %30.20e %30.20e %30.20e \n",new_elem_rhs[0],new_elem_rhs[3],new_elem_rhs[6]);
        printf("fe_sw2_wet_dry_wrapper :: new_elem_rhs.x_eq before redistribution :: %30.20e %30.20e %30.20e \n",new_elem_rhs[1],new_elem_rhs[4],new_elem_rhs[7]);
        printf("fe_sw2_wet_dry_wrapper :: new_elem_rhs.y_eq before redistribution :: %30.20e %30.20e %30.20e \n",new_elem_rhs[2],new_elem_rhs[5],new_elem_rhs[8]);
    }
#endif
    
    double factor = 1;//(total_vtx >> 1) - (total_vtx % 2);
    //this will need to possibly be fixed, or new routine added
    if (redistribute_flag == ON) {    /* neither for avg term nor for PG */
        RE_DISTRIBUTE_1D(x, new_x, new_elem_rhs, only_vtx);
    }
    
#ifdef _DEBUG
    if (DEBUG) {
        printf("fe_sw2_wet_dry_wrapper :: factor before redistribution: %30.20e \n",factor);
        printf("fe_sw2_wet_dry_wrapper :: new_elem_rhs.c_eq after redistribution :: %30.20e %30.20e %30.20e \n",new_elem_rhs[0],new_elem_rhs[3],new_elem_rhs[6]);
        printf("fe_sw2_wet_dry_wrapper :: new_elem_rhs.x_eq after redistribution :: %30.20e %30.20e %30.20e \n",new_elem_rhs[1],new_elem_rhs[4],new_elem_rhs[7]);
        printf("fe_sw2_wet_dry_wrapper :: new_elem_rhs.y_eq after redistribution :: %30.20e %30.20e %30.20e \n",new_elem_rhs[2],new_elem_rhs[5],new_elem_rhs[8]);
    }
#endif
    
    for (i = 0; i < NDONTRI; i++) {
        elem_rhs[i*3] += factor * new_elem_rhs[i*3];
        elem_rhs[i*3+1] += factor * new_elem_rhs[i*3+1];
        elem_rhs[i*3+2] += factor * new_elem_rhs[i*3+2];
    }
    factor *= djac / elem_djac;
    
#ifdef _DEBUG
    if (DEBUG) {
        printf("fe_sw2_wet_dry_wrapper :: elem_rhs.c_eq :: %30.20e %30.20e %30.20e \n",elem_rhs[0],elem_rhs[3],elem_rhs[6]);
        printf("fe_sw2_wet_dry_wrapper :: elem_rhs.x_eq :: %30.20e %30.20e %30.20e \n",elem_rhs[1],elem_rhs[4],elem_rhs[7]);
        printf("fe_sw2_wet_dry_wrapper :: elem_rhs.y_eq :: %30.20e %30.20e %30.20e \n",elem_rhs[2],elem_rhs[5],elem_rhs[8]);
        printf("fe_sw2_wet_dry_wrapper :: factor returned: %30.20e \n",factor + ((factor < 0) ? 1.0 : 0.0) );
    }
#endif
    
    if (new_f_wet_dry != NULL) new_f_wet_dry = (double *)  tl_free(sizeof(double),  nnodes, new_f_wet_dry);
    if (new_v_wet_dry != NULL) new_v_wet_dry = (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, new_v_wet_dry);
    //if (new_grad_phi  != NULL) new_grad_phi =  (SVECT2D *) tl_free(sizeof(SVECT2D), nnodes, new_grad_phi);
    
    return factor + ((factor < 0) ? 1.0 : 0.0);
}