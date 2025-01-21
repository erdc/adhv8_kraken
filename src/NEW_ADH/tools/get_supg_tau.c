/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_get_supg_tau.c This file collections functions responsible for
 *          calculating the SUPG coefficient, tau and element length comprising this.       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "adh.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the SUPG tau parameter using the spectral radius for the shallow water equations
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * \returns (double) tau
 *
 * @param[in] nnodes (int) the number of nodes on the element
 * @param[in] nodes (SVECT *) the elemental node locations
 * @param[in] g (double) gravity
 * @param[in] elem_avg_depth (double) the elementally averaged depth
 * @param[in] elem_avg_u (double) the elementally averaged x-velocity
 * @param[in] elem_avg_v (double) the elementally averaged y-velocity
 * @param[in] elem_avg_w (double) the elementally averaged z-velocity
 * @param[in] grad_shp_x (double *) the x-gradient in cartesian space of the basis functions
 * @param[in] grad_shp_y (double *) the y-gradient in cartesian space of the basis functions
 * @param[in] grad_shp_z (double *) the z-gradient in cartesian space of the basis functions
 * @param[in] volume_or_area (double) the element volume/area
 * @param[in] alpha (double) a user-input SUPG scalar ranging from 0 to 0.5
 * @param[in] ndim (int) the dimensionality of the element
 * @param[in] tau_method_flag (int) a flag to determine which tau method to use
 * @param[in] le_method_flag (int) a flag to determine which element length method to use
 *
 * \details Elemental lengths are generally calculated using the Tezduyar and Park \cite TezduyarPark1986
 *          or through the elemental area/volumes.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double fe_get_supg_tau_sw(int nnodes, SVECT *nodes, double g, double elem_avg_depth, double elem_avg_u, double elem_avg_v,
                                 double elem_avg_w,double *grad_shp_x, double *grad_shp_y, double *grad_shp_z, double djac,
                                 double alpha, int ndim, int tau_method_flag, int le_method_flag) {
    
    // estimate elemental length
    double le = get_element_length(nnodes, nodes, elem_avg_u,elem_avg_v,elem_avg_w,grad_shp_x,grad_shp_y,grad_shp_z,djac,ndim,le_method_flag); // UNITS: [L]
    
    // shallow water wave celerity
    double celerity_squared = g * elem_avg_depth;
    if (elem_avg_depth <= 0.0) celerity_squared = 0.;
    
    // get eigenvalue matrix spectral radius
    double spectral_radius = 0;
    if (ndim == 2) {
        spectral_radius = sqrt(pow(elem_avg_u,2) + pow(elem_avg_v,2) + celerity_squared + SMALL);
    } else {
        spectral_radius = sqrt(pow(elem_avg_u,2) + pow(elem_avg_v,2) + pow(elem_avg_w,2) + celerity_squared + SMALL);
        
        // to exactly match prism branch
        double elem_avg_v_mag = sqrt(pow(elem_avg_u,2) + pow(elem_avg_v,2) + pow(elem_avg_w,2));
        spectral_radius = sqrt(g * elem_avg_depth + pow(elem_avg_v_mag,2));
    }
    
    // get tau
    double tau = alpha * le / spectral_radius;
    
    //printf("le: %20.10e \t alpha: %20.10e \t spectral_radius: %20.10e \t celerity_squared: %20.10e \t elem_avg_u: %20.10e \t elem_avg_v: %20.10e \t elem_avg_w: %20.10e \n",le,alpha,spectral_radius,celerity_squared,elem_avg_u,elem_avg_v,elem_avg_w);
    
#if _DEBUG
    if (tau < -SMALL) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        printf("alpha: %20.10f \t le: %20.10f \t spectral_radius: %20.10f \t tau: %20.10f\n",alpha,le,spectral_radius,tau);
        tl_error("tau < -SMALL");
    }
#endif
    
    return tau;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns an element length estimate for SUPG calculations
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * \returns (double) the element length
 *
 * @param[in] nnodes (int) the number of nodes on the element
 * @param[in] nodes (SVECT *) the elemental node locations
 * @param[in] elem_avg_u (double) the elementally averaged x-velocity 
 * @param[in] elem_avg_v (double) the elementally averaged y-velocity
 * @param[in] elem_avg_w (double) the elementally averaged z-velocity
 * @param[in] grad_shp_x (double *) the x-gradient in cartesian space of the basis functions
 * @param[in] grad_shp_y (double *) the y-gradient in cartesian space of the basis functions
 * @param[in] grad_shp_z (double *) the z-gradient in cartesian space of the basis functions
 * @param[in] volume_or_area (double) the element volume/area
 * @param[in] ndim (int) the dimensionality of the element
 * @param[in] method_flag (int) a flag to determine which method to use
 *
 * \details Elemental lengths are generally calculated using the Tezduyar and Park \cite TezduyarPark1986
 *          or through the elemental area/volumes.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double get_element_length(int nnodes, SVECT *node, double elem_avg_u, double elem_avg_v, double elem_avg_w,
                          double *grad_shp_x, double *grad_shp_y, double *grad_shp_z, double volume_or_area,
                          int ndim, int method_flag) {
   
    int i;
    double le = 0.;
    
#if _DEBUG
    assert(ndim > 0 && ndim < 4);
#endif
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (ndim == 1) {
#if _DEBUG
        //assert(method_flag == 1 || method_flag == 2);
#endif
#if _DEBUG
        assert(grad_shp_x != NULL);
#endif
        // Tezduyar and Park element length
        // h = (1/2) * (sum[ (u/u_mag) dot grad_phi_i ]) = (1/2) * (u_mag / sum[ u dot grad_phi_i])
        // for 3D continuity calculations, set elem_avg_u = 1, to force upwind direction
        
        double sum = 0.;
        for (i=0; i<nnodes; i++) {
            sum += fabs(elem_avg_u * grad_shp_x[i]);
        }
        le = one_2 * fabs(elem_avg_u) / (sum + 2. * NOT_QUITE_SMALL);
#if _DEBUG
        assert(le > 0);
#endif
        return le;
        
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    } else if (ndim == 2) {
        
#if _DEBUG
        assert(method_flag == 1 || method_flag == 2);
#endif
        
        if (method_flag == 1) {
            
            // Tezduyar and Park element length
            // h = (1/2) * (sum[ (u/u_mag) dot grad_phi_i ]) = (1/2) * (u_mag / sum[ u dot grad_phi_i])
#if _DEBUG
            assert(grad_shp_x != NULL);
            assert(grad_shp_y != NULL);
#endif
            
            double flow_norm = sqrt(elem_avg_u * elem_avg_u + elem_avg_v * elem_avg_v);
            double sum = 0.;
            for (i=0; i<nnodes; i++) {
                sum += fabs(elem_avg_u * grad_shp_x[i] + elem_avg_v * grad_shp_y[i]);
            }
            le = one_2 * flow_norm / (sum + 2. * NOT_QUITE_SMALL);
            
        } else if (method_flag == 2) {
            // Area based element length
            // CJT \:: here we assume an equilateral triangle (length *not* in flow direction)
            // A = (1/2) * le_x * le_y = (1/2) * le^2 --> le = sqrt( 2*A)
            le = sqrt(volume_or_area); //sqrt(2 * volume_or_area); // CJT :: the latter is more correct, but does not match old trunk
        }
#if _DEBUG
        //assert(le > 0);
#endif
        return le;
        
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    } else if (ndim == 3) {
        
#if _DEBUG
        //assert(method_flag == 1 || method_flag == 2);
#endif
        
        if (nnodes == NDONTET) { // Tetrahedral element
            
            // Tezduyar and Park element length
            // h = (1/2) * (sum[ (u/u_mag) dot grad_phi_i ]) = (1/2) * (u_mag / sum[ u dot grad_phi_i])
#if _DEBUG
            assert(grad_shp_x != NULL);
            assert(grad_shp_y != NULL);
            assert(grad_shp_z != NULL);
#endif
            double flow_norm = sqrt(elem_avg_u * elem_avg_u + elem_avg_v * elem_avg_v + elem_avg_w * elem_avg_w) + 2. * NOT_QUITE_SMALL;
            double sum = 0.;
            for (i=0; i<nnodes; i++) {
                sum += fabs(elem_avg_u * grad_shp_x[i] + elem_avg_v * grad_shp_y[i] + elem_avg_w * grad_shp_z[i]);
            }
            le = one_2 * flow_norm / (sum + 2. * NOT_QUITE_SMALL);
            
            // cap the element length with the edge length of a regular tet's volume
            double len = pow(6 * volume_or_area, one_3);
            if (le > len) le = len;
#if _DEBUG
            if (le < 1e-6) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                printf("le: %20.10e \t flow_norm: %20.10e \t volume: %20.10e \t elem_avg_u: %20.10e \t elem_avg_v: %20.10e \t elem_avg_w: %20.10e\n",le,flow_norm,volume_or_area,elem_avg_u,elem_avg_v,elem_avg_w);
                tl_error(">> SUPG length calculation is <= 0\n");
            }
#endif
        } else { // triangular prism element
            
#if _DEBUG
            assert(node[0].x == node[3].x); assert(node[1].x == node[4].x); assert(node[2].x == node[5].x);
            assert(node[0].y == node[3].y); assert(node[1].y == node[4].y); assert(node[2].y == node[5].y);
#endif
            
            // cjt :: for prism cell measure, take the average of cross element diagonal lengths || there are a billion ways to do this ...
            // cjt :: this does not include displacement for now
            
            //            5
            //            o
            //           /|\
            //          / | \
            //       3 o-----o 4
            //         |  |  |
            //         |  |  |
            //         |  o  |
            //         | /2\ |
            //         |/   \|
            //       0 o-----o 1
            
            double L_15 = sqrt( pow(node[0].x - node[4].x,2) + pow(node[0].y - node[4].y,2) + pow(node[0].z - node[4].z,2)  );
            double L_26 = sqrt( pow(node[1].x - node[5].x,2) + pow(node[1].y - node[5].y,2) + pow(node[1].z - node[5].z,2)  );
            double L_34 = sqrt( pow(node[2].x - node[3].x,2) + pow(node[2].y - node[3].y,2) + pow(node[2].z - node[3].z,2)  );
#ifdef _DEBUG
            assert(L_15 > 1e-6);
            assert(L_26 > 1e-6);
            assert(L_34 > 1e-6);
#endif
            le = one_3 * (L_15 + L_26 + L_34);
            
            // cjt :: assume prism has equilateral triangles, so that volume = (1/2) * le_x * le_y * le_z = (1/2) * le^2 * le_z
            // cjt :: assume most of the time, flow is across the element laterally, so element length is flow direction is this way
            // cjt :: V = (1/2) * le^2 * le_z --> le = sqrt( 2 * V / le_z );
            double le_z = fabs(node[0].z - node[3].z);
            le = pow(2. * volume_or_area / le_z, one_2);
            
            // cap the element length with the edge length of a regular prisms volume
            double len = pow(2. * volume_or_area, one_3);
            if (le > len) le = len;
            
            //le = pow(5. * volume_or_area, one_3); // old adh code
        }
#if _DEBUG
        if (le < 0) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            printf("volume_or_area: %20.10e\n",volume_or_area);
            tl_error(">> le <= 0 \n");
        }
#endif
        return le;
    }
    
    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
    tl_error(">> bad ndim sent\n");

    return (UNSET_INT);
}