#include "global_header.h"

static double TOLERANCE = 1E-8; // for point inside checks

double orient2d(const double pa[], const double pb[], const double pc[]) {
    const double acx = pa[0] - pc[0];
    const double bcx = pb[0] - pc[0];
    const double acy = pa[1] - pc[1];
    const double bcy = pb[1] - pc[1];
    return acx * bcy - acy * bcx;
}

// LINEAR ELEMENTS -------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------

double interpolate1D(double x, double xL, double xR, double yL, double yR) {
    double denom = (xR - xL);
    if (fabs(denom) > 1e-6) {
        return yL + (x - xL) * (yR - yL) / denom;
    }
    return -9999998.0;
}

#define UNPACK_SVECT_TO_1D_ARRAY(_in) { _in.x, _in.y }

int compute_interpolate_weights_2D_triangle(SGRID *grid, int ie, double x, double y, double *weights, int *lnd, int *iedge) {
    
    
    int i;
    const double node[3][2] = {
        UNPACK_SVECT_TO_1D_ARRAY(grid->node[grid->elem2d[ie].nodes[0]]),
        UNPACK_SVECT_TO_1D_ARRAY(grid->node[grid->elem2d[ie].nodes[1]]),
        UNPACK_SVECT_TO_1D_ARRAY(grid->node[grid->elem2d[ie].nodes[2]])
    };
    const double p0[2] = {x, y};
    
    *lnd = UNSET_INT;    // ID of local node if point falls on element node
    *iedge = UNSET_INT;  // ID of local edge if the point falls on element edge

    // i-th barycentric coordinate <==> orient2d(node_0 of opposite edge, node_1 of opposite edge, p0) + normalization (below)
    for (i=0; i<3; i++)
        weights[i] = orient2d(node[(i+1)%3], node[(i+2)%3], p0);

    // normalization for barycentric coordinates
    const double s = weights[0] + weights[1] + weights[2];
    for (i=0; i<3; i++)
        weights[i] /= s;

    //printf("weights: %20.10f %20.10f %20.10f\n",weights[0],weights[1],weights[2]);

    // is point on element?
    if (weights[0] > -TOLERANCE &&  weights[0] < 1+TOLERANCE &&
        weights[1] > -TOLERANCE &&  weights[1] < 1+TOLERANCE &&
        weights[2] > -TOLERANCE &&  weights[2] < 1+TOLERANCE) {

        // is point on a node?
        for (i=0; i<3; i++) {
            if (fabs(weights[i] - 1)<TOLERANCE) {
                *lnd = i;
                return 1; // point found in element
            }
        }

        // is point on edge?
        for (i=0; i<3; i++) {
            if (fabs(weights[grid->nd_on_TriEdge[i][0]] + weights[grid->nd_on_TriEdge[i][1]] - 1) < TOLERANCE) {
                *iedge = i;
                return 1; // point found in element
            }
        }

        // otherwise node is in interior
        return 1; // point found in element
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
//    // works for large grids!
//    //printf("xL: %20.10f \t yL: %20.10f\n",grid->xL,grid->yL);
//    x /= grid->xL;
//    y /= grid->yL;
//
//    double nd0_x = node[0].x/grid->xL; double nd0_y = node[0].y/grid->yL;
//    double nd1_x = node[1].x/grid->xL; double nd1_y = node[1].y/grid->yL;
//    double nd2_x = node[2].x/grid->xL; double nd2_y = node[2].y/grid->yL;
//
//
//    double x12 = nd0_x - nd1_x, x23 = nd1_x - nd2_x, x31 = nd2_x - nd0_x;
//    double y12 = nd0_y - nd1_y, y23 = nd1_y - nd2_y, y31 = nd2_y - nd0_y;
//    double djac = nd0_x*y23 + nd1_x*y31 + nd2_x*y12;
//
//    weights[0] = (y23*x - x23*y + nd1_x*nd2_y - nd2_x*nd1_y)/djac;
//    weights[1] = (y31*x - x31*y + nd2_x*nd0_y - nd0_x*nd2_y)/djac;
//
////    weights[0] = y23*x/djac - x23*y/djac + nd1_x*nd2_y/djac - nd2_x*nd1_y/djac;
////    weights[1] = y31*x/djac - x31*y/djac + nd2_x*nd0_y/djac - nd0_x*nd2_y/djac;
////    TOLERANCE /= djac;
//
//    weights[2] = 1 - weights[0] - weights[1];
//    //printf("weights: %20.10f %20.10f %20.10f || djac: %20.10f \n",weights[0],weights[1],weights[2],djac);
//
//    // is point on element?
//    if (weights[0] > -TOLERANCE &&  weights[0] < 1+TOLERANCE &&
//        weights[1] > -TOLERANCE &&  weights[1] < 1+TOLERANCE &&
//        weights[2] > -TOLERANCE &&  weights[2] < 1+TOLERANCE) {
//
//        // is point on a node?
//        for (i=0; i<3; i++) {
//            if (fabs(weights[i] - 1)<TOLERANCE) {
//                *lnd = i;
//                return 1; // point found in element
//            }
//        }
//
//        // is point on edge?
//        for (i=0; i<3; i++) {
//            if (fabs(weights[grid->nd_on_TriEdge[i][0]] + weights[grid->nd_on_TriEdge[i][1]] - 1) < TOLERANCE) {
//                *iedge = i;
//                return 1; // point found in element
//            }
//        }
//
//        // otherwise node is in interior
//        return 1; // point found in element
//     }
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
//        printf("xL: %20.10f \t yL: %20.10f\n",grid->xL,grid->yL);
//        x /= grid->xL;
//        y /= grid->yL;
//
//        double nd0_x = node[0].x/grid->xL; double nd0_y = node[0].y/grid->yL;
//        double nd1_x = node[1].x/grid->xL; double nd1_y = node[1].y/grid->yL;
//        double nd2_x = node[2].x/grid->xL; double nd2_y = node[2].y/grid->yL;
//
//
//        double x12 = nd0_x - nd1_x, x23 = nd1_x - nd2_x, x31 = nd2_x - nd0_x;
//        double y12 = nd0_y - nd1_y, y23 = nd1_y - nd2_y, y31 = nd2_y - nd0_y;
//        double djac = nd0_x*y23 + nd1_x*y31 + nd2_x*y12;
//
//        weights[0] = (y23*x - x23*y + nd1_x*nd2_y - nd2_x*nd1_y);
//        weights[1] = (y31*x - x31*y + nd2_x*nd0_y - nd0_x*nd2_y);
//        weights[2] = djac - weights[0] - weights[1];
//
//        printf("weights: %20.10f %20.10f %20.10f\n",weights[0],weights[1],weights[2]);
//
//        // is point on element?
//        if (weights[0] > -TOLERANCE * djac &&  weights[0] < (1+TOLERANCE) * djac &&
//            weights[1] > -TOLERANCE * djac &&  weights[1] < (1+TOLERANCE) * djac &&
//            weights[2] > -TOLERANCE * djac &&  weights[2] < (1+TOLERANCE) * djac) {
//
//            // is point on a node?
//            for (i=0; i<3; i++) {
//                if (fabs(weights[i] - djac)<TOLERANCE) {
//                    *lnd = i;
//                    return 1; // point found in element
//                }
//            }
//
//            // is point on edge?
//            for (i=0; i<3; i++) {
//                if (fabs(weights[grid->nd_on_TriEdge[i][0]] + weights[grid->nd_on_TriEdge[i][1]] - djac) < TOLERANCE) {
//                    *iedge = i;
//                    return 1; // point found in element
//                }
//            }
//
//            // otherwise node is in interior
//            return 1; // point found in element
//         }
    
    return 0; // point not found in element, on element nodes or on element edges
}

int compute_interpolate_weights_3D_tetrahedron(SGRID *grid, int ie, double x, double y, double z, double *weights, int *lnd, int *iedge, int *iface) {
    
    int i, k1, k2, k3, k4;
    double djac = 0., a[4], b[4], c[4], d[4];
    SVECT node[4];
    for (i=0; i<4; i++) node[i] = grid->node[grid->elem3d[ie].nodes[i]]; // shallow copy
    
    *lnd = UNSET_INT;    // ID of local node if point falls on element node
    *iedge = UNSET_INT;  // ID of local edge if the point falls on element edge
    
    for (i=0; i<4; i++) {
        if (i == 0) {
            k1 = 1; k2 = 2; k3 = 3;
        } else if (i == 1) {
            k1 = 0; k2 = 2; k3 = 3;
        } else if (i == 2) {
            k1 = 0; k2 = 1; k3 = 3;
        } else {
            k1 = 0; k2 = 1; k3 = 2;
        }
        a[i] = pow((-1.0),(i+1+1)) * (node[k1].x*node[k2].y*node[k3].z +
                                      node[k1].y*node[k3].x*node[k2].z +
                                      node[k2].x*node[k3].y*node[k1].z -
                                      node[k3].x*node[k2].y*node[k1].z -
                                      node[k3].y*node[k1].x*node[k2].z -
                                      node[k2].x*node[k1].y*node[k3].z);
        b[i] = pow((-1.0),(i+1))   * (node[k1].y*node[k2].z +
                                      node[k2].y*node[k3].z +
                                      node[k3].y*node[k1].z -
                                      node[k3].y*node[k2].z -
                                      node[k2].y*node[k1].z -
                                      node[k1].y*node[k3].z );
        c[i] = pow((-1.0),(i+1+1)) * (node[k1].x*node[k2].z +
                                      node[k2].x*node[k3].z +
                                      node[k3].x*node[k1].z -
                                      node[k3].x*node[k2].z -
                                      node[k2].x*node[k1].z -
                                      node[k1].x*node[k3].z );
        d[i] = pow((-1.0),(i+1))   * (node[k1].x*node[k2].y +
                                      node[k2].x*node[k3].y +
                                      node[k3].x*node[k1].y -
                                      node[k3].x*node[k2].y -
                                      node[k2].x*node[k1].y -
                                      node[k1].x*node[k3].y);
        
        djac = djac + a[i];
        weights[i] = a[i] + b[i]*x + c[i]*y + d[i]*z;
    }
    for (i=0; i<4; i++) {
        weights[i]=weights[i]/djac;
    }
    
    // is point on element?
    if (weights[0] > -TOLERANCE &&  weights[0] < 1+TOLERANCE &&
        weights[1] > -TOLERANCE &&  weights[1] < 1+TOLERANCE &&
        weights[2] > -TOLERANCE &&  weights[2] < 1+TOLERANCE &&
        weights[3] > -TOLERANCE &&  weights[3] < 1+TOLERANCE ) {
        
        // is point on a node?
        for (i=0; i<4; i++) {
            if (fabs(weights[i] - 1)<TOLERANCE) {
                *lnd = i;
                return 1; // point found in element
            }
        }
        
        // is point on edge?
        for (i=0; i<6; i++) {
            if (fabs(weights[grid->nd_on_TetEdge[i][0]] + weights[grid->nd_on_TetEdge[i][1]] - 1) < TOLERANCE) {
                *iedge = i;
                return 1; // point found in element
            }
        }
        
        // is point on face?
        for (i=0; i<4; i++) {
            if (fabs(weights[nd_on_fc[i][0]] + weights[nd_on_fc[i][1]] + weights[nd_on_fc[i][2]] - 1) < TOLERANCE) {
                *iface = i;
                return 1; // point found in element
            }
        }
        
        // otherwise node is in interior
        return 1; // point found in element
     }
    
    return 0; // point not found in element, on element nodes or on element edges
}

// BILINEAR ELEMENTS -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------
