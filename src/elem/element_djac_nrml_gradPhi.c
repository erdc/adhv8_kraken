
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  elem2d3d_linear.c This file collects functions that calculate elemental normals, djacs, shape and shape gradients. */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define SCALE_FOR_MASTER_TET_VOLUME (1.0/6.0) // tetrahedral djac to volume scale
#define SCALE_FOR_MASTER_TRI_VOLUME (1.0/2.0) // tetrahedral djac to volume scale

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* 1D ELEMENTS */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

inline void get_elem1d_linear_djac_gradPhi(SGRID *grid, SELEM_1D *elem1d) {
    
    /* sets the nodes on the element */
    int nd1 = elem1d->nodes[0];
    int nd2 = elem1d->nodes[1];
    
    /* computes the jacobian */
    elem1d->djac = DIST_2D(grid->node[nd1], grid->node[nd2]);
    
    /* computes the gradient of the shape functions */
    elem1d->grad_shp[0] = -1.0 / elem1d->djac;
    elem1d->grad_shp[1] =  1.0 / elem1d->djac;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* 2D ELEMENTS */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// local basis function values
//  2 o
//    | \
//    x   x
//    |    \
//  0 o--x--o 1
//
inline void get_triangle_local_shape(double xhat, double yhat, double zhat, double *lshape) {
    lshape[0] = 1-yhat-xhat;
    lshape[1] = xhat;
    lshape[2] = yhat;
}

// local basis function values
//  3 o----o 2
//    |    |
//    |    |
//  0 o----o 1
//
inline void get_quadrilateral_local_shape(double xhat, double yhat, double zhat, double *lshape) {
    lshape[0] = one_4*(1-xhat)*(1-yhat); // [-1,-1]
    lshape[1] = one_4*(1+xhat)*(1-yhat); // [ 1,-1]
    lshape[2] = one_4*(1+xhat)*(1+yhat); // [ 1, 1]
    lshape[3] = one_4*(1-xhat)*(1+yhat); // [-1, 1]
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// local quadratic basis function values
//  2 o
//    |\
//  5 x  x 4
//    |   \
//  0 o--x--o 1
//       3
inline void get_triangle_local_shape_quad(double xhat, double yhat, double zhat, double *lshape) {
    lshape[0] = (1-yhat-xhat)*(1-2*xhat-2*yhat);
    lshape[1] = xhat*(2*xhat-1);
    lshape[2] = yhat*(2*yhat-1);
    lshape[3] = 4*xhat*(1-xhat-yhat);
    lshape[4] = 4*xhat*yhat;
    lshape[5] = 4*yhat*(1-xhat-yhat);
}

// local quadratic basis function values
//       6
//  3 o--x--o 2
//    |     |
//  7 x     x 5
//    |     |
//  0 o--x--o 1
//       4
inline void get_quadrilateral_local_shape_quad(double xhat, double yhat, double zhat, double *lshape) {
    
    lshape[0] = one_4*(1-xhat)*(1-yhat)*(-xhat-yhat-1);    // [-1,-1]
    lshape[1] = one_4*(1+xhat)*(1-yhat)*( xhat-yhat-1);    // [ 1,-1]
    lshape[2] = one_4*(1+xhat)*(1+yhat)*( xhat+yhat-1);    // [ 1, 1]
    lshape[3] = one_4*(1-xhat)*(1+yhat)*(-xhat+yhat-1);    // [-1, 1]
    lshape[4] = one_2*(1-xhat*xhat)*(1-yhat);              // [ 0,-1]
    lshape[5] = one_2*(1-yhat*yhat)*(1+xhat);              // [ 1, 0]
    lshape[6] = one_2*(1-xhat*xhat)*(1+yhat);              // [ 0, 1]
    lshape[7] = one_2*(1-yhat*yhat)*(1-xhat);              // [-1, 0]
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// local basis function parent space gradients
inline void get_triangle_local_shape_gradients(SVECT *lgrad_shp) {
    lgrad_shp[0].x =  -1.;
    lgrad_shp[0].y =  -1.;
    lgrad_shp[1].x =   1.;
    lgrad_shp[1].y =   0.;
    lgrad_shp[2].x =   0.;
    lgrad_shp[2].y =   1.;
}

inline void get_quadrilateral_local_shape_gradients(double xhat, double yhat, double zhat, SVECT *lgrad_shp) {
    lgrad_shp[0].x =  one_4*(yhat - 1);
    lgrad_shp[0].y =  one_4*(xhat - 1);
    lgrad_shp[1].x = -one_4*(yhat - 1);
    lgrad_shp[1].y = -one_4*(xhat + 1);
    lgrad_shp[2].x =  one_4*(yhat + 1);
    lgrad_shp[2].y =  one_4*(xhat + 1);
    lgrad_shp[3].x = -one_4*(yhat + 1);
    lgrad_shp[3].y = -one_4*(xhat - 1);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// get the 2D elements horizontally projected area :: if element is a triangle, assume djac2d has been calculated
inline double get_elem2d_area2d(SELEM_2D *elem2d, SNODE *node) {
    if (elem2d->nnodes == NDONTRI) {
        return elem2d->djac;  // djac2d is actually the area in AdH
    } else {
        double d1 = node[elem2d->nodes[1]].x - node[elem2d->nodes[3]].x;
        double d2 = node[elem2d->nodes[0]].x - node[elem2d->nodes[2]].x;
        double t1 = d1 * node[elem2d->nodes[0]].y;
        double t2 = d2 * node[elem2d->nodes[1]].y;
        double t3 = d1 * node[elem2d->nodes[2]].y;
        double t4 = d2 * node[elem2d->nodes[3]].y;
        return (0.5*(-t1 + t2 + t3 - t4));
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// get centroid of triangle
inline SVECT get_triangle_centroid(SNODE nd1, SNODE nd2, SNODE nd3) {
    SVECT center;
    center.x = one_3 * (nd1.x + nd2.x + nd3.x);
    center.y = one_3 * (nd1.y + nd2.y + nd3.y);
    center.z = one_3 * (nd1.z + nd2.z + nd3.z);
    return center;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Returns a double which established which way the nodes on the triangule are numbers according to a
// given reference vector
// Return variable < 0 :: counter clock-wise node numbering w.r.t. the reference vector
inline double get_triangle_orientation(SNODE nd1, SNODE nd2, SNODE nd3, SVECT ref_vec) {
    
    SVECT side1, side2;		/* two sides - the cross product is the normal */
    double magnitude;		/* the magnitude of the normal */
    
    /* initialize the variables */
    svect_init(&side1);
    svect_init(&side2);
    
    /* computes two side vectors */
    side1.x = nd2.x - nd1.x;
    side1.y = nd2.y - nd1.y;
    side1.z = nd2.z - nd1.z;
    side2.x = nd3.x - nd1.x;
    side2.y = nd3.y - nd1.y;
    side2.z = nd3.z - nd1.z;
    
    SVECT cross = svect_cross(side1,side2);
    return svect_dotp(cross,ref_vec);
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// these are all constants on a triangle and can be calculated once and stored if grid points do not move
// note :: nd_SNODE is a global node vector
inline void get_triangle_linear_djac_nrml_gradPhi(SELEM_2D *elem2d, SNODE *nd_SNODE, SVECT *nd_SVECT) {
    
    SVECT side1, side2;		/* two sides - the cross product is the normal */
    double magnitude;		/* the magnitude of the normal */
    
#ifdef _DEBUG
    if (nd_SVECT == NULL) assert(nd_SNODE != NULL);
    if (nd_SNODE == NULL) assert(nd_SVECT != NULL);
#endif
    
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    if (nd_SNODE != NULL) { /* Gajanan gkc : Changing the to allow passing mod->grid->node as nd_SNODE in the function. */
        x1 = nd_SNODE[elem2d->nodes[0]].x; x2 = nd_SNODE[elem2d->nodes[1]].x; x3 = nd_SNODE[elem2d->nodes[2]].x;
        y1 = nd_SNODE[elem2d->nodes[0]].y; y2 = nd_SNODE[elem2d->nodes[1]].y; y3 = nd_SNODE[elem2d->nodes[2]].y;
        z1 = nd_SNODE[elem2d->nodes[0]].z; z2 = nd_SNODE[elem2d->nodes[1]].z; z3 = nd_SNODE[elem2d->nodes[2]].z;
    } else {
        x1 = nd_SVECT[0].x; x2 = nd_SVECT[1].x; x3 = nd_SVECT[2].x;
        y1 = nd_SVECT[0].y; y2 = nd_SVECT[1].y; y3 = nd_SVECT[2].y;
        z1 = nd_SVECT[0].z; z2 = nd_SVECT[1].z; z3 = nd_SVECT[2].z;
    }
    
    //printf("x1: %20.10f y1: %20.10f z1: %20.10f\n",x1,y1,z1);
    //printf("x2: %20.10f y2: %20.10f z2: %20.10f\n",x2,y2,z2);
    //printf("x3: %20.10f y3: %20.10f z3: %20.10f\n",x3,y3,z3);
    
    /* initialize the variables */
    svect_init(&elem2d->nrml);
    svect_init(&side1);
    svect_init(&side2);
    
    /* computes two side vectors */
    side1.x = x2 - x1;
    side1.y = y2 - y1;
    side1.z = z2 - z1;
    side2.x = x3 - x1;
    side2.y = y3 - y1;
    side2.z = z3 - z1;
    
    /* computes 2d element un-normalized normal and normal mag */
    elem2d->nrml.x = side1.y * side2.z - side1.z * side2.y;
    elem2d->nrml.y = side1.z * side2.x - side1.x * side2.z;
    elem2d->nrml.z = side1.x * side2.y - side1.y * side2.x;
    magnitude = svect_mag(elem2d->nrml);
    
    /* normalizes the normal */
    svect_scale_replace(&(elem2d->nrml), 1./magnitude);
    
    /* sets the jacobian */
    elem2d->djac3d = SCALE_FOR_MASTER_TRI_VOLUME * magnitude;
    assert(fabs(elem2d->djac3d) > SMALL);
    
    /* if 3d jacobian for original element has not been set yet, do it */
    if (fabs(elem2d->djac3d_fixed) < 1e-7) {
        elem2d->djac3d_fixed = elem2d->djac3d;
    }
    
    /* computes the 2d jacobian */
    elem2d->djac = (side1.x * side2.y) - (side2.x * side1.y);
    
    /* sets the 2d jacobian and the 2d shape functions */
    // in 3d, this jacobian is (-)!!!!! causes problem in bedload 

    if(fabs(elem2d->djac) > SMALL) {
        elem2d->grad_shp[0].x = (side1.y - side2.y);
        elem2d->grad_shp[0].y = (side2.x - side1.x);
        elem2d->grad_shp[1].x =  side2.y;
        elem2d->grad_shp[1].y = -side2.x;
        elem2d->grad_shp[2].x = -side1.y;
        elem2d->grad_shp[2].y =  side1.x;
        svect2d_scale_replace_array(elem2d->grad_shp, (1./elem2d->djac), elem2d->nnodes);
    } else {
        elem2d->grad_shp[0].x = UNSET_FLT;
        elem2d->grad_shp[0].y = UNSET_FLT;
        elem2d->grad_shp[1].x = UNSET_FLT;
        elem2d->grad_shp[1].y = UNSET_FLT;
        elem2d->grad_shp[2].x = UNSET_FLT;
        elem2d->grad_shp[2].y = UNSET_FLT;
    }
    
    /* scales the 2d jacobian to make it the area */
    elem2d->djac *= SCALE_FOR_MASTER_TRI_VOLUME;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// returns a normal vector to the 2D element plane
inline SVECT get_elem2d_normals(SVECT *nd) {
    SVECT side1, side2;		    /* two sides - the cross product is the normal */
    double magnitude = 0.;		/* the magnitude of the normal */
    SVECT nrml;
    svect_init(&nrml);
    svect_init(&side1);
    svect_init(&side2);
    
    /* computes two side vectors */
    side1.x = nd[1].x - nd[0].x;
    side1.y = nd[1].y - nd[0].y;
    side1.z = nd[1].z - nd[0].z;
    side2.x = nd[2].x - nd[0].x;
    side2.y = nd[2].y - nd[0].y;
    side2.z = nd[2].z - nd[0].z;
    
    /* computes their cross product */
    nrml.x = side1.y * side2.z - side1.z * side2.y;
    nrml.y = side1.z * side2.x - side1.x * side2.z;
    nrml.z = side1.x * side2.y - side1.y * side2.x;
    magnitude = svect_mag(nrml);
    
    /* normalizes the normal */
    svect_scale_replace(&nrml, 1./magnitude);
    
    return nrml;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// returns the 2D projected quadrilateral Jacobian
inline double get_quadrilateral_linear_djac2d(double xhat, double yhat, SVECT *nd) {
    double x12, x23, x34, x14, x24, x13, tx, ty, tc, term;
    x12 = (nd[0].x - nd[1].x);
    x23 = (nd[1].x - nd[2].x);
    x34 = (nd[2].x - nd[3].x);
    x14 = (nd[0].x - nd[3].x);
    x24 = (nd[1].x - nd[3].x);
    x13 = (nd[0].x - nd[2].x);
    tx = ( x34*nd[0].y - x34*nd[1].y - x12*nd[2].y + x12*nd[3].y);
    ty = ( x23*nd[0].y - x14*nd[1].y + x14*nd[2].y - x23*nd[3].y);
    tc = (-x24*nd[0].y + x13*nd[1].y + x24*nd[2].y - x13*nd[3].y);
    term = tx*xhat + ty*yhat + tc;
    return (one_8 * term);
}

// cjt :: when called on sidewall, we just use the regular djac2d with coordinate exchange, since our grid is constrained to columns
//inline double elem2d_get_quad_djac2d3d(double xhat, double yhat, double zhat, double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double z1, double z2, double z3) {
//    double djac2d = elem2d_get_quad_djac2d_2(xhat, yhat, x1, x2, x3, x4, y1, y2, y3, y4);
//    //printf("djac2d: %20.10f \t ((x2 - x3)*y1 - (x1 - x3)*y2 + (x1 - x2)*y3): %20.10f \n",djac2d,((x2 - x3)*y1 - (x1 - x3)*y2 + (x1 - x2)*y3));
//    double xx1 = x1*x1, xx2 = x2*x2, xx3 = x3*x3, yy1 = y1*y1, yy2 = y2*y2, yy3 = y3*y3, zz1 = z1*z1, zz2 = z2*z2, zz3 = z3*z3;
//    double c1 = one_8*sqrt((xx2 - 2*x2*x3 + xx3)*yy1 - 2*(x1*x2 - (x1 + x2)*x3 + xx3)*y1*y2 + (xx1 - 2*x1*x3 + xx3)*yy2 + (xx1 - 2*x1*x2 + xx2)*yy3 + (xx2 - 2*x2*x3 + xx3 + yy2 - 2*y2*y3 + yy3)*zz1 - 2*(x1*x2 - (x1 + x2)*x3 + xx3 + y1*y2 - (y1 + y2)*y3 + yy3)*z1*z2 + (xx1 - 2*x1*x3 + xx3 + yy1 - 2*y1*y3 + yy3)*zz2 + (xx1 - 2*x1*x2 + xx2 + yy1 - 2*y1*y2 + yy2)*zz3 + 2*((x1*x2 - xx2 - (x1 - x2)*x3)*y1 - (xx1 - x1*x2 - (x1 - x2)*x3)*y2)*y3 + 2*((x1*x2 - xx2 - (x1 - x2)*x3 + y1*y2 - yy2 - (y1 - y2)*y3)*z1 - (xx1 - x1*x2 - (x1 - x2)*x3 + yy1 - y1*y2 - (y1 - y2)*y3)*z2)*z3)/((x2 - x3)*y1 - (x1 - x3)*y2 + (x1 - x2)*y3);
//    return (-c1 * djac2d * 8);
//}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// returns the 2D projected quadrilateral Jacobian and cartesian space shape function gradients
inline double get_quadrilateral_linear_djac_gradPhi(double xhat, double yhat, SVECT *nd, SVECT *grad_shp) {
    
    double djac2d = get_quadrilateral_linear_djac2d(xhat, yhat, nd);
    double con = (1./(djac2d * 8.));
    grad_shp[ 0 ].x = con * (       -xhat*nd[2].y + (xhat - 1)*nd[3].y - (nd[1].y - nd[2].y)*yhat + nd[1].y );
    grad_shp[ 1 ].x = con * (  (xhat + 1)*nd[2].y -       xhat*nd[3].y + (nd[0].y - nd[3].y)*yhat - nd[0].y );
    grad_shp[ 2 ].x = con * (        xhat*nd[0].y - (xhat + 1)*nd[1].y - (nd[0].y - nd[3].y)*yhat + nd[3].y );
    grad_shp[ 3 ].x = con * ( -(xhat - 1)*nd[0].y +       xhat*nd[1].y + (nd[1].y - nd[2].y)*yhat - nd[2].y );
    
    grad_shp[ 0 ].y = con * (  (nd[2].x - nd[3].x)*xhat + (nd[1].x - nd[2].x)*yhat - nd[1].x + nd[3].x );
    grad_shp[ 1 ].y = con * ( -(nd[2].x - nd[3].x)*xhat - (nd[0].x - nd[3].x)*yhat + nd[0].x - nd[2].x );
    grad_shp[ 2 ].y = con * ( -(nd[0].x - nd[1].x)*xhat + (nd[0].x - nd[3].x)*yhat + nd[1].x - nd[3].x );
    grad_shp[ 3 ].y = con * (  (nd[0].x - nd[1].x)*xhat - (nd[1].x - nd[2].x)*yhat - nd[0].x + nd[2].x );
    
    return djac2d;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// 3D ELEMENTS
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// local basis function values
inline void get_tet_local_shape(double xhat, double yhat, double zhat, double *lshape) {
    lshape[0] = 1-xhat-yhat-zhat;
    lshape[1] = xhat;
    lshape[2] = yhat;
    lshape[3] = zhat;
}

inline void get_triprism_local_shape(double xhat, double yhat, double zhat, double *lshape) {
    double t1 = 1-zhat, t2 = 1+zhat, t3 = 1-xhat-yhat;
    lshape[0] = one_2*t1*t3;   // {0,0,-1}
    lshape[1] = one_2*t1*xhat; // {1,0,-1}
    lshape[2] = one_2*t1*yhat; // {0,1,-1}
    lshape[3] = one_2*t2*t3;   // {0,0, 1}
    lshape[4] = one_2*t2*xhat; // {1,0, 1}
    lshape[5] = one_2*t2*yhat; // {0,1, 1}
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// local quadratic basis function values
inline void get_tet_local_shape_quad(double xhat, double yhat, double zhat, double *lshape_quad) {
    double lshape[4]; sarray_init_dbl(lshape, 4);
    get_tet_local_shape(xhat, yhat, zhat, lshape);
    lshape_quad[0] = lshape[0] * (2*lshape[0] - 1);
    lshape_quad[1] = lshape[1] * (2*lshape[1] - 1);
    lshape_quad[2] = lshape[2] * (2*lshape[2] - 1);
    lshape_quad[3] = lshape[3] * (2*lshape[3] - 1);
    lshape_quad[4] = 4 * lshape[0] * lshape[1];
    lshape_quad[5] = 4 * lshape[0] * lshape[2];
    lshape_quad[6] = 4 * lshape[0] * lshape[3];
    lshape_quad[7] = 4 * lshape[1] * lshape[2];
    lshape_quad[8] = 4 * lshape[1] * lshape[3];
    lshape_quad[9] = 4 * lshape[2] * lshape[3];
}

inline void get_triprism_local_shape_quad(double xhat, double yhat, double zhat, double *lshape_quad) {
    double t1 = 1-zhat, t2 = 1+zhat, t3 = 1-xhat-yhat;
    lshape_quad[0]  =  one_2*t3*t1*(-2*xhat-2*yhat-zhat);  // {0,0,-1}
    lshape_quad[3]  =  one_2*t3*t2*(-2*xhat-2*yhat+zhat);  // {0,0, 1}
    lshape_quad[1]  =  one_2*xhat*t1*(2*xhat - zhat - 2);  // {1,0,-1}
    lshape_quad[4]  =  one_2*xhat*t2*(2*xhat + zhat - 2);  // {1,0, 1}
    lshape_quad[2]  =  one_2*yhat*t1*(2*yhat - zhat - 2);  // {0,1,-1}
    lshape_quad[5]  =  one_2*yhat*t2*(2*yhat + zhat - 2);  // {0,1, 1}
    lshape_quad[6]  =  2*xhat*t3*t1;                       // {(1/2),0,-1}
    lshape_quad[12]  =  2*xhat*t3*t2;                      // {(1/2),0, 1}
    lshape_quad[7]  =  2*xhat*yhat*t1;                     // {(1/2),(1/2),-1}
    lshape_quad[13]  =  2*xhat*yhat*t2;                    // {(1/2),(1/2), 1}
    lshape_quad[8] =  2*yhat*t3*t1;                        // {0,(1/2),-1}
    lshape_quad[14] =  2*yhat*t3*t2;                       // {0,(1/2), 1}
    lshape_quad[9] =  t3*t2*t1;                            // {0,0,0}
    lshape_quad[10] =  xhat*t2*t1;                         // {1,0,0}
    lshape_quad[11] =  yhat*t2*t1;                         // {0,1,0}
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// local basis function gradients in parent space
inline void get_tet_local_shape_gradients(SVECT *lgrad_shp) {
    lgrad_shp[0].x =  -1.;
    lgrad_shp[0].y =  -1.;
    lgrad_shp[0].z =  -1.;
    
    lgrad_shp[1].x =  1.;
    lgrad_shp[1].y =  0.;
    lgrad_shp[1].z =  0.;
    
    lgrad_shp[2].x =  0.;
    lgrad_shp[2].y =  1.;
    lgrad_shp[2].z =  0.;
    
    lgrad_shp[3].x =  0.;
    lgrad_shp[3].y =  0.;
    lgrad_shp[3].z =  1.;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void get_triprism_local_shape_gradients(double xhat, double yhat, double zhat, SVECT *lgrad_shp) {
    double t1 = zhat-1, t3 = xhat + yhat - 1;
    lgrad_shp[0].x =  one_2*t1;
    lgrad_shp[0].y =  one_2*t1;
    lgrad_shp[0].z =  one_2*t3;
    lgrad_shp[1].x =  one_2*t1;
    lgrad_shp[1].y =  one_2*t1;
    lgrad_shp[1].z =  one_2*t3;
    lgrad_shp[2].x =  one_2*t1;
    lgrad_shp[2].y =  one_2*t1;
    lgrad_shp[2].z =  one_2*t3;
    lgrad_shp[3].x =  one_2*t1;
    lgrad_shp[3].y =  one_2*t1;
    lgrad_shp[3].z =  one_2*t3;
    lgrad_shp[4].x =  one_2*t1;
    lgrad_shp[4].y =  one_2*t1;
    lgrad_shp[4].z =  one_2*t3;
    lgrad_shp[5].x =  one_2*t1;
    lgrad_shp[5].y =  one_2*t1;
    lgrad_shp[5].z =  one_2*t3;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// calculates determinant given {xhat, yhat, zhat} for a constrained prism
inline double get_triprism_djac(double xhat, double yhat, double zhat, SVECT *nd)  {
#if _DEBUG
    assert(nd[0].x == nd[3].x); assert(nd[1].x == nd[4].x); assert(nd[2].x == nd[5].x);
    assert(nd[0].y == nd[3].y); assert(nd[1].y == nd[4].y); assert(nd[2].y == nd[5].y);
#endif
    
    return (-one_2*(nd[1].x*nd[0].y - nd[2].x*nd[0].y - nd[0].x*nd[1].y + nd[2].x*nd[1].y + nd[0].x*nd[2].y - nd[1].x*nd[2].y)*
    (xhat*nd[0].z + yhat*nd[0].z - xhat*nd[1].z - yhat*nd[2].z - xhat*nd[3].z - yhat*nd[3].z + xhat*nd[4].z + yhat*nd[5].z - nd[0].z + nd[3].z));
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// calculates determinant and cartesian shape function gradients for a constrained prism (x4=x1, x5=x2, ....)
// 415 operations
inline double get_triprism_linear_djac_gradPhi(double xhat, double yhat, double zhat, SVECT *nd, SVECT *grad_shp) {
    
#if _DEBUG
    assert(nd[0].x == nd[3].x); assert(nd[1].x == nd[4].x); assert(nd[2].x == nd[5].x);
    assert(nd[0].y == nd[3].y); assert(nd[1].y == nd[4].y); assert(nd[2].y == nd[5].y);
#endif
    
    double x1 = nd[0].x, x2 = nd[1].x, x3 = nd[2].x, x4 = nd[3].x, x5 = nd[4].x, x6 = nd[5].x;
    double y1 = nd[0].y, y2 = nd[1].y, y3 = nd[2].y, y4 = nd[3].y, y5 = nd[4].y, y6 = nd[5].y;
    double z1 = nd[0].z, z2 = nd[1].z, z3 = nd[2].z, z4 = nd[3].z, z5 = nd[4].z, z6 = nd[5].z;
    
    double x12 = x1 - x2, x13 = x1 - x3, x23 = x2 - x3;
    double y12 = y1 - y2, y13 = y1 - y3, y23 = y2 - y3;
    
    double phi1 = xhat + yhat - 1;
    double phi2 = xhat - 1;
    double phi3 = yhat - 1;
    double d1 = (x23*y1 - x13*y2 + x12*y3), d2 = (x23*xhat*y1 - x13*xhat*y2 + x12*xhat*y3), d3 = (x23*xhat - x23);
    double d4 = (x12*xhat - x12), d5 = (x13*xhat - x13), d6 = (d3*y1 - d5*y2 + d4*y3 + d1*yhat);
    double denom = (d1*yhat*z3 - d1*yhat*z6 - d6*z1 + d2*z2 + d6*z4 - d2*z5);
    double denom2 = (phi1*z1 - xhat*z2 - yhat*z3 - phi1*z4 + xhat*z5 + yhat*z6);
    
    double djac3d = one_2 * denom;
    
    // adh prism code
    djac3d = (-one_2*(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3)*(xhat*z1 + yhat*z1 - xhat*z2 - yhat*z3 - xhat*z4 - yhat*z4 + xhat*z5 + yhat*z6 - z1 + z4));
    
    double t1 = (phi2*y1 - xhat*y2 + y13*yhat + y3);
    double t2 = (phi2*y1 - phi2*y2 + y13*yhat);
    double t3 = (phi2*y2 - phi2*y3 + y23*yhat);
    double t4 = (phi2*y1 + xhat*y2 - (2*xhat - 1)*y3 + y13*yhat);
    double t5 = (phi2*y1 - phi2*y2 + (y1 - 2*y2 + y3)*yhat);
    double t6 = (x12*xhat + x13*yhat - x13);
    double t7 = (x12*xhat + x13*yhat - x12);
    double t8 = (x23*xhat + x23*yhat - x23);
    double t9 = ((x1 + x2 - 2*x3)*xhat + x13*yhat - x13);
    double t10 = (x12*xhat + (x1 - 2*x2 + x3)*yhat - x12);
    double t11 = (xhat*y1 - xhat*y2 - y13*yhat);
    double t12 = (x12*xhat - x13*yhat);
    double t13 = (x12*xhat + x13*yhat);
    double t14 = (xhat*y1 - xhat*y2 + y13*yhat);
    
    double con = one_4 / djac3d;
    
    grad_shp[0].x =   con*(t1*z2 - t2*z3 - 2*t3*z4 + t4*z5 -  t5*z6 - (t1*z2 - t2*z3 - t1*z5 + t2*z6)*zhat);
    grad_shp[0].y =  -con*(t6*z2 - t7*z3 - 2*t8*z4 + t9*z5 - t10*z6 - (t6*z2 - t7*z3 - t6*z5 + t7*z6)*zhat);
    grad_shp[0].z =  phi1/denom2;
    
    grad_shp[1].x =  -con*(t1*z1 - t14*z3 - t4*z4 + 2*(xhat*y1 - xhat*y3)*z5 - t11*z6 - (t1*z1 - t14*z3 - t1*z4 + t14*z6)*zhat);
    grad_shp[1].y =   con*(2*x13*xhat*z5 + t6*z1 - t13*z3 - t9*z4 - t12*z6 - (t6*z1 - t13*z3 - t6*z4 + t13*z6)*zhat);
    grad_shp[1].z =  -xhat/denom2;
    
    grad_shp[2].x =   con*(2*y12*yhat*z6 + t2*z1 - t14*z2 -  t5*z4 + t11*z5 - (t2*z1 - t14*z2 - t2*z4 + t14*z5)*zhat);
    grad_shp[2].y =  -con*(2*x12*yhat*z6 + t7*z1 - t13*z2 - t10*z4 + t12*z5 - (t7*z1 - t13*z2 - t7*z4 + t13*z5)*zhat);
    grad_shp[2].z =  -yhat/denom2;
    
    grad_shp[3].x =   con*(2*t3*z1 - t4*z2 +  t5*z3 - t1*z5 + t2*z6 + (t1*z2 - t2*z3 - t1*z5 + t2*z6)*zhat);
    grad_shp[3].y =  -con*(2*t8*z1 - t9*z2 + t10*z3 - t6*z5 + t7*z6 + (t6*z2 - t7*z3 - t6*z5 + t7*z6)*zhat);
    grad_shp[3].z =  -phi1/denom2;
    
    grad_shp[4].x =  -con*(t4*z1 - 2*(xhat*y1 - xhat*y3)*z2 + t11*z3 - t1*z4 + t14*z6 + (t1*z1 - t14*z3 - t1*z4 + t14*z6)*zhat);
    grad_shp[4].y =  -con*(2*x13*xhat*z2 - t9*z1 - t12*z3 + t6*z4 - t13*z6 - (t6*z1 - t13*z3 - t6*z4 + t13*z6)*zhat);
    grad_shp[4].z =  xhat/denom2;

    grad_shp[5].x =  -con*(2*y12*yhat*z3 -  t5*z1 + t11*z2 + t2*z4 - t14*z5 - (t2*z1 - t14*z2 - t2*z4 + t14*z5)*zhat);
    grad_shp[5].y =   con*(2*x12*yhat*z3 - t10*z1 + t12*z2 + t7*z4 - t13*z5 - (t7*z1 - t13*z2 - t7*z4 + t13*z5)*zhat);
    grad_shp[5].z =  yhat/denom2;
    
    // old prism code way
    grad_shp[0].x =  one_2*(t1*z2 - t2*z3 - 2*t3*z4 + t4*z5 - t5*z6 - (t1*z2 - t2*z3 - t1*z5 + t2*z6)*zhat)/denom;
    grad_shp[0].y =  -one_2*(t6*z2 - t7*z3 - 2*t8*z4 + t9*z5 - t10*z6 - (t6*z2 - t7*z3 - t6*z5 + t7*z6)*zhat)/denom;
    grad_shp[0].z =  phi1/denom2;
    
    grad_shp[1].x =  -one_2*(t1*z1 - t14*z3 - t4*z4 + 2*(xhat*y1 - xhat*y3)*z5 - t11*z6 - (t1*z1 - t14*z3 - t1*z4 + t14*z6)*zhat)/denom;
    grad_shp[1].y =  one_2*(2*x13*xhat*z5 + t6*z1 - t13*z3 - t9*z4 - t12*z6 - (t6*z1 - t13*z3 - t6*z4 + t13*z6)*zhat)/denom;
    grad_shp[1].z =  -xhat/denom2;
    
    grad_shp[2].x =  one_2*(2*y12*yhat*z6 + t2*z1 - t14*z2 - t5*z4 + t11*z5 - (t2*z1 - t14*z2 - t2*z4 + t14*z5)*zhat)/denom;
    grad_shp[2].y =  -one_2*(2*x12*yhat*z6 + t7*z1 - t13*z2 - t10*z4 + t12*z5 - (t7*z1 - t13*z2 - t7*z4 + t13*z5)*zhat)/denom;
    grad_shp[2].z =  -yhat/denom2;
    
    grad_shp[3].x =  one_2*(2*t3*z1 - t4*z2 + t5*z3 - t1*z5 + t2*z6 + (t1*z2 - t2*z3 - t1*z5 + t2*z6)*zhat)/denom;
    grad_shp[3].y =  -one_2*(2*t8*z1 - t9*z2 + t10*z3 - t6*z5 + t7*z6 + (t6*z2 - t7*z3 - t6*z5 + t7*z6)*zhat)/denom;
    grad_shp[3].z =  -phi1/denom2;
    
    grad_shp[4].x =  -one_2*(t4*z1 - 2*(xhat*y1 - xhat*y3)*z2 + t11*z3 - t1*z4 + t14*z6 + (t1*z1 - t14*z3 - t1*z4 + t14*z6)*zhat)/denom;
    grad_shp[4].y =  -one_2*(2*x13*xhat*z2 - t9*z1 - t12*z3 + t6*z4 - t13*z6 - (t6*z1 - t13*z3 - t6*z4 + t13*z6)*zhat)/denom;
    grad_shp[4].z =  xhat/denom2;
    
    grad_shp[5].x =  -one_2*(2*y12*yhat*z3 - t5*z1 + t11*z2 + t2*z4 - t14*z5 - (t2*z1 - t14*z2 - t2*z4 + t14*z5)*zhat)/denom;
    grad_shp[5].y =  one_2*(2*x12*yhat*z3 - t10*z1 + t12*z2 + t7*z4 - t13*z5 - (t7*z1 - t13*z2 - t7*z4 + t13*z5)*zhat)/denom;
    grad_shp[5].z =  yhat/denom2;
    
    return djac3d;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// calculates 3D tetrahedral djac (area) and cartesian space shape function gradients, all constant on element
inline void get_tet_linear_djac_gradPhi(SELEM_3D *elem3d, SNODE *nd_SNODE, SVECT *nd_SVECT) {

#ifdef _DEBUG
    if (nd_SVECT == NULL) assert(nd_SNODE != NULL);
    if (nd_SNODE == NULL) assert(nd_SVECT != NULL);
#endif
    
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    if (nd_SNODE != NULL) {
        x1 = nd_SNODE[0].x; x2 = nd_SNODE[1].x; x3 = nd_SNODE[2].x; x4 = nd_SNODE[3].x;
        y1 = nd_SNODE[0].y; y2 = nd_SNODE[1].y; y3 = nd_SNODE[2].y; y4 = nd_SNODE[3].y;
        z1 = nd_SNODE[0].z; z2 = nd_SNODE[1].z; z3 = nd_SNODE[2].z; z4 = nd_SNODE[3].z;
    } else {
        x1 = nd_SVECT[0].x; x2 = nd_SVECT[1].x; x3 = nd_SVECT[2].x; x4 = nd_SVECT[3].x;
        y1 = nd_SVECT[0].y; y2 = nd_SVECT[1].y; y3 = nd_SVECT[2].y; y4 = nd_SVECT[3].y;
        z1 = nd_SVECT[0].z; z2 = nd_SVECT[1].z; z3 = nd_SVECT[2].z; z4 = nd_SVECT[3].z;
    }
    
    /* computes three side vectors */
    SVECT side1, side2, side3;
    side1.x = x2 - x1;
    side1.y = y2 - y1;
    side1.z = z2 - z1;
    side2.x = x3 - x1;
    side2.y = y3 - y1;
    side2.z = z3 - z1;
    side3.x = x4 - x1;
    side3.y = y4 - y1;
    side3.z = z4 - z1;
    
    /* determinant of Jacobian map */
    elem3d->djac  = side1.x * (side2.y * side3.z - side3.y * side2.z);
    elem3d->djac -= side2.x * (side1.y * side3.z - side3.y * side1.z);
    elem3d->djac += side3.x * (side1.y * side2.z - side2.y * side1.z);
    
    /* elements of the Jacobian inverse map matrix */
    double dksidx_0, dksidx_1, dksidx_2, dksidx_3, dksidx_4, dksidx_5, dksidx_6, dksidx_7, dksidx_8;
    dksidx_0 = (side2.y * side3.z - side3.y * side2.z);
    dksidx_1 = (side3.x * side2.z - side2.x * side3.z);
    dksidx_2 = (side2.x * side3.y - side3.x * side2.y);
    dksidx_3 = (side3.y * side1.z - side1.y * side3.z);
    dksidx_4 = (side1.x * side3.z - side3.x * side1.z);
    dksidx_5 = (side3.x * side1.y - side1.x * side3.y);
    dksidx_6 = (side1.y * side2.z - side2.y * side1.z);
    dksidx_7 = (side2.x * side1.z - side1.x * side2.z);
    dksidx_8 = (side1.x * side2.y - side2.x * side1.y);
    
    /* gradients */
    elem3d->grad_shp[0].x = -dksidx_0 - dksidx_3 - dksidx_6;
    elem3d->grad_shp[0].y = -dksidx_1 - dksidx_4 - dksidx_7;
    elem3d->grad_shp[0].z = -dksidx_2 - dksidx_5 - dksidx_8;
    elem3d->grad_shp[1].x = dksidx_0;
    elem3d->grad_shp[1].y = dksidx_1;
    elem3d->grad_shp[1].z = dksidx_2;
    elem3d->grad_shp[2].x = dksidx_3;
    elem3d->grad_shp[2].y = dksidx_4;
    elem3d->grad_shp[2].z = dksidx_5;
    elem3d->grad_shp[3].x = dksidx_6;
    elem3d->grad_shp[3].y = dksidx_7;
    elem3d->grad_shp[3].z = dksidx_8;
    
    /* one list jacobian determinate divide */
    svect_scale_replace_array(elem3d->grad_shp, (1./elem3d->djac), elem3d->nnodes);
    
    // convert djac to area
    elem3d->djac *= SCALE_FOR_MASTER_TET_VOLUME;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// calculates 3D tetrahedral djac (area) and cartesian space shape function gradients, all constant on element
inline double get_tet_linear_djac(SNODE *nd_SNODE, SVECT *nd_SVECT) {
    
#ifdef _DEBUG
    if (nd_SVECT == NULL) assert(nd_SNODE != NULL);
    if (nd_SNODE == NULL) assert(nd_SVECT != NULL);
#endif
    
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    if (nd_SNODE != NULL) {
        x1 = nd_SNODE[0].x; x2 = nd_SNODE[1].x; x3 = nd_SNODE[2].x; x4 = nd_SNODE[3].x;
        y1 = nd_SNODE[0].y; y2 = nd_SNODE[1].y; y3 = nd_SNODE[2].y; y4 = nd_SNODE[3].y;
        z1 = nd_SNODE[0].z; z2 = nd_SNODE[1].z; z3 = nd_SNODE[2].z; z4 = nd_SNODE[3].z;
    } else {
        x1 = nd_SVECT[0].x; x2 = nd_SVECT[1].x; x3 = nd_SVECT[2].x; x4 = nd_SVECT[3].x;
        y1 = nd_SVECT[0].y; y2 = nd_SVECT[1].y; y3 = nd_SVECT[2].y; y4 = nd_SVECT[3].y;
        z1 = nd_SVECT[0].z; z2 = nd_SVECT[1].z; z3 = nd_SVECT[2].z; z4 = nd_SVECT[3].z;
    }
    
    /* computes three side vectors */
    SVECT side1, side2, side3;
    side1.x = x2 - x1;
    side1.y = y2 - y1;
    side1.z = z2 - z1;
    side2.x = x3 - x1;
    side2.y = y3 - y1;
    side2.z = z3 - z1;
    side3.x = x4 - x1;
    side3.y = y4 - y1;
    side3.z = z4 - z1;
    
    /* determinant of Jacobian map */
    double djac = 0.;
    djac  = side1.x * (side2.y * side3.z - side3.y * side2.z);
    djac -= side2.x * (side1.y * side3.z - side3.y * side1.z);
    djac += side3.x * (side1.y * side2.z - side2.y * side1.z);
    return djac;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// calculates 3D tetrahedral djac (area) and cartesian space shape function gradients, all constant on element
inline double get_tet_linear_djac_gradPhi2(SNODE *nd_SNODE, SVECT *nd_SVECT, SVECT *grad_shp) {
    
#ifdef _DEBUG
    if (nd_SVECT == NULL) assert(nd_SNODE != NULL);
    if (nd_SNODE == NULL) assert(nd_SVECT != NULL);
#endif
    
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    if (nd_SNODE != NULL) {
        x1 = nd_SNODE[0].x; x2 = nd_SNODE[1].x; x3 = nd_SNODE[2].x; x4 = nd_SNODE[3].x;
        y1 = nd_SNODE[0].y; y2 = nd_SNODE[1].y; y3 = nd_SNODE[2].y; y4 = nd_SNODE[3].y;
        z1 = nd_SNODE[0].z; z2 = nd_SNODE[1].z; z3 = nd_SNODE[2].z; z4 = nd_SNODE[3].z;
    } else {
        x1 = nd_SVECT[0].x; x2 = nd_SVECT[1].x; x3 = nd_SVECT[2].x; x4 = nd_SVECT[3].x;
        y1 = nd_SVECT[0].y; y2 = nd_SVECT[1].y; y3 = nd_SVECT[2].y; y4 = nd_SVECT[3].y;
        z1 = nd_SVECT[0].z; z2 = nd_SVECT[1].z; z3 = nd_SVECT[2].z; z4 = nd_SVECT[3].z;
    }
    
    /* computes three side vectors */
    SVECT side1, side2, side3;
    side1.x = x2 - x1;
    side1.y = y2 - y1;
    side1.z = z2 - z1;
    side2.x = x3 - x1;
    side2.y = y3 - y1;
    side2.z = z3 - z1;
    side3.x = x4 - x1;
    side3.y = y4 - y1;
    side3.z = z4 - z1;
    
    /* determinant of Jacobian map */
    double djac = 0.;
    djac  = side1.x * (side2.y * side3.z - side3.y * side2.z);
    djac -= side2.x * (side1.y * side3.z - side3.y * side1.z);
    djac += side3.x * (side1.y * side2.z - side2.y * side1.z);
    
    /* elements of the Jacobian inverse map matrix */
    double dksidx_0, dksidx_1, dksidx_2, dksidx_3, dksidx_4, dksidx_5, dksidx_6, dksidx_7, dksidx_8;
    dksidx_0 = (side2.y * side3.z - side3.y * side2.z);
    dksidx_1 = (side3.x * side2.z - side2.x * side3.z);
    dksidx_2 = (side2.x * side3.y - side3.x * side2.y);
    dksidx_3 = (side3.y * side1.z - side1.y * side3.z);
    dksidx_4 = (side1.x * side3.z - side3.x * side1.z);
    dksidx_5 = (side3.x * side1.y - side1.x * side3.y);
    dksidx_6 = (side1.y * side2.z - side2.y * side1.z);
    dksidx_7 = (side2.x * side1.z - side1.x * side2.z);
    dksidx_8 = (side1.x * side2.y - side2.x * side1.y);
    
    /* gradients */
    grad_shp[0].x = -dksidx_0 - dksidx_3 - dksidx_6;
    grad_shp[0].y = -dksidx_1 - dksidx_4 - dksidx_7;
    grad_shp[0].z = -dksidx_2 - dksidx_5 - dksidx_8;
    grad_shp[1].x = dksidx_0;
    grad_shp[1].y = dksidx_1;
    grad_shp[1].z = dksidx_2;
    grad_shp[2].x = dksidx_3;
    grad_shp[2].y = dksidx_4;
    grad_shp[2].z = dksidx_5;
    grad_shp[3].x = dksidx_6;
    grad_shp[3].y = dksidx_7;
    grad_shp[3].z = dksidx_8;
    
    /* one list jacobian determinate divide */
    svect_scale_replace_array(grad_shp, (1./djac), NDONTET);
    
    // convert djac to area
    djac *= SCALE_FOR_MASTER_TET_VOLUME;
    
    return djac;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// returns a 3D element volume
inline double get_elem3d_volume(SVECT *node, int nnodes) {
    if (nnodes == NDONTET) {
        double d1 = node[1].x - node[2].x, d2 = node[0].x - node[2].x, d3 = node[0].x - node[1].x;
        double d4 = node[2].x - node[3].x, d5 = node[1].x - node[3].x, d6 = node[0].x - node[3].x;
        return (one_6 * ((d4*node[1].y - d5*node[2].y + d1*node[3].y)*node[0].z -
                         (d4*node[0].y - d6*node[2].y + d2*node[3].y)*node[1].z +
                         (d5*node[0].y - d6*node[1].y + d3*node[3].y)*node[2].z -
                         (d1*node[0].y - d2*node[1].y + d3*node[2].y)*node[3].z));
    } else if (nnodes == NDONPRISM) {
        return get_triprism_volume(node);
    } else {
        tl_error("bad call to elem2d_get_volume");
    }
    return (-1);
}

inline double get_triprism_volume(SVECT *node) {
    double d1 = node[1].x - node[2].x, d2 = node[0].x - node[2].x, d3 = node[0].x - node[1].x;
    double prism_term = d1*node[0].y - d2*node[1].y + d3*node[2].y;
    double t1 = node[0].z + node[1].z + node[2].z;
    t1 -= (node[3].z + node[4].z + node[5].z);
    return (one_6 * prism_term * t1);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
































