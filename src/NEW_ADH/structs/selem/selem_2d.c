#include "adh.h"
#define SCALE_FOR_MASTER_TRI_VOLUME (1.0/2.0) // tetrahedral djac to volume scale

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// local basis function values
//  2 o
//    | \
//    x   x
//    |    \
//  0 o--x--o 1
//
void selem2d_get_triangle_local_shape(double xhat, double yhat, double zhat, double *lshape) {
    lshape[0] = 1-yhat-xhat;
    lshape[1] = xhat;
    lshape[2] = yhat;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// local basis function values
//  3 o----o 2
//    |    |
//    |    |
//  0 o----o 1
//
void selem2d_get_quadrilateral_local_shape(double xhat, double yhat, double zhat, double *lshape) {
    lshape[0] = one_4*(1-xhat)*(1-yhat); // [-1,-1]
    lshape[1] = one_4*(1+xhat)*(1-yhat); // [ 1,-1]
    lshape[2] = one_4*(1+xhat)*(1+yhat); // [ 1, 1]
    lshape[3] = one_4*(1-xhat)*(1+yhat); // [-1, 1]
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void selem2d_alloc(SELEM_2D *elem2d, int nnodes_on_elem) {
    assert(nnodes_on_elem == 3 || nnodes_on_elem == 4);
    elem2d->nnodes = nnodes_on_elem;
    elem2d->grad_shp = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes_on_elem);
    elem2d->nodes = (int *) tl_alloc(sizeof(int), nnodes_on_elem);
    elem2d->levels = (int *) tl_alloc(sizeof(int), nnodes_on_elem);
    int i;
    for (i=0; i<nnodes_on_elem; i++) {
        svect2d_init(&(elem2d->grad_shp[i]));
        elem2d->nodes[i] = UNSET_INT;
        elem2d->levels[i] = 0;
    }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void selem2d_load(SELEM_2D *elem2d, int gid, int lid, int elem_nnodes, int *local_node_ids, int bflag, SVECT *nds, int mat) {
    
    int i;
    
    elem2d->id = lid;
    elem2d->gid = gid;
    elem2d->id_orig = lid; // hmmm, what if we call this later?
    elem2d->bflag = bflag;
    elem2d->mat = mat;
    elem2d->nnodes = elem_nnodes;
    selem2d_alloc(elem2d, elem_nnodes);
    for (i=0;i<elem_nnodes;i++) elem2d->nodes[i] = local_node_ids[i];
   
    if (elem_nnodes == 3) {
        
        elem2d->nedges = 3;
        elem2d->nnodes_quad = 6;
        
        // calculate constant jacobians, normal and basis gradients
        selem2d_get_triangle_linear_djac_nrml_gradPhi(elem2d, NULL, nds);
        
        // node ID direction check
        // cjt -- need to check bed and sidewall ordering too!
        if (bflag == SURFACE || bflag == BODY) {
            if (elem2d->djac < SMALL6) {
                fprintf(stderr, "ERROR :: Improperly numbered triangle gid = %d || Nodes %d %d %d Jacobian is: %20.10f \n",
                        gid+1,elem2d->nodes[0]+1,elem2d->nodes[1]+1,elem2d->nodes[2]+1,elem2d->djac);
                tl_error("ERROR: Improperly numbered triangle");
            }
        } else {
            // 2d element 3D bed jacobian is negative since normal is flipped, so we must flip it back
            elem2d->djac = fabs(elem2d->djac);
            elem2d->djac3d = fabs(elem2d->djac3d);
            elem2d->djac3d_fixed = fabs(elem2d->djac3d_fixed);
        }
        
        // save t=0 grid 2D jacobian in 3D space
        if (elem2d->djac3d_fixed < SMALL6) {
            elem2d->djac3d_fixed = elem2d->djac3d;
        }
        
    } else {
        
        elem2d->nedges = 4;
        elem2d->nnodes_quad = 8;
        elem2d->nrml = selem2d_get_elem2d_normals(nds);
    }
    
    if (bflag != SIDEWALL) {
        elem2d->area = selem2d_get_elem2d_area2d(elem2d,nds);
    } else {
        // CJT -- CURRENTLY ONLY FOR TRIANGLES!
        elem2d->area = elem2d->djac3d;
    }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void selem2d_free(SELEM_2D *elem2d) {
    assert(elem2d->nnodes > 0);
    elem2d->grad_shp = (SVECT2D *) tl_free(sizeof(SVECT2D), elem2d->nnodes, elem2d->grad_shp);
    elem2d->nodes = (int *) tl_free(sizeof(int), elem2d->nnodes, elem2d->nodes);
    elem2d->levels = (int *) tl_free(sizeof(int), elem2d->nnodes, elem2d->levels);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void selem2d_alloc_array(SELEM_2D **elem2d, int nelems2d) {
    assert(nelems2d > 0);
    (*elem2d) = (SELEM_2D *) tl_alloc(sizeof(SELEM_2D), nelems2d);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void selem2d_free_array(SELEM_2D *elem2d, int nelems2d) {
    assert(nelems2d > 0);
    int ie;
    for (ie=0; ie<nelems2d; ie++) {
        selem2d_free(&(elem2d[ie]));
    }
    elem2d = (SELEM_2D *) tl_free(sizeof(SELEM_2D), nelems2d, elem2d);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void selem2d_init(SELEM_2D *elem2d) {
    elem2d->id = UNSET_INT;
    elem2d->gid = UNSET_INT;
    elem2d->id_orig = UNSET_INT;
    elem2d->id_3d  = UNSET_INT;
    elem2d->djac = 0.0;
    elem2d->djac3d = 0.0;
    elem2d->djac3d_fixed = 0.0;
    elem2d->interface=0;
    elem2d->flux_elem_tot = 0;
    svect_init(&(elem2d->nrml));
    elem2d->string = UNSET_INT;
    elem2d->mat = UNSET_INT;
    elem2d->nvars = 0;
    elem2d->vars = NULL;
    elem2d->bflag = UNSET_INT;
    elem2d->nedges = UNSET_INT;
    elem2d->nnodes_quad = UNSET_INT;
    elem2d->flux_elem = NULL;
    elem2d->resident_pe = UNSET_INT;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void selem2d_init_array(SELEM_2D *elem2d, int nelems2d) {
    int ie=0;
    for (ie=0; ie<nelems2d; ie++) {
        selem2d_init(&(elem2d[ie]));
    }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void selem2d_init_alloc_array(SELEM_2D **elem2d, int nelems2d) {
    selem2d_alloc_array(elem2d, nelems2d);
    selem2d_init_array((*elem2d), nelems2d);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void selem2d_copy(SELEM_2D *to, SELEM_2D from) {
    int i=0;
    to->id = from.id;
    to->gid = from.gid;
    to->id_orig = from.id_orig;
    to->id_3d  = from.id_3d;
    to->djac = from.djac;
    to->djac3d = from.djac3d;
    to->djac3d_fixed = from.djac3d_fixed;
    to->interface = from.interface;
    to->flux_elem_tot = from.flux_elem_tot;
    svect_copy_array(&(to->nrml), &(from.nrml), 1);
    to->string = from.string;
    to->mat = from.mat;
    to->bflag = from.bflag;
    to->nnodes = from.nnodes;
    to->nnodes_quad = from.nnodes_quad;
    for (i=0; i<from.nnodes; i++) {
        to->grad_shp[i].x = from.grad_shp[i].x;
        to->grad_shp[i].y = from.grad_shp[i].y;
        to->nodes[i] = from.nodes[i];
        to->levels[i] = from.levels[i];
    }
    to->nedges = from.nedges;
    to->flux_elem = from.flux_elem;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void selem2d_printScreen(SELEM_2D *elem2d) {
    int i;
    printf("\n--------------------------------------------\n");
    printf("2D ELEMENT: local ID: %d :: original ID: %d :: global ID: %d \n",elem2d->id,elem2d->id_orig,elem2d->gid);
    printf("nnodes: %d \t nnodes_quad: %d\n",elem2d->nnodes, elem2d->nnodes_quad);
    printf("djac 2d: %20.10f  djac 3d: %20.10f  djac 2d fixed: %20.10f \n",elem2d->djac,elem2d->djac3d,elem2d->djac3d_fixed);
    printf("element normal: "); svect_printScreen(elem2d->nrml,"nrml");
    printf("string: %d\n",elem2d->string);
    printf("material id: %d\n",elem2d->mat);
    printf("boundary flag: %d\n",elem2d->bflag);
    printf("owning 3d element id: %d\n",elem2d->id_3d);
    printf("original id: %d\n",elem2d->id_orig);
    
    printf("interface: %d\n",elem2d->interface);
    printf("flux_elem_tot: %d\n",elem2d->flux_elem_tot);
    printf("original id: %d\n",elem2d->id_orig);
    
    for (i=0; i<elem2d->nnodes; i++) {
        printf("2d element local node id: %d :: global id: %d grad_shp.x: %15.10e grad_shp.y: %15.10e levels: %d\n",i,elem2d->nodes[i],elem2d->grad_shp[i].x,elem2d->grad_shp[i].y,elem2d->levels[i]);
    }
    printf("nedges: %d\n",elem2d->nedges);
    printf("--------------------------------------------\n");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SVECT selem2d_get_triangle_centroid(SNODE nd1, SNODE nd2, SNODE nd3) {
    SVECT center;
    center.x = one_3 * (nd1.x + nd2.x + nd3.x);
    center.y = one_3 * (nd1.y + nd2.y + nd3.y);
    center.z = one_3 * (nd1.z + nd2.z + nd3.z);
    return center;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Returns a double which established which way the nodes on the triangule are numbers according to a
// given reference vector
// Return variable < 0 :: counter clock-wise node numbering w.r.t. the reference vector
double selem2d_get_triangle_orientation(SNODE nd1, SNODE nd2, SNODE nd3, SVECT ref_vec) {
    
    SVECT side1, side2;        /* two sides - the cross product is the normal */
    double magnitude;        /* the magnitude of the normal */
    
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

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// these are all constants on a triangle and can be calculated once and stored if grid points do not move
// note :: nd_SNODE is a global node vector
void selem2d_get_triangle_linear_djac_nrml_gradPhi(SELEM_2D *elem2d, SNODE *nd_SNODE, SVECT *nd_SVECT) {
    
    SVECT side1, side2;        /* two sides - the cross product is the normal */
    double magnitude;        /* the magnitude of the normal */
    
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

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// returns a normal vector to the 2D element plane
SVECT selem2d_get_elem2d_normals(SVECT *nd) {
    SVECT side1, side2;            /* two sides - the cross product is the normal */
    double magnitude = 0.;        /* the magnitude of the normal */
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

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// returns the 2D projected quadrilateral Jacobian
double selem2d_get_quadrilateral_linear_djac2d(double xhat, double yhat, SVECT *nd) {
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

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// returns the 2D projected quadrilateral Jacobian and cartesian space shape function gradients
double selem2d_get_quadrilateral_linear_djac_gradPhi(double xhat, double yhat, SVECT *nd, SVECT *grad_shp) {
    
    double djac2d = selem2d_get_quadrilateral_linear_djac2d(xhat, yhat, nd);
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
// local quadratic basis function values
//  2 o
//    |\
//  5 x  x 4
//    |   \
//  0 o--x--o 1
//       3
void selem2d_get_triangle_local_shape_quad(double xhat, double yhat, double zhat, double *lshape) {
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
void selem2d_get_quadrilateral_local_shape_quad(double xhat, double yhat, double zhat, double *lshape) {
    
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
void selem2d_get_triangle_local_shape_gradients(SVECT *lgrad_shp) {
    lgrad_shp[0].x =  -1.;
    lgrad_shp[0].y =  -1.;
    lgrad_shp[1].x =   1.;
    lgrad_shp[1].y =   0.;
    lgrad_shp[2].x =   0.;
    lgrad_shp[2].y =   1.;
}

void selem2d_get_quadrilateral_local_shape_gradients(double xhat, double yhat, double zhat, SVECT *lgrad_shp) {
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
double selem2d_get_elem2d_area2d(SELEM_2D *elem2d, SVECT *nds) {
    if (elem2d->nnodes == NDONTRI) {
        return elem2d->djac;  // djac2d is actually the area in AdH
    } else {
        double d1 = nds[1].x - nds[3].x;
        double d2 = nds[0].x - nds[2].x;
        double t1 = d1 * nds[0].y;
        double t2 = d2 * nds[1].y;
        double t3 = d1 * nds[2].y;
        double t4 = d2 * nds[3].y;
        return (0.5*(-t1 + t2 + t3 - t4));
    }
}
