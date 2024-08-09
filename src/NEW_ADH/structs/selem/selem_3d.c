#include "adh.h"
#define SCALE_FOR_MASTER_TET_VOLUME (1.0/6.0) // tetrahedral djac to volume scale

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void selem3d_alloc(SELEM_3D *elem3d, int nnodes_on_elem) {
    assert(nnodes_on_elem == 4 || nnodes_on_elem == 6);
    elem3d->nnodes = nnodes_on_elem;
    elem3d->grad_shp = (SVECT *) tl_alloc(sizeof(SVECT), nnodes_on_elem);
    elem3d->nodes = (int *) tl_alloc(sizeof(int), nnodes_on_elem);
    elem3d->levels = (int *) tl_alloc(sizeof(int), nnodes_on_elem);
    int i;
    for (i=0; i<nnodes_on_elem; i++) {
        elem3d->grad_shp[i].x = 0.0;
        elem3d->grad_shp[i].y = 0.0;
        elem3d->grad_shp[i].z = 0.0;
        elem3d->nodes[i] = UNSET_INT;
        elem3d->levels[i] = 0;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void selem3d_load(SELEM_3D *elem3d, int gid, int lid, int elem_nnodes, int *local_node_ids, int column, SVECT *nds, int mat) {
    
    int i;
    
    elem3d->id = lid;
    elem3d->gid = gid;
    elem3d->id_orig = lid; // hmmm, what if we call this later?
    elem3d->icol = column;
    elem3d->mat = mat;
    elem3d->nnodes = elem_nnodes;
    selem3d_alloc(elem3d, elem_nnodes);
    for (i=0;i<elem_nnodes;i++) elem3d->nodes[i] = local_node_ids[i];
   
    if (elem_nnodes == 4) {
        elem3d->nedges = 6;
        elem3d->nnodes_quad = 10;
        
        // calculate constant jacobians, normal and basis gradients
        selem3d_get_tet_linear_djac_gradPhi(elem3d, NULL, nds);
        
        // node ID direction check
        if (elem3d->djac < SMALL6) {
            fprintf(stderr, "ERROR :: Improperly numbered tetrahedron gid = %d || Nodes %d %d %d %d || Jacobian is: %20.10f \n",
                    gid+1,elem3d->nodes[0],elem3d->nodes[1],elem3d->nodes[2],elem3d->nodes[3],elem3d->djac);
            tl_error("ERROR: Improperly numbered tetrahedron");
        }
        
    } else {

        elem3d->nedges = 9;
        elem3d->nnodes_quad = 16;
    }
    
    elem3d->volume = selem3d_get_elem3d_volume(nds,elem_nnodes);
    //printf("elem3d->volume: %f\n",elem3d->volume);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void selem3d_free(SELEM_3D *elem3d) {
    assert(elem3d->nnodes > 0);
    elem3d->grad_shp = (SVECT *) tl_free(sizeof(SVECT), elem3d->nnodes, elem3d->grad_shp);
    elem3d->nodes = (int *) tl_free(sizeof(int), elem3d->nnodes, elem3d->nodes);
    elem3d->levels = (int *) tl_free(sizeof(int), elem3d->nnodes, elem3d->levels);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void selem3d_alloc_array(SELEM_3D **elem3d, int nelems3d) {
    assert(nelems3d > 0);
    (*elem3d) = (SELEM_3D *) tl_alloc(sizeof(SELEM_3D), nelems3d);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void selem3d_free_array(SELEM_3D *elem3d, int nelems3d) {
    assert(nelems3d > 0);
    int ie;
    for (ie=0; ie<nelems3d; ie++) {
        selem3d_free(&(elem3d[ie]));
    }
    elem3d = (SELEM_3D *) tl_free(sizeof(SELEM_3D), nelems3d, elem3d);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void selem3d_init(SELEM_3D *elem3d) {
    elem3d->id = UNSET_INT;
    elem3d->gid = UNSET_INT;
    elem3d->id_orig = UNSET_INT;
    elem3d->djac = UNSET_FLT;
    elem3d->string = UNSET_INT;
    elem3d->mat = UNSET_INT;
    elem3d->nvars = 0;
    elem3d->vars = NULL;
    elem3d->icol = UNSET_INT;
    elem3d->elem2d_sur = UNSET_INT;
    elem3d->elem2d_bed = UNSET_INT;
    elem3d->error = 0.;
    elem3d->flux_ptr = UNSET_INT;
    elem3d->nnodes = UNSET_INT;
    elem3d->nnodes_quad = UNSET_INT;
    elem3d->nedges = UNSET_INT;
    elem3d->interface=0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void selem3d_init_array(SELEM_3D *elem3d, int nelems3d) {
    int ie=0;
    for (ie=0; ie<nelems3d; ie++) {
        selem3d_init(&(elem3d[ie]));
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void selem3d_init_alloc_array(SELEM_3D **elem3d, int nelems3d) {
    selem3d_alloc_array(elem3d, nelems3d);
    selem3d_init_array((*elem3d), nelems3d);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void selem3d_copy(SELEM_3D *to, SELEM_3D from) {
    to->id = from.id;
    to->gid = from.gid;
    to->id_orig = from.id_orig;
    to->nnodes = from.nnodes;
    to->nnodes_quad = from.nnodes_quad;
    to->djac = from.djac;
    to->string = from.string;
    to->mat = from.mat;
    to->icol = from.icol;
    to->elem2d_sur = from.elem2d_sur;
    to->elem2d_bed = from.elem2d_bed;
    to->error = from.error;
    to->flux_ptr = from.flux_ptr;
    to->interface = from.interface;
    int i;
    for (i=0; i<from.nnodes; i++) {
        to->grad_shp[i].x = from.grad_shp[i].x;
        to->grad_shp[i].y = from.grad_shp[i].y;
        to->grad_shp[i].z = from.grad_shp[i].z;
        to->nodes[i] = from.nodes[i];
        to->levels[i] = from.levels[i];
    }
    to->nedges = from.nedges;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void selem3d_printScreen(SELEM_3D *elem3d) {
    int i;
    printf("\n");
    printf("3D ELEMENT: local ID: %d :: original ID: %d :: global ID: %d ---------\n",elem3d->id,elem3d->id_orig,elem3d->gid);
    printf("nnodes: %d \t nnodes_quad: %d \n",elem3d->nnodes, elem3d->nnodes_quad);
    printf("djac: %20.10e \n",elem3d->djac);
    printf("string: %d\n",elem3d->string);
    printf("material id: %d\n",elem3d->mat);
    printf("column id: %d\n",elem3d->icol);
    printf("interface: %d\n",elem3d->interface);
    printf("2d surface element: %d\n", elem3d->elem2d_sur);
    printf("2d bottom element: %d\n", elem3d->elem2d_bed);
    printf("max residual error for all equations: %20.10e\n",elem3d->error);
    for (i=0; i<elem3d->nnodes; i++) {
        printf("element local node id: %d :: global id: %d grad_shp.x: %15.10e grad_shp.y: %15.10e grad_shp.z: %15.10e level: %d\n",i,elem3d->nodes[i],elem3d->grad_shp[i].x, elem3d->grad_shp[i].y, elem3d->grad_shp[i].z, elem3d->levels[i]);
    }
    printf("nedges: %d\n",elem3d->nedges);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// local basis function values
void selem3d_get_tet_local_shape(double xhat, double yhat, double zhat, double *lshape) {
    lshape[0] = 1-xhat-yhat-zhat;
    lshape[1] = xhat;
    lshape[2] = yhat;
    lshape[3] = zhat;
}

void selem3d_get_triprism_local_shape(double xhat, double yhat, double zhat, double *lshape) {
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
void selem3d_get_tet_local_shape_quad(double xhat, double yhat, double zhat, double *lshape_quad) {
    double lshape[4]; sarray_init_dbl(lshape, 4);
    selem3d_get_tet_local_shape(xhat, yhat, zhat, lshape);
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

void selem3d_get_triprism_local_shape_quad(double xhat, double yhat, double zhat, double *lshape_quad) {
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
void selem3d_get_tet_local_shape_gradients(SVECT *lgrad_shp) {
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
void selem3d_get_triprism_local_shape_gradients(double xhat, double yhat, double zhat, SVECT *lgrad_shp) {
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
double selem3d_get_triprism_djac(double xhat, double yhat, double zhat, SVECT *nd)  {
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
double selem3d_get_triprism_linear_djac_gradPhi(double xhat, double yhat, double zhat, SVECT *nd, SVECT *grad_shp) {
    
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
void selem3d_get_tet_linear_djac_gradPhi(SELEM_3D *elem3d, SNODE *nd_SNODE, SVECT *nd_SVECT) {

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
double selem3d_get_tet_linear_djac(SNODE *nd_SNODE, SVECT *nd_SVECT) {
    
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
double selem3d_get_tet_linear_djac_gradPhi2(SNODE *nd_SNODE, SVECT *nd_SVECT, SVECT *grad_shp) {
    
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
double selem3d_get_elem3d_volume(SVECT *node, int nnodes) {
    if (nnodes == NDONTET) {
        double d1 = node[1].x - node[2].x, d2 = node[0].x - node[2].x, d3 = node[0].x - node[1].x;
        double d4 = node[2].x - node[3].x, d5 = node[1].x - node[3].x, d6 = node[0].x - node[3].x;
        return (one_6 * ((d4*node[1].y - d5*node[2].y + d1*node[3].y)*node[0].z -
                         (d4*node[0].y - d6*node[2].y + d2*node[3].y)*node[1].z +
                         (d5*node[0].y - d6*node[1].y + d3*node[3].y)*node[2].z -
                         (d1*node[0].y - d2*node[1].y + d3*node[2].y)*node[3].z));
    } else if (nnodes == NDONPRISM) {
        return selem3d_get_triprism_volume(node);
    } else {
        tl_error("bad call to elem2d_get_volume");
    }
    return (-1);
}

double selem3d_get_triprism_volume(SVECT *node) {
    double d1 = node[1].x - node[2].x, d2 = node[0].x - node[2].x, d3 = node[0].x - node[1].x;
    double prism_term = d1*node[0].y - d2*node[1].y + d3*node[2].y;
    double t1 = node[0].z + node[1].z + node[2].z;
    t1 -= (node[3].z + node[4].z + node[5].z);
    return (one_6 * prism_term * t1);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
