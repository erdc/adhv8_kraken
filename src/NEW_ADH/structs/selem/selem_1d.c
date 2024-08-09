#include "adh.h"

/***********************************************************/
/***********************************************************/
void selem1d_alloc(SELEM_1D *elem1d) {
    elem1d->nnodes = 2;
    elem1d->grad_shp = (double *) tl_alloc(sizeof(double), 2);
    elem1d->nodes =  (int *) tl_alloc(sizeof(int), 2);
    elem1d->levels = (int *) tl_alloc(sizeof(int), 2);
    int i;
    for (i=0; i<2; i++) {
        elem1d->grad_shp[i] = UNSET_FLT;
        elem1d->nodes[i] = UNSET_INT;
        elem1d->levels[i] = 0;
    }
}

/***********************************************************/
/***********************************************************/
void selem1d_alloc_array(SELEM_1D **elem1d, int nelems1d) {
    assert(nelems1d > 0);
    (*elem1d) = (SELEM_1D *) tl_alloc(sizeof(SELEM_1D), nelems1d);
}

/***********************************************************/
/***********************************************************/
void selem1d_load(SELEM_1D *elem1d, int gid, int lid, int elem_nnodes, int *local_node_ids, int bflag, SVECT *nds, int mat) {
    
    assert(elem_nnodes == 2);
    
    elem1d->id = lid;
    elem1d->gid = gid;
    elem1d->id_orig = lid; // hmmm, what if we call this later?
    elem1d->bflag = bflag;
    elem1d->mat = mat;
    elem1d->nnodes = elem_nnodes;
    selem1d_alloc(elem1d);
    elem1d->nodes[0] = local_node_ids[0];
    elem1d->nodes[1] = local_node_ids[1];
    
    /* computes the jacobian */
    elem1d->djac = DIST_2D(nds[0], nds[1]);
    elem1d->length = elem1d->djac;
    if (elem1d->djac < SMALL6) {
        fprintf(stderr,"ERROR :: Improperly numbered line segment GID = %d || Nodes %d %d Jacobian is: %20.10f \n",
                gid,elem1d->nodes[0],elem1d->nodes[1],elem1d->djac);
        tl_error("ERROR: Improperly numbered line segment given");
    }
    
    /* computes the gradient of the shape functions */
    elem1d->grad_shp[0] = -1.0 / elem1d->djac;
    elem1d->grad_shp[1] =  1.0 / elem1d->djac;
    
    /* computes 1D element normal in the 2D plane */
    elem1d->nrml.x =   (nds[1].y - nds[0].y) / elem1d->djac;
    elem1d->nrml.y =  -(nds[1].x - nds[0].x) / elem1d->djac;
}

/***********************************************************/
/***********************************************************/
void selem1d_realloc_array(SELEM_1D **elem1d, int nelems1d_old, int nelems1d_new) {
    (*elem1d) = (SELEM_1D *) tl_realloc(sizeof(SELEM_1D), nelems1d_new, nelems1d_old, (*elem1d) );
    int ie;
    for (ie=nelems1d_old; ie<nelems1d_new; ie++) {
        (*elem1d)[ie].nnodes = 2;
        (*elem1d)[ie].grad_shp = (double *) tl_alloc(sizeof(double), (*elem1d)[ie].nnodes);
        (*elem1d)[ie].nodes = (int *) tl_alloc(sizeof(int), (*elem1d)[ie].nnodes);
        (*elem1d)[ie].levels = (int *) tl_alloc(sizeof(int), (*elem1d)[ie].nnodes);
    }
}

/***********************************************************/
/***********************************************************/
void selem1d_free_array(SELEM_1D *elem1d, int nelems1d) {
    int ie;
    for (ie=0; ie<nelems1d; ie++) {
        elem1d[ie].grad_shp = (double *) tl_free(sizeof(double), elem1d[ie].nnodes, elem1d[ie].grad_shp);
        elem1d[ie].nodes = (int *) tl_free(sizeof(int), elem1d[ie].nnodes, elem1d[ie].nodes);
        elem1d[ie].levels = (int *) tl_free(sizeof(int), elem1d[ie].nnodes, elem1d[ie].levels);
    }
    elem1d = (SELEM_1D *) tl_free(sizeof(SELEM_1D), nelems1d, elem1d);
}

/***********************************************************/
/***********************************************************/
void selem1d_init(SELEM_1D *elem1d) {
    int i = 0;
    elem1d->id = UNSET_INT;
    elem1d->gid = UNSET_INT;
    elem1d->id_orig = UNSET_INT;
    elem1d->elem2d = UNSET_INT;
    elem1d->djac = UNSET_FLT;
    elem1d->nrml.x = UNSET_FLT;
    elem1d->nrml.y = UNSET_FLT;
    elem1d->string = UNSET_INT;
    elem1d->mat = UNSET_INT;
    elem1d->nvars = 0;
    elem1d->vars = NULL;
    for (i = 0; i < elem1d->nnodes; i++) {
        elem1d->grad_shp[i] = UNSET_FLT;
        elem1d->nodes[i] = UNSET_INT;
        elem1d->levels[i] = 0;
    }
    elem1d->flux_elem = NULL;
}

/***********************************************************/
/***********************************************************/
void selem1d_init_array(SELEM_1D *elem1d, int nelems1d) {
    int ie=0;
    for (ie=0; ie<nelems1d; ie++) {
        selem1d_init(&(elem1d[ie]));
    }
}

/***********************************************************/
/***********************************************************/
void selem1d_init_alloc_array(SELEM_1D **elem1d, int nelems1d) {
    selem1d_alloc_array(elem1d, nelems1d);
    selem1d_init_array((*elem1d), nelems1d);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
void selem1d_copy(SELEM_1D *to, SELEM_1D from) {
    to->id = from.id;
    to->gid = from.gid;
    to->id_orig = from.id_orig;
    to->elem2d = from.elem2d;
    to->djac = from.djac;
    to->nrml.x = from.nrml.x;
    to->nrml.y = from.nrml.y;
    to->string = from.string;
    to->mat = from.mat;
    to->nnodes = from.nnodes;
    int i;
    for (i = 0; i < from.nnodes; i++) {
        to->grad_shp[i] = from.grad_shp[i];
        to->nodes[i] = from.nodes[i];
        to->levels[i] = from.levels[i];
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
void selem1d_printScreen(SELEM_1D *elem1d) {
    int i = -1;
    
    printf("\n1D ELEMENT ------------------------------------------\n");
    printf("local id: %d \t global id: %d \t 2d element: %d \t string: %d \t material: %d\n", elem1d->id, elem1d->gid, elem1d->elem2d, elem1d->string, elem1d->mat);
    printf("nrml.x: %10.5f nrml.y: %10.5f\n", elem1d->nrml.x, elem1d->nrml.y);
    for (i = 0; i < elem1d->nnodes; i++) {
        printf("local id %d: \t global id: %d \t adaption level: %d \t gradient: %10.5f\n", i, elem1d->nodes[i], elem1d->levels[i], elem1d->grad_shp[i]);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
void selem1d_outward_nrml(SELEM_1D *elem1d, SNODE *node, SELEM_2D elem2d){
    
    double distance;            /* used as a distance meas and also the dot product */
    double x, x1, x2;           /* x measures */
    double y, y1, y2;           /* y measures */
    double delx;                /* x difference */
    double dely;                /* y difference */
    int ie;                     /* loop counter over the elements */
    
    /* the outward normal in the z plane is needed to do the overland flow calculations */
    /* right now this is only run in the groundwater flow type */
    
    if (elem1d->string != UNSET_INT) {
        /* here we calculate the z plane outward normal to the 1d element */
        x1 = node[elem1d->nodes[0]].x;
        x2 = node[elem1d->nodes[1]].x;
        y1 = node[elem1d->nodes[0]].y;
        y2 = node[elem1d->nodes[1]].y;
        delx = x2 - x1;
        dely = y2 - y1;
        distance = sqrt(delx * delx + dely * dely);
        elem1d->nrml.x = dely / distance;
        elem1d->nrml.y = -delx / distance;
        
        /* we only know that we are normal to the 1d element */
        /* we must now make sure that the normal points outward */
        /* to do this we find the 2d centroid of the associated 2d element */
        /* then we contruct a vector from the first node in the element to this centroid */
        /* the dot product of this with our normal vector should be negative */
        /* if not reverse the sign on the normal */
        
        x = one_3 * (node[elem2d.nodes[0]].x + node[elem2d.nodes[1]].x + node[elem2d.nodes[2]].x);
        y = one_3 * (node[elem2d.nodes[0]].y + node[elem2d.nodes[1]].y + node[elem2d.nodes[2]].y);
        delx = x - x1;
        dely = y - y1;
        distance = delx * elem1d->nrml.x + dely * elem1d->nrml.y;
        
        if (distance > 0.) {
            elem1d->nrml.x *= -1.;
            elem1d->nrml.y *= -1.;
        }
        
    }
}

/***********************************************************/
/***********************************************************/

void selem1d_get_elem1d_linear_djac_gradPhi(SELEM_1D *elem1d, SVECT nd1, SVECT nd2) {
    
    
    /* computes the jacobian */
    elem1d->djac = DIST_2D(nd1, nd2);
    
    /* computes the gradient of the shape functions */
    elem1d->grad_shp[0] = -1.0 / elem1d->djac;
    elem1d->grad_shp[1] =  1.0 / elem1d->djac;
}

/***********************************************************/


