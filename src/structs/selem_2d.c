#include "global_header.h"

/***********************************************************/
/***********************************************************/
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

/***********************************************************/
/***********************************************************/
void selem2d_free(SELEM_2D *elem2d) {
    assert(elem2d->nnodes > 0);
    elem2d->grad_shp = (SVECT2D *) tl_free(sizeof(SVECT2D), elem2d->nnodes, elem2d->grad_shp);
    elem2d->nodes = (int *) tl_free(sizeof(int), elem2d->nnodes, elem2d->nodes);
    elem2d->levels = (int *) tl_free(sizeof(int), elem2d->nnodes, elem2d->levels);
}

/***********************************************************/
/***********************************************************/
void selem2d_alloc_array(SELEM_2D **elem2d, int nelems2d) {
    assert(nelems2d > 0);
    (*elem2d) = (SELEM_2D *) tl_alloc(sizeof(SELEM_2D), nelems2d);
}

/***********************************************************/
/***********************************************************/
void selem2d_free_array(SELEM_2D *elem2d, int nelems2d) {
    assert(nelems2d > 0);
    int ie;
    for (ie=0; ie<nelems2d; ie++) {
        selem2d_free(&(elem2d[ie]));
    }
    elem2d = (SELEM_2D *) tl_free(sizeof(SELEM_2D), nelems2d, elem2d);
}

/***********************************************************/
/***********************************************************/
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
    elem2d->resident_pe = UNSET_INT;
    svect_init(&(elem2d->nrml));
    elem2d->string = UNSET_INT;
    elem2d->mat = UNSET_INT;
    elem2d->my_pe = UNSET_INT;
    elem2d->bflag = UNSET_INT;
    elem2d->nedges = UNSET_INT;
    elem2d->nnodes_quad = UNSET_INT;
    elem2d->edges = NULL;
    elem2d->flux_elem = NULL;
}

/***********************************************************/
/***********************************************************/
void selem2d_init_array(SELEM_2D *elem2d, int nelems2d) {
    int ie=0;
    for (ie=0; ie<nelems2d; ie++) {
        selem2d_init(&(elem2d[ie]));
    }
}

/***********************************************************/
/***********************************************************/
void selem2d_init_alloc_array(SELEM_2D **elem2d, int nelems2d) {
    selem2d_alloc_array(elem2d, nelems2d);
    selem2d_init_array((*elem2d), nelems2d);
}

/***********************************************************/
/***********************************************************/
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
    to->resident_pe = from.resident_pe ;
    svect_copy_array(&(to->nrml), &(from.nrml), 1);
    to->string = from.string;
    to->mat = from.mat;
    to->my_pe = from.my_pe;
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
    to->edges = from.edges;
    to->flux_elem = from.flux_elem;
}

/***********************************************************/
/***********************************************************/
void selem2d_printScreen(SELEM_2D *elem2d) {
    int i;
    printf("\n--------------------------------------------\n");
    printf("2D ELEMENT: local ID: %d :: original ID: %d :: global ID: %d \n",elem2d->id,elem2d->id_orig,elem2d->gid);
    printf("nnodes: %d \t nnodes_quad: %d\n",elem2d->nnodes, elem2d->nnodes_quad);
    printf("djac 2d: %20.10f  djac 3d: %20.10f  djac 2d fixed: %20.10f \n",elem2d->djac,elem2d->djac3d,elem2d->djac3d_fixed);
    printf("element normal: "); svect_printScreen(elem2d->nrml,"nrml");
    printf("string: %d\n",elem2d->string);
    printf("material id: %d\n",elem2d->mat);
    printf("owning processor: %d\n",elem2d->my_pe);
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
    if (elem2d->nedges > 0) {
        for (i=0;i<elem2d->nedges;i++) {
          printf("edge: %d \t nd1: %d \t nd2: %d \n",i,elem2d->edges[i][0],elem2d->edges[i][1]);
        }
    }
    printf("--------------------------------------------\n");
}

/***********************************************************/
/***********************************************************/
