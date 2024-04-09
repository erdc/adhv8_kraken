//#include "global_header.h"
#include "local_header.h"

/***********************************************************/
/***********************************************************/
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

/***********************************************************/
/***********************************************************/
void selem3d_free(SELEM_3D *elem3d) {
    assert(elem3d->nnodes > 0);
    elem3d->grad_shp = (SVECT *) tl_free(sizeof(SVECT), elem3d->nnodes, elem3d->grad_shp);
    elem3d->nodes = (int *) tl_free(sizeof(int), elem3d->nnodes, elem3d->nodes);
    elem3d->levels = (int *) tl_free(sizeof(int), elem3d->nnodes, elem3d->levels);
}

/***********************************************************/
/***********************************************************/
void selem3d_alloc_array(SELEM_3D **elem3d, int nelems3d) {
    assert(nelems3d > 0);
    (*elem3d) = (SELEM_3D *) tl_alloc(sizeof(SELEM_3D), nelems3d);
}

/***********************************************************/
/***********************************************************/
void selem3d_free_array(SELEM_3D *elem3d, int nelems3d) {
    assert(nelems3d > 0);
    int ie;
    for (ie=0; ie<nelems3d; ie++) {
        selem3d_free(&(elem3d[ie]));
    }
    elem3d = (SELEM_3D *) tl_free(sizeof(SELEM_3D), nelems3d, elem3d);
}

/***********************************************************/
/***********************************************************/
void selem3d_init(SELEM_3D *elem3d) {
    elem3d->id = UNSET_INT;
    elem3d->gid=UNSET_INT;
    elem3d->id_orig = UNSET_INT;
    elem3d->djac = UNSET_FLT;
    elem3d->string = UNSET_INT;
    elem3d->mat = UNSET_INT;
    elem3d->my_pe = UNSET_INT;
    elem3d->icol = UNSET_INT;
    elem3d->elem2d_sur = UNSET_INT;
    elem3d->elem2d_bed = UNSET_INT;
    elem3d->error = 0.;
    elem3d->flux_ptr = UNSET_INT;
    elem3d->nnodes_quad = UNSET_INT;
    elem3d->nedges = UNSET_INT;
    elem3d->edges = NULL;
    elem3d->interface=0;
}

/***********************************************************/
/***********************************************************/
void selem3d_init_array(SELEM_3D *elem3d, int nelems3d) {
    int ie=0;
    for (ie=0; ie<nelems3d; ie++) {
        selem3d_init(&(elem3d[ie]));
    }
}

/***********************************************************/
/***********************************************************/
void selem3d_init_alloc_array(SELEM_3D **elem3d, int nelems3d) {
    selem3d_alloc_array(elem3d, nelems3d);
    selem3d_init_array((*elem3d), nelems3d);
}

/***********************************************************/
/***********************************************************/
void selem3d_copy(SELEM_3D *to, SELEM_3D from) {
    to->id = from.id;
    to->gid = from.gid;
    to->id_orig = from.id_orig;
    to->nnodes = from.nnodes;
    to->nnodes_quad = from.nnodes_quad;
    to->djac = from.djac;
    to->string = from.string;
    to->mat = from.mat;
    to->my_pe = from.my_pe;
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
    to->edges = from.edges;
}

/***********************************************************/
/***********************************************************/
void selem3d_printScreen(SELEM_3D *elem3d) {
    int i;
    printf("\n");
    printf("3D ELEMENT: local ID: %d :: original ID: %d :: global ID: %d ---------\n",elem3d->id,elem3d->id_orig,elem3d->gid);
    printf("nnodes: %d \t nnodes_quad: %d \n",elem3d->nnodes, elem3d->nnodes_quad);
    printf("djac: %20.10e \n",elem3d->djac);
    printf("string: %d\n",elem3d->string);
    printf("material id: %d\n",elem3d->mat);
    printf("owning processor: %d\n",elem3d->my_pe);
    printf("column id: %d\n",elem3d->icol);
    printf("interface: %d\n",elem3d->interface);
    printf("2d surface element: %d\n", elem3d->elem2d_sur);
    printf("2d bottom element: %d\n", elem3d->elem2d_bed);
    printf("max residual error for all equations: %20.10e\n",elem3d->error);
    for (i=0; i<elem3d->nnodes; i++) {
        printf("element local node id: %d :: global id: %d grad_shp.x: %15.10e grad_shp.y: %15.10e grad_shp.z: %15.10e level: %d\n",i,elem3d->nodes[i],elem3d->grad_shp[i].x, elem3d->grad_shp[i].y, elem3d->grad_shp[i].z, elem3d->levels[i]);
    }
    printf("nedges: %d\n",elem3d->nedges);
    for (i=0;i<elem3d->nedges;i++) {
        printf("edge: %d \t nd1: %d \t nd2: %d \n",i,elem3d->edges[i][0],elem3d->edges[i][1]);
    }
}
