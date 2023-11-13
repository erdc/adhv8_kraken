#include "global_header.h"

/***********************************************************/
/***********************************************************/
void selem2d_alloc(SELEM_2D *elem2d, int nnodes_on_elem) {
    int i;
    
    assert(nnodes_on_elem == 3 || nnodes_on_elem == 4);
    
    elem2d->nnodes = nnodes_on_elem;
    elem2d->nodes = (int *) calloc(nnodes_on_elem, sizeof(int));
    for (i=0; i<nnodes_on_elem; i++) {
        elem2d->nodes[i] = UNSET_INT;
    }
    
    //elem2d->nadjc_elems = nadjc_elems;
    //elem2d->adjc_elems = (int *) calloc(nadjc_elems, sizeof(int));
    //for (i=0; i<nadjc_elems; i++) {
    //    elem2d->adjc_elems[i] = UNSET_INT;
    //}
}

/***********************************************************/
/***********************************************************/
void selem2d_free(SELEM_2D *elem2d) {
    
    free(elem2d->nodes);
    free(elem2d->adjc_elems);
}

/***********************************************************/
/***********************************************************/
void selem2d_alloc_array(SELEM_2D **elem2d, int nelems2d) {
    assert(nelems2d > 0);
    (*elem2d) = (SELEM_2D *) calloc(nelems2d, sizeof(SELEM_2D));
}

/***********************************************************/
/***********************************************************/
void selem2d_free_array(SELEM_2D *elem2d, int nelems2d) {
    assert(nelems2d > 0);
    int ie;
    for (ie=0; ie<nelems2d; ie++) {
        selem2d_free(&(elem2d[ie]));
    }
    free(elem2d);
}

/***********************************************************/
/***********************************************************/
void selem2d_init(SELEM_2D *elem2d) {
    elem2d->id = -1;
    elem2d->djac = 0.0;
    elem2d->djac3d = 0.0;
    svect_init(&(elem2d->nrml));
    elem2d->mat = -1;
    elem2d->bflag = UNSET_INT;
    elem2d->nnodes = -1;
    elem2d->nodes = NULL;
    elem2d->nadjc_elems = -1;
    elem2d->adjc_elems = NULL;
    elem2d->edge_flag[0] = UNSET_INT;
    elem2d->edge_flag[1] = UNSET_INT;
    elem2d->edge_flag[2] = UNSET_INT;
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
    to->djac = from.djac;
    to->djac3d = from.djac3d;
    svect_copy_array(&(to->nrml), &(from.nrml), 1);
    to->mat = from.mat;
    to->bflag = from.bflag;
    for (i=0; i<from.nnodes; i++) {
        to->nodes[i] = from.nodes[i];
    }
}

/***********************************************************/
/***********************************************************/
void selem2d_printScreen(SELEM_2D *elem2d) {
    int i;
    printf("\n--------------------------------------------\n");
    printf("2D ELEMENT: ID: %d \n",elem2d->id);
    printf("nnodes: %d \n",elem2d->nnodes);
    printf("djac 2d: %20.10f  djac 3d: %20.10f  \n",elem2d->djac,elem2d->djac3d);
    printf("element normal: "); svect_printScreen(elem2d->nrml,"nrml");
    printf("material id: %d\n",elem2d->mat);
    printf("boundary flag: %d\n",elem2d->bflag);
    
    for (i=0; i<elem2d->nnodes; i++) {
        printf("2d element %d local node id: %d \n",i,elem2d->nodes[i]);
    }
    printf("--------------------------------------------\n");
}

/***********************************************************/
/***********************************************************/
