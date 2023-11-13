#include "global_header.h"

/***********************************************************/
/***********************************************************/
void selem3d_alloc(SELEM_3D *elem3d, int nnodes_on_elem) {
    int i;
    
    assert(nnodes_on_elem == 4 || nnodes_on_elem == 6);
    
    elem3d->nnodes = nnodes_on_elem;
    elem3d->nodes = (int *) calloc(nnodes_on_elem, sizeof(int));
    for (i=0; i<nnodes_on_elem; i++) {
        elem3d->nodes[i] = UNSET_INT;
    }
    
    //elem3d->nadjc_elems = nadjc_elems;
    //elem3d->adjc_elems = (int *) calloc(nadjc_elems, sizeof(int));
    //for (i=0; i<nadjc_elems; i++) {
    //    elem3d->adjc_elems[i] = UNSET_INT;
    //}
}

/***********************************************************/
/***********************************************************/
void selem3d_free(SELEM_3D *elem3d) {
    
    free(elem3d->nodes);
    free(elem3d->adjc_elems);
}

/***********************************************************/
/***********************************************************/
void selem3d_alloc_array(SELEM_3D **elem3d, int nelems3d) {
    assert(nelems3d > 0);
    (*elem3d) = (SELEM_3D *) calloc(nelems3d, sizeof(SELEM_3D));
}

/***********************************************************/
/***********************************************************/
void selem3d_free_array(SELEM_3D *elem3d, int nelems3d) {
    assert(nelems3d > 0);
    int ie;
    for (ie=0; ie<nelems3d; ie++) {
        selem3d_free(&(elem3d[ie]));
    }
    free(elem3d);
}

/***********************************************************/
/***********************************************************/
void selem3d_init(SELEM_3D *elem3d) {
    elem3d->id = -1;
    elem3d->djac = 0.0;
    elem3d->mat = -1;
    elem3d->nnodes = -1;
    elem3d->nodes = NULL;
    elem3d->nadjc_elems = -1;
    elem3d->adjc_elems = NULL;
    elem3d->elem2d_surface = -1;

    int i;
    for (i=0; i<4; i++) {
      elem3d->face_flag[i] = UNSET_INT;         // the element ID across this face, UNSET_INT for external face
      elem3d->elem2d_face[i] = UNSET_INT;       // the 2d element ID for each face of a tet
    }
    for (i=0; i<6; i++) {
      elem3d->edge_flag[i] = UNSET_INT;
    };
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
    int i=0;
    to->id = from.id;
    to->djac = from.djac;
    to->mat = from.mat;
    to->elem2d_surface = from.elem2d_surface;
    for (i=0; i<from.nnodes; i++) {
        to->nodes[i] = from.nodes[i];
    }
}

/***********************************************************/
/***********************************************************/
void selem3d_printScreen(SELEM_3D *elem3d) {
    int i;
    printf("\n--------------------------------------------\n");
    printf("3D ELEMENT ID: %d || 2D ELEMENT SURFACE ID: %d\n",elem3d->id,elem3d->elem2d_surface);
    printf("nnodes: %d \n",elem3d->nnodes);
    printf("djac: %20.10f  \n",elem3d->djac);
    printf("material id: %d\n",elem3d->mat);
    
    for (i=0; i<elem3d->nnodes; i++) {
        printf("3d element %d local node id: %d \n",i,elem3d->nodes[i]);
    }
    printf("--------------------------------------------\n");
}

/***********************************************************/
/***********************************************************/
