
#ifndef H_SELEM_1D_
#define H_SELEM_1D_

// dependencies :: SVECT2D, SNODE, SELEM_2D

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
  int elem;
  int sub;
  int super;
} FLUX_ELEM;

typedef struct {

    int id;             /* current 1d element number */
    int gid;
    int id_orig;        /* original 1d element number (every 1d element should have a value) */
    int nnodes;         /* total # of nodes on element */
    int elem2d;         /* this is the 2d element that has an edge that matches this 1d element */
    double djac;        /* the jacobian of the element */
    SVECT2D nrml;       /* this is the 2d outward normal to the element this is connected to in the x-y plane */
    int string;         /* boundary conditions to to this 1d element if there is one */
    int mat;            /* material id */
    double *grad_shp;   /* the gradients of the shape functions */
    int *nodes;         /* the nodes in the element */
    int *levels;        /* the node levels in the element */
    FLUX_ELEM *flux_elem;
    int flux_elem_tot;

} SELEM_1D;

/*********************************************************/
/* struct methods -------------------------------------- */

void selem1d_alloc_array(SELEM_1D **elem1d, int nelems1d, int nnodes);
void selem1d_realloc_array(SELEM_1D **elem1d, int nelems1d_old, int nelem1d_new, int nnodes);
void selem1d_free_array(SELEM_1D *elem1d, int nelems1d);
void selem1d_init(SELEM_1D *elem1d);
void selem1d_init_array(SELEM_1D *elem1d, int nelems1d);
void selem1d_init_alloc_array(SELEM_1D **elem1d, int nelems1d, int nnodes);
void selem1d_copy(SELEM_1D *to, SELEM_1D from);
void selem1d_printScreen(SELEM_1D *elem1d);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
