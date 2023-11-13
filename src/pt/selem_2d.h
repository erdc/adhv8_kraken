#ifndef H_SELEM_2D_
#define H_SELEM_2D_

// dependencies :: SVECT, SVECT2D

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
    int id;                   /* current 2d element number */
    int id_3d;                /* the 3d element that has a face coincident with this 2d element */
    int nnodes;               /* total # of nodes on element */
    double djac;              /* jacobian given only 2D coordinates */
    double djac3d;            /* the 2D jacobian surface in 3D */
    SVECT nrml;               /* the normal to the face */
    int mat;                  /* 2d element material type */
    int bflag;                /* element boundary flag :: 0 = surface, 1 = bottom, 2 = sidewall (3d) */
    int *nodes;               /* the nodes in the element */
    int nadjc_elems;          /* the number of adjacent elements */
    int *adjc_elems;          /* and array of 2d element IDs adjacent to the element */
    int edge_flag[3];         // the element ID across this edge, UNSET_INT for external edge
} SELEM_2D;

/*********************************************************/
/* struct methods -------------------------------------- */


void selem2d_alloc(SELEM_2D *elem2d, int nnodes_on_elem);
void selem2d_free(SELEM_2D *elem2d);
void selem2d_alloc_array(SELEM_2D **elem2d, int nelems2d);
void selem2d_free_array(SELEM_2D *elem2d, int nelems2d);
void selem2d_init(SELEM_2D *elem2d);
void selem2d_init_array(SELEM_2D *elem2d, int nelems2d);
void selem2d_init_alloc_array(SELEM_2D **elem2d, int nelems2d);
void selem2d_copy(SELEM_2D *to, SELEM_2D from);
void selem2d_printScreen(SELEM_2D *elem2d);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
