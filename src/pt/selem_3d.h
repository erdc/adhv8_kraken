#ifndef H_SELEM_3D_
#define H_SELEM_3D_

// dependencies :: SVECT, SVECT2D

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
    int id;                   /* current 2d element number */
    int nnodes;               /* total # of nodes on element */
    double djac;              /* jacobian given only 2D coordinates */
    int mat;                  /* 2d element material type */
    int *nodes;               /* the nodes in the element */
    int nadjc_elems;          /* the number of adjacent elements */
    int *adjc_elems;          /* and array of 2d element IDs adjacent to the element */
    int face_flag[4];         // the element ID across this face, UNSET_INT for external face
    int elem2d_face[4];       // the 2d element ID for each face of a tet
    int edge_flag[6];
    int nedge2elem[6];        // holds the count of the number of edge-adjacent elemnts to each edge of this element
    int **edge2elem;          // maps all adjacent elements to each edge of the 3d element
    int elem2d_surface;       // for ** columnar grids **, this it the 2d element surface ID
} SELEM_3D;

/* struct methods -------------------------------------- */

void selem3d_alloc(SELEM_3D *elem3d, int nnodes_on_elem);
void selem3d_free(SELEM_3D *elem3d);
void selem3d_alloc_array(SELEM_3D **elem3d, int nelems3d);
void selem3d_free_array(SELEM_3D *elem3d, int nelems3d);
void selem3d_init(SELEM_3D *elem3d);
void selem3d_init_array(SELEM_3D *elem3d, int nelems3d);
void selem3d_init_alloc_array(SELEM_3D **elem3d, int nelems3d);
void selem3d_copy(SELEM_3D *to, SELEM_3D from);
void selem3d_printScreen(SELEM_3D *elem3d);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
