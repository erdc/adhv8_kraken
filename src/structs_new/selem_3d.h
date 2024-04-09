#ifndef H_SELEM_3D_
#define H_SELEM_3D_

// dependencies: SVECT

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {

    int id;                   /* current 3d element number */
    int gid;                  /* global id should == id_orig for non-mpi */
    int id_orig;              /* original 2d element number (every 2d element should have a value) */
    int nnodes;               /* total # of nodes on element */
    int nnodes_quad;          /* the total number of quadratic nodes on element */
    double djac;              /* the jacobian of the element */
    int string;               /* string number associated with this element */
    int mat;                  /* 2d element material type */
    int my_pe;                /* owning processor */
    int icol;                 /* the column this tet belongs to */
    int elem2d_sur;           /* the 2d surface element of the column this tet belongs to */
    int elem2d_bed;           /* the 2d bed element of the column this tet belongs to */
    double error;             /* max elemental error in continuity equations for the current ts */
    int flux_ptr;             /* 2d element pointer for flux line output */
    SVECT *grad_shp;          /* the gradients of the shape functions */
    int *nodes;               /* the nodes in the element */
    int *levels;              /* the node levels in the element */
    int nedges;               /* the number of edges on the element */
    int **edges;              /* the 3D element local nodes IDs for each edge*/
    
    /*interface flag to prevent adaption */
    int interface;

} SELEM_3D;

/*********************************************************/
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
