#ifndef H_SELEM_2D_
#define H_SELEM_2D_

// dependencies :: SVECT, SVECT2D

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
    
    int id;                   /* current 2d element number */
    int gid;                  /*global element id */
    int id_orig;              /* original 2d element number (every 2d element should have a value) */
    int id_3d;                /* the 3d element that has a face coincident with this 2d element */
    int nnodes;               /* total # of nodes on element */
    int nnodes_quad;          /* the total # of quadratic nodes on element */
    double djac;              /* jacobian given only 2D coordinates */
    double djac3d;            /* the 2D jacobian surface in 3D */
    double djac3d_fixed;      /* the area of the original element in 3D */
    SVECT nrml;               /* the normal to the face */
    int string;               /* string # associated with this 2d element*/
    int mat;                  /* 2d element physics material type */
    int nvars;          /* the total number of independent variables on this element */
    int *vars;          /* an array of independent variables on this element */
    int bflag;                /* element boundary flag :: -1 = body, 0 = surface, 1 = bottom, 2 = sidewall (3d) */
    SVECT2D *grad_shp;        /* the gradients of the shape functions */
    int *nodes;               /* the nodes in the element */
    int *levels;              /* the node levels in the element */
    int nedges;               /* the number of edges on the element */
    double area;
    int resident_pe;          /* arbitrary since AdH decomp is node-based. Defines as lowest node PE */
    
    /*interface flag to prevent adaption */
    int interface;
    
    FLUX_ELEM *flux_elem;      /* supermodel flux coupling information */
    int flux_elem_tot;
    
} SELEM_2D;

/*********************************************************/
/* struct methods -------------------------------------- */

void selem2d_get_triangle_local_shape(double xhat, double yhat, double zhat, double *lshape);
void selem2d_get_quadrilateral_local_shape(double xhat, double yhat, double zhat, double *lshape);
void selem2d_alloc(SELEM_2D *elem2d, int nnodes_on_elem);
void selem2d_load(SELEM_2D *elem2d, int gid, int lid, int elem_nnodes, int *local_node_ids, int bflag, SVECT *nds, int mat);
void selem2d_free(SELEM_2D *elem2d);
void selem2d_alloc_array(SELEM_2D **elem2d, int nelems2d);
void selem2d_free_array(SELEM_2D *elem2d, int nelems2d);
void selem2d_init(SELEM_2D *elem2d);
void selem2d_init_array(SELEM_2D *elem2d, int nelems2d);
void selem2d_init_alloc_array(SELEM_2D **elem2d, int nelems2d);
void selem2d_copy(SELEM_2D *to, SELEM_2D from);
void selem2d_printScreen(SELEM_2D *elem2d);
SVECT selem2d_get_triangle_centroid(SNODE nd1, SNODE nd2, SNODE nd3);
double selem2d_get_triangle_orientation(SNODE nd1, SNODE nd2, SNODE nd3, SVECT ref_vec);
void selem2d_get_triangle_linear_djac_nrml_gradPhi(SELEM_2D *elem2d, SNODE *nd_SNODE, SVECT *nd_SVECT);
SVECT selem2d_get_elem2d_normals(SVECT *nd);
double selem2d_get_quadrilateral_linear_djac2d(double xhat, double yhat, SVECT *nd);
double selem2d_get_quadrilateral_linear_djac_gradPhi(double xhat, double yhat, SVECT *nd, SVECT *grad_shp);
void selem2d_get_triangle_local_shape_quad(double xhat, double yhat, double zhat, double *lshape);
void selem2d_get_quadrilateral_local_shape_quad(double xhat, double yhat, double zhat, double *lshape);
void selem2d_get_triangle_local_shape_gradients(SVECT *lgrad_shp);
void selem2d_get_quadrilateral_local_shape_gradients(double xhat, double yhat, double zhat, SVECT *lgrad_shp);
double selem2d_get_elem2d_area2d(SELEM_2D *elem2d, SVECT *nds);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
