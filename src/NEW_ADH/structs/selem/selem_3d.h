#ifndef H_SELEM_3D_
#define H_SELEM_3D_

// dependencies: SVECT

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

typedef struct {

    int id;                   /* current 3d element number */
    int gid;                  /* global id should == id_orig for non-mpi */
    int id_orig;              /* original 2d element number (every 2d element should have a value) */
    int nnodes;               /* total # of nodes on element */
    int nnodes_quad;          /* the total number of quadratic nodes on element */
    double djac;              /* the jacobian of the element */
    int string;               /* string number associated with this element */
    int mat;                  /* 3d element physics material type */
    int nvars;          /* the total number of independent variables on this element */
    int *vars;          /* an array of independent variables on this element */
    int icol;                 /* the column this tet belongs to */
    int elem2d_sur;           /* the 2d surface element of the column this tet belongs to */
    int elem2d_bed;           /* the 2d bed element of the column this tet belongs to */
    double error;             /* max elemental error in continuity equations for the current ts */
    int flux_ptr;             /* 2d element pointer for flux line output */
    SVECT *grad_shp;          /* the gradients of the shape functions */
    int *nodes;               /* the nodes in the element */
    int *levels;              /* the node levels in the element */
    int nedges;               /* the number of edges on the element */
    double volume;
    int resident_pe;
    
    /*interface flag to prevent adaption */
    int interface;

} SELEM_3D;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// struct methods

void selem3d_alloc(SELEM_3D *elem3d, int nnodes_on_elem);
void selem3d_load(SELEM_3D *elem3d, int gid, int lid, int elem_nnodes, int *local_node_ids, int bflag, SVECT *nds, int mat);
void selem3d_free(SELEM_3D *elem3d);
void selem3d_alloc_array(SELEM_3D **elem3d, int nelems3d);
void selem3d_free_array(SELEM_3D *elem3d, int nelems3d);
void selem3d_init(SELEM_3D *elem3d);
void selem3d_init_array(SELEM_3D *elem3d, int nelems3d);
void selem3d_init_alloc_array(SELEM_3D **elem3d, int nelems3d);
void selem3d_copy(SELEM_3D *to, SELEM_3D from);
void selem3d_printScreen(SELEM_3D *elem3d);
void selem3d_get_tet_local_shape(double xhat, double yhat, double zhat, double *lshape);
void selem3d_get_triprism_local_shape(double xhat, double yhat, double zhat, double *lshape);
void selem3d_get_tet_local_shape_quad(double xhat, double yhat, double zhat, double *lshape_quad);
void selem3d_get_triprism_local_shape_quad(double xhat, double yhat, double zhat, double *lshape_quad);
void selem3d_get_tet_local_shape_gradients(SVECT *lgrad_shp);
void selem3d_get_triprism_local_shape_gradients(double xhat, double yhat, double zhat, SVECT *lgrad_shp);
double selem3d_get_triprism_djac(double xhat, double yhat, double zhat, SVECT *nd);
double selem3d_get_triprism_linear_djac_gradPhi(double xhat, double yhat, double zhat, SVECT *nd, SVECT *grad_shp);
void selem3d_get_tet_linear_djac_gradPhi(SELEM_3D *elem3d, SNODE *nd_SNODE, SVECT *nd_SVECT);
double selem3d_get_tet_linear_djac(SNODE *nd_SNODE, SVECT *nd_SVECT);
double selem3d_get_tet_linear_djac_gradPhi2(SNODE *nd_SNODE, SVECT *nd_SVECT, SVECT *grad_shp);
double selem3d_get_elem3d_volume(SVECT *node, int nnodes);
double selem3d_get_triprism_volume(SVECT *node);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
