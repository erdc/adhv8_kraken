#ifndef H_SQUAD_
#define H_SQUAD_

// Dependencies :: SVECT

// quad[order].nqp
// quad[order].pt[ipoint].xhat

/***********************************************************/
/***********************************************************/
/***********************************************************/
// stores values for a quadrature point

typedef struct {

    double xhat; // quadrature x-coord on parent element
    double yhat; // quadrature y-coord on parent element
    double zhat; // quadrature z-coord on parent element
    double w;    // quadrature weight
    
    double *lshape;         // local shape function value
    double *lshape_quad;    // local quadratic shape function values
    SVECT *lgrad_shp;       // local shape function gradient
    SVECT *grad_shp;        // cartesian shape function gradient

    double djac;            // jacobian determinate
    double djac2d;
    double djac2d_3d;       // used for 2d quadrature
    double djac3d_fixed;    // used for 2d quadrature
    
    // shallow water storage
    SVECT vel;
    SVECT vel_old;
    SVECT vel_12dt;
    SVECT vel_32dt;
    SVECT vel_rel;
    SVECT grad_u;
    SVECT grad_v;
    SVECT grad_w;
    SVECT grad_p;
    SVECT grad_density;
    double prs;
    double density;
    double relvel;
    double djac_old;
    double djac_12dt;
    double djac_32dt;
    SVECT grad_shp_quad[NDONPRISM];

} SQUAD_PT;

/***********************************************************/
/***********************************************************/
// stores values for a quadrature point
typedef struct {
    int n;              // # of quadrature points for given order
    SQUAD_PT *pt;       // list of quadrature points for given order
} SQUAD;

/*********************************************************/
/* struct methods -------------------------------------- */

//int squad_alloc_init(SQUAD **, SQUAD **);
int squad_segment_alloc_init(SQUAD **quad_seg);
int squad_triangle_alloc_init(SQUAD **quad_tri);
int squad_rectangle_alloc_init(SQUAD **quad_rect);
int squad_tetrahedron_alloc_init(SQUAD **quad_tet);
int squad_triprism_alloc_init(SQUAD **quad_prism);
void squad_free(SQUAD *quad_tri, int);
void squad_pt_alloc_array(SQUAD_PT **qp, int nqp, int ndof);
void squad_pt_free_array(SQUAD_PT *qp, int nqp, int ndof);
void squad_pt_init(SQUAD_PT *qp, int ndof);
void squad_pt_init_array(SQUAD_PT *qp, int nqp, int ndof);
void squad_printScreen(SQUAD *qp, int flag, int ndof, char * geometry);
SVECT squad_get_svect(SQUAD_PT *qp, SVECT *v, int nnodes);
SVECT2D squad_get_svect2d(SQUAD_PT *qp, SVECT2D *v, int nnodes);
double squad_get_function(SQUAD_PT *qp, double *f, int nnodes);
double squad_get_quadFunction(SQUAD_PT *qp, double *f, int nnodes);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
