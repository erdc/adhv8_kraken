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

//int SQUAD_alloc_init(SQUAD **, SQUAD **);
int SQUAD_triangle_alloc_init(SQUAD **quad_tri);
int SQUAD_rectangle_alloc_init(SQUAD **quad_rect);
int SQUAD_tetrahedron_alloc_init(SQUAD **quad_tet);
int SQUAD_triprism_alloc_init(SQUAD **quad_prism);
void SQUAD_free(SQUAD *quad_tri, int);
void SQUAD_PT_alloc_array(SQUAD_PT **qp, int nqp, int ndof);
void SQUAD_PT_free_array(SQUAD_PT *qp, int nqp, int ndof);
void SQUAD_PT_init(SQUAD_PT *qp, int ndof);
void SQUAD_PT_init_array(SQUAD_PT *qp, int nqp, int ndof);
void SQUAD_printScreen(SQUAD *qp, int flag, int ndof, char * geometry);
SVECT SQUAD_get_svect(SQUAD_PT *qp, SVECT *v, int nnodes);
SVECT2D SQUAD_get_svect2d(SQUAD_PT *qp, SVECT2D *v, int nnodes);
double SQUAD_get_function(SQUAD_PT *qp, double *f, int nnodes);
double SQUAD_get_quadFunction(SQUAD_PT *qp, double *f, int nnodes);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
