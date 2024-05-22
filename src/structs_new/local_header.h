// type of grid
#define UNSTRUCTURED 2
#define COLUMNAR 3

#define NEWIO_HASHSIZE 10069

#define BODY -1
#define SURFACE 0
#define BED 1
#define SIDEWALL 2
#define COLUMN 3

#define NDONSEG 2       /* nodes on line segment */
#define NDONTRI 3       /* nodes on triangle */
#define NDONTRIQUAD 6
#define NDONQUAD 4       /* nodes on quadrilateral */
#define NDONQUADQUAD 8
#define NDONTET 4         /* nodes on tetrahedron */
#define NDONPRISM 6       /* nodes on prism */
#define NDONPYRAMID 5     /* nodes on a pyramid */
#define NDONTETQUAD 10    /* nodes on a tetrahedron with midpoints for quadratic baasis */
#define NDONPRISMQUAD 15  /* nodes on a triangular prism with midpoints for quadratic baasis */

#define MAX_NNODES_ON_ELEM1D 2  /* max number of nodes on all possible 2d elements (line here) */
#define MAX_NNODES_ON_ELEM2D 4  /* max number of nodes on all possible 2d elements (quads here) */
#define MAX_NNODES_ON_ELEM3D 6  /* max number of nodes on all possible 3d elements (triprisms here) */
#define MAX_NNODES_ON_ELEM3D_QUAD 15  /* max number of quadratic nodes on all possible 3d elements (triprisms here) */

#define MAX_QUAD_ORDER 4 /* the max integrand order for quadrature */
#define MAX_QUAD_POINTS 18 /* max number of quadrature points for an element */

#define DIST_2D(a,b) sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y))
#define SMALL6 1.0e-6

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>
#include <stddef.h>
#include <limits.h>
#include <ctype.h>

#include "debug.h"
#include "header_tl_alloc.h"

#include "assert.h"

#include "stdbool.h"
#include "svect2d.h"        // dependencies :: smpi.h with _MESSG active
#include "svect.h"          // dependencies :: smpi.h with _MESSG active
#include "squad.h"
#include "slist_items.h"
#include "smpi.h"

#include "snode.h"
#include "selem_1d.h"       // dependencies :: svect2d.h
#include "selem_2d.h"       // dependencies :: svect.h, svect2d.h
#include "selem_3d.h"       // dependencies :: svect

#include "elem_physics.h"
#include "sgrid.h"          // dependencies :: selem_1d.h, selem_2d.h, selem_3d.h,snode.h, slist_items.h, smpi.h
#include "super_model.h"
#include "design_model.h"

 #include "sarray.h"

#define UNSET_INT -3
#define UNSET_FLT -9999999.9
#define one_3  0.333333333333333333333333
#define ON 1            /* on */
#define OFF 0           /* off */
#define YES 1           /* yes */
#define NO -3           /* no */
#define NORMAL -1
#define NOT_QUITE_SMALL FLT_EPSILON
#define SMALL DBL_EPSILON

#define one_2  0.5e+0
#define one_3  0.333333333333333333333333
#define one_4  2.5e-1
#define one_5  0.200000000000000000000000
#define one_6  0.166666666666666666666667
#define one_8  0.125000000000000000000000
#define one_9  0.111111111111111111111111
#define one_12 0.083333333333333333333333
#define one_18 0.055555555555555555555556
#define one_20 0.050000000000000000000000
#define one_24 0.041666666666666666666667
#define one_30 0.033333333333333333333333
#define one_36 0.027777777777777800000000
#define one_45 0.022222222222222222222222
#define one_48 0.020833333333333332000000
#define one_60 0.016666666666666666666667
#define one_72 0.013888888888888888888888
#define one_120 0.00833333333333333000000
#define one_144 0.00694444444444444444444
#define one_180 0.00555555555555555555556
#define one_360 0.00277777777777778000000
#define one_720 0.00138888888888889000000
#define one_1440 0.0006944444444444440000


// ADH CANONICAL EDGE TO LOCAL ELEMENT NODE ID MAPS
#define GET_ND_ON_TRI_EDGE(nd_on_edge) \
    nd_on_edge[0][0] = 0; \
    nd_on_edge[0][1] = 1; \
    nd_on_edge[1][0] = 1; \
    nd_on_edge[1][1] = 2; \
    nd_on_edge[2][0] = 2; \
    nd_on_edge[2][1] = 0;

#define GET_ND_ON_QUAD_EDGE(nd_on_edge) \
    nd_on_edge[0][0] = 0; \
    nd_on_edge[0][1] = 1; \
    nd_on_edge[1][0] = 1; \
    nd_on_edge[1][1] = 2; \
    nd_on_edge[2][0] = 2; \
    nd_on_edge[2][1] = 3; \
    nd_on_edge[3][0] = 3; \
    nd_on_edge[3][1] = 0;

#define GET_ND_ON_PRISM_EDGE(nd_on_edge) \
    nd_on_edge[0][0] = 0; \
    nd_on_edge[0][1] = 1; \
    nd_on_edge[1][0] = 1; \
    nd_on_edge[1][1] = 2; \
    nd_on_edge[2][0] = 2; \
    nd_on_edge[2][1] = 0; \
    nd_on_edge[6][0] = 3; \
    nd_on_edge[6][1] = 4; \
    nd_on_edge[7][0] = 4; \
    nd_on_edge[7][1] = 5; \
    nd_on_edge[8][0] = 5; \
    nd_on_edge[8][1] = 3; \
    nd_on_edge[3][0] = 0; \
    nd_on_edge[3][1] = 3; \
    nd_on_edge[4][0] = 1; \
    nd_on_edge[4][1] = 4; \
    nd_on_edge[5][0] = 2; \
    nd_on_edge[5][1] = 5;

#define GET_ND_ON_Tet_EDGE(nd_on_edge) \
    nd_on_edge[0][0] = 0; \
    nd_on_edge[0][1] = 1; \
    nd_on_edge[1][0] = 0; \
    nd_on_edge[1][1] = 2; \
    nd_on_edge[2][0] = 0; \
    nd_on_edge[2][1] = 3; \
    nd_on_edge[3][0] = 1; \
    nd_on_edge[3][1] = 2; \
    nd_on_edge[4][0] = 1; \
    nd_on_edge[4][1] = 3; \
    nd_on_edge[5][0] = 2; \
    nd_on_edge[5][1] = 3;

#define VT_3D_VCOPY(vect, target) \
        target.x = vect.x; \
        target.y = vect.y; \
        target.z = vect.z;

int messg_comm_rank(MPI_Comm);
int messg_comm_size(MPI_Comm);
void messg_barrier(MPI_Comm);
double messg_dmax(double, MPI_Comm);
double messg_dmin(double, MPI_Comm);
double messg_dsum(double, MPI_Comm);
int    messg_imax(int, MPI_Comm);
int messg_isum(int, MPI_Comm);
int    messg_imin(int, MPI_Comm);
void   messg_max_norm_loc(double *, int *, SVECT *, MPI_Comm, int);
void messg_err(int);
void messg_buffer_alloc(int, size_t, MESSG_BUFFER *);
void messg_buffer_free( MESSG_BUFFER *);
void messg_buffer_init( MESSG_BUFFER *, int);
void messg_pack_alloc(MESSG_BUFFER *, MESSG_BUFFER *, MPI_Comm);
void messg_asend(MESSG_BUFFER *, int, SMPI *);
int messg_incoming(int *, SMPI *);
int messg_precv(MESSG_BUFFER *, int, SMPI *);
void messg_arecv(MESSG_BUFFER *, int, SMPI *);
void messg_unpack(MESSG_BUFFER *, MESSG_BUFFER *, MPI_Comm);
void messg_wait(SMPI *);
void messg_pack(MESSG_BUFFER *, MESSG_BUFFER *, MPI_Comm);
void messg_gather_int(int root,
                      int *my_x, /* my part of the array */
                      int *x,    /* the resulting array */
                      int size  /* the size of x */
    );
void messg_gather_dbl(int root,
                      double *my_x, /* my part of the array */
                      double *x,    /* the resulting array */
                      int size  /* the size of x */
    );
int sort_key(const void *, const void *);
