/* standard header files */

#define NOT_QUITE_SMALL FLT_EPSILON /* the smallest single precision number */
#define UNSET_INT -3    /* an integer variable that has not been set */

#ifdef _ADH_HDF5
#include "adh_hdf5.h"
#endif

#include <stdio.h>
//#include <unistd.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>
#include <stddef.h>
#include <limits.h>
#include <ctype.h>
#include "assert.h"
#include "stdbool.h"
#include "svect.h"          // dependencies :: smpi.h with _MESSG active
#include "selem_2d.h"       // dependencies :: svect.h, svect2d.h
#include "selem_3d.h"       // dependencies :: svect
#include "sgrid.h"
#include "sarray.h"
#include "sparticle.h"
#include "soyster.h"
#include "spolygon.h"
#include "utilities.h"
#include "smodel.h"
#include "sensory_point.h"
#include "sensor_species.h"

#ifdef _MPI
#include <mpi.h>
#endif


// particle types
#define GENERAL 0
#define OYSTER 1

#define ON 1      /* on */
#define OFF 0     /* off */
#define YES 1     /* yes */
#define NO -3     /* no */
#define TRUE 1    /* true */
#define FALSE 0   /* false */

// GLOBAL VARIABLES
bool RUN_VERBOSE;

/* the nodes on an edge of a 2d element */
static const int nd_on_2dedge[3][2] = { {0, 1}, {0, 2}, {1, 2} };

/* the nodes on an edge of a 3d element */
static const int nd_on_3dedge[6][2] = { {0, 1}, {0, 2}, {0, 3}, {1, 2}, {3, 2}, {1, 3} };

/* the nodes on a face */
static const int nd_on_fc[4][3] = { {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1} };

/* face opposite a node - redundant with this numbering but maintained for readability */
static const int fc_opp_nd[4] = { 0, 1, 2, 3 };

/* node opposite a face - also redundant */
static const int nd_opp_fc[4] = { 0, 1, 2, 3 };

/* edge opposite a node (2d element) - redundant with this numbering but maintained for readability */
static const int edge_opp_nd[3] = { 0, 1, 2 };

/* node opposite an edge (2d element) - also redundant */
static const int nd_opp_edge[3] = { 0, 1, 2 };

/* finds an edge given the local */
static const int find_edge[4][4] = { {UNSET_INT, 0, 1, 2}, {0, UNSET_INT, 3, 5}, {1, 3, UNSET_INT, 4}, {2, 5, 4, UNSET_INT} };

#define PI 3.14159265359
#define rad2deg 57.29577951308230
#define deg2rad 0.01745329251994330
#define DEFFLEQEPSILON 0.001
#define FLOAT_EQE(x,v,e)((((v)-(e))<(x))&&((x)<((v)+(e))))

#define dot(u,v)   ((u).x * (v).x + (u).y * (v).y + (u).z * (v).z)
#define norm(v)    sqrt(dot(v,v))  // norm = length of  vector
#define d(u,v)     norm(u-v)        // distance = norm of difference
#define abs(x)     ((x) >= 0 ? (x) : -(x))   //  absolute value

int pt_test_searchEngine();
int pt_verification();
int pt_run(SGRID *grid, int np, SVECT *p, int *p_elemID, double t0, double tf, double dt, FILE *fp_vel, FILE *fp_out, FILE *fp_error, char *vel_file, SVECT (*get_analytic_p)(SVECT p0, double t), HDF5 *hdf5, double *abs_error_tmax);

void project_velocity_2D(SVECT v1, SVECT v2, SVECT nd1, SVECT nd2, SVECT *v1_proj, SVECT *v2_proj);
void project_point_velocity_2D(SVECT *v, SVECT nd1, SVECT nd2);
void project_point_velocity_3D(SVECT *v, SVECT nd1, SVECT nd2, SVECT nd3);
void project_velocity_3D(SVECT v1, SVECT v2, SVECT v3, SVECT nd1, SVECT nd2, SVECT nd3, SVECT *v1_proj, SVECT *v2_proj, SVECT *v3_proj);

double interpolate1D(double x, double xL, double xR, double yL, double yR);
int compute_interpolate_weights_2D_triangle(SGRID *grid, int ie, double x, double y, double *weights, int *lnd, int *iedge);
int compute_interpolate_weights_3D_tetrahedron(SGRID *grid, int ie, double x, double y, double z, double *weights, int *lnd, int *iedge, int *iface);

SVECT euler_step(SMODEL *mod, int ip, double t, double dt, SVECT p_vel, int p_vel_elem, int *new_elem, char *string);
void euler_update(SMODEL *mod, int ip, double t1, double *dpl1);
void rk2_update(SMODEL *mod, int ip, double t1, double t2, double *dpl1, double *dp2);
void time_update(SMODEL *mod, int rk_order);
double read_adh_velocity_frame(FILE *fp, SGRID *grid, SVECT *vel, int time_units, bool normalize);
double read_adh_wq_frame(FILE *fp, SGRID *grid, double *oxy, double *sal, double *sun, int time_units);
double read_adh_dpl_frame(FILE *fp, SGRID *grid, double *dpl, int time_units);;

bool pointOnLine(SVECT p1, SVECT p2, SVECT p);
bool pointOnLine2D(SVECT p1, SVECT p2, SVECT p);
int intersect2D_SegSeg(SVECT p1,SVECT p2,SVECT p3,SVECT p4,SVECT *pi);
int intersect3D_SegSeg(SVECT p1,SVECT p2,SVECT p3,SVECT p4,SVECT *pi);
int intersect2D_LineLine(SVECT line1_p1, SVECT line1_p2, SVECT line2_p1, SVECT line2_p2, double *x, double *y);
int intersect3D_SegTriangle(SVECT tri_node0, SVECT tri_node1, SVECT tri_node2, SVECT p0, SVECT p1, SVECT *I);
SVECT shakeAndBake(SGRID *grid, int ielem, int lnode);

SVECT get_errors_test(SVECT p, SVECT p0, double t, SVECT (*get_analytic_p)(SVECT p0, double t), SVECT *abs_error, SVECT *rel_error);
void test_geometry();

int elementSearch2D(SGRID *grid, SVECT p1, SVECT p2, int *ielem, bool cycle, SVECT *pi);
int elementSearch3D(SGRID *grid, SVECT p1, SVECT p2, int *ielem, bool cycle, SVECT *pi);
int elementSearch(SGRID *grid, SVECT p1, SVECT p2, int *ielem, SVECT *pi, double *weights, int *call_ID);
int checkElementAfterDisplacement(SGRID *grid, SVECT *p, int *ielem);

// velocity field creation functions
SVECT get_vel_2d_uniform_oscillating(SVECT p, double t);
SVECT get_vel_2d_uniform_oscillating_off_grid(SVECT p, double t);
SVECT get_vel_2d_rotation(SVECT p, double t);
SVECT get_vel_3d_rotation(SVECT p, double t);
SVECT get_vel_3d_off_grid_node(SVECT p, double t);
SVECT get_vel_3d_off_grid_edge(SVECT p, double t);
SVECT get_vel_3d_off_grid_face(SVECT p, double t);
SVECT get_vel_2d_closed_boundary(SVECT p, double t);
SVECT get_vel_3d_oyster_fall(SVECT p, double t);
SVECT get_vel_3d_oyster_swim(SVECT p, double t);
SVECT get_vel_3d_oyster_tvel(SVECT p, double t);
SVECT get_vel_3d_surface_dpl(SVECT p, double t);

// analytic position functions
SVECT get_analytic_p_2d_uniform_oscillating(SVECT p0, double t);
SVECT get_analytic_p_2d_uniform_oscillating_off_grid(SVECT p0, double t);
SVECT get_analytic_p_2d_rotation(SVECT p0, double t);
SVECT get_analytic_p_3d_rotation(SVECT p0, double t);
SVECT get_analytic_p_3d_off_grid_node(SVECT p0, double t);
SVECT get_analytic_p_3d_off_grid_edge(SVECT p0, double t);
SVECT get_analytic_p_3d_off_grid_face(SVECT p0, double t);
SVECT get_analytic_p_2d_closed_boundary(SVECT p0, double t);
SVECT get_analytic_p_3d_oyster_fall(SVECT p0, double t);
SVECT get_analytic_p_3d_oyster_swim(SVECT p0, double t);
SVECT get_analytic_p_3d_oyster_tvel(SVECT p0, double t);
SVECT get_wq_3d_oyster_oxygen(SVECT p, double t);
SVECT get_wq_3d_oyster_salinity(SVECT p, double t);
SVECT get_wq_3d_oyster_sunlight(SVECT p, double t);
SVECT get_analytic_p_3d_surface_dpl(SVECT p0, double t);

double find_random_in_range(double min, double max);
double rand_val(int seed);
double find_random_from_normal(double mean, double std_dev);
double RandomFromvonMises(double mu, double kappa, int Seed);


#ifdef _ADH_HDF5
void xdmf_init(SGRID *, int, int, const char *);
void ps_print_geo_xdmf(SGRID *);
void ps_print_xdmf(SGRID *grid, SVECT *vel, double *dpl, double time);
void xdmf_finalize(HDF5*, const char *, int, int, int);
void xdmf_init_particles(HDF5 *hdf5, int np, SPARTICLE *p, int npes, int myid, const char *proj_name);
void ps_print_particle_t0_xdmf(HDF5 *hdf5, int np, SPARTICLE *p);
void ps_print_particle_xdmf(HDF5 *hdf5, SVECT *dpl, int np, double time);
void xdmf_particles_print_ts(HDF5 *hdf5,char *data_name,SVECT *p,int np,double time,int ntime,int data_type);
void xdmf_particle_finalize(HDF5 *hdf5, const char *proj_name, int npes, int myid, int flag);
#endif
