#ifndef H_FE_PROTOTYPES_
#define H_FE_PROTOTYPES_

#include "./sw2/fe_sw2_prototypes.h"
#include "./sw3/fe_sw3_prototypes.h"
#include "./ns3/fe_ns3_prototypes.h"
#include "./sw_hybrid/fe_sw_hybrid_prototypes.h"
#include "./transport/fe_transport_prototypes.h"
#include "./sediment/fe_sediment_prototypes.h"
#include "./diffusive_wave/fe_diffusive_prototypes.h"

#ifdef _ADH_GROUNDWATER
#include "./gw/fe_gw_prototypes.h"
#ifdef _DWGW_COUPLING
#include "./dwgw_hybrid/fe_dwgw_hybrid_prototypes.h"
#endif
#endif
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int fe_main(SSUPER_MODEL *sm);


int fe_newton(SSUPER_MODEL *sm, int isuperModel, int my_nnodes, int nnodes, int macro_nnodes,
#ifdef _MESSG
              SMPI *supersmpi,                            /* either sm->supersmpi or sm->supersmpi_wvel */
#endif
        void (*init_fnctn) (SSUPER_MODEL *, int), 
        void (*update_fnctn) (SSUPER_MODEL *, int), 
        void (*residual_fnctn) (SSUPER_MODEL *, int), 
        void (*load_fnctn) (SSUPER_MODEL *, int), 
        void  (*inc_fnctn) (SSUPER_MODEL *, int));

double fe_sw2_wet_dry_factor(SVECT *x, double *h, double djac);
double fe_sw2_wet_dry_wrapper(DOF_3 *elem_rhs, SVECT *x, double *h, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, int REDISTRIBUTE_FLAG, int DEBUG, double *, void (*fe_sw_func) ());
void calculate_supermodel_fluxes(int isuper_model, int isubmodel, SMODEL *mod);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// fe_assign_elem_mat_db
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void fe_assign_mom_db_dof3(int inode, int nnodes, DOF_3 *elem_mat_u, DOF_3 *elem_mat_v, DOF_3 *elem_mat_h);
void fe_assign_dacont_db_dof3(int inode, int nnodes, DOF_3 *elem_mat_u, DOF_3 *elem_mat_v, DOF_3 *elem_mat_h);
void fe_assign_mom_db_dof4(int inode, int nnodes, DOF_4 *elem_mat_u, DOF_4 *elem_mat_v, DOF_4 *elem_mat_w, DOF_4 *elem_mat_p);
void fe_assign_prs_db_dof4(int inode, int nnodes, DOF_4 *elem_mat_u, DOF_4 *elem_mat_v, DOF_4 *elem_mat_w, DOF_4 *elem_mat_p);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// fe_assemble_matrix
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void fe_global_matrix_assemble_scalar(SGRID *grid, int nsys, int nnodes_on_elem, int *GnodeID, int *fmap, double *elem_mat//, double *diagonal, SPARSE_VECT *matrix
#ifdef _PETSC
        , Mat A, int Istart, int Iend, const int *ownership_range
#else
        , double *diagonal, SPARSE_VECT *matrix
#endif
        );
void fe_global_matrix_assemble_sw2(SGRID *grid, int nsys, int nnodes_on_elem, int *GnodeID, int *fmap, DOF_3 *elem_mat_u, DOF_3 *elem_mat_v, DOF_3 *elem_mat_h,
#ifdef _PETSC
        Mat A, int Istart, int Iend, const int *ownership_range
#else
        double *diagonal, SPARSE_VECT *matrix
#endif
        );
void fe_global_matrix_assemble_sw3(int dim, SGRID *grid, int nsys, int nnodes_on_elem, int ie, int *fmap, DOF_3 *elem_mat_u, DOF_3 *elem_mat_v, DOF_3 *elem_mat_dpl,
#ifdef _PETSC
        Mat A
#else
        double *diagonal, SPARSE_VECT *matrix
#endif
        );
void fe_global_matrix_assemble_ns3(int dim, SGRID *grid, int nsys, int nodes_on_elem, int ie, int *fmap, DOF_4 *elem_mat_u, DOF_4 *elem_mat_v, DOF_4 *elem_mat_w, DOF_4 *elem_mat_p, double *diagonal, SPARSE_VECT *matrix);
#ifdef _PETSC
void fe_global_matrix_assemble_wvel(SGRID *grid, int nsys, int nodes_on_elem, int *GnodeID, int *fmap, double *elem_mat, Mat A, int Istart, int Iend, const int *ownership_range);
#endif

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// fe_elem1d_integrations
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double integrate_line_h_f         (double djac, double c, double *h, double *f);
void  integrate_line_phi_h_v_dot_n(double djac, double c, double *h, SVECT2D *v, SVECT2D nrml, double *integral);
void  integrate_line_phi_h_h      (double djac, double c, double h1, double h2, double *integral);
void  integrate_line_phi_h_h_g    (double djac, double c, double h1, double h2, double g1, double g2, double *integral);
void  integrate_line_phi          (double djac, double c, double *integral);
void  integrate_line_phi_f        (double djac, double c, double *f, double *integral);
void  integrate_line_phi_f_g      (double djac, double c, double *f, double *g, double *integral);
void  integrate_line_phi_f_g_h    (double djac, double c, double f1, double f2, double g1, double g2, double h1, double h2, double *integral);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// fe_elem2d_integrations
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void extractNodesQuad(SVECT *v, double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4);
void extractVect2DQuad(SVECT2D *v, double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4);
void extractFunctionQuad(double *f, double *f1, double *f2, double *f3, double *f4);
double  integrate_triangle_area      (double djac, double c);
double  integrate_triangle_f         (double djac, double c, double *f);
double  integrate_triangle_f_f       (double djac, double c, double *f);
void    integrate_triangle_phi       (double djac, double c, double *integral);
void    integrate_triangle_phi_f     (double djac, double c, double *f, double *integral);
void    integrate_triangle_phi_f_lump(double djac, double c, double *f, double *integral);
void    integrate_triangle_phi_f_g   (double djac, double c, double *f, double *g, double *integral);
void    integrate_triangle_phi_fquad (double djac, double c, double *f, double *integral);
void    integrate_triangle_phi_h_df  (double djac, double c, double *h, SVECT2D df, double *integral_x, double *integral_y);
void    integrate_triangle_phi_h_g_df(double djac, double c, double *h, double *g, SVECT2D df, double *integral_x, double *integral_y) ;
void    integrate_triangle_dphi_f_f         (SVECT2D *grad_shp, double djac, double c, double *f, double *integral_x, double *integral_y);
void    integrate_triangle_dphi_f_g         (SVECT2D *grad_shp, double djac, double c, double *f, double *g, double *integral_x, double *integral_y);
void    integrate_triangle_dphi_f_g_h       (SVECT2D *grad_shp, double djac, double c, double *f, double *g, double *h, double *integral_x, double *integral_y);
void    integrate_triangle_dphi_f_f_h       (SVECT2D *grad_shp, double djac, double c, double *f, double *h, double *integral_x, double *integral_y);
void    integrate_triangle_gradPhi_dot_vbar (SVECT2D *grad_shp, double djac, double c, SVECT2D vbar, double *integral);
void    integrate_triangle_gradPhi_dot_v    (SVECT2D *grad_shp, double djac, double c, SVECT2D *v, double *integral);
double  integrate_triangle_grad_dot_f_vbar  (SVECT2D *grad_shp, double djac, double c, double *f, SVECT2D vbar);
void    integrate_triangle_gradPhi_dot_f_v  (SVECT2D *grad_shp, double djac, double c, double *f, SVECT2D *v, double *integral);
void    integrate_triangle_gradPhi_dot_f_g_v(SVECT2D *grad_shp, double djac, double c, double *f, double *g, SVECT2D *v, double *integral);
double   integrate_quadrilateral_area             (SVECT *nd, double c);
double   integrate_quadrilateral_f                (SVECT *nd, double c, double *f);
SVECT2D  integrate_quadrilateral_gradF            (SVECT *nd, double c, double *f);
double   integrate_quadrilateral_f_f              (SVECT *nd, double c, double *f);
void     integrate_quadrilateral_phi              (SVECT *nd, double c, double *integral);
void     integrate_quadrilateral_phi_f            (SVECT *nd, double c, double *f, double *integral);
void     integrate_quadrilateral_phi_f_lump       (SVECT *nd, double c, double *f, double *integral);
void     integrate_quadrilateral_phi_f_g          (SVECT *nd, double c, double *f, double *g, double *integral) ;
void     integrate_quadrilateral_dphi_f_f         (SVECT *nd, double c, double *f, double *integral_x, double *integral_y);
void     integrate_quadrilateral_gradPhi_dot_vbar (SVECT *nd, double c, SVECT2D vbar, double *integral);
void     integrate_quadrilateral_gradPhi_dot_v    (SVECT *nd, double c, SVECT2D *v, double *integral);
double   integrate_quadrilateral_grad_dot_f_vbar  (SVECT *nd, double c, double *f, SVECT2D vbar);
void     integrate_quadrilateral_gradPhi_dot_f_g_v(SVECT *nd, double c, double *f, double *g, SVECT2D *v, double *integral);
void     integrate_quadrilateral_gradPhi_dot_f_v  (SVECT *nd, double c, double *f, SVECT2D *v, double *integral);
SVECT2D  integrate_quadrilateral_dv               (SVECT *nd, double c, SVECT2D *v);
void     integrate_quadrilateral_dphi_f_h_h       (SVECT *nd, double c, double *f, double *h, double *integral_x, double *integral_y);
void     integrate_quadrilateral_phi_h_df         (SVECT *nd, double c, double *df, double *h, double *integral_x, double *integral_y);
void     integrate_quadrilateral_phi_h_g_df       (SVECT *nd, double c, double *df, double *h, double *g, double *integral_x, double *integral_y);
double   integrate_quadrilateral_grad_dot_h_v     (SVECT *nd, double c, SVECT2D *v, double *h);
double   integrate_quadrilateral_grad_dot_h_vbar  (SVECT *nd, double c, SVECT2D vbar, double *h);
double   integrate_quadrilateral_grad_dot_h_g_v   (SVECT *nd, double c, SVECT2D *v, double *h, double *g);
double   integrate_quadrilateral_grad_dot_h_g_vbar(SVECT *nd, double c, SVECT2D vbar, double *h, double *g);
void     integrate_quadrilateral_gradPhi_dot_Df   (SVECT *nd, SQUAD *quad, double c, STENSOR2D T, double *f, double *integral);
void     integrate_quadrilateral_gradPhi_dot_Dv   (SVECT *nd, SQUAD *quad, double c, STENSOR2DAI T, SVECT2D *v, double *integral_x, double *integral_y);
double  integrate_quadZ_f      (SVECT *nd, double c, double *f);
void    integrate_quadZ_phi    (SVECT *nd, double c, double *integral);
void    integrate_quadZ_phi_f  (SVECT *nd, double c, double *f, double *integral);
void    integrate_quadZ_phi_f_g(SVECT *nd, double c, double *f, double *g, double *integral);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// fe_elem3d_integrations
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double  integrate_tetrahedron_c      (double djac, double c);
void    integrate_tetrahedron_fi     (double djac, double c,  double *f, double *integral);
double  integrate_tetrahedron_f      (double djac, double c,  double *f);
double  integrate_tetrahedron_f_g    (double djac, double c,  double *f, double *g);
void    integrate_tetrahedron_phi_f  (double djac, double c,  double *f, double *integral);
void    integrate_tetrahedron_phi_f_g(double djac, double c,  double *f, double *g, double *integral);
void    integrate_tetrahedron_ci_f_g (double djac, double *c, double *f, double *g, double *integral);
void    integrate_tetrahedron_gradPhi_dot_v  (SVECT *grad_phi, double djac, double c, SVECT *vector, double *integral);
void    integrate_tetrahedron_gradPhi_dot_vcon(SVECT *grad_phi, double djac, double c, SVECT vector, double *integral);
void    integrate_tetrahedron_gradPhi_dot_f_v(SVECT *grad_phi, double djac, double c, double *f, SVECT *vector, double *integral);
double  integrate_triPrism_area              (SVECT *nd, double c);
double  integrate_triPrism_f                 (SVECT *nd, double c, double *f);
void    integrate_triPrism_phi               (SVECT *nd, double c, double *integral);
void    integrate_triPrism_phi_f             (SVECT *nd, double c, double *f, double *integral);
double  integrate_triPrism_df                (SVECT *nd, double c, double *f, int direction);
SVECT   integrate_triPrism_df_full           (SVECT *nd, double c, double *f);
void    integrate_triPrism_dphi              (SVECT *nd, double c, double *integral, int direction);
void    integrate_triPrism_dphi_f            (SVECT *nd, double c, double *f, double *integral, int direction);
void    integrate_triPrism_dphi_f_g          (SVECT *nd, double c, double *f, double *g, double *integral, int direction);
void    integrate_triPrism_gradPhi_dot_f_v   (SVECT *nd, double c, double *f, SVECT *v, double *integral);
void    integrate_triPrism_gradPhi_dot_v     (SVECT *nd, double c, SVECT *v, double *integral);
void    integrate_triPrism_gradPhi_dot_DgradF(SVECT *nd, SQUAD *quad, double c, STENSOR T, double *f, double *integral);
void    integrate_triPrism_gradPhi_dot_DgradV(SVECT *nd, SQUAD *quad, double c, SVECT tx, SVECT ty, double *u, double *v, double *integral_x, double *integral_y);
void    integrate_triPrism_dphi_f_quadrature (SVECT *nd, SQUAD *quad, double c, double *f, double *integral_x, double *integral_y, int quad_order);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// fe_elem_utilities
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void rhs_contrib_3dof(char *string, int nnodes, int ie, int *nodes, DOF_3 *elem_rhst, DOF_3 *elem_rhs);
void rhs_4dof(char *string, int nnodes, int ie, int *nodes, DOF_4 *elem_rhs, int printFieldWidth, int printPrecision);
void rhs_3dof(char *string, int nnodes, int ie, int *nodes, DOF_3 *elem_rhs);
void rhs_contrib_1dof(char *string, int nnodes, int ie, int *nodes, double *elem_rhst, double *elem_rhs);
void rhs_1dof(char *string, int nnodes, int ie, int *nodes, double *elem_rhs);
void printScreen_debug_int(char *descript, int *f, int n);
void printScreen_debug_dbl(char *descript, double *f, int n);
void printScreen_debug2_dbl(char *descript, double *f, int n, int *global_nd_ids);
void printScreen_debug_vec(char *descript, SVECT *f, int n);
void printScreen_debug_svect(char *descript, SVECT *v, int n, int *global_nd_ids);
void printScreen_debug_svec2d(char *descript, SVECT2D *v, int n, int *global_nd_ids);
void printScreen_debug_vec2d(char *descript, SVECT2D *f, int n);
void elem_get_tposition(double *local1, double *local2, double *local3, double scale, double time1, double time2, int n);
double elem_get_tposition_scalar(double local2, double local3, double scale, double time1, double time2);
void elem_get_tposition_vect(SVECT *local1, SVECT *local2, SVECT *local3, double scale, double time1, double time2, int n);
double get_second_order(double local1, double local2) ;
void get_second_order_array(double *local, double *local1, double *local2, int n);
void elem_get_vect_second_order(SVECT *local1, SVECT *local2, SVECT *local3, int n);
void elem_matrix_deriv(int indx, int nnodes, DOF_3 *local1, DOF_3 *local2, DOF_3 *mat, double diff_ep);
void elem_matrix_deriv2(int indx, int nnodes, DOF_3 *local1, DOF_3 *local2, DOF_3 *mat, double diff_ep_inv);
void elem_matrix_deriv_4dof(int indx, int nnodes, DOF_4 *local1, DOF_4 *local2, DOF_4 *mat, double diff_ep);
void elem_matrix_deriv_1dof(int indx, int nnodes, double *local1, double *local2, double *mat, double diff_ep);
void grad_phi_dot_v(SVECT *grad_phi, SVECT *v, SVECT *grad_x, SVECT *grad_y, SVECT *grad_z, int ndof);
void grad_phi_f(SVECT *grad_phi, double *f, SVECT *grad, int ndof);
void grad2d_phi_dot_v(SVECT2D *grad_phi, SVECT2D *v, SVECT2D *grad_x, SVECT2D *grad_y, int ndof);
void grad2d_phi_f(SVECT2D *grad_phi, double *f, SVECT2D *grad, int ndof);
SVECT tensor_vector_product(STENSOR tens, SVECT vect);
SVECT2D tensor2d_vector2d_product(STENSOR2D tens, SVECT2D vect);
SVECT2D tensor2d_vector2d_AI_product(STENSOR2D tens, SVECT2D vect, SVECT2D vdir);
void get_elem3d_node_coordinates(SGRID *grid, int ie, double *elem_dpl, SVECT *elem_nodes);
void get_elem2d_node_coordinates(SGRID *grid, int ie, double *elem_dpl, SVECT *elem_nodes);
void global_to_local_dbl(double *global, double *local, int *nodes, int nnodes);
void global_to_local_int(int *global, int *local, int *nodes, int nnodes);
void global_to_local_svect2d(SVECT2D *global, SVECT2D *local, int *nodes, int nnodes) ;
void global_to_local_svect(SVECT *global, SVECT *local, int *nodes, int nnodes);
void elem_get_midpt_locations(SNODE *nodes, double *z, SVECT *p, int nnodes, int nnodes_quad, int **nd_on_Edge);
void elem_get_midpt_locations2(SNODE *nodes, double *dpl, SVECT *p, int nnodes, int nnodes_quad, int **nd_on_Edge);
void elem_get_relative_velocities(SVECT *elem_vel, double *elem_displacement, double *elem_old_displacement, double *elem_older_displacement, double dt, int nnodes, SVECT *elem_rel_vel);
void  dumpVector2D(SVECT2D *vel, int nnodes, double *u, double *v);
void  dumpVector(SVECT *vel, int nnodes, double *u, double *v, double *w);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// fe_get_supg_tau
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double fe_get_supg_tau_sw(int nnodes, SVECT *node, double g, double elem_avg_depth, double elem_avg_u, double elem_avg_v,
                          double elem_avg_w,double *grad_shp_x, double *grad_shp_y, double *grad_shp_z, double djac,
                          double alpha, int ndim, int tau_method_flag, int le_method_flag);
double fe_get_supg_tau(int nnodes, SVECT *node, double diffusion, double elem_avg_u, double elem_avg_v,
                       double elem_avg_w,double *grad_shp_x, double *grad_shp_y, double *grad_shp_z, double djac,
                       double alpha, int ndim, int tau_method_flag, int le_method_flag);
double get_element_length(int nnodes, SVECT *node, double elem_avg_u, double elem_avg_v, double elem_avg_w,
                          double *grad_shp_x, double *grad_shp_y, double *grad_shp_z, double volume_or_area,
                          int ndim, int method_flag);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
