#include <stdbool.h>

/* cross contamination prototypes */
/* cjt :: this can probably be remedied by taking a further look */
/* cjt :: this is a good guide for code modularization */

/* series */
int sseries_read_id(SIO io, char **data, int);
int sseries_set_type(SMODEL *, char **, int);
void sseries_set_meteor_stations(SSERIES *, SGRID *, int);

/* grid */
#ifdef _MESSG
void sgrid_init(SGRID *, int, MPI_Comm);
void sgrid_alloc_init(SGRID **, SIO *, int, SFILE_OUTPUT, MPI_Comm);
#else
void sgrid_init(SGRID *, int);
void sgrid_alloc_init(SGRID **, SIO *, int, SFILE_OUTPUT);
#endif
void sgrid_print_adapted_ts(SGRID *, SFILE, char *, int, int, int **, int, int *);

void sfile_output_init(SFILE_OUTPUT *file_output);


/* shallow water 2d */
void ssw_2d_print_ts(SSW_2D *, SIO *, SGRID *, double, int, SFILE, char *, int, int, int **, int, int *, int, SFILE_OUTPUT);
void ssw_2d_print_adapted_ts(SSW_2D *, SGRID *, double, int, SFILE, char *, int, int);
void ssw_2d_open_output(SMODEL *);
void ssw_2d_free(SSW_2D *, SGRID *, SFLAGS);
void ssw_2d_calculate_elem_error(SSW_2D *, SGRID *, SMAT *, double);
void ssw_diffusive_calculate_elem_error(SSW_2D *, SGRID *, SMAT *, double);

/* navier stokes 2d */
void sns_2d_print_ts(SNS_2D *, SIO *, SGRID *, double, int, SFILE, char *, int, int, int **, int, int *, int, SFILE_OUTPUT);
void sns_2d_print_adapted_ts(SNS_2D *, SGRID *, double, int, SFILE, char *, int, int);
void sns_2d_open_output(SMODEL *);
void sns_2d_free(SNS_2D *, SGRID *, SFLAGS);
void sns_2d_calculate_elem_error(SNS_2D *, SGRID *, SMAT *, double);
void sns_diffusive_calculate_elem_error(SNS_2D *, SGRID *, SMAT *, double);

/* shallow water 3d */
void ssw_3d_alloc_init(SSW_3D **, SGRID *, SIO *, SFLAGS);
void ssw_3d_init(SSW_3D *, SGRID *);
void ssw_3d_free(SSW_3D *, SGRID *, SFLAGS);
void ssw_3d_print_ts(SSW_3D *, SIO *, SGRID *, double, int, SFILE, char *, int, int, int **, int, int *, int **, int, int *, int, SFILE_OUTPUT);
void ssw_3d_open_output(SMODEL *);
void ssw_3d_open_input(SMODEL *);
void ssw_3d_calculate_elem_error(SSW_3D *, SGRID *, SMAT *, double);
void ssw_3d_node_avg(SSW_3D *, int, int, int, SGRID *);
void ssw_3d_node_avg_sur(SSW_3D *, int, int, int, SGRID *);
void ssw_3d_node_avg_bed(SSW_3D *, int, int, int, SGRID *);

/* navier stokes */
void sns_init(SNS **,  SGRID *, SIO *, SFLAGS);
void sns_print_ts_MPI(SIO *, SNS *, SGRID *, double, SFLAGS, STR_VALUE *);
void sns_print_adapted_ts(SNS *, SGRID *, double, int, SFLAGS, SFILE, char *, int, int);
void sns_print_ts_winds(SGRID *, SWIND *, FILE *, double, int **, int, int *);
void sns_print_ts_waves(SGRID *, SWAVE *, FILE *, double, int **, int, int *);
void sns_print_ts_winds_MPI(SGRID *, SWIND *, FILE *, double, int);
void sns_print_ts_waves_MPI(SGRID *, SWAVE *, FILE *, double, int);
void sns_print_ts(SNS *ns, SIO *io, SGRID *grid, double time, int outfact, SFLAGS flag, SFILE sup, char *proj_name, int it1, int it2, int **ndata, int my_nnode_max, int *my_nnode_ext, int **ndata_sur, int my_nnode_max_sur, int *my_nnode_ext_sur, int flag1, SFILE_OUTPUT);
void sns_print_ts_winds(SGRID *grid, SWIND *winds, FILE *fout, double time, int **ndata, int my_nnode_max, int *my_nnode_ext);
void sns_print_ts_waves(SGRID *grid, SWAVE *waves, FILE *fout, double time, int **ndata, int my_nnode_max, int *my_nnode_ext);
void sns_read_hot(SMODEL *);
void sns_open_output(SMODEL *);
void sns_free(SNS *, SGRID *, SFLAGS);
void sns_alloc_init(SNS **, SGRID *, SIO *, SFLAGS);

/* navier-stokes 3d */
void sns_3d_alloc_init(SNS_3D **, SGRID *, SIO *, SFLAGS);
void sns_3d_init(SNS_3D *, SGRID *);
void sns_3d_free(SNS_3D *, SGRID *, SFLAGS);
void sns_3d_print_ts(SNS_3D *, SIO *, SGRID *, double, int, SFILE, char *, int, int, int **, int, int *, int **, int, int *, int, SFILE_OUTPUT);
void sns_3d_open_output(SMODEL *);
void sns_3d_open_input(SMODEL *);
void sns_3d_calculate_elem_error(SNS_3D *, SGRID *, SMAT *, double);
void sns_3d_node_avg(SNS_3D *, int, int, int, SGRID *);
void sns_3d_node_avg_sur(SNS_3D *, int, int, int, SGRID *);
void sns_3d_node_avg_bed(SNS_3D *, int, int, int, SGRID *);

/* shallow water */
void ssw_init(SSW **,  SGRID *, SIO *, SFLAGS);
void ssw_print_ts(SSW *, SIO *, SGRID *, double, int, SFLAGS, SFILE, char *, int,int,int **, int, int *, int **, int, int *, int, SFILE_OUTPUT);
void ssw_print_ts_MPI(SIO *, SSW *, SGRID *, double, SFLAGS, STR_VALUE *);
void ssw_print_adapted_ts(SSW *, SGRID *, double, int, SFLAGS, SFILE, char *, int, int);
void ssw_print_ts_winds(SGRID *, SWIND *, FILE *, double, int **, int, int *);
void ssw_print_ts_waves(SGRID *, SWAVE *, FILE *, double, int **, int, int *);
void ssw_print_ts_winds_MPI(SGRID *, SWIND *, FILE *, double, int);
void ssw_print_ts_waves_MPI(SGRID *, SWAVE *, FILE *, double, int);
void ssw_read_hot(SMODEL *);
void ssw_open_output(SMODEL *);
void ssw_free(SSW *, SGRID *, SFLAGS);
void ssw_alloc_init(SSW **, SGRID *, SIO *, SFLAGS);

/* ground water */
#ifdef _ADH_GROUNDWATER
void sgw_alloc_init(SGW **, SGRID *, SIO *, SFLAGS );
void sgw_init(SGW *, int, int, int, int);
void sgw_free(SGW *, SGRID*, SFLAGS);
void sgw_read_hot(SMODEL *);
void sgw_print_ts(SMODEL *, SGW *, SIO *, SGRID *, double, int, SFLAGS, SFILE, char *, int,int,int **, int, int *, int **, int, int *, int);
void sgw_open_output(SMODEL *);
void sgw_open_input(SMODEL *);
void sgw_nodal_sat_avg(SGW* gw, SGRID* grid, int max_num_nodes, double *sat_avg, int* cnt_arr);
void sgw_elem3d_flux_to_nodal_vel(SGRID* grid, SMAT *mat, SGW* gw, SVECT *vel_avg);
void sgw_set_psk_series(SMODEL *);
void sgw_3d_calculate_elem_error(SGW *gw, SGRID *grid, SMAT *mat, double dt);
void sgw_3d_renumber(SGW *gw, int max_nnode, int *new_numbers, int *order_tmp, double *dtmp, SVECT2D *v2tmp, SVECT *vtmp);
void sgw_3d_realloc_init(SGW *gw, int nnodes_old, int nnodes_new);
void sgw_3d_elem_realloc_init(SGW *gw, int nelems3d_new, int nelems3d_old, int nalloc_inc);
void sgw_3d_node_avg(SGW *gw, int node_new, int node1, int node2, SGRID *grid);
void sgw_3d_elem_avg(SGW *gw, int orig_elem, int new_elem, int new_node_in_orig_elem, int new_node_in_new_elem, SGRID *grid);
double sgw_eval_sat(SSERIES*, double , double *, double *);
double sgw_eval_kr(SSERIES*, double , double *, double *);
double sgw_eval_kr_elem(SSERIES*, int, double *);
void sgw_evaluate_element_saturations(SMODEL *mod);
void smat_gw_init(SMODEL *, SMAT_GW *);
void smat_gw_check(SMODEL *, SMAT_GW, int);
void sgw_eval_flux(int, SSERIES *, const STENSOR, double, double*, double *, SVECT *, SVECT *);
void sgw_evaluate_element_fluxes(SMODEL *);
#endif 
/* materials */
void smat_init(SMODEL *, SMAT **, int);
void smat_free(SMODEL *, SMAT *, int);
void smat_checkall(SMODEL *);
void smat_wnsm_init(SMODEL *, SMAT_NSM *);
void smat_trn_init(SMODEL *, SMAT_TRN *);
void smat_trn_check(SMODEL *mod, int EEVF, int EVSF, int DIFF_FLAG, int imat, int id);
void smat_sw_init(SMODEL *, SMAT_SW *);
void smat_sw_check(SMODEL *, SMAT_SW, int);
void smat_ns_init(SMODEL *, SMAT_NS *);
void smat_ns_check(SMODEL *, SMAT_NS, int);

/* strings */
void sstr_value_set_total_area(STR_VALUE *, SGRID *);

/* sediment */
#ifdef _SEDIMENT
void ssediment_node_avg(SSED *, int, int, int, SGRID *);
void ssediment_prep(SMODEL *, SSED *, int);
void ssediment_open_output(SMODEL *);
void ssediment_print_ts(SGRID *, SIO *, SSED *, double, double, int, SFILE, char *, int, int, int **, int, int *, int **, int, int *, int);
void ssediment_print_adapted_ts(SGRID *, SSED *, double, double, int, SFILE, char *, int, int);
void ssediment_read_hot(SMODEL *);
void ssediment_calculate_2d_error(SSED *, SSW_2D *, SGRID *, SMAT *, double);
void ssediment_calculate_3d_error(SSED *, SSW_3D *, SGRID *, SMAT *, double);
#endif

/* scon */
void scon_print_ts(int, SCON *, SIO *, SGRID *, double, int, SFLAGS, SFILE, char *, int, int, int **, int, int *, int);
void scon_print_adapted_ts(int, double, SCON *, SGRID *, double, SFLAGS, SFILE, char *, int, int);
void scon_open_output(SMODEL *);
void scon_calculate_elem2d_error(int, SCON *, SVECT2D *, double *, SGRID *, SMAT *, double);
void scon_calculate_elem3d_error(int, SCON *, SVECT *, double *, double *, double *, SGRID *, SMAT *, double);

/* swindlib */
#ifdef WINDLIB
void swindlib_firstread(SMODEL *mod, SWINDLIB *wl, int nnodes);
void swindlib_update(SMODEL *mod, SWINDLIB *wl, double t_curr, int nnodes, double **values);
#endif

void svect2d_print_array_MPI(SGRID *, FILE *, FILE *, SVECT2D *, int **, int, int *, int);
void svect_print_array_MPI(SGRID *, FILE *, FILE *, SVECT *, int **, int, int *, int);

/*mpi*/
#ifdef _MESSG
void smpi_init(SMPI *, MPI_Comm);
void smpi_realloc(SMPI *,MPI_Comm);
#else
void smpi_init(SMPI *);
void smpi_realloc(SMPI *);
#endif
void smpi_free(SMPI *);
void smpi_defaults(SMPI *);
