#include "fe_prototypes.h"

#ifdef _ADH_HDF5
void xdmf_init(SMODEL *mod, int npes, int myid);
void xdmf_finalize(HDF5 *hdf5, SIO *io, int npes, int myid, int flag);
void xdmf_finalize_supermodel(SSUPER_MODEL *sm, int nsupermodels, int myid);
void ps_print_xdmf(SMODEL *, double, double);
void ps_print_geo_xdmf(SMODEL *);
void ps_print_surf_geo_xdmf(SMODEL *);
#else
/* Functions from xdmf_utils.c */
void ps_print_geo_xdmf(void);
void ps_print_surf_geo_xdmf(void);
void ps_print_ts_xdmf(void);
void ps_print_surf_ts_xdmf(void);
void xdmf_write_dataset(void);
void xdmf_read_dataset(void);
void xdmf_read_npoints(void);
void xdmf_read_extents(void);
void xdmf_read_ndims(void);
void xdmf_write_attribute(void);
void xdmf_read_attribute(void);
#endif
// main

#ifdef _MESSG
int adh_init_func_(SSUPER_MODEL **, SSUPER_INTERFACE **, int *, int *, int *argc, char ***argv, MPI_Comm *comm);
#else
int adh_init_func_(SSUPER_MODEL **, SSUPER_INTERFACE **, int *, int *, int *argc, char ***argv);
#endif
int adh_run_func_(SSUPER_MODEL *, SSUPER_INTERFACE *, int, int, double *);
int adh_finalize_func_(SSUPER_MODEL **, SSUPER_INTERFACE **, int, int);

// elem
void edge_hash_add_entry(int, int, EDGE_LIST_ITEM **);
int  edge_hash_index(int, int);
EDGE_LIST_ITEM *edge_hash_lookup(int, int, EDGE_LIST_ITEM **);
ELEM2D_LIST_ITEM *elem1d_find_elem2d(int, int,  ELEM2D_LIST_ITEM **);
void elem1d_find_elem2d_init(SGRID *);
void elem1d_hash_add_entry(int, int, ELEM1D_LIST_ITEM **, int);
ELEM1D_LIST_ITEM *elem1d_hash_lookup(int, int, ELEM1D_LIST_ITEM **);
int  elem1d_new(SGRID *, int);
void elem1d_outward_nrml(SGRID *);
void elem1d_renumber(SGRID *);
void elem1d_split(SMODEL *mod, int, EDGE_LIST_ITEM **);
#ifdef _MESSG
void elem1d_packi(void *, int, SGRID *);
void elem1d_packd(void *, int, SGRID *);
void elem1d_unpacki( void *, void *, NODE_LIST_ITEM **, ELEM1D_LIST_ITEM **, SMODEL *);
int  elem1d_pack_icnt(void);
int  elem2d_pack_icnt(void);
void elem2d_unpacki( void *, void *, NODE_LIST_ITEM **, ELEM2D_LIST_ITEM **, SMODEL *);
void elem3d_unpacki( void *, void *, NODE_LIST_ITEM **, ELEM3D_LIST_ITEM **, SMODEL *);
int  elem3d_pack_icnt(void);
void elem2d_packi(void *, int, SGRID *);
void elem2d_packd(void *, int, SGRID *);
void elem3d_packi(void *, int, SGRID *);
void elem3d_packd(void *, int, SGRID *);
#else
void elem1d_packi(void);
void elem1d_packi(void);
void elem1d_pack_icnt(void);
void elem2d_pack_icnt(void);
void elem3d_pack_icnt(void);
void elem2d_packi(void );
void elem2d_packd(void );
void elem3d_packi(void );
void elem3d_packd(void );
#endif
void elem2d_hash_add_entry(int, int, int, ELEM2D_LIST_ITEM **, int);
ELEM2D_LIST_ITEM *elem2d_hash_lookup(int, int, int, ELEM2D_LIST_ITEM **);
int elem2d_new(SGRID *, int, int); /* Gajanan gkc adding int nnodes_on_elem */
SVECT  elem2d_normal(SNODE *);
void   elem2d_renumber(SGRID *);
void   elem2d_split(SMODEL *, int, EDGE_LIST_ITEM **);
void   elem3d_split(SMODEL *, int, EDGE_LIST_ITEM **, int);
void   elem3d_hash_add_entry(int, int, int, int, ELEM3D_LIST_ITEM **, int);
ELEM3D_LIST_ITEM *elem3d_hash_lookup(int, int, int, int, ELEM3D_LIST_ITEM **);
void   elem3d_get_mg_lin_grad(SGRID *, SELEM_3D *, double *);
void   elem3d_lin_grad(SGRID *, SELEM_3D *);
void   elem3d_get_djacs_gradPhi(SELEM_3D *elem3d, SNODE *node);
int    elem3d_new(SGRID *, int, int); /* Gajanan gkc adding int nnodes_on_elem */
void   elem3d_renumber(SGRID *);
double elem3d_get_volume(SVECT *node, int nnodes);
int elem_hash_index(int);
int elem3d_level(SELEM_3D);
int elem2d_level(SELEM_2D);
int elem1d_level(SELEM_1D);
#ifdef _MESSG
int elem1d_merge(SGRID *, int, int *, NODE_LIST_ITEM **, ELEM1D_LIST_ITEM **, int *);
int elem2d_merge(SGRID *, int, int *, NODE_LIST_ITEM **, ELEM2D_LIST_ITEM **, int *);
int elem3d_merge(SGRID *, int, int *, NODE_LIST_ITEM **, ELEM3D_LIST_ITEM **, int *);
#else
int elem1d_merge(SGRID *, int, int *, ELEM1D_LIST_ITEM **, int *);
int elem2d_merge(SGRID *, int, int *, ELEM2D_LIST_ITEM **, int *);
int elem3d_merge(SGRID *, int, int *, ELEM3D_LIST_ITEM **, int *);
#endif

 double get_elem2d_area2d(SELEM_2D *elem2d, SNODE *node);
 double get_elem3d_volume(SVECT *node, int nnodes);
 SVECT  get_elem2d_normals(SVECT *nd);
 void get_elem1d_linear_djac_gradPhi(SGRID *grid, SELEM_1D *elem1d);
 SVECT get_triangle_centroid(SNODE nd1, SNODE nd2, SNODE nd3);
 double get_triangle_orientation(SNODE nd1, SNODE nd2, SNODE nd3, SVECT ref_vec);
 void get_triangle_local_shape(double xhat, double yhat, double zhat, double *lshape);
 void get_triangle_local_shape_quad(double xhat, double yhat, double zhat, double *lshape);
 void get_triangle_local_shape_gradients(SVECT *lgrad_shp);
 void get_triangle_linear_djac_nrml_gradPhi(SELEM_2D *elem2d, SNODE *nd_snode, SVECT *nd_svect);
 void get_quadrilateral_local_shape(double xhat, double yhat, double zhat, double *lshape);
 void get_quadrilateral_local_shape_quad(double xhat, double yhat, double zhat, double *lshape);
 void get_quadrilateral_local_shape_gradients(double xhat, double yhat, double zhat, SVECT *lgrad_shp);
 double get_quadrilateral_linear_djac2d(double xhat, double yhat, SVECT *nd);
 double get_quadrilateral_linear_djac_gradPhi(double xhat, double yhat, SVECT *nd, SVECT *grad_shp);
 void get_tet_local_shape(double xhat, double yhat, double zhat, double *lshape);
 void get_tet_local_shape_quad(double xhat, double yhat, double zhat, double *lshape_quad);
 void get_tet_local_shape_gradients(SVECT *lgrad_shp);
 void get_tet_linear_djac_gradPhi(SELEM_3D *elem3d, SNODE *nd_SNODE, SVECT *nd_SVECT);
 double get_tet_linear_djac_gradPhi2(SNODE *nd_SNODE, SVECT *nd_SVECT, SVECT *grad_shp); 
 double get_tet_linear_djac(SNODE *nd_SNODE, SVECT *nd_SVECT);
 void get_triprism_local_shape(double xhat, double yhat, double zhat, double *lshape);
 void get_triprism_local_shape_quad(double xhat, double yhat, double zhat, double *lshape_quad);
 void get_triprism_local_shape_gradients(double xhat, double yhat, double zhat, SVECT *lgrad_shp);
 double get_triprism_djac(double xhat, double yhat, double zhat, SVECT *nd);
 double get_triprism_linear_djac_gradPhi(double xhat, double yhat, double zhat, SVECT *nd, SVECT *grad_shp);
double get_triprism_volume(SVECT *node);

// grid
int adpt_get_node(SMODEL *, EDGE_LIST_ITEM *, EDGE_LIST_ITEM **);
void adpt_main(SMODEL *);
void adpt_fix_global_ids(SGRID *);
#ifdef _MESSG
void adpt_rank_edges(EDGE_LIST_ITEM **, SGRID *, int * ,int *);
#else
void adpt_rank_edges(EDGE_LIST_ITEM **, SGRID *);
#endif
long adpt_ref(SMODEL *, long);
void adpt_set_flags(SGRID *, int *);
void adpt_unref(SMODEL *);
void adpt_set_node_number(SGRID *);
void adpt_unset_node_number(SGRID *);
void findBoundaryFaces(SGRID * );
void error_main(SMODEL *);
#ifdef _MESSG
void findrank(EDGE_RANK *, int, int *, SMPI *);
void adpt_fix_adj(SGRID *grid, int *, NODE_LIST_ITEM **);
void adpt_merge_elem(SGRID *, int *, ELEM_REF_LIST_ITEM **, ELEM_REF_LIST_ITEM **, ELEM_REF_LIST_ITEM **, NODE_LIST_ITEM **, ELEM3D_LIST_ITEM **, ELEM2D_LIST_ITEM **, ELEM1D_LIST_ITEM **);
#else
void adpt_fix_adj(SGRID *grid, int *);
void adpt_merge_elem(SGRID *, int *, ELEM3D_LIST_ITEM **, ELEM2D_LIST_ITEM **, ELEM1D_LIST_ITEM **);
#endif

// tools
double projected_node_distance(SGRID *grid1, int n1, SGRID *grid2, int n2);
int is_my_node(int global_id, SGRID *g);
int is_my_elem(int node1, int node2, SGRID *g);
int doesFileExist(const char *fname);
double get_coriolis_angular_speed(double coriolis_factor);
void createChildMesh(SGRID *grid, SIO *io);
double erf(double);
double erfc(double);
int tc_end(SMODEL *);
void tc_init(SMODEL *);
void tc_timeunits(FILE *, int);
double tc_eval_series(SSERIES, int, double, int);
#ifdef _ADH_GROUNDWATER
double tc_eval_series_slope(SSERIES, int, double, int);
#endif
double tc_conversion_factor(int, int);
double tc_trap_area(double, double, double, double);
void tc_scale(double *
#ifdef _MESSG
    , MPI_Comm
#endif
    );
SVECT2D tl_bendway_correction(SVECT2D *, SVECT2D *, double *, double *, double, double, double, double, double, int, int);
double tl_find_grid_mass_elem2d(double density, STR_VALUE *str, SSERIES *series_head, double *depth, SGRID *grid, SFLAGS flags);
double tl_find_grid_mass_elem3d(double density, SGRID *, double *);
double tl_find_grid_mass_error_elem2d(double density, double *depth, SVECT2D *vel, SGRID *grid, SFLAGS flags, double initial_grid_mass, SSERIES *series_head, STR_VALUE *str, double dt, double *total_time_mass_flux_T);
double tl_find_3d_grid_mass_error(STR_VALUE *str, SSERIES *series_head, double initial_grid_mass, double density, SGRID *grid, SVECT *vel, double *displacement, double *        old_displacement, double *, double dt, double *new_grid_mass, double *);
void tl_3d_to_2d(void);
void tl_list_setup(void);
void *tl_list_alloc(int);
void tl_list_free_all(int);
int tl_long_edge(SNODE *, int, int, int, int, EDGE_LIST_ITEM **);
void create_dbovl_elem1d(SMODEL *);
void printScreen_int_array(char *, int *, int, int, char *);
void printScreen_dble_array(char *, double *, int, int, char *);
void printScreen_resid(char *, double *, int, int, int, char *);
void printScreen_matrix(char *, double *, SPARSE_VECT *, int, int, int, char *);
void printFile_matrix(FILE *fp, char *descript, double *diagonal, SPARSE_VECT * matrix, int nnode, int nsys_sq, int linenumber, char *filename);
void printScreen_matrixProfile(char *, Profile_Info *, int, char *);
void tl_model_units(SMODEL *);
void tl_find_edge_mdpt_pressure(SGRID *, SSW_3D *, double, double);
void tl_calculate_pressure(SGRID *, SSW_3D *, double, double);
void tl_calculate_depthavgvel(SGRID *, SSW_3D *);
void tl_density_calculator_metric(double,  double *, double, double *, double, int, double *, int);
double tl_find_avg_column_depth(SGRID *, int *, double *);
void tl_find_common_int(int *);
void tl_vertical_adapt(SGRID *, double *);
void tl_vertical_adapt_ALE(SGRID *, SSW_3D *);
void tl_get_dpl_perturbation(SGRID *, double, double *, double *);
void tl_nodal_grad_avg_2d(SGRID *, STR_VALUE *, double *, SVECT2D *);
void tl_nodal_grad_avg_3d(SGRID *, STR_VALUE *, double *, SVECT2D *);
int tl_find_nodes_on_sliceX(SGRID *, int, int, int **);
int tl_find_nodes_on_sliceY(SGRID *, int, int, int **);
int tl_find_surface_nodes_on_sliceX(SGRID *, int, int, int **);
int tl_find_surface_nodes_on_sliceY(SGRID *, int, int, int **);
void root_print(char *);
int isPointInTriangle(double x1, double y1, double x2, double y2, double x3, double y3, double x, double y);
double factorial(int n);
double besselj(double p, double z);


// initio
void read_flux_geo(int super, int sub,SIO *io, SGRID *grid);
void read_superfile(char *infile, int *nmodels, char ***modelnames, int ***coupling);
void pdata(char *, char *, FILE *, int, char *, char *, char *, char *, int, int);
void print_build_info(void);
void print_runtime_info(SIO);
void print_trailer(FILE *);
void print_xms_trailers(SMODEL *);
void print_header(SMODEL *, FILE *, int);
void print_header_sw2d(SMODEL *, FILE *, int);
void print_header_sw3d(SMODEL *, FILE *, int);
void print_header_sediment(SMODEL *, FILE *, int, int);
void open_output_file(SFILE *, char *, int);
void open_input_file(SFILE_IN *, char *, int);
char * build_filename(char *, const int, const char *, const char *);
char * build_filename2(char *, const int, const char *, const char *, const int, const char *, const int);
char * build_filename3(char *, const int, const char *, const char *, const int, const char *, const int, const char *, const int);
char * convert_to_uppercase(char *);
int get_index_as_text(const int, char *);
FILE * io_fopen_prep_ds(SMODEL *, const char *, const int, const int);
FILE * io_fopen_prep_ds_index(SMODEL *, const char *, const int, const int, const int);
int io_read_error(SIO, const char *, const int);
int io_save_line(SIO *, FILE *, const char *, const char *);
CARD parse_card(char *, char **);
double read_dbl_field(SIO, char **);
double read_dbl_field_custom(SIO, char **, int *, char *, int, int);
int read_int_field(SIO, char **);
int read_int_field_custom(SIO, char **, int *, char *, int, int);
char * read_text_field(SIO, char **, char *, int );
char * read_text_field_custom(SIO, char **, char *, int, int *, char *, int, int);
double read_query_dbl_field(SIO, char **, int *);
int read_query_int_field(SIO, char **, int *);
char * read_query_text_field(SIO, char **, char *, int, int *);
int strip_comments(char *);
int strip_white(char **);
int get_string_id(SIO, char **, int);
int get_series_id(SIO, char **, int);
int get_material_id(SIO, char **, int);
int get_transport_id(SIO, char **, int);
int get_node_id(SIO, char **, int);
int get_bed_layer(SIO, char **, int);
void read_physics(SMODEL *);
int read_elems_nodes(SIO *, int *, int *, int *, int *, int *, int *);
void read_elems_nodes_mpi(SIO *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int, int);
int free_node_list(NODE_LIST_ITEM *);
int read_geo(SIO *, SGRID *);
#ifdef _MESSG
int read_geo_mpi(SIO *, SGRID *);
#endif
void read_faces(SIO *io, SGRID *grid, int JUST_COUNT);
void read_bc_prep(SMODEL *);
void read_bc(SMODEL *);
void read_bc_DB(SMODEL *, char *);
void read_bc_FR(SMODEL *, char *);
void read_bc_IP(SMODEL *, char *);
void read_bc_MP(SMODEL *, char *);
void read_bc_NB(SMODEL *, char *);
void read_bc_OFF(SMODEL *, char *);
void read_bc_OP(SMODEL *, char *);
void read_bc_PC(SMODEL *, char *);
void read_bc_STRINGS(SMODEL *, CARD, char *, int *, int *);
void read_bc_TC(SMODEL *, char *);
void read_bc_DEBUG(SMODEL *, char *);
void read_bc_NOTERM(SMODEL *, char *);
void read_bc_SERIES(SMODEL *, int, char *);
void read_bc_AWRITE(SMODEL *, char *);
void read_bc_FILE_OUTPUT(SMODEL *, char *);
void read_bc_SCREEN_OUTPUT(SMODEL *, char *);
void read_bc_CN(SMODEL *, char *);
void read_bc_nsm(SMODEL *, char *);
void read_bc_SEDLIB_prep(SMODEL *, char *);
void read_bc_SEDLIB(SMODEL *, char *);
void read_bc_structures(SMODEL *, char *);
void read_bc_WINDLIB(SMODEL *, char *);
void read_bc_NSM(SMODEL *, char *);
int read_data_set(SIO, SFILE_IN *, SGRID *, double *, char *);
void print_double_array(FILE *, double *, int, double);
int assign_series_type(SMODEL *, char **, int);
void write_testcase_error(SSUPER_MODEL *);
void testcase_prep(SSUPER_MODEL *);
void testcase_clean();
void read_bc_TEST(SMODEL *, char *);
void read_bc_NSM(SMODEL *mod, char *data);
int read_interface_geo(SIO *io, SGRID *grid);
#ifdef _ADH_GROUNDWATER
void read_bc_GW(SMODEL *, char*);
void read_bc_GW_MP(SMODEL *,char*);
#endif

// tl_columns
void build_columns(SGRID *, int);
void build_column_hash(SGRID **);
void build_column_list(SGRID **);
void build_vertical_list(SGRID **);
void build_midpt_list(SGRID **);
void classify_2d_elements(SGRID **);
void tl_find_surface_elements(SGRID **);
int tl_column_init(SGRID *);
SVECT2D tl_calc_3dcentroid(SGRID *, int);
SVECT2D tl_calc_2dcentroid(SGRID *, int);
int column_hash_function(SGRID *, SVECT2D);
void Add_column_hash_entry(SGRID *, CENT_LIST_ITEM **, int);
void Add_column_hash2d_entry(SGRID *, CENT_LIST_ITEM **, int);
void Add_column_entry(ID_LIST_ITEM **, int);
void hash_histogram(SGRID *, CENT_LIST_ITEM **);
int find_column_from_centroid(SGRID *, int, CENT_LIST_ITEM **);
int find_column_from_center(SGRID *, int, CENT_LIST_ITEM **);
int find_column_from_point(SGRID *, SVECT2D, CENT_LIST_ITEM **);
int find_2d_hashvalue(SGRID *, int);
int find_3d_hashvalue(SGRID *, int);
int pre_intcompare(const void *, const void *);
int pre_dblcompare(const void *, const void *);
void Add_sorted_vertical_entry(SGRID *, ID_LIST_ITEM **, int);
int length_node_list(ID_LIST_ITEM *);
void Add_vertical_node_hash_entry(SGRID *, CENT_LIST_ITEM **, int);
int find_sidewall_column(SGRID *, int, CENT_LIST_ITEM **);
void build_sidewall_list(SGRID **);
void Add_sidewall_entry(ID_LIST_ITEM **, int);
void push_vertical_hash_entry(SGRID *, CENT_LIST_ITEM **, SVECT2D);
void push_column_hash_entry(SGRID *, CENT_LIST_ITEM **, SVECT2D);
void free_column_hash(SGRID *, CENT_LIST_ITEM **);
void free_id_list(ID_LIST_ITEM **, int);
int find_vertical_segment(SGRID *, int, CENT_LIST_ITEM **);
void find_bottom_elements2d(SGRID **);
void find_bottom_nodes(SGRID **);
void build_vertical_hash(SGRID **);
void free_vertical_hash(SGRID *, CENT_LIST_ITEM **);
int length_column_list(CENT_LIST_ITEM *);
void free_midpt_hash(SGRID *, CENT_LIST_ITEM **);
void build_midpt_hash(SGRID **);
void Add_midpt_hash_entry(SGRID *, CENT_LIST_ITEM **, int, int);
void push_midpt_hash_entry(SGRID *, CENT_LIST_ITEM **, SVECT2D);
int find_midpt_segment(SGRID *, SVECT, CENT_LIST_ITEM **);
void free_midpt_list(MIDPT_LIST_ITEM **, int);
void Add_sorted_midpt_entry(SGRID *, MIDPT_LIST_ITEM **, SVECT, int, int);
int set_value_midpt_list(SGRID *, MIDPT_LIST_ITEM **, int, int, double, int);
int get_value_midpt_list(SGRID *, MIDPT_LIST_ITEM **, int, int, double *);
void build_midpt_connectivity(SGRID *);
void set_midpt_elem(SGRID *, MIDPT_LIST_ITEM **, int, SVECT, int, int, int);
int above_or_below(SGRID *, int, int, int);
MIDPT_LIST_ITEM *get_midpt_entry(SGRID *, MIDPT_LIST_ITEM **, int, int);
void free_column_hash2d(SGRID *, CENT_LIST_ITEM **);
void global_free_columns(SGRID *);
void column_adpt_merge_elem(SGRID *grid, int, int, int *, ELEM_REF_LIST_ITEM **, NODE_LIST_ITEM **,  ELEM3D_LIST_ITEM **, ELEM2D_LIST_ITEM **);
long column_adpt_ref(SMODEL *mod, long cycle);
void column_adpt_unref(SMODEL *mod);
void column_adpt_set_flags(SGRID *grid, int *surf_unref_flags, int *node_unref_flags);
int column_elem2d_merge(
#ifdef _MESSG
                        SGRID *grid,
                        int ielem, /* the right element to be merged */
                        int *new_local, /* the node of the element that is being removed */
                        int *new_elem_num_pntr,    /* the pointer to the merged element number */
                        NODE_LIST_ITEM ** node_hashtab,    /* the node hashtable */
                        ELEM2D_LIST_ITEM ** elem2d_hashtab /* the 2D element hash table */
#else
                        SGRID *grid,
                        int ielem, /* the right element to be merged */
                        int *new_local, /* the node of the element that is being removed */
                        int *new_elem_num_pntr,    /* the pointer to the merged element number */
                        ELEM2D_LIST_ITEM ** elem2d_hashtab /* the 2D element hash table */
#endif
); /* Gajanan GKC changing to add Quads. */
int column_elem2d_split(SMODEL *mod,
                        int ielem,  /* the element to be split */
                        int *split_edge, /* the edge to split */
                        EDGE_LIST_ITEM ** edge_hashtab  /* the hash table of edges */
); /* Gajanan GKC changing to add Quads. */
void column_elem3d_merge(SGRID *grid,
                         int ielem, /* the right element to be merged */
                         int *new_local, /* the node of the element that is being removed */
                         int *node_unref_flags, NODE_LIST_ITEM ** node_hashtab, /* the node hashtable */
                         ELEM3D_LIST_ITEM ** elem3d_hashtab /* the 3D element hash table */
);  /* Gajanan GKC changing to add prisms*/
void column_elem3d_split(
                         SMODEL *mod, int icol,  /* the column we are currently in */
                         ID_LIST_ITEM ** column_list,   /* the linked list for the columns */
                         EDGE_LIST_ITEM ** edge_hashtab /* the hash table of edges */
);






/*mesh*/
void mesh_pack_node(int, int *, MESSG_BUFFER *, MESSG_BUFFER *, MESSG_BUFFER *, MESSG_BUFFER *, SMODEL *);
void mesh_unpack_node(MESSG_BUFFER *, MESSG_BUFFER *, MESSG_BUFFER *, MESSG_BUFFER *, SMODEL *, NODE_LIST_ITEM **, int);
void mesh_pack_elem1d(int *, MESSG_BUFFER *, MESSG_BUFFER *, SGRID *);
void mesh_pack_elem2d(int *, MESSG_BUFFER *, MESSG_BUFFER *, SGRID *);
void mesh_pack_elem3d(int *, MESSG_BUFFER *, MESSG_BUFFER *, SGRID *);
void mesh_unpack_elem1d(NODE_LIST_ITEM **, ELEM1D_LIST_ITEM **, MESSG_BUFFER *, MESSG_BUFFER *, SMODEL *);
void mesh_unpack_elem2d(NODE_LIST_ITEM **, ELEM2D_LIST_ITEM **,MESSG_BUFFER *, MESSG_BUFFER *, SMODEL *);
void mesh_unpack_elem3d(NODE_LIST_ITEM **, ELEM3D_LIST_ITEM **,MESSG_BUFFER *, MESSG_BUFFER *, SMODEL *);

/* messg */
#ifdef _MESSG
void comm_update_double(double *, int, SMPI *);
void comm_update_double_surf(double *, int, SMPI *);
void comm_update_int(int *, int, SMPI *);
void comm_update_int_surf(int *, int, SMPI *);
void comm_update_VECT2D(SVECT2D *, SMPI *);
void comm_update_VECT(SVECT *, SMPI *);
void comm_update_snode(SGRID *);
void comm_update_ns3(SNS_3D *, SGRID *);
#ifdef _ADH_GROUNDWATER
void comm_update_gw(SGW *gw, SGRID *grid);
#endif
void comm_update_sw3(SSW_3D *, SGRID *);
void comm_update_sw3_surface(SSW_3D *, SGRID *);
void comm_update_sw2(SSW_2D *, SGRID *);
void comm_update_con(SCON *, SGRID *, int);
void comm_update_swaves(SGRID *grid, SWAVE *waves, int cstorm_flag);
#ifdef _SEDIMENT
void comm_update_sed(SSED *, SGRID *, int, int, int);
#endif
void comm_new_node_nums(int *, int **, int , int, SMODEL *);
void comm_update_GN(int, SGRID *);
void comm_node_data(int *, int **, SNODE *, int, SMODEL *);

void comm_set_keys(SGRID *);
void comm_set_keys_supermodel(SSUPER_MODEL *);
void comm_set_keys_supermodel_wvel(SSUPER_MODEL *);
void comm_check(SGRID *);
void comm_edge_setup(EDGE_LIST_ITEM **, SGRID *, int *, int *);
void comm_flag_edges(EDGE_LIST_ITEM **, SGRID *, int *, int *, SMODEL *);
void comm_update_edges(EDGE_LIST_ITEM **, SGRID *, int *, int *, SMODEL *);
void comm_adj_nodes(int *, NODE_LIST_ITEM **, SMODEL *);
void comm_elem_levels(ELEM_REF_LIST_ITEM **, ELEM_REF_LIST_ITEM **, ELEM_REF_LIST_ITEM **, NODE_LIST_ITEM **, ELEM3D_LIST_ITEM **, ELEM2D_LIST_ITEM **, ELEM1D_LIST_ITEM **, SGRID *);
void partition_main(SMODEL *,int);
void partition_cleanup(SGRID *, int);
void partition_form(SGRID *);
void partition_form_final(SGRID *g);
void partition_adpt(SGRID *);
void partition_adpt_final(SGRID *);
void partition_transfer(SMODEL *);
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
int sort_key(const void *, const void *);

#else
int messg_comm_rank(void);
int messg_comm_size(void);
double messg_dmax(double);
double messg_dmin(double);
double messg_dsum(double);
int    messg_imax(int);
int    messg_imin(int);
void   messg_max_norm_loc(double *, int *, SVECT *);

#endif

/* meteor */
void station_node_contrib(SSERIES *, double *, double *, int);
void meteor_update_series(int, double, double, SSERIES *, double **);
void meteor_update(SMODEL *, double);
void winds_proc_teeter(double, double, double, double, double *, double *);
void winds_proc_wu(double, double, double, double *, double *);

/* NSM WQ */
void nsmwq(SMODEL *);
/* Structures */
void weir_flow(SNODE *, int, int, double *, SVECT2D *, STR_VALUE *, SWEIR_C *, SFLAGS, double, SGRID *);
void weir_compound_flux(int, SNODE *, int, double *, SVECT2D *, int, STR_VALUE *, SWEIR_C *, double, SGRID *);
void flap_compound_flux(int, SNODE *, int, double *, SVECT2D *, int, STR_VALUE *, SFLAP_C *, double, SGRID *);
void flap_flow(SNODE *, int, int, double *, SVECT2D *, STR_VALUE *, SFLAP_C *, SFLAGS, double, SGRID *);
void sluice_flow(SNODE *, int, int, double *, SVECT2D *, STR_VALUE *, SSERIES *, SSLUICE_C *, SFLAGS, double, SGRID *);
void sluice_compound_flux(int, SNODE *, int, double *, SVECT2D *, int, STR_VALUE *, SSERIES *, SSLUICE_C *, double, SGRID *);
/* turbulence */
double tur_MY_2 (double, double, double, double, double, double, double, int, int);
double tur_smag (double, double, double, double, double, double);
double tur_ws (double, double, double);

/* sedlib wrapping */
SVECT2D sed_vorticity_velocity_components (double, SVECT2D, SCON, int);
void sedlib_link_bed(SMODEL *, int, int, int, SVECT2D);
void sedlib_link_transport(SMODEL *, int, int, int, SVECT2D);

/* node */
void node2elem2d_map(SGRID *);
void node_avg(SMODEL *, int, int, int);
int node_cmp(SNODE *, SNODE *);
int node_new(SMODEL *);
int node_new_surface(SMODEL *);
int node_new_bed(SMODEL *);
void node_order(SGRID *, int, int *);
void node_renumber(SMODEL *, int);
void node_renumber_surface(SMODEL *);
void node_renumber_double(int, double *, double *, int *, int *);
void node_renumber_int(int, int *, int *, int *, int *);
void node_renumber_snode(int, SNODE *, int *, int *);
void node_renumber_tensor2d(int, STENSOR2D *, STENSOR2D *, int *);
void node_renumber_vect(int, SVECT *, SVECT *, int *, int *);
void node_renumber_vect2d(int, SVECT2D *, SVECT2D *, int *, int *);
#ifdef _MESSG
int node_get_local(SNODE,  NODE_LIST_ITEM **, SMODEL *);
int node_hash_lookup(SNODE, NODE_LIST_ITEM **, int);
void node_hash_add_entry(SNODE, NODE_LIST_ITEM **, int, int);
int node_hash_index(SNODE, int);
int node_pack_icnt(void);
int node_pack_dcnt(SMODEL *);
int node_pack_dcnt_surface(SMODEL *);
int node_pack_dcnt_bed(SMODEL *);
void node_packi(void *, int, SGRID *);
void node_packd(void *, int, SSW *,
#ifdef _ADH_GROUNDWATER
        SGW *,
#endif
        SGRID *, int, SCON *
#ifdef _SEDIMENT
    , SSED *,
    int,
    int,
    int
#endif
    );
void node_packd_surface(void *, int, SSW *,SGRID *

    );
void node_packd_bed(void *, int, SSW *,SGRID *
#ifdef _SEDIMENT
    , SSED *,
    int,
    int,
    int
#endif
    );
void node_unpacki(void *, int, SGRID *);
void node_unpackd(void *, int, SSW *,
#ifdef _ADH_GROUNDWATER
        SGW *,
#endif
        SGRID *, int, SCON *
#ifdef _SEDIMENT
    , SSED *,
    int,
    int,
    int
#endif
    );
void node_unpackd_surface(void *, int, SSW *
    );
void node_unpackd_bed(void *, int, SSW *
#ifdef _SEDIMENT
    , SSED *,
    int,
    int,
    int
#endif
    );
#else
void node_get_local(void);
void node_hash_add_entry(void);
void node_hash_index(void);
void node_hash_lookup(void);
void node_renumber_gn(void);
void node_pack_icnt(void);
void node_pack_dcnt(void);
void node_pack_dcnt_surface(void);
void node_packi(void);
void node_packd(void);
void node_packd_surface(void);
#endif

/* tools - simple math */
void init_adh_matrix(int nnodes, int max_nsys_sq, SPARSE_VECT *matrix, double *diagonal);
void check_matrix_diagonal_for_nonexistant_entries(int nnodes, int max_nsys_sq, int nsys, double *diagonal);
void printScreen_matrix(char * descript, double *diagonal, SPARSE_VECT * matrix, int nnode, int nsys_sq, int linenumber, char *filename);
void printScreen_matrixProfile(char * descript, Profile_Info *profile, int linenumber, char *filename);

void solv_init_int(int, int *);
void solv_init_dbl(int, double *);
int intcompare(const void *, const void *);
int dblcompare(const void *, const void *);
int is_double_small(double );
int compare_double(double , double );
void check_for_small_dbl(double *);
void Is_DoubleArray_Inf_or_NaN(double *, int, char *, int);
void Is_Double_Inf_or_NaN(double, char *, int);

int Is_Double_Inf_or_NaN_noExit(double X, char *filename, int linenumber);
int Is_DoubleArray_Inf_or_NaN_noExit(double *X, int arraybounds, char *filename, int linenumber);
int Is_vectorArray_Inf_or_NaN_noExit(SVECT *X, int arraybounds, char *filename, int linenumber);
int Is_Vector_Inf_or_NaN_noExit(SVECT X, char *filename, int linenumber);
double solv_l2_norm(int, double *);
double solv_l2_norm_scaled(int, double *, double);
int solv_isnan(double);
int solv_isinf(double);
#ifdef _MESSG
double solv_infty_norm(int, double *,MPI_Comm);
double solv_infty_norm(int, double *,MPI_Comm);
#else
double solv_infty_norm(int, double *);
double solv_infty_norm(int, double *);
#endif

void elem2d_to_node_int(SGRID *, int *, int *, int, int);
void elem2d_to_node_double(SGRID *, STR_VALUE *, double *, double *);
void elem3d_to_node_int(SGRID *, int *, int *, int, int);
void elem3d_to_node_double(SGRID *, STR_VALUE *, double *, double *);

void elem2dbed_to_node_int(SGRID *, int *, int *, int, int);
void elem2dbed_to_node_double(SGRID *, STR_VALUE *, double *, double *);

/* prototypes for blas_wrapper.c */
void solv_copy(int, double *, double *);
void solv_daxpy(int, double, double *, double *);
double solv_dot(int, double *, double *);
void solv_scal(int, double,double *);


/* file io */
FILE * io_fopen(const char *, const char *, const int);
void create_fileio(int *);
void close_fileio(int *);
int write_fileio_dbl(int *, double *, int *);
int read_fileio_dbl(int *, double *, int *);

/* Function Declarations if no header file already exists */
#ifdef _LAPACK
void dgetrf_(int *, int *, double *, int *, int *, int *);
void dgetrs_(char *, int *, int *, double *, int *, int *, double *, int *, int *, int);
void dgecon_(char *, int *, double *, int *, double *, double *, double *, int *, int *);
#endif

#ifdef _BLAS
void daxpy_(int *, double *, double *,int *, double *,  int *);
void dcopy_(int *, double *, int *, double *, int *);
double ddot_(int *, double *, int *, double *, int *);
void dscal_(int *, double *, double *, int *);
void dgemv_(char *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void dgemm_(char *, char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void dgglse_(int *, int *, int *, double *, int *, double *, int *, double *, double *, double *, double *, int *, int *);
#endif
