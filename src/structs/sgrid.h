#ifndef H_SGRID_
#define H_SGRID_

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
    
    int ndim;               /* grid dimension */
    int isGridMixedElement; /* for mixed element grids, such as 3D prisms */
    int haveTets;
    int havePrisms;
    int haveTris;
    int haveQuads;
    int type;               // Grid descriptor :: options :: UNSTRUCTURED, COLUMNAR
    
    int **nd_on_TriEdge;
    int **nd_on_QuadEdge;
    int **nd_on_PrismEdge;
    int **nd_on_TetEdge;
    
    int macro_nnodes;       /* total number of nodes in the global mesh */
    int macro_nelems1d;       /* total number of elements in the global mesh */
    int macro_nelems2d;       /* total number of elements in the global mesh */
    int macro_nelems3d;       /* total number of elements in the global mesh */
    int orig_macro_nnodes;       /* total number of original nodes in the global mesh */
    int orig_macro_nelems1d;       /* total number of original elements in the global mesh */
    int orig_macro_nelems2d;       /* total number of original elements in the global mesh */
    int orig_macro_nelems3d;       /* total number of original elements in the global mesh */
    int macro_nelems2d_bed;   /* total number of 2d bed/surface elements across all PEs on the global mesh */
    int initial_nnodes;     /* total number of initial (before adaption) number of ghost + residential nodes */
    int initial_nelems;     /* total number of initial (before adaption) number of ghost + residential 2d or 3d elements */
    int initial_nnodes_bed; /* total number of initial (before adaption) number of ghost + residential bed nnodes */
    int initial_nelems_bed; /* total number of initial (before adaption) number of ghost + residential bed elements */
    int orig_initial_nnodes;/* when writing adapted grids, which temporarily change initial_nnodes, this stores them */
    int orig_initial_nelems;/* when writing adapted grids, which temporarily change initial_nelems, this stores them */
    int nnodes_prev;         /* partiioning and build_columns needs this */
    int nnodes;             /* number of ghost + residential nodes */
    int nelems3d;           /* number of ghost + resdiential 3d elements */
    int nelems2d;           /* number of ghost + residential 2d elements */
    int nelems1d;           /* number of ghost + residential 1d elements */
    int nelems1d2d;
    int nedges;
    int nmat;               /* the number of materials */
    int interface;

    /* quadrature variables */
    int quadrature_order;   /* max order of quadrature integration */
    int nqp_2d;             /* the total number of 2d quadrature points */
    int nqp_3d;             /* the total number of 3d quadrature points */
    SQUAD *quad_tri;        /* stores all 2d quadrature data, including points, weights, basis values, etc. */
    SQUAD *quad_rect;
    SQUAD *quad_tet;        /* stores all 2d quadrature data, including points, weights, basis values, etc. */
    SQUAD *quad_prism;
    
    /* store last nnodes and elements (before adaption/partitioning) for reallocation */
    int nnodes_old;
    int nelems2d_old;
    int nelems3d_old;
    
    int nnodes_on_elem1d;  /* the number of local nodes on 1d element */
    int nnodes_on_elem2d;  /* the number of local nodes on 2d element */
    int nnodes_on_elem3d;  /* the number of local nodes on 3d element */
    int nnodes_on_elem1d_quad;  /* the number of local nodes for a quadratic 1d element */
    int nnodes_on_elem2d_quad;  /* the number of local nodes for a quadratic 2d element */
    int nnodes_on_elem3d_quad;  /* the number of local nodes for a quadratic 3d element */
    
    int nnodes_bed_prev;
    int nnodes_sur_prev;
    int nnodes_sur_orig;
    int nnodes_matrix;      /* the last number of nodes used to calc dof of the fe matrix */
    int max_nnodes;         /* the maximum number of nodes allocated (nnode + nalloc_inc) */
    int max_nelems1d;       /* the maximum number of 1d elements allocated (nelems1d + nalloc_inc) */
    int max_nelems2d;       /* the maximum number of 2d elements allocated (nelems2d + nalloc_inc) */
    int max_nelems3d;       /* the maximum number of 3d elements allocated (nelems3d + nalloc_inc) */
    int my_nnodes;          /* the number of residential nodes*/
    
    /* grid elements */
    SELEM_3D *elem3d;        /* an array of 3d elements on this grid */
    SELEM_2D *elem2d;        /* an array of 2d elements on this grid */
    SELEM_1D *elem1d;        /* an array of 1d elements on this grid */
    
    /* maximum elemental solution error */
    double *elem_error;
    
    /* mesh volume */
    double mesh_volume;
    
    /* wet dry flag */
    int *wd_flag;         /* flag if 2d element is wet or dry */
    
    /* Eddy viscosity and diffusivity */
    double *hyd_eddy; //GSAVANT
    double *trn_diff; //GSAVANT
    
    /* grid nodes */
    SNODE *node;
    
    /* grid bounds */
    double x_min, x_max, y_min, y_max, xL, yL;
    
    /* 3d to 2d projected node ID mappings */
    int *nodeID_3d_to_2d_sur;        /* maps 3d grid node IDs to overlaying 2d grid node IDs if on surface */
    int *nodeID_3d_to_2d_bed;        /* maps 3d grid node IDs to overlaying 2d grid node IDs if on bed */
    int *nodeID_2d_to_3d_sur;        /* maps 2d grid node IDs to overlaying 3d surace grid node IDs */
    int *nodeID_2d_to_3d_bed;        /* maps 2d grid node IDs to overlaying 3d bed grid node IDs */
    int *old_global_surf;
    int *old_global_bed;
    
    /* processor id counts this is the global id # for starting and ending node numbers of a given processor */
    int *start_node_ID; /* sized to number of processors */
    int *end_node_ID; /* sized to number of procesors*/
    
    /* columns ***************************************************/
    int nnodes_sur;
    int nnodes_bed;
    int my_nnodes_sur;
    int my_nnodes_bed;
    int macro_nnodes_sur;
    int macro_nnodes_bed;
    int orig_macro_nnodes_sur;
    int orig_macro_nnodes_bed;
    int max_nnodes_sur;
    int max_nnodes_bed;
    int old_max_nnodes_sur;
    int old_max_nnodes_bed;
    int num_midpts;
    int num_vert_segments;
    int nelems2d_sur;      /* = ncolumns, keep for sanity checking */
    int nelems2d_bed;      /* = ncolumns, keep for sanity checking */
    int nelems2d_sidewall;
    
    int *elem3d_sur;       /* 3d surface element id for each column */
    int *elem2d_sur;       /* 2d surface element id for each column */
    int *elem3d_bed;
    int *elem2d_bed;
    int *elem3d_sidewall;
    int *elem2d_sidewall;
    
    int isize_sidewall_elems;
    int isize_sur_elems;
    int isize_bed_elems;
    int isize3, isize4, isize5;
    int scale_factor;
    int shift_factor;
    
    /* linked lists in columns */
    int ncolumns;
    int isize_ncolumns;
    ID_LIST_ITEM **column_list;         /* linked list of 3d elements that lie in columns */
    ID_LIST_ITEM **column_list2d;       /* linked list of 2d faces on the columns */
    ID_LIST_ITEM **vertical_list;       /* Linked list of nodes that lie on vertical segment */
    ID_LIST_ITEM **sidewall_list;       /* Linked list of 2D sidewall elements */
    MIDPT_LIST_ITEM **midpt_list;
    
    /* hash tables */
    int hash_size;
    int is_allocated_column_hash;
    CENT_LIST_ITEM **column_hash;       /* hash table for columns of 3d elements */
    int is_allocated_column_hash2d;
    CENT_LIST_ITEM **column_hash2d;     /* hash table for columns of 2d elements */
    int is_allocated_midpt_hash;
    CENT_LIST_ITEM **midpt_hash;        /* */
    int is_allocated_vertical_hash;
    CENT_LIST_ITEM **vertical_hash;     /* hash table for xy-coords of points; used to build list of nodes that lie on a vertical line */
    
    SMPI *smpi;      /* MPI struct */
    SMPI *part_smpi; /* MPI struct for processors participating in re-partioning */
    int *part_map;   /* maps partition processor ids to main grid MPI ids */
    SMPI *supersmpi; /* MPI struct of the supermodel - only an alias! */ /* gkc for 2D-3D coupling / multimodel runs */
    SMPI *supersmpi_wvel; /* MPI struct of the supermodel - only an alias! */ /* gkc for 2D-3D coupling / multimodel runs */ 

    
} SGRID;

/***********************************************************/
/* struct methods ---------------------------------------- */

void sgrid_free(SGRID *);
void sgrid_printScreen(SGRID *grid, char *filename);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif


