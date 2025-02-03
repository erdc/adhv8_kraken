#ifndef H_SGRID_
#define H_SGRID_

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++//

typedef struct {

    // winds
    SMETEOR_FILE *wind_file;
    SWIND *winds;

    // waves
    SMETEOR_FILE *wave_file;
    SWAVE *waves;

    // rain
    SMETEOR_FILE *rain_file;
    double *rain;
    
} NODAL_ATTRIBUTES;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++//

typedef struct {
    char filename[100];  // grid filename
    char type[15];       // grid descriptor :: options :: UNSTRUCTURED, COLUMNAR
    int ndim;            // highest grid dimension
    bool columnar;       // does grid have columnar parts
    bool haveSegs;       // does grid have line segments
    bool haveTets;       // does grid have tetrahedrons
    bool havePrisms;     // does grid have trianguler prisms
    bool haveTris;       // does grid have triangles
    bool haveQuads;      // does grid have quadrilaterals
    
    // Local totals (used for serial and HPC)
    int nnodes;          // number of ghost + residential nodes
    int nnodes_sur;      // number of ghost + residential surface nodes
    int nnodes_bed;      // number of ghost + residential bed nodes
    int nelems3d;        // number of ghost + residential 3d elements
    int nelems2d;        // number of ghost + residential 2d elements
    int nelems1d;        // number of ghost + residential 1d elements
    int nTets;
    int nPrisms;
    int nTris;
    int nQuads;


    SNODE *node;         // number of ghost + residential nodes
    SELEM_3D *elem3d;    // an array of 3d elements on this grid
    SELEM_2D *elem2d;    // an array of 2d elements on this grid
    SELEM_1D *elem1d;    // an array of 1d elements on this grid
    
    // HPC totals
    int macro_nnodes;          // total # of nodes in the global mesh
    int macro_nnodes_sur;      // total # of surface nodes
    int macro_nnodes_bed;      // total # of bed nodes
    int macro_nelems3d;        // total # of 3d elements in the global mesh
    int macro_nelems2d;        // total # of 2d elements in the global mesh
    int macro_nelems1d;        // total # of 1d elements in the global mesh
    int macro_nelems2d_bed;    // total # of 2d bed/surface elements across all PEs on the global mesh
    int macro_nTets;           // total # of tetrahedron elements in global mesh
    int macro_nPrisms;         // total # of prism elements in global mesh
    int macro_nQuads;          // total # of quadrilateral elements in global mesh
    int macro_nTris;           // total # of triangle elements in global mesh
    int orig_macro_nnodes;     // total number of original nodes in the global mesh
    int orig_macro_nnodes_sur; // total number of original surface nodes in the global mesh
    int orig_macro_nnodes_bed; // total number of original bed nodes in the global mesh
    int orig_macro_nelems1d;   // total number of original elements in the global mesh
    int orig_macro_nelems2d;   // total number of original elements in the global mesh
    int orig_macro_nelems3d;   // total number of original elements in the global mesh
    int my_nnodes;          // total # of residential only nodes in the global mesh
    int my_nnodes_sur;      // total # of residential only surface nodes in the global mesh
    int my_nnodes_bed;      // total # of residential only bed nodes in the global mesh
    int my_nelems3d;        // total # of residential only 3d elements in the global mesh
    int my_nelems2d;        // total # of residential only 2d elements in the global mesh
    int my_nelems1d;        // total # of residential only 1d elements in the global mesh
    int my_nTets;           // total # of residential only tetrahedron elements in global mesh
    int my_nPrisms;         // total # of residential only prism elements in global mesh
    int my_nQuads;          // total # of residential only quadrilateral elements in global mesh
    int my_nTris;           // total # of residential only triangle elements in global mesh

    // Refinement-oriented counts
    int max_nnodes;          // number of ghost + residential nodes + nalloc_inc
    int max_nnodes_sur;      // number of ghost + residential surface nodes + nalloc_inc
    int max_nnodes_bed;      // number of ghost + residential bed nodes + nalloc_inc
    int max_nelems3d;        // number of ghost + resdiential 3d elements + nalloc_inc
    int max_nelems2d;        // number of ghost + residential 2d elements + nalloc_inc
    int max_nelems1d;        // number of ghost + residential 1d elements + nalloc_inc
    int nnodes_old;          // store last nnodes (before adaption/partitioning) for reallocation
    int nnodes_sur_old;      // store last surface nnodes (before adaption/partitioning) for reallocation
    int nnodes_bed_old;      // store last bed nnodes (before adaption/partitioning) for reallocation
    int nelems3d_old;        // store last nelem3d's (before adaption/partitioning) for reallocation
    int nelems2d_old;        // store last nelem2d's (before adaption/partitioning) for reallocation
    int nelems1d_old;        // store last nelem1d's (before adaption/partitioning) for reallocation
                             
    int initial_nnodes;      // total # of initial (before adaption) number of ghost + residential nodes
    int initial_nelems;      // total # of initial (before adaption) number of ghost + residential 2d/ 3d elements
    int initial_nnodes_bed;  // total # of initial (before adaption) number of ghost + residential bed nnodes
    int initial_nelems_bed;  // total # of initial (before adaption) number of ghost + residential bed elements
    int orig_initial_nnodes; // when writing adapted grids, which temporarily change initial_nnodes, this stores them
    int orig_initial_nelems; // when writing adapted grids, which temporarily change initial_nelems, this stores them

    // Mesh extents and length, area, volume
    double mesh_length, my_mesh_length;
    double mesh_area, my_mesh_area;
    double mesh_area_surface, my_mesh_area_surface;
    double mesh_area_bed, my_mesh_area_bed;
    double mesh_area_sidewalls, my_mesh_area_sidewalls;
    double mesh_volume, my_mesh_volume;
    double x_min, x_max, y_min, y_max, z_min, z_max; // grid bounds
    double my_x_min, my_x_max, my_y_min, my_y_max, my_z_min, my_z_max; // grid bounds
    
    // Quadrature variables
    int quadrature_order;   // max order of quadrature integration
    int nqp_1d;             // the total number of 1d quadrature points
    int nqp_2d;             // the total number of 2d quadrature points
    int nqp_3d;             // the total number of 3d quadrature points
    SQUAD *quad_seg;        // stores all 1d quadrature data, including points, weights, basis values, etc.
    SQUAD *quad_tri;        // stores all 2d quadrature data, including points, weights, basis values, etc.
    SQUAD *quad_rect;       // stores all 2d quadrature data, including points, weights, basis values, etc.
    SQUAD *quad_tet;        // stores all 3d quadrature data, including points, weights, basis values, etc.
    SQUAD *quad_prism;      // stores all 3d quadrature data, including points, weights, basis values, etc.
    
    // 3d to 2d projected node ID mappings 
    int *nodeID_3d_to_2d_sur;        // maps 3d grid node IDs to overlaying 2d grid node IDs if on surface
    int *nodeID_3d_to_2d_bed;        // maps 3d grid node IDs to overlaying 2d grid node IDs if on bed
    int *nodeID_2d_to_3d_sur;        // maps 2d grid node IDs to overlaying 3d surace grid node IDs
    int *nodeID_2d_to_3d_bed;        // maps 2d grid node IDs to overlaying 3d bed grid node IDs

    // MPI info for partitioning grid
    SMPI *smpi;
    SMPI *part_smpi; // MPI struct for processors participating in re-partioning
    int *part_map;   // maps partition processor ids to main grid MPI ids
    
    // Nodal Attributes
    NODAL_ATTRIBUTES nodal_attribute;

    //inverse permutation for output (serial only)
    int *inv_per_node;

    //Mark added for now, need to ask
    int **nd_on_QuadEdge;
    
} SGRID;

/***********************************************************/
/* struct methods ---------------------------------------- */

void sgrid_read_node(SGRID *g, char *line, int *start_node_id, int *end_node_id, NODE_LIST_ITEM *ghost_nodes);
int sgrid_read_elem(SGRID *g, char *line, int *start_node_id, int *end_node_id, int *num_ghosts, NODE_LIST_ITEM **ghost_nodes, int nnodes_on_elem, int elem_dim, bool JUST_COUNT);
void sgrid_free(SGRID *);
void sgrid_printScreen(SGRID *);
#ifdef _MESSG
void sgrid_read(SGRID **pgrid, char *filename, MPI_Comm model_comm);
#else
void sgrid_read(SGRID **pgrid, char *filename);
#endif
void sgrid_write_hdf5(SGRID *g);
void init_hdf5_file(SGRID *g);
void sgrid_write_xdmf(SGRID *g);
void sgrid_write_nodal_pe(SGRID *g);
void sgrid_write_elemental_pe(SGRID *g);
void sgrid_write_xdmf_nodal_pe(SGRID *g);
void sgrid_write_xdmf_elemental_pe(SGRID *g);
void sgrid_read_nodal_attribute(SGRID *g);
int sgrid_reorder(SGRID *grid, int option);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
