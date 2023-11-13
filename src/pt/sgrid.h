#ifndef H_SGRID_
#define H_SGRID_

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {

    int ndim;               // grid dimension
    int isGridMixedElement;
    int haveTets;
    int havePrisms;
    int haveTris;
    int haveQuads;
    int type;               // Grid descriptor :: options :: UNSTRUCTURED, COLUMNAR

    int nnodes;             // number residential nodes
    int nelems3d;           // number of 3d elements
    int nelems2d;           // number of 2d elements
    int nelems1d;           // number of 1d elements
    int nedges;
    int nmat;               // the number of material

    SELEM_3D *elem3d;       // an array of 3d elements on this grid
    SELEM_2D *elem2d;       // an array of 2d elements on this grid
    //SELEM_1D *elem1d;     // an array of 1d elements on this grid
    SVECT *node;            // an array of grid nodes at some time t
    SVECT *node_t0;         // an array of grid nodes at the initial model time
    SVECT *node_sur;       // an array of surface node ids for all 3d nodes on a quasi-structured grid

    int *node_flag;         // a flag to determind if node is on boundary
    
    double xmin, xmax, ymin, ymax, xL, yL, xLinv, yLinv; // grid bounds. Also, 1/xL, 1/yL, which are used for normalization, if requested.
    double mesh_volume;                        // grid volume

    int *nc_nelems;
    int **nc_elems;

    int nd_on_TriEdge[3][2];
    int nd_on_TetEdge[6][2];

    int vertical_dpl_flag;  // a flag to determine if the grid displaces vertically
    
#ifdef _ADH_HDF5
    HDF5 hdf5;
#endif
    
} SGRID;

/***********************************************************/
/* struct methods ---------------------------------------- */

void sgrid_read_adh(SGRID **grid, char *file_base);
void sgrid_free(SGRID *grid);
void sgrid_printScreen(SGRID *grid);
void sgrid_create2d(SGRID *grid, double xmin, double xmax, double ymin, double ymax, int npx, int npy, double theta, double dz, double za, double zb, double zc);
void sgrid_create3d(SGRID *grid, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int npx, int npy, int npz, double theta);
void sgrid_write_binary(SGRID *grid, FILE *fptr);
void sgrid_read_binary(SGRID *grid, FILE *fptr);
void sgrid_dpl_update(SGRID *grid, double *dpl, double tL_dpl, double tR_dpl, double *dpl_tL, double *dpl_tR, double t);

/***********************************************************/
/* normalization methods --------------------------------- */

void svect_normalize_inplace(bool normalize, SGRID * grid, SVECT * v);
void svect_denormalize_inplace(bool normalize, SGRID * grid, SVECT * v);
SVECT svect_normalize(bool normalize, SGRID * grid, SVECT in);
SVECT svect_denormalize(bool normalize, SGRID * grid, SVECT in);
void sgrid_normalize(bool normalize, SGRID * g);


/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif



