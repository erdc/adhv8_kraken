#ifndef H_SGRID_
#define H_SGRID_

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
    char filename[100];      // grid filename
    char type[15];          // grid descriptor :: options :: UNSTRUCTURED, COLUMNAR
    int nnodes;          // number of ghost + residential nodes
    int nelems3d;        // number of ghost + resdiential 3d elements
    int nelems2d;        // number of ghost + residential 2d elements
    int nelems1d;        // number of ghost + residential 1d elements
    SNODE *node;         // number of ghost + residential nodes
    SELEM_3D *elem3d;    // an array of 3d elements on this grid
    SELEM_2D *elem2d;    // an array of 2d elements on this grid
    SELEM_1D *elem1d;    // an array of 1d elements on this grid
    double total_volume; // grid volume
    double x_min, x_max, y_min, y_max, z_min, z_max; // grid bounds
    
    // refinement-oriented counts
    int max_nnodes;          // number of ghost + residential nodes + nalloc_inc
    int max_nelems3d;        // number of ghost + resdiential 3d elements + nalloc_inc
    int max_nelems2d;        // number of ghost + residential 2d elements + nalloc_inc
    int max_nelems1d;        // number of ghost + residential 1d elements + nalloc_inc
    
} SGRID;

/***********************************************************/
/* struct methods ---------------------------------------- */

void sgrid_alloc_init(SGRID **pgrid, char *filename, char *type, int nnodes, int nelems1d, int nelems2d, int nelems3d);
void sgrid_free(SGRID *);
void sgrid_read(SGRID **pgrid, FILE *fp, char *filename);
void sgrid_printScreen(SGRID *);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
