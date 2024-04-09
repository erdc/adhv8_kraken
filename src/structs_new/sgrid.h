#ifndef H_SGRID_
#define H_SGRID_

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
    
    char *filename;      // grid filename
    char *type;          // grid descriptor :: options :: UNSTRUCTURED, COLUMNAR
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
    
    // physics on the grid
    int *nSubMods1d;                // [nelems1d] the total number of physics modules on each 1D element
    int *nSubMods2d;                // [nelems2d] the total number of physics modules on each 2D element
    int *nSubMods3d;                // [nelems3d] the total number of physics modules on each 3D element
    ELEM_PHYSICS **elem1d_physics;  // [nelems1d][nsubmods_1d] the fe routines for each type of physics on each 1D element
    ELEM_PHYSICS **elem2d_physics;  // [nelems2d][nsubmods_2d] the fe routines for each type of physics on each 2D element
    ELEM_PHYSICS **elem3d_physics;  // [nelems3d][nsubmods_3d] the fe routines for each type of physics on each 3D element
    
} SGRID;

/***********************************************************/
/* struct methods ---------------------------------------- */

void sgrid_alloc_init(SGRID **pgrid, char *filename, char *type, int nnodes, int nelems1d, int nelems2d, int nelems3d, int *nsubMods1d, int *nsubMods2d, int *nsubMods3d);
void sgrid_free(SGRID *);
void sgrid_read(SGRID **pgrid, FILE *fp, char *filename);
void sgrid_printScreen(SGRID *);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
