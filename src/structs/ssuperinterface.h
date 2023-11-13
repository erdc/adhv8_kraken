/* This is the interface structure in AdH
 * to allow for 2D-3D flux coupling
 * between multiple supermodels
 */

#ifndef H_SSUPER_INTERFACE_
#define H_SSUPER_INTERFACE_


/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
    int surfnode[2][2];   // surfnode[i][j]: corresponds to submodel[i],node[j]
    int size[2];          // Number of faces in each coupled column. Can vary with each edge
    int *couplededges[2]; // Length of couplededges[][] to be set to = size[]
    double *coupled_normal_flux[2]; // Length to be set to = size[]
    double total_flux[2];
    /* For wvel calculations, if required */
    double total_master_flux;
    double total_slave_flux;

} SSUP_IFCE_BOUNDARY_LIST;


/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
    int sm_id[2];            // SuperModel numbering to start from 0,1,...
    int model_id[2];         // Model numbering to start from 0,1,...
    int bdrystrings[2];      // Edge or Face string ID's that are flux coupled
    int NumEdges;

    SSUP_IFCE_BOUNDARY_LIST *bdrylist; // Length to be set to NumEdges;

} SSUPER_INTERFACE;

void ssuperinterface_alloc_init(SSUPER_INTERFACE *, int, int*, int*, int *);
void ssuperinterface_boundarylist_alloc_init(SSUP_IFCE_BOUNDARY_LIST *, int *);
void ssuperinterface_free(SSUPER_INTERFACE *);
//void generate_2d2d_superinterface(SSUPER_INTERFACE *, SSUPER_MODEL *, int, int);
//void generate_2d3d_superinterface(SSUPER_INTERFACE *, SSUPER_MODEL *, int, int);

#endif
