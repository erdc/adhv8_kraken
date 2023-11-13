/* This is the interface structure in AdH
 * to allow for 2D-3D algebraic coupling
 * between multiple submodels within a 
 * supermodel.
 */

#ifndef H_SINTERFACE_
#define H_SINTERFACE_


/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
    int surfnode[2];      // Surface node numbers of the coupled node columns
    int size[2];          // Number of nodes in each coupled columns. Can vary with each node!
    int *couplednodes[2]; // Length of couplednodes[][] to be set to = size[]

    /* For wvel calculations, if required */
    int wvel_residual_length;
    double *save_wvel_residual;

} SINTERFACE_NODELIST;


/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
    int model_id[2];    // Model numbering to start from 0,1,...
    int NumNodeColumns;

    SINTERFACE_NODELIST *nodelist; // Length to be set to NumNodeColumns;

} SINTERFACE;

void sinterface_alloc_init(SINTERFACE *, int, int, int);
void sinterface_free(SINTERFACE *);

#endif
