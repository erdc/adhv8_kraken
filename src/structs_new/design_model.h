// An AdH DESIGN MODEL

#ifndef H_DESIGN_MODEL_
#define H_DESIGN_MODEL_

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Structure
typedef struct {

    int nSuperModels;          // total # of superModels
    SUPER_MODEL *superModel;   // Array of superModels to be flux coupled
    
    // grid and physics on the grid
    SGRID *grid;               // 1 designer grid or multiple supermodel grids. Decide later.
    int *nSubMods1d;                // [nelems1d] the total number of physics modules on each 1D element
    int *nSubMods2d;                // [nelems2d] the total number of physics modules on each 2D element
    int *nSubMods3d;                // [nelems3d] the total number of physics modules on each 3D element
    ELEM_PHYSICS **elem1d_physics;  // [nelems1d][nsubmods_1d] the fe routines for each type of physics on each 1D element
    ELEM_PHYSICS **elem2d_physics;  // [nelems2d][nsubmods_2d] the fe routines for each type of physics on each 2D element
    ELEM_PHYSICS **elem3d_physics;  // [nelems3d][nsubmods_3d] the fe routines for each type of physics on each 3D element
    
    // int *nFluxInterfaces; // total # of flux interfaces between supermodesl
    // SINTERFACE_FLUX  will hold all data for interfacing between supermodels
    
} DESIGN_MODEL;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void design_model_alloc_init(DESIGN_MODEL **dmod, int nSuperModels);
void design_model_free(DESIGN_MODEL *dmod);
void design_model_printScreen(DESIGN_MODEL *dmod);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
