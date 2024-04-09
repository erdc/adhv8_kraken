// An AdH DESIGN MODEL

#ifndef H_DESIGN_MODEL_
#define H_DESIGN_MODEL_

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Structure
typedef struct {

    int nSuperModels;          // total # of superModels
    SUPER_MODEL *superModel;   // Array of superModels to be flux coupled
    SGRID *grid;               // 1 designer grid or multiple supermodel grids. Decide later.
    
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
