// An AdH DESIGN MODEL

#ifndef H_SMODEL_DESIGN_
#define H_SMODEL_DESIGN_

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Structure
typedef struct {

    int nSuperModels;          // total # of superModels
    SMODEL_SUPER *superModel;   // Array of superModels to be flux coupled
    
    // grid and physics on the grid
    SGRID *grid;               // 1 designer grid or multiple supermodel grids. Decide later.

    
    // int *nFluxInterfaces; // total # of flux interfaces between supermodesl
    // SINTERFACE_FLUX  will hold all data for interfacing between supermodels
    
} SMODEL_DESIGN;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void design_model_alloc_init(SMODEL_DESIGN **dmod, int nSuperModels);
void design_model_free(SMODEL_DESIGN *dmod);
void design_model_printScreen(SMODEL_DESIGN *dmod);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
