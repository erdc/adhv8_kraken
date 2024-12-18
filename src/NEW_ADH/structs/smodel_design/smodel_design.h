// An AdH DESIGN MODEL
#ifndef H_SMODEL_DESIGN_
#define H_SMODEL_DESIGN_

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Structure
typedef struct {

    //entire design model has dt, supermodels are multiples of this?
    double dt, old_dt, dt_err;
    double dt_prev;
    double inc_nonlin;
    double tol_nonlin;
    double t_init;
    double t_prev;
    double t_final;
    int t_adpt_flag;

    int nseries;              // the number of series in this model 
    int itrns;

    int nSuperModels;          // total # of superModels
    SMODEL_SUPER *superModel;   // Array of superModels to be flux coupled

    //int *superType; //[nSuperModels], for each supermodel, which is Mono and which are single models
    int nMono; //number of monolothic-coupled supermodels 
    int nSimple; //number of single physics supermodels (nMono+nSimple = nSuperModels)
    int nUnique; //number of unique supermodels, leading to unique sparsity structure
                 //a Unique supermodel is one that is either monolithically coupled or
                 // unique sparsity
    int *unique_id; //integer array of size nUnique, saves index of first Unique super model for each sparsity pattern

    
    // grid and physics on the grid
    SGRID *grid;               // 1 designer grid or multiple supermodel grids. Decide later.

    
    //array of linear systems, some supermodels may share the same
    int *lin_sys_id; //array of [nSuperModels] that gives the index of the linear system it belongs to
    SLIN_SYS *lin_sys; //array of [nUnique] systems

    //Mark added local to process, this should be same as local_range[1]-local_range[0]
    //are these redundant in any way?
    int *my_ndofs; //array of [nUnique]
    int *my_ndofs_old; //array of [nUnique]
    int *ndofs; // array of [nUnique] local number of degrees of freedom
    int *ndofs_old; //array of [nUnique] local numer of solution variables the processor is in charge of
    int *macro_ndofs; //array of [nUnique]
    int *macro_ndofs_old; //array of [nUnique]



    // int *nFluxInterfaces; // total # of flux interfaces between supermodesl
    // SINTERFACE_FLUX  will hold all data for interfacing between supermodels
    
} SMODEL_DESIGN;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
//void design_model_alloc_init(SMODEL_DESIGN **dmod, int nSuperModels);
void smodel_design_alloc(SMODEL_DESIGN *dmod, int nSuperModels, int nMono, int nSimple,
    int nUnique);
void smodel_design_free(SMODEL_DESIGN *dm);
void smodel_design_printScreen(SMODEL_DESIGN *dmod);
void smodel_design_no_read_simple(SMODEL_DESIGN *dm, double dt_in, double t_init, double t_final,
    int nphysics_mat_1d, int nphysics_mat_2d, int nphysics_mat_3d, char elemVarCode[4] ,
    SGRID *grid);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
