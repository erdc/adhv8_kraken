// An AdH SuperModel
#ifndef H_SUPER_MODEL_
#define H_SUPER_MODEL_

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {

    double dt;
    double dt_prev;
    double inc_nonlin;
    double tol_nonlin;
    
    //FE_MATRIX *matrix;  // stores matrix

} SUPER_MODEL;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void super_model_alloc_init(SUPER_MODEL **smod, int nSuperModels);
void super_model_free(SUPER_MODEL *smod, int nSuperModels);
void super_model_read(SUPER_MODEL *smod, FILE *fp);
void super_model_printScreen(SUPER_MODEL *smod);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
