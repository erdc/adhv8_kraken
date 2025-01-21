// Stores physics on each element
#ifndef H_SMODEL_
#define H_SMODEL_


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {
    
    void *fe_inc;
    //void *fe_init;
    int fe_init;
    void *fe_update;
    void *fe_solve;
    // this could be body or boundary, depending on what element
    //void *fe_boundary_resid;
    //void *fe_body_resid;
    //change to just integer that is a code
    int fe_resid;
    //int (*fe_resid)(); 
    void *fe_load;
    
    int nvar;
    //also needs physics vars
    int* physics_vars;//int physics_vars[ndof];
    
} SMODEL;

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
//void smodel_alloc_init_array(SMODEL ***elemPhys,int nelems,int *nSubMods);
void smodel_alloc_init_array(SMODEL **elemPhys, int nSubMods, int *nSubMod_nvar);
void smodel_alloc_init(SMODEL *physics,int nvar);
void smodel_free_array(SMODEL *Models, int nModels);
void smodel_free(SMODEL *model);
//int smodel_fe_resid()

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
