// Stores physics on each element
#ifndef H_SMODEL_
#define H_SMODEL_

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {
    
    void *fe_inc;
    void *fe_init;
    void *fe_update;
    void *fe_solve;
    int (*fe_resid)(); // this could be body or boundary, depending on what element
    //void *fe_boundary_resid;
    //void *fe_body_resid;
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
void smodel_free(SMODEL **elemPhys,int nelems,int *nSubMods);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
