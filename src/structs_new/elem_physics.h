// Stores physics on each element
#ifndef H_ELEM_PHYSICS_
#define H_ELEM_PHYSICS_

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {
    
    void *fe_inc;
    //void *fe_init;
    void *fe_update;
    void *fe_solve;
    void *fe_resid; // this could be body or boundary, depending on what element
    //void *fe_boundary_resid;
    //void *fe_body_resid;
    void *fe_load;
    
    int ndof;
} ELEM_PHYSICS;

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void elem_physics_alloc_init(ELEM_PHYSICS **elemPhys,int nelems,int *nSubMods);
void elem_physics_free(ELEM_PHYSICS **elemPhys,int nelems,int *nSubMods);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
