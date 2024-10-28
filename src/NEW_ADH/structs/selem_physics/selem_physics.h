// Stores physics on each element
#ifndef H_SELEM_PHYSICS_
#define H_SELEM_PHYSICS_

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
    int ndof;
    
} SELEM_PHYSICS;

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void selem_physics_alloc_init(SELEM_PHYSICS ***elemPhys,int nelems,int *nSubMods);
void selem_physics_free(SELEM_PHYSICS **elemPhys,int nelems,int *nSubMods);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
