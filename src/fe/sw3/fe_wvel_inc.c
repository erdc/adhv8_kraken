/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Increments the SW 3D WVEL solutions after a Newton iterate.
 * \author    Charlie Berger, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Gary Brown, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout]  mod (SMODEL *) pointer to the model struct
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_wvel_inc(SSUPER_MODEL *sm, int imod) {
    int i;
    SMODEL *mod = &(sm->submodel[imod]);
#ifdef _PETSC
    PetscScalar const *values;
    PetscInt ierr;
    
    ierr = VecGetArrayRead(sm->sol, &(values));

    for(i = 0; i < mod->grid->my_nnodes; i++) {
        mod->sw->d3->vel[i].z += values[3*mod->fmap_wvel[i]];
    }

    ierr = VecRestoreArrayRead(sm->sol,&(values));
#else
    for(i = 0; i < mod->grid->nnodes; i++) {
        mod->sw->d3->vel[i].z += sm->sol[mod->fmap_wvel[i]];
    }
#endif
}
