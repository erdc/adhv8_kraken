/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Increments the 2D and 3D transport solutions after a Newton iterate.
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Gary Brown, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout]  mod (SMODEL *) pointer to the model struct
 *
 * \note CJT \:: transport here can be either sediment suspended load our a general constituent
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_transport_inc(SSUPER_MODEL *sm, int imod) {

    SMODEL *mod = &(sm->submodel[imod]);
    int i=0, index=0;
    
    double *c;
    if (mod->is_sediment_running) {
#ifdef _SEDIMENT
        c = mod->sed->susload[mod->ised].c;
#endif
    } else {
        c = mod->con[mod->itrns].concentration;
    }

#ifdef _PETSC
    // Convert PETSc sol vector to double array
    PetscScalar const *values;
    PetscInt ierr;

    ierr = VecGetArrayRead(sm->sol,&(values));
    for (i = 0; i < mod->grid->my_nnodes; i++) {
        // Multiply the index by max_nsys because
        // the matrix and vectors are allocated
        // based on max_nsys which may be greater
        // than 1 when solving multiple types of
        // equations within on model (e.g. sw2 and
        // transport).
        // We need to skip over the extra values
        // to get the ones we care about here.
        c[i] += values[mod->fmap[i]*mod->max_nsys];
    }
    ierr = VecRestoreArrayRead(sm->sol,&(values));
#else
    for (i = 0; i < mod->grid->nnodes; i++) {
        c[i] += sm->sol[mod->fmap[i]];
    }
#endif
}
