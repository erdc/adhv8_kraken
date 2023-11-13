/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Increments the SW 2D solutions after a Newton iterate.
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
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_sw2_inc(SSUPER_MODEL *sm, int isubModel) {

    SMODEL *mod = &(sm->submodel[isubModel]);
#ifdef _DEBUG
    assert(mod->flag.SW2_FLOW == ON);
#endif
    
    int i, i3;
    double factor = 1.;

#ifdef _PETSC
    PetscScalar const *values;
    PetscInt ierr;

    ierr = VecGetArrayRead(sm->sol,&(values));
    for (i = 0; i < mod->grid->my_nnodes; i++) {
        i3 = mod->fmap[i] * 3;
        mod->sw->d2->vel[i].x += factor * values[i3];
        mod->sw->d2->vel[i].y += factor * values[i3 + 1];
        mod->sw->d2->head[i]  += factor * values[i3 + 2];
    }
    ierr = VecRestoreArrayRead(sm->sol,&(values));
#else
    for (i = 0; i < mod->grid->nnodes; i++) {
        i3 = mod->fmap[i] * 3;
        mod->sw->d2->vel[i].x += factor * sm->sol[i3];
        mod->sw->d2->vel[i].y += factor * sm->sol[i3 + 1];
        mod->sw->d2->head[i]  += factor * sm->sol[i3 + 2];

       //printf("i increments: %d head: %30.20e u: %30.20e v: %30.20e\n",i,mod->sol[i3 + 2],mod->sol[i3 + 0],mod->sol[i3 + 1]);
    }
#endif
    ///exit(-1);
}
