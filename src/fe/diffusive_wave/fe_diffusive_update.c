/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     MPI updates the sw2d solutions after a Newton iterate.
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

void fe_diffusive_update(SSUPER_MODEL *sm, int imod) {

#ifdef _MESSG
    SMODEL *mod = &(sm->submodel[imod]);
    comm_update_double(mod->sw->d2->head, 1, mod->grid->smpi);
    comm_update_double(mod->sw->d2->vel,  2, mod->grid->smpi);
#endif
}
