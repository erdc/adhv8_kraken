/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     MPI updates the 3D NS solutions after a Newton iterate.
 * \author    Corey Trahan, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    David Smith
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

void fe_ns3_update(SSUPER_MODEL *sm, int isubModel) {
#ifdef _MESSG
    SMODEL *mod = &(sm->submodel[isubModel]);
    comm_update_double(mod->ns->d3->prs, 1, mod->grid->smpi);
    comm_update_VECT(mod->ns->d3->vel, mod->grid->smpi);
#endif
}
