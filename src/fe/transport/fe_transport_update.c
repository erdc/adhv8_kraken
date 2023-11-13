/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     MPI updates the 2D and 3D transport solutions after a Newton iterate.
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

void fe_transport_update(SSUPER_MODEL *sm, int imod) {

    SMODEL *mod = &(sm->submodel[imod]);
    
#ifdef _MESSG
    if (mod->is_sediment_running) {
#ifdef _SEDIMENT
        comm_update_double(mod->sed->susload[mod->ised].c, 1, mod->grid->smpi);
#endif
    } else {
        comm_update_double(mod->con[mod->itrns].concentration, 1, mod->grid->smpi);
    }
#endif
}
