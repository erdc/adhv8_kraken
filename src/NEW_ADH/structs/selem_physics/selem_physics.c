/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  sgrid.c This file collects methods of the SELEM\_PHYSICS structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and array of size nelem of physice FE routines
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elemPhys           (SELEM_PHYSICS **)  the struct double pointer
 * @param[in]  nelems           (int) the total number of elements (1D, 2D, or 3D)
 * @param[in]  nSubModels                (int*) the total number of submodels on each element
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void selem_physics_alloc_init(SELEM_PHYSICS **elemPhys,int nelems,int *nSubMods) { // triple pointer?
    int ie, ifun;
    
    elemPhys = (SELEM_PHYSICS **) tl_alloc(sizeof(SELEM_PHYSICS *), nelems);
    for (ie=0; ie<nelems; ie++) {
        assert(nSubMods[ie] > 0);
        elemPhys[ie] = (SELEM_PHYSICS *) tl_alloc(sizeof(SELEM_PHYSICS), nSubMods[ie]);
        for (ifun=0; ifun<nSubMods[ie]; ifun++) {
            elemPhys[ie][ifun].fe_inc = NULL;
            elemPhys[ie][ifun].fe_init = NULL;
            elemPhys[ie][ifun].fe_update = NULL;
            elemPhys[ie][ifun].fe_solve = NULL;
            elemPhys[ie][ifun].fe_resid = NULL;
            elemPhys[ie][ifun].fe_load = NULL;
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees an ELEM\_PHYSICS struct
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elemPhys           (SELEM_PHYSICS **)  the struct double pointer
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void selem_physics_free(SELEM_PHYSICS **elemPhys,int nelems,int *nSubMods) {
    int ie;
    
    for (ie=0; ie<nelems; ie++) {
        elemPhys[ie] = (SELEM_PHYSICS *) tl_free(sizeof(SELEM_PHYSICS), nSubMods[ie], elemPhys[ie]);
    }
    elemPhys = (SELEM_PHYSICS **) tl_free(sizeof(SELEM_PHYSICS *), nelems, elemPhys);
    
}
