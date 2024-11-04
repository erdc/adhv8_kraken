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
void selem_physics_alloc_init(SELEM_PHYSICS ***elemPhys,int nelems,int *nSubMods) { // triple pointer?
    int ie, ifun;
    
    assert(nelems>=0);
    if (nelems == 0){
        elemPhys = NULL;
        return;
    }else{
        (*elemPhys) = (SELEM_PHYSICS **) tl_alloc(sizeof(SELEM_PHYSICS *), nelems);

        SELEM_PHYSICS *physics; // for alias
        for (ie=0; ie<nelems; ie++) {
            assert(nSubMods[ie] > 0);
            (*elemPhys[ie]) = (SELEM_PHYSICS *) tl_alloc(sizeof(SELEM_PHYSICS), nSubMods[ie]);
            for (ifun=0; ifun<nSubMods[ie]; ifun++) {
                physics = &(*elemPhys)[ie][ifun]; // alias
                physics->fe_inc = NULL;
                physics->fe_init = NULL;
                physics->fe_update = NULL;
                physics->fe_solve = NULL;
                physics->fe_resid = NULL;
                physics->fe_load = NULL;
            }
        }
        return;
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
