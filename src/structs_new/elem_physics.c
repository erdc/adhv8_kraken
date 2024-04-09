/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  sgrid.c This file collects methods of the SELEM\_PHYSICS structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//#include "global_header.h"
#include "local_header.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and array of size nelem of physice FE routines
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elemPhys           (ELEM_PHYSICS **)  the struct double pointer
 * @param[in]  nelems           (int) the total number of elements (1D, 2D, or 3D)
 * @param[in]  nSubModels                (int*) the total number of submodels on each element
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void elem_physics_alloc_init(ELEM_PHYSICS **elemPhys,int nelems,int *nSubMods) { // triple pointer?
    int ie, ifun;
    
    elemPhys = (ELEM_PHYSICS **) tl_alloc(sizeof(ELEM_PHYSICS *), nelems);
    for (ie=0; ie<nelems; ie++) {
        assert(nSubMods[ie] > 0);
        elemPhys[ie] = (ELEM_PHYSICS *) tl_alloc(sizeof(ELEM_PHYSICS), nSubMods[ie]);
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
 * @param[inout] elemPhys           (ELEM_PHYSICS **)  the struct double pointer
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void elem_physics_free(ELEM_PHYSICS **elemPhys,int nelems,int *nSubMods) {
    int ie;
    
    for (ie=0; ie<nelems; ie++) {
        elemPhys[ie] = (ELEM_PHYSICS *) tl_free(sizeof(ELEM_PHYSICS), nSubMods[ie], elemPhys[ie]);
    }
    elemPhys = (ELEM_PHYSICS **) tl_free(sizeof(ELEM_PHYSICS *), nelems, elemPhys);
    
}
