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
 * @param[inout] elemPhys           (SMODEL **)  the struct double pointer
 * @param[in]  nelems           (int) the total number of elements (1D, 2D, or 3D)
 * @param[in]  nSubModels                (int*) the total number of submodels on each element
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//void smodel_alloc_init_double_pointer(SMODEL ***elemPhys,int nelems,int *nSubMods) { // triple pointer?
//    int ie, ifun;
//    
//    assert(nelems>=0);
//    if (nelems == 0){
//        elemPhys = NULL;
//        return;
//    }else{
//        //al
//        (*elemPhys) = (SMODEL **) tl_alloc(sizeof(SMODEL *), nelems);//

//        SMODEL *physics; // for alias
//        for (ie=0; ie<nelems; ie++) {
//            assert(nSubMods[ie] > 0);
//            (*elemPhys[ie]) = (SMODEL *) tl_alloc(sizeof(SMODEL), nSubMods[ie]);
//            for (ifun=0; ifun<nSubMods[ie]; ifun++) {
//                physics = &(*elemPhys)[ie][ifun]; // alias
//                physics->fe_inc = NULL;
//                physics->fe_init = NULL;
//                physics->fe_update = NULL;
//                physics->fe_solve = NULL;
//                physics->fe_resid = NULL;
//                physics->fe_load = NULL;
//            }
//        }
//        return;
//    }   
//}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and array of size nelem of physice FE routines
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elemPhys           (SMODEL **)  the struct double pointer
 * @param[in]  nelems           (int) the total number of elements (1D, 2D, or 3D)
 * @param[in]  nSubModels                (int*) the total number of submodels on each element
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_alloc_init_array(SMODEL **elemPhys, int nSubMods, int *nSubMod_nvar) { // triple pointer?
    int ie;
    
    assert(nSubMods>0);
    //allocate an array of elemPhys
    (*elemPhys) = (SMODEL *) tl_alloc(sizeof(SMODEL), nSubMods);

    SMODEL *physics; // for alias
    for (ie=0; ie<nSubMods; ie++) {
        physics = &(*elemPhys)[ie]; // alias
        smodel_alloc_init(physics, nSubMod_nvar[ie]);
    }
        return;
}   




/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and array of size nelem of physice FE routines
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elemPhys           (SMODEL **)  the struct double pointer
 * @param[in]  nelems           (int) the total number of elements (1D, 2D, or 3D)
 * @param[in]  nSubModels                (int*) the total number of submodels on each element
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_alloc_init(SMODEL *physics,int nvar) { // triple pointer?
    assert(nvar>0);
    physics->fe_inc = NULL;
    physics->fe_init = UNSET_INT;
    physics->fe_update = NULL;
    physics->fe_solve = NULL;
    physics->fe_resid = UNSET_INT;
    physics->fe_load = NULL;
    physics->nvar = nvar;
    //allocate array of physics variables and set to default
    (physics->physics_vars) = (int*) tl_alloc(sizeof(int), nvar);
    for(int i=0; i<nvar; i++){
        physics->physics_vars[i] = 0;
    }
    return;
}   
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and array of size nelem of physice FE routines
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elemPhys           (SMODEL **)  the struct double pointer
 * @param[in]  nelems           (int) the total number of elements (1D, 2D, or 3D)
 * @param[in]  nSubModels                (int*) the total number of submodels on each element
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_free_array(SMODEL *Models, int nModels) { // triple pointer?
    int ie;
    for (ie=0; ie<nModels; ie++) {
        smodel_free(&(Models[ie]));
    }
    //free the array of elemPhys
    Models = (SMODEL *) tl_free(sizeof(SMODEL), nModels, Models);
    return;
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
 * @param[inout] elemPhys           (SMODEL **)  the struct double pointer
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_free(SMODEL *model) {
 
    (model->physics_vars) = (int*) tl_free(sizeof(int), model->nvar, model->physics_vars);
    
}
