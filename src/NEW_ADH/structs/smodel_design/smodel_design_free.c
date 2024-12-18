#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees an AdH Designer Model
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] dmod           (SDMODEL *)  a pointer to an AdH design-level model
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_design_free(SMODEL_DESIGN *dm) {
    
    //free more complex structures
    if(dm->lin_sys!=NULL){slin_sys_free_array(dm->lin_sys,dm->nUnique);}
    if(dm->superModel != NULL) {smodel_super_free_array(dm->superModel,dm->nSuperModels);}
    if(dm->grid !=NULL){sgrid_free(dm->grid);}
    //free basic arrays, complex models have pointers to these so must go last
    if(dm->unique_id!=NULL){dm->unique_id = (int *) tl_free(sizeof(int), dm->nUnique, dm->unique_id);}
    if(dm->lin_sys_id!=NULL){dm->lin_sys_id = (int *) tl_free(sizeof(int), dm->nSuperModels, dm->lin_sys_id);}
    if(dm->my_ndofs!=NULL){dm->my_ndofs = (int *) tl_free(sizeof(int), dm->nUnique, dm->my_ndofs);}
    if(dm->my_ndofs_old!=NULL){dm->my_ndofs_old = (int *) tl_free(sizeof(int), dm->nUnique, dm->my_ndofs_old);}
    if(dm->ndofs!=NULL){dm->ndofs = (int *) tl_free(sizeof(int), dm->nUnique, dm->ndofs);}
    if(dm->ndofs_old!=NULL){dm->ndofs_old = (int *) tl_free(sizeof(int), dm->nUnique, dm->ndofs_old);}
    if(dm->macro_ndofs!=NULL){dm->macro_ndofs = (int *) tl_free(sizeof(int), dm->nUnique, dm->macro_ndofs);}
    if(dm->macro_ndofs_old!=NULL){dm->macro_ndofs_old = (int *) tl_free(sizeof(int), dm->nUnique, dm->macro_ndofs_old);}


    

}
