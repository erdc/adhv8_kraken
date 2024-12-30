/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  sgrid.c This file collects methods of the SMODEL_DESIGN structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes an AdH Designer Model
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] dmod           (SDMODEL **)  a pointer to an AdH design-level model
 * @param[in]  nSuperModels            (int) the total number of supermodels in the design model
 * @param[in]  nSubModels                (int*) the total number of submodels in each supermodel
 * @param[in]  nFluxInterfaces      (int*) the total number of flux interfaces between supermodels
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_design_alloc(SMODEL_DESIGN *dmod, int nSuperModels, int nMono, int nSimple,
    int nUnique) {
    // allocate the design grid (may choose to do different superModel grids later)
    // Need to add SFILE_IN gridFile to SGRID STRUCTURE
    // Need to add SFILE_IN facesFile to SGRID STRUCTURE
    // MPI and Communicate need sit somewhere.  Think about this
    //sgrid_read(&((*dmod)->grid),grid_file,filename);
    
    
    // initialize the design model and all its super and sub models
    dmod->nSuperModels = nSuperModels;
    dmod->nMono = nMono;
    dmod->nSimple = nSimple;
    dmod->nUnique = nUnique;

    dmod->ndofs = (int*) tl_alloc(sizeof(int), nUnique);
    dmod->ndofs_old = (int*) tl_alloc(sizeof(int), nUnique);
    dmod->my_ndofs = (int*) tl_alloc(sizeof(int), nUnique);
    dmod->my_ndofs_old = (int*) tl_alloc(sizeof(int), nUnique);
    dmod->macro_ndofs = (int*) tl_alloc(sizeof(int), nUnique);
    dmod->macro_ndofs_old = (int*) tl_alloc(sizeof(int), nUnique);
    dmod->lin_sys_id = (int*) tl_alloc(sizeof(int), nSuperModels);
    dmod->unique_id = (int*) tl_alloc(sizeof(int), nUnique);

    printf("callimg smodel_super_alloc\n");
    //allocate arra of nSuperModels super models
    smodel_super_alloc_init_array(&(dmod->superModel),nSuperModels); //,nFluxInterfaces);
    printf("callimg slin_sys_alloc_array\n");
    //allocate array of nUnique linear systems
    slin_sys_alloc_array(&(dmod->lin_sys),nUnique); 
    //array of nMono dof maps maybe? but this comes from supermodels so maybe not
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Prints an AdH Designer Model to screen
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] dmod           (SDMODEL *)  a pointer to an AdH design-level model
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_design_printScreen(SMODEL_DESIGN *dmod) {
    int iSuperMod = 0;
    
    printf("++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("nSuperModels: %d\n",dmod->nSuperModels);
    for (iSuperMod=0; iSuperMod<dmod->nSuperModels; iSuperMod++) {
        printf("--------------------------------------------\n");
        printf("SuperModel #%d \n",iSuperMod);
        printf("--------------------------------------------\n");
        smodel_super_printScreen(&(dmod->superModel[iSuperMod]));
    }
    //sgrid_printScreen(dmod->grid);
    printf("++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("++++++++++++++++++++++++++++++++++++++++++++\n");
}
