/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smodel_super.c This file collects methods of the SUPER_MODEL structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes an array of AdH Super Models
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] smod           (SUPER_MODEL **)  a double pointer to an array of AdH supermodels
 * @param[in]  nSuperModels            (int) the total number of supermodels in the design model
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_super_alloc_init(SMODEL_SUPER **smod, int nSuperModels) {
    
    // assertions
    assert(nSuperModels > 0);  // must have a least one superModel defined to run
    
    // allocate all the supermodels
    (*smod) = (SMODEL_SUPER *) tl_alloc(sizeof(SMODEL_SUPER), nSuperModels);
    SMODEL_SUPER *sm = (*smod); // alias
    
    // initialize
    int isup;
    for (isup=0; isup<nSuperModels; isup++) {
        sm[isup].dt = 0.0;
        sm[isup].dt_prev = 0.0;
        sm[isup].inc_nonlin = 0.0;
        sm[isup].tol_nonlin = 0.0;
        //sm[isup].matrix_petsc = NULL;
        //sm[isup].matrix_adh = NULL;
    }
 
    sm->physics_mat_code = NULL;
    //sm->dof_map_local = NULL;
    //sm-dof_map_global = NULL;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees an array of AdH supermodels
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] smod           (SMODEL_SUPER *)  an array of AdH supermodels
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_super_free(SMODEL_SUPER *smod, int nSuperModels) {
    smod = (SMODEL_SUPER *) tl_free(sizeof(SMODEL_SUPER), nSuperModels, smod);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Read SuperModel data from a file
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] smod                (SMODEL_SUPER *)  an AdH superModel
 * @param[in]  FILE                    (FILE *) the SuperModel input file
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_super_read(SMODEL_SUPER *smod, FILE *fp) {
    // read it
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Prints SuperModel data to screen
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] supmod           (SMODEL_SUPER *)  an array of AdH supermodels
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_super_printScreen(SMODEL_SUPER *smod) {
    
    printf("dt: %10.5f || dt_prev: %10.5f\n",smod->dt,smod->dt_prev);
    printf("inc_nonlin: %10.5e || tol_nonlin: %10.5e\n",smod->inc_nonlin,smod->tol_nonlin);
    
}
