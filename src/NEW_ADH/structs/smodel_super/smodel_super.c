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

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Initialize super model without reading a file for testing
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] smod                (SMODEL_SUPER *)  an AdH superModel
 * @param[in]  FILE                    (FILE *) the SuperModel input file
 * \note This supermodel is already assumed to have a grid pointer within it that is populated
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_super_no_read_simple(SMODEL_SUPER *sm, double dt_in, double t_init, double t_final,
    int nphysics_mat_1d, int nphysics_mat_2d, int nphysics_mat_3d, char elemVarCode[4] ) {
    // assign scalars
    sm->dt = dt_in;
    sm->old_dt = dt_in;
    sm->dt_err = dt_in;
    sm->dt_prev = dt_in;
    sm->inc_nonlin = 1e-3;
    sm->tol_nonlin = 1e-5;
    sm->t_init = t_init;
    sm->t_prev = t_init;
    sm->t_final = t_final;
    sm->t_adpt_flag = 0;
    sm->nseries = 0;              // the number of series in this model 
    sm->itrns = 0;

    //start building up materials

    //this is how many different unique physics materials in each dimension there are
    sm->nphysics_mat_1d = nphysics_mat_1d;
    sm->nphysics_mat_2d = nphysics_mat_2d;
    sm->nphysics_mat_3d = nphysics_mat_3d;

    //now we allocate the physics_mat arrays
    smat_physics_alloc_init(&sm->elem1d_physics_mat, sm->nphysics_mat_1d, NULL);
    smat_physics_alloc_init(&sm->elem2d_physics_mat, sm->nphysics_mat_2d, NULL);
    smat_physics_alloc_init(&sm->elem3d_physics_mat, sm->nphysics_mat_3d, NULL);

    //now we fill in the physics mat objects with the provided codes to specify the equations to be solved
    //for now we just assume the elemVarCode is for all the elements, either 1,2, or 3d
    //check inputs are valid, only 1 of these should be nonzero and equal to 1
    assert( nphysics_mat_1d == 0 || nphysics_mat_2d == 0);
    assert( nphysics_mat_1d == 0 || nphysics_mat_3d == 0);
    assert( nphysics_mat_2d == 0 || nphysics_mat_3d == 0);
    assert( nphysics_mat_1d == 1 || nphysics_mat_2d == 1 || nphysics_mat_3d == 1);

    //fill in elemVarCode for the physics mat that is nonempty
    int i,j;
    for(i=0;i<nphysics_mat_1d;i++){
        for(j=0;j<4;j++){
            strcpy(&sm->elem1d_physics_mat[i].elemVarCode[j],&elemVarCode[j]);
        }
    }
    for(i=0;i<nphysics_mat_2d;i++){
        for(j=0;j<4;j++){
            strcpy(&sm->elem2d_physics_mat[i].elemVarCode[j],&elemVarCode[j]);
        }
    }
    for(i=0;i<nphysics_mat_3d;i++){
        for(j=0;j<4;j++){
            strcpy(&sm->elem3d_physics_mat[i].elemVarCode[j],&elemVarCode[j]);
        }
    }



}