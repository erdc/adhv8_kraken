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
 *             Only designed to take one material for an entire mesh
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveland
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
    int i,j;
    printf("Initializing smodel super without file read\n");
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



    //now we fill in the physics mat objects with the provided codes to specify the equations to be solved
    //for now we just assume the elemVarCode is for all the elements, either 1,2, or 3d
    //check inputs are valid, only 1 of these should be nonzero and equal to 1
    assert( nphysics_mat_1d == 0 || nphysics_mat_2d == 0);
    assert( nphysics_mat_1d == 0 || nphysics_mat_3d == 0);
    assert( nphysics_mat_2d == 0 || nphysics_mat_3d == 0);
    assert( nphysics_mat_1d == 1 || nphysics_mat_2d == 1 || nphysics_mat_3d == 1);

    //this is how many different unique physics materials in each dimension there are
    // all but one are 0 and the other one is 1
    sm->nphysics_mat_1d = nphysics_mat_1d;
    sm->nphysics_mat_2d = nphysics_mat_2d;
    sm->nphysics_mat_3d = nphysics_mat_3d;


    //use this simple set up to mark each cell with a single 0
    //this should be easy to determine from grid read
    //we need go put this in other routine maybe?
    if (sm->grid->nelems3d >0){
        sm->elem3d_physics_mat_id = (int*) tl_alloc(sizeof(int), sm->grid->nelems3d);
    }
    if(sm->grid->nelems2d>0){
        sm->elem2d_physics_mat_id = (int*) tl_alloc(sizeof(int), sm->grid->nelems2d);
    }
    if(sm->grid->nelems1d>0){
        sm->elem1d_physics_mat_id = (int*) tl_alloc(sizeof(int), sm->grid->nelems1d);
    }

    //now assign the physics single mat id for this simple set up
    //in general this would come from some sort of grid read
    for(i=0;i<sm->grid->nelems3d;i++){
        sm->elem3d_physics_mat_id[i] = 0;
    }
    for(i=0;i<sm->grid->nelems2d;i++){
        sm->elem2d_physics_mat_id[i] = 0;
    }
    for(i=0;i<sm->grid->nelems1d;i++){
        sm->elem1d_physics_mat_id[i] = 0;
    }

    //now we allocate the physics_mat arrays
    printf("Number of physics mats in each dimension (1d,2d,3d): (%d,%d,%d)\n",sm->nphysics_mat_1d,sm->nphysics_mat_2d,sm->nphysics_mat_3d);
    int ntrns[1];
    int nSubModels[1];
    //utilize fact only one of these models has a physics material defined
    nSubModels[0] = 1;
    ntrns[0] = 0;
    //initialize pointers and allocate memory for physics materials
    smat_physics_alloc_init(&sm->elem1d_physics_mat, sm->nphysics_mat_1d, ntrns);
    smat_physics_alloc_init(&sm->elem2d_physics_mat, sm->nphysics_mat_2d, ntrns);
    smat_physics_alloc_init(&sm->elem3d_physics_mat, sm->nphysics_mat_3d, ntrns);
    printf("Initialized phyisics materials\n");


    //fill in elemVarCode for the physics mat that is nonempty
    //this assumes only one elemVarCode was given, need to generalize and put in a routine
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

    //now fill in the number of physics modules on each physics material
    //given the elemVarCode
    //must be single quotes!!!
    //Maybe store these in global variable instead of hardcoded comparison on LHS
    //replace this with a call to routine which uses elemVarCode and establises
    // physics routines
    for(i=0;i<nphysics_mat_1d;i++){

        if (sm->elem1d_physics_mat[i].elemVarCode[0] == '1'){
            printf("SW1 Activated for 1d physics mat id %d\n",i);
        }
        if (sm->elem1d_physics_mat[i].elemVarCode[1] == '0'){
            printf("GW NOT ACTIVATED for 1d physics mat id %d\n",i);
        }
        if (sm->elem1d_physics_mat[i].elemVarCode[2] == '0'){
            printf("Transport NOT ACTIVATED for 1d physics mat id %d\n",i);   
        }
    }
    for(i=0;i<nphysics_mat_2d;i++){
        if (sm->elem2d_physics_mat[i].elemVarCode[0] == '2'){
            printf("SW2 Activated for 2d physics mat id %d\n",i);
        }else if(sm->elem2d_physics_mat[i].elemVarCode[0] == '6'){
            printf("DW Activated for 2d physics mat id %d\n",i);
        }

        if (sm->elem2d_physics_mat[i].elemVarCode[1] == '0'){
            printf("GW NOT ACTIVATED for 2d physics mat id %d\n",i);
        }
        if (sm->elem2d_physics_mat[i].elemVarCode[2] == '0'){
            printf("Transport NOT ACTIVATED for 2d physics mat id %d\n",i);   
        }

    }
    for(i=0;i<nphysics_mat_3d;i++){
        if (sm->elem3d_physics_mat[i].elemVarCode[0] == '3'){
            printf("SW3 Activated for 3d physics mat id %d\n",i);
        }else if(sm->elem2d_physics_mat[i].elemVarCode[0] == '4'){
            printf("NS Activated for 3d physics mat id %d\n",i);
        }

        if (sm->elem3d_physics_mat[i].elemVarCode[1] == '0'){
            printf("GW NOT ACTIVATED for 3d physics mat id %d\n",i);
        }else if(sm->elem3d_physics_mat[i].elemVarCode[1] == '1'){
            printf("GW ACTIVATED for 3d physics mat id %d\n",i);
        }

        if (sm->elem3d_physics_mat[i].elemVarCode[2] == '0'){
            printf("Transport NOT ACTIVATED for 3d physics mat id %d\n",i);   
        }
    }





    //Use infor in physics_mat to establish correct pointers to residual routines
    //i think it may make sense to include the elem_physics structure as
    //part of the physics mat
    printf("Initializing elem physics routines\n");
    selem_physics_alloc_init(&sm->elem1d_physics, sm->nphysics_mat_1d, nSubModels);
    selem_physics_alloc_init(&sm->elem2d_physics, sm->nphysics_mat_2d, nSubModels);
    selem_physics_alloc_init(&sm->elem3d_physics, sm->nphysics_mat_3d, nSubModels);
    printf("Initialized elem physics routines\n");
    
    //now use elemVarCode info to pick out the correct pointers to routines
    //maybe store this in some sort of dictionary sort of thing?
    //for now just hard code this one example, will go back and fix later once
    //we have the structs and methods figured out
    //how does this work?
    //sm->elem2d_physics[0][0].fe_resid = fe_sw2_body_resid;
    //hard code to poisson routine for now
    sm->elem2d_physics[0][0].fe_resid =poisson_residual;
    sm->elem2d_physics[0][0].nvar = 1;
    sm->elem2d_physics[0][0].physics_vars = (int*) tl_alloc(sizeof(int), sm->elem2d_physics[0][0].nvar);
    sm->elem2d_physics[0][0].physics_vars[0] = PERTURB_U;
    //sm->elem2d_physics[0][0].physics_vars[1] = PERTURB_U;
    //sm->elem2d_physics[0][0].physics_vars[2] = PERTURB_V;




    //hard code nodal vars too and dof_map_local, easy because it is all same
    //need to form routines to form correct nodal mats
    //this will give us proper mappings
    sm->nphysics_mat_node = 1;
    sm->node_physics_mat_id = (int*) tl_alloc(sizeof(int), sm->grid->nnodes);
    for(i=0;i<sm->grid->nnodes;i++){
        sm->node_physics_mat_id[i] = 0;
    }
    smat_physics_alloc_init(&sm->node_physics_mat, sm->nphysics_mat_node, ntrns);
    //node code is same as elemVarCode
    for(i=0;i<sm->nphysics_mat_node;i++){
        for(j=0;j<4;j++){
            strcpy(&sm->node_physics_mat[i].elemVarCode[j],&elemVarCode[j]);
        }
    }


    //also need to set the variables, not just equations
    //elemental materials first
    sm->elem2d_physics_mat[0].nvar=3;
    sm->elem2d_physics_mat[0].vars = (int*) tl_alloc(sizeof(int), sm->elem2d_physics_mat->nvar);
    sm->elem2d_physics_mat[0].vars[0] = PERTURB_V;
    sm->elem2d_physics_mat[0].vars[1] = PERTURB_U;
    sm->elem2d_physics_mat[0].vars[2] = PERTURB_H;
    sm->elem2d_physics_mat[0].nSubmodels = 1;
    //same for nodes
    sm->node_physics_mat[0].nvar=3;
    sm->node_physics_mat[0].vars = (int*) tl_alloc(sizeof(int), sm->node_physics_mat->nvar);
    sm->node_physics_mat[0].vars[0] = PERTURB_V;
    sm->node_physics_mat[0].vars[1] = PERTURB_U;
    sm->node_physics_mat[0].vars[2] = PERTURB_H;

    //initalize residual vector
    //need a routine to get ndofs for all of this
    //hard coded for now
    sm->ndofs_old = 0;
    sm->ndofs = sm->grid->nnodes*3;
    sm->my_ndofs = sm->grid->nnodes*3;
    sm->local_range[1] = sm->grid->nnodes*3;
    sm->local_range[0] = 0;
    sm->local_range_old[1] = 0;
    sm->local_range_old[0] = 0;
    sm->residual=NULL;
    sm->sol=NULL;
    sm->indptr_diag=NULL;
    sm->indptr_off_diag=NULL;
    sm->residual = (double*) tl_alloc(sizeof(double), sm->ndofs);
    sm->sol = (double*) tl_alloc(sizeof(double), sm->ndofs);
    sm->indptr_diag = (int*) tl_alloc(sizeof(int), sm->my_ndofs+1);
    //when in init need to change this to check if we have more than one processor or not
    sm->indptr_off_diag = NULL;
    
    //make a trivial fmap_local
    //maybe have special case when nmat is 1 to avoid storing this and use special case mapping

    //in practice this fmap will be built element by element, not nodally
    //for CG only nodal is OK
    sm->dof_map_local = (int*) tl_alloc(sizeof(int), sm->grid->nnodes);
    for(i=0;i<sm->grid->nnodes;i++){
        sm->dof_map_local[i] = i*3;
    }

    //allocate_adh_system(sm);

    printf("Assigned sw2 to residual structure\n");

    //call residual
    //double elemrhs[3]; 
    //elemrhs[0] = 0.0; elemrhs[1] = 0.0; elemrhs[2] = 0.0;
    //sm->elem2d_physics[0][0].fe_resid(sm, elemrhs, 0, 0.0,0,0,0, 0);
    //printf("Elemental resid: {%f,%f,%f}\n",elemrhs[0],elemrhs[1],elemrhs[2]);


}