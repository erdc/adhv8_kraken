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
void smodel_super_alloc_init_array(SMODEL_SUPER **smod, int nSuperModels) {
    
    // assertions
    assert(nSuperModels > 0);  // must have a least one superModel defined to run
    
    // allocate array of supermodels
    (*smod) = (SMODEL_SUPER *) tl_alloc(sizeof(SMODEL_SUPER), nSuperModels);
    SMODEL_SUPER *sm; // alias
    
    // initialize
    int isup;
    for (isup=0; isup<nSuperModels; isup++) {
        sm = &(*smod)[imat]; // alias
        smodel_super_alloc_init(sm);
        //sm[isup].matrix_petsc = NULL;
        //sm[isup].matrix_adh = NULL;
    }
    //sm->dof_map_local = NULL;
    //sm-dof_map_global = NULL;
}


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
void smodel_super_alloc_init(SMODEL_SUPER *sm) {
    
    // assertions
    assert(nSuperModels > 0);  // must have a least one superModel defined to run
    
    //SIO?
    sm->io = NULL;  /* file names associated with this model application */
    //SFLAGS flag;
    sm->o_flag = 0;
    sm->grid = NULL; // just a pointer to the grid in the design model

    // initialize with default values
    sm->isSimple = 0;
    sm->dt = 0.0;
    sm->dt_old = 0.0;
    sm->dt_err = 0.0;
    sm->dt_prev = 0.0;
    sm->inc_nonlin = 0.0;
    sm->tol_nonlin = 0.0;
    sm->t_init = 0.0;
    sm->t_prev = 0.0;
    sm->t_final = 0.0;
    sm->t_adpt_flag = 0;
    sm->nseries = 0;              // the number of series in this model 
    sm->itrns = 0;

    //Mark adding other things that used to be part of solver
    //should some of these be part of LINSYS?
    sm->max_nonlin_linesearch_cuts = 0;
    sm->max_nonlin_it = 0;
    sm->it_count_nonlin = 0;
    sm->force_nonlin_it = 0;
    sm->nonlinear_it_total = 0;
    sm->LINEAR_PROBLEM = 0;
    sm->it_count_nonlin_failed = 0;


    // now all held within slin_sys structure
    sm->lin_sys = NULL; //pointer to the design model's linear system

    //IDK if physics_mat_code is still used? dont believe we need dof_map_global
    sm->physics_mat_code = NULL; // an code for each physics material that
                           // determinds which equations will be
                           // solved on that material
    sm->dof_map_local = NULL;    // a local map from the local node ID to the
                           // local equation number for building the
                           // FE residual and matrix 
                           //onl allocate if isSimple is not 1

    sm->dof_map_global = NULL;   // for ghost nodes - a map from the local
                           // node ID to the global equation number
                           // for building the FE matrix 

    sm->meshcode = 0; //need a way if supermodel is defined on entire mesh (0), surface(1), or floor(2)
 

    //Mark adding sizes for convenience
    sm->nphysics_mat_1d = 0;
    sm->nphysics_mat_2d = 0;
    sm->nphysics_mat_3d = 0;

    //Mark, elements need integer to store the physics mat id
    //do we want to assume grid is given so we can actuall allocate these?
    sm->elem1d_physics_mat_id = NULL; //[nelem1d]
    sm->elem2d_physics_mat_id = NULL; //[nelem2d]
    sm->elem3d_physics_mat_id = NULL; //[nelem3d]
    
    //same with elem physics mats, these will come from super file
    sm->elem1d_physics_mat = NULL;
    sm->elem2d_physics_mat = NULL;
    sm->elem3d_physics_mat = NULL;

     //Mark proposes swapping nodal vars above to node-based material, this will cut down on memory
    //but may be challenging to form. This won't be set by user but implicitly built
    //at run time
    sm->nphysics_mat_node = 0;
    sm->node_physics_mat_id = NULL; //[nnode] ? local vs what idk
    sm->node_physics_mat = NULL; //[nphysics_mat_node]


    //SHOULD THESE BE IN SLINSYS OR SUPERMODEL??
    /* boundary conditions mask */
    //maybe we can get rid of this through weak enforcement
    sm->bc_mask = NULL;
    //instead of mask, allocate array of ints that is dof of each dirichlet dof
    //int *dirchlet_dofs;
    sm->dirichlet_data = NULL;
    //Mark added local to process, this should be same as local_range[1]-local_range[0]
    //are these redundant in any way to what we have in SLIN_SYS?
    sm->my_ndofs = NULL; //pointers to design model, not arrays
    sm->my_ndofs_old = NULL;
    sm->ndofs = NULL; // local number of degrees of freedom
    sm->ndofs_old = NULL; //local numer of solution variables the processor is in charge of
    sm->macro_ndofs = NULL;
    sm->macro_ndofs_old = NULL;
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
    int nphysics_mat_1d, int nphysics_mat_2d, int nphysics_mat_3d, char elemVarCode[4], int isSimple,
    SGRID *grid, SLIN_SYS *sys) {
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
    sm->isSimple = isSimple;

    // need to build the lin sys and the grid
    //start building up materials
    sm->grid = grid;
    sm->lin_sys = sys;



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
    int nvars[1];

    int **subMod_nvars;

    //utilize fact only one of these models has a physics material defined
    nSubModels[0] = 1;
    ntrns[0] = 0;
    nvars[0] = 3;
    //this double pointer will be a set of nmat arrays, each of size nSubModels[i]
    //tells us for each submodel how many variables there are
    (subMod_nvars) = (int **) tl_alloc(sizeof(int *), sm->nphysics_mat_2d);
    for (i=0;i<sm->nphysics_mat_2d;i++){
        (subMod_nvars[i]) = (int*) tl_alloc(sizeof(int), nSubModels[i]);
        for(j=0;j<subMod_nvars[i];j++){
            subMod_nvars[i][j] = 1;
        }
    }


    //initialize pointers and allocate memory for physics materials
    smat_physics_alloc_init_array(&sm->elem1d_physics_mat, sm->nphysics_mat_1d, ntrns, nvars, nSubModels, subMod_nvars);
    smat_physics_alloc_init_array(&sm->elem2d_physics_mat, sm->nphysics_mat_2d, ntrns, nvars, nSubModels, subMod_nvars);
    smat_physics_alloc_init_array(&sm->elem3d_physics_mat, sm->nphysics_mat_3d, ntrns, nvars, nSubModels, subMod_nvars);

//    smat_physics_alloc_init(&sm->elem1d_physics_mat, sm->nphysics_mat_1d, ntrns);
//    smat_physics_alloc_init(&sm->elem2d_physics_mat, sm->nphysics_mat_2d, ntrns);
//    smat_physics_alloc_init(&sm->elem3d_physics_mat, sm->nphysics_mat_3d, ntrns);
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
    
    //now use elemVarCode info to pick out the correct pointers to routines
    //maybe store this in some sort of dictionary sort of thing?
    //for now just hard code this one example, will go back and fix later once
    //we have the structs and methods figured out
    //how does this work?
    //sm->elem2d_physics[0][0].fe_resid = fe_sw2_body_resid;
    //hard code to poisson routine for now
    sm->elem2d_physics_mat[0]->elem_physics[0][0].fe_resid = poisson_residual;
    sm->elem2d_physics_mat[0]->elem_physics[0][0].nvar = 1;
    sm->elem2d_physics_mat[0]->elem_physics[0][0].physics_vars[0] = PERTURB_U;
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
    //smat_physics_alloc_init(&sm->node_physics_mat, sm->nphysics_mat_node, ntrns);
    nSubModels[0] = 0;
    smat_physics_alloc_init_array(&sm->node_physics_mat, sm->nphysics_mat_node, ntrns, nvars, nSubModels, subMod_nvars);
    //node code is same as elemVarCode
    for(i=0;i<sm->nphysics_mat_node;i++){
        for(j=0;j<4;j++){
            strcpy(&sm->node_physics_mat[i].elemVarCode[j],&elemVarCode[j]);
        }
    }


    //also need to set the variables, not just equations
    //elemental materials first
    //should already be set
    assert(sm->elem2d_physics_mat[0].nvar==3);
    sm->elem2d_physics_mat[0].vars[0] = PERTURB_V;
    sm->elem2d_physics_mat[0].vars[1] = PERTURB_U;
    sm->elem2d_physics_mat[0].vars[2] = PERTURB_H;
    sm->elem2d_physics_mat[0].nSubmodels = 1;
    //same for nodes
    //sm->node_physics_mat[0].nvar = 3;
    assert(sm->node_physics_mat[0].nvar == 3;);
    //sm->node_physics_mat[0].vars = (int*) tl_alloc(sizeof(int), sm->node_physics_mat->nvar);
    sm->node_physics_mat[0].vars[0] = PERTURB_V;
    sm->node_physics_mat[0].vars[1] = PERTURB_U;
    sm->node_physics_mat[0].vars[2] = PERTURB_H;

    //initalize residual vector
    //need a routine to get ndofs for all of this
    //hard coded for now
    
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