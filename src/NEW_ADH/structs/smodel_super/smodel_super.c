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
        sm = &(*smod)[isup]; // alias
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
    
    //SIO?
    sm->io = NULL;  /* file names associated with this model application */
    //SFLAGS flag;
    sm->o_flag = 0;
    sm->grid = NULL; // just a pointer to the grid in the design model

    //Mark adding some default vals
    sm->tau_temporal = 0; //can change from super model to super model
    sm->gravity=9.81; //idk where this should be? probably design, this can just be pointer
    sm->density=1.0;
    sm->vorticity_id = UNSET_INT;

    // initialize with default values
    sm->isSimple = 0;
    sm->dt = NULL;
    sm->old_dt = NULL;
    sm->dt_err = NULL;
    sm->dt_prev = NULL;
    sm->inc_nonlin = 0;
    sm->tol_nonlin = 0;
    sm->t_init = NULL;
    sm->t_prev = NULL;
    sm->t_final = NULL;
    sm->t_adpt_flag = NULL;
    sm->nsubsteps = 1;
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
    sm->sol = NULL;
    sm->sol_old = NULL;
    sm->sol_older = NULL;
    sm->dof_map_local = NULL;    // a local map from the local node ID to the
                           // local equation number for building the
                           // FE residual and matrix 
                           //onl allocate if isSimple is not 1

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
    //sm->nphysics_mat_node = 0;
    //sm->node_physics_mat_id = NULL; //[nnode] ? local vs what idk
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

    sm->forward_step = 0;
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
    
    printf("dt: %10.5f || dt_prev: %10.5f\n",*(smod->dt),*(smod->dt_prev));
    printf("inc_nonlin: %10.5e || tol_nonlin: %10.5e\n",smod->inc_nonlin,smod->tol_nonlin);
    
}

