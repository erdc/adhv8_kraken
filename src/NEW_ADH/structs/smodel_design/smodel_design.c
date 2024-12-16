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
void smodel_design_alloc_init(SMODEL_DESIGN **dmod, int nSuperModels, int nMono, int nSimple,
    int nUnique) {
    int i;
    
    // allocate array/pointer of design models, for now hard set at 1
    (*dmod) = (SMODEL_DESIGN *) tl_alloc(sizeof(SMODEL_DESIGN), 1);

    
    // allocate the design grid (may choose to do different superModel grids later)
    // Need to add SFILE_IN gridFile to SGRID STRUCTURE
    // Need to add SFILE_IN facesFile to SGRID STRUCTURE
    // MPI and Communicate need sit somewhere.  Think about this
    //sgrid_read(&((*dmod)->grid),grid_file,filename);
    
    
    // initialize the design model and all its super and sub models
    (*dmod)->nSuperModels = nSuperModels;
    (*dmod)->nMono = nMono;
    (*dmod)->nSimple = nSimple;
    (*dmod)->nUnique = nUnique;

    (*dmod)->ndofs = (int*) tl_alloc(sizeof(int), nUnique);
    (*dmod)->ndofs_old = (int*) tl_alloc(sizeof(int), nUnique);
    (*dmod)->my_ndofs = (int*) tl_alloc(sizeof(int), nUnique);
    (*dmod)->my_ndofs_old = (int*) tl_alloc(sizeof(int), nUnique);
    (*dmod)->macro_ndofs = (int*) tl_alloc(sizeof(int), nUnique);
    (*dmod)->macro_ndofs_old = (int*) tl_alloc(sizeof(int), nUnique);

    //allocate arra of nSuperModels super models
    smodel_super_alloc_init_array(&((*dmod)->superModel),nSuperModels); //,nFluxInterfaces);

    //allocate array of nUnique linear systems
    slin_sys_alloc_init_array(&((*dmod)->lin_sys),nUnique); 


    //array of nMono dof maps maybe? but this comes from supermodels so maybe not
    
    
}
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
void smodel_design_free(SMODEL_DESIGN *dmod) {
    
    // aliases
    SMODEL_SUPER *sm_p = dmod->superModel;
    SGRID *grid_p = dmod->grid;
    
    smodel_super_free(sm_p,dmod->nSuperModels);
    sgrid_free(grid_p);
    //nFluxInterfaces = (int *) tl_free(sizeof(int), sm_p->nSuperModels, nFluxInterfaces);
    dmod = (SMODEL_DESIGN *) tl_free(sizeof(SMODEL_DESIGN), 1, dmod);
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
void smodel_design_no_read_simple(SMODEL_DESIGN **dm_ptr, double dt_in, double t_init, double t_final,
    int nphysics_mat_1d, int nphysics_mat_2d, int nphysics_mat_3d, char elemVarCode[4] ,
    SGRID *grid) {
    
    SMODEL_DESIGN *dm;
    int i,j;
    int isSimple=0;
    printf("Initializing design model without file read\n");
    smodel_design_alloc_init(dm_ptr, 1, 1, 0,1);
    printf("Design model alloc init completed\n");

    //creat alias for shorter code
    dm = *dm_ptr;
    //assign the pointer to the grid
    dm->grid = grid;

    printf("DM GRID NELEMS2D %d\n",dm->grid->nelems2d);

    //simple case where design model has one super model
    //these should be the values
    //dm->nSuperModels = 1;
    //even though it is simple supermodel, lets call it Mono to activate things
    //dm->nMono = 1;
    //dm->nUnique = 1;
    //allocate the proper number of super models
    //dm->superModel = (SMODEL_SUPER*) tl_alloc(sizeof(SMODEL_SUPER), dm->nSuperModels);
    //allocate the proper number of linear systems
    //dm->lin_sys = (SMODEL_SUPER*) tl_alloc(sizeof(SMODEL_SUPER), dm->nSuperModels);
    // assign scalars
    dm->dt = dt_in;
    dm->old_dt = dt_in;
    dm->dt_err = dt_in;
    dm->dt_prev = dt_in;
    dm->inc_nonlin = 1e-3;
    dm->tol_nonlin = 1e-5;
    dm->t_init = t_init;
    dm->t_prev = t_init;
    dm->t_final = t_final;
    dm->t_adpt_flag = 0;
    dm->nseries = 0;              // the number of series in this model 
    dm->itrns = 0;
    printf("set some basic params in design model %d\n", dm->nSuperModels);
    //fill out materials in each super model without reading superfile

    for(i=0;i<dm->nSuperModels;i++){
        printf("Calling super model init\n");
        printf("Some grid stuff %d\n",dm->grid->nelems2d);
        smodel_super_no_read_simple(&(dm->superModel[i]), dt_in, t_init, t_final,
        nphysics_mat_1d, nphysics_mat_2d, nphysics_mat_3d, elemVarCode, isSimple,
        dm->grid, &(dm->lin_sys[i]));
        printf("Super model init completed\n");
        
        dm->superModel[i].tol_nonlin = 1e-5;
        dm->superModel[i].inc_nonlin = 1e-3;
        dm->superModel[i].max_nonlin_linesearch_cuts = 5;
        dm->superModel[i].it_count_nonlin_failed = 0;
        dm->superModel[i].max_nonlin_it = 20;
        dm->superModel[i].LINEAR_PROBLEM = NO;
        dm->superModel[i].force_nonlin_it = NO;
        dm->superModel[i].force_nonlin_it = NO;
        dm->superModel[i].nonlinear_it_total = 0; 


    }

    for(i=0;i<dm->nUnique;i++){
        //define some pointers in each supermodel
        //assign some hard coded values first, then set the pointers
        //in general these would require info taken from each superModel
        dm->ndofs[i] = dm->grid->nnodes*3;
        dm->ndofs_old[i] = 0;
        dm->my_ndofs[i] = dm->grid->nnodes*3;
        dm->my_ndofs_old[i] = 0;
        dm->macro_ndofs[i] = dm->grid->nnodes*3;
        dm->macro_ndofs_old[i] = 0;



        
        dm->superModel[i].my_ndofs = &(dm->my_ndofs[i]); //pointers to design model, not arrays
        dm->superModel[i].my_ndofs_old = &(dm->my_ndofs_old[i]);
        dm->superModel[i].ndofs = &(dm->ndofs[i]);
        dm->superModel[i].ndofs_old = &(dm->ndofs_old[i]);
        dm->superModel[i].macro_ndofs = &(dm->macro_ndofs[i]);
        dm->superModel[i].macro_ndofs_old = &(dm->macro_ndofs_old[i]);
        dm->superModel[i].bc_mask = (int*) tl_alloc(sizeof(int), *(dm->superModel[i].ndofs));
        dm->superModel[i].dirichlet_data = (double*) tl_alloc(sizeof(double), *(dm->superModel[i].ndofs));


        printf("Assigning some pointers, ndofs = %d , %d\n",dm->ndofs[i],  *(dm->superModel[i].ndofs));
        //also set the lin sys variables
        dm->lin_sys[i].local_range[1] = dm->grid->nnodes*3;
        dm->lin_sys[i].local_range[0] = 0;
        dm->lin_sys[i].local_range_old[1] = 0;
        dm->lin_sys[i].local_range_old[0] = 0;
        dm->lin_sys[i].nghost = 0;
        //should ths just be pointers back to dm->ndofs[]?
        dm->lin_sys[i].size = dm->grid->nnodes*3;
        dm->lin_sys[i].local_size = dm->grid->nnodes*3;
        dm->lin_sys[i].global_size = dm->grid->nnodes*3;

        //when in init need to change this to check if we have more than one processor or not
        dm->lin_sys[i].residual=NULL;
        dm->lin_sys[i].sol=NULL;
        dm->lin_sys[i].indptr_diag=NULL;
        dm->lin_sys[i].indptr_off_diag=NULL;

        //should this go somewhere else? Like still in super model since we have pointer now?
        //but this should only happen over each unique model
        dm->lin_sys[i].residual = (double*) tl_alloc(sizeof(double), dm->lin_sys[i].size);
        dm->lin_sys[i].sol = (double*) tl_alloc(sizeof(double), dm->lin_sys[i].size);
        dm->lin_sys[i].indptr_diag = (int*) tl_alloc(sizeof(int), dm->lin_sys[i].local_size+1);
        dm->lin_sys[i].sol_old = (double*) tl_alloc(sizeof(double), dm->lin_sys[i].size);
        dm->lin_sys[i].sol_older = (double*) tl_alloc(sizeof(double), dm->lin_sys[i].size);
        //actual solution of linear system is an increment within Newton iteration
        dm->lin_sys[i].dsol = (double*) tl_alloc(sizeof(double), dm->lin_sys[i].size);
        dm->lin_sys[i].scale_vect = (double*) tl_alloc(sizeof(double), dm->lin_sys[i].size);

        //init to 0's here?
        sarray_init_dbl(dm->lin_sys[i].residual,dm->lin_sys[i].size);
        sarray_init_dbl(dm->lin_sys[i].sol,dm->lin_sys[i].size);
        sarray_init_dbl(dm->lin_sys[i].sol_old,dm->lin_sys[i].size);
        sarray_init_dbl(dm->lin_sys[i].sol_older,dm->lin_sys[i].size);
        sarray_init_dbl(dm->lin_sys[i].dsol,dm->lin_sys[i].size);
        sarray_init_dbl(dm->lin_sys[i].scale_vect,dm->lin_sys[i].size);

        //now with info from the super models, assign appropriate vals to the design model
        //for each unique model we build up the sparsity and solution variables
        //allocate CSR sparsity structure for each unique model  
        //missing a step to find proper index of supermodel that is unique
        printf("Calling split CSR creator\n");
        create_sparsity_split_CSR(dm->superModel[i].lin_sys, &(dm->superModel[i]), dm->grid);
    }

    printf("DM ndofs[0] = %d\n",dm->ndofs[0]);
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
