#include "global_header.h"
//void couple_edge_to_edge(SSUPER_MODEL *sm);

/***********************************************************/
/***********************************************************/
static int DEBUG = OFF;
static int DEBUG_FMAP = OFF;

/***********************************************************/
/***********************************************************/

void ssuperModel_alloc_init(SSUPER_MODEL **superModel_ptr, int nsupermodels, int *nsubmodels, int *nInterfaces) {
    
    int i,j;
    assert(nsupermodels>0); // global variable
    (* superModel_ptr) = (SSUPER_MODEL *) tl_alloc(sizeof(SSUPER_MODEL), nsupermodels);
    SSUPER_MODEL *sm = (* superModel_ptr);
    //superModel = (SSUPER_MODEL *) tl_alloc(sizeof(SSUPER_MODEL), nsupermodels);
    for (i=0; i<nsupermodels; i++) {
        assert(nsubmodels[i] > 0);
        if(nsupermodels==1) {
            sm[i].proc_flag = 1; // no division of processors necessary
        } else{
            sm[i].proc_flag = 0; // processors get divided in superfile read
        }
        sm[i].nsubmodels = nsubmodels[i];
        sm[i].submodel = (SMODEL *) tl_alloc(sizeof(SMODEL), nsubmodels[i]);
        for (j=0;j<nsubmodels[i];j++){
            if(nsubmodels[i]==1){
                sm[i].submodel[j].proc_flag=1; // no division of processors necessary
            }else{
                sm[i].submodel[j].proc_flag=0; // processors get divided in superfile read
            }
        }
        // gkc below
        sm[i].NumInterfaces = nInterfaces[i];
        if (nInterfaces[i] != 0) {
            sm[i].interface = (SINTERFACE *) tl_alloc(sizeof(SINTERFACE), nInterfaces[i]);
        }
        else{
            sm[i].interface = (SINTERFACE *) NULL;
        }
        
        
        sm[i].supersmpi = (SMPI *) tl_alloc(sizeof(SMPI), 1);
        sm[i].supersmpi_wvel = (SMPI *) tl_alloc(sizeof(SMPI), 1);
        
        //gkc above
        
        sm[i].nsys = 0;
        sm[i].nsys_sq = 0;
        sm[i].max_nsys = 0;
        sm[i].max_nsys_sq = 0;
        
        sm[i].my_nnodes = 0;
        sm[i].nnodes = 0;
        sm[i].macro_nnodes = 0;
        sm[i].wvel_my_nnodes = 0;
        sm[i].wvel_nnodes = 0;
        sm[i].wvel_macro_nnodes = 0;
        
        sm[i].macro_nnodes_original = 0;
        sm[i].macro_nnodes = 0;
        sm[i].total_macro_nnodes = 0;
        sm[i].total_nnodes = 0;
        sm[i].nnodes_matrix = 0;
        sm[i].total_nnodes_matrix = 0;
        
        sm[i].nblock = 1;
        sm[i].inc_nonlin = 0.0;
        sm[i].tol_nonlin = 0.0;
        sm[i].max_nonlin_it = 0;
        sm[i].nalloc_inc = 0;
        sm[i].nonlinear_it_total = 0;
        sm[i].nonlinear_it_total_hvel = 0;
        sm[i].nonlinear_it_total_wvel = 0;
        sm[i].ntransport = 0;
        
        solv_initialize(&(sm[i].solver_info));
        sm[i].solver_info.prec_value = 1;
        sm[i].solver_info.node_block = NULL;
        
        /* WVEL */
        sm[i].save_wvel_firstcall = 1;
        
        sm[i].fmap = NULL;
        sm[i].fmap_nodes = NULL;
        sm[i].fmap_wvel = NULL;
        sm[i].fmap_wvel_nodes = NULL;
#ifdef _PETSC
        int ierr;
        sm[i].A = PETSC_NULL;
        sm[i].P = PETSC_NULL;
        sm[i].sol = PETSC_NULL;
        sm[i].residual = PETSC_NULL;
        sm[i].ksp = PETSC_NULL;
        sm[i].bc_mask = NULL;
        sm[i].old_bc_mask_size = 0;
        sm[i].scale_vect = NULL;
        sm[i].PREALLOC_FLAG = ON;
        // Get flag for preallocating; default is to preallocate
        ierr = PetscOptionsGetInt(NULL,NULL,"-prealloc_flag",&(sm[i].PREALLOC_FLAG),NULL);
#else
        sm[i].matrix = NULL;
        sm[i].bc_mask = NULL;
        sm[i].residual = NULL;
        sm[i].sol = NULL;
        sm[i].scale_vect = NULL;
        sm[i].diagonal = NULL;
#endif
        sm[i].proc_map = NULL;
        sm[i].proc_map_wvel = NULL;
    }
}

/***********************************************************/
/***********************************************************/

void ssuperModel_forward_map(SSUPER_MODEL *sm) {
    int i, j, k, sub_pe,ii;
    SINTERFACE *ifce;
    SINTERFACE_NODELIST *ndlist;
    int** old_fmap;
    
    if (sm->fmap_nodes == NULL){
        sm->fmap_nodes = (int *) tl_alloc(sizeof(int ) , sm->nsubmodels);
    }
    if (sm->fmap == NULL){
        for (i=0; i<sm->nsubmodels; i++) sm->fmap_nodes[i]=0;
        sm->fmap = (int **) tl_alloc(sizeof(int *) , sm->nsubmodels);
    }
    //old_fmap = (int **) tl_alloc(sizeof(int *) , sm->nsubmodels);
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            if (sm->submodel[i].grid->max_nnodes > sm->fmap_nodes[i]){
                sm->submodel[i].fmap = NULL; /* to prevent memory leak? I think! */
                sm->fmap[i] = (int *) tl_realloc(sizeof(int) , sm->submodel[i].grid->max_nnodes, sm->fmap_nodes[i], sm->fmap[i]); // realloc for adaptivity
                // for(j=0;j<sm->submodel[i].grid->nnodes; j++){
                //     sm->fmap[i][j] = 0;
                // }
                sm->fmap_nodes[i] = sm->submodel[i].grid->max_nnodes;
            }
        }
    }
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            for(j=0;j< sm->submodel[i].grid->nnodes;j++){
                sm->fmap[i][j] = 0.0;
            }
        }
    }
    
    int **fmap = sm->fmap; /* alias */
    int nnodes_cum_sum = 0, my_nnodes_cum_sum=0;
    //sm->total_nnodes_matrix = 0;
    sm->total_macro_nnodes = 0;
#ifdef _MESSG
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            // collecting all resident nodes at the beginning of the matrix
            for (j=0; j<sm->submodel[i].grid->my_nnodes; j++){
                fmap[i][j] = my_nnodes_cum_sum + j; // Slightly complicated mapping!
            }
            my_nnodes_cum_sum += sm->submodel[i].grid->my_nnodes;
            //sm->total_nnodes_matrix += sm->submodel[i].grid->nnodes_matrix;
            sm->total_macro_nnodes += sm->submodel[i].grid->macro_nnodes;
        }
    }
    
    // // Sending ghost nodes to the end
    // nnodes_cum_sum = my_nnodes_cum_sum; // will be updated in the following loop
    // for (i=0; i<sm->nsubmodels; i++){
    //     for (j=sm->submodel[i].grid->my_nnodes; j<sm->submodel[i].grid->nnodes; j++){
    //         fmap[i][j] = nnodes_cum_sum + (j-sm->submodel[i].grid->my_nnodes); // Slightly complicated mapping!
    //     }
    //     nnodes_cum_sum += (sm->submodel[i].grid->nnodes - sm->submodel[i].grid->my_nnodes);
    // }
    
    int i_pe; /* counter for PE's */
    nnodes_cum_sum = my_nnodes_cum_sum; // will be updated in the following loop
    int *lastpoint = (int *) tl_alloc(sizeof(int), sm->nsubmodels);
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            lastpoint[i] = sm->submodel[i].grid->my_nnodes;
        }else{
            lastpoint[i] = 0;
        }
    }
    for (i_pe=0; i_pe<sm->supersmpi->npes; i_pe++){
        for (i=0; i<sm->nsubmodels; i++){
            if(sm->submodel[i].proc_flag==1){
                sub_pe = sm->proc_map[i_pe][i];
                int jj = sm->submodel[i].grid->my_nnodes;
                if(sub_pe != UNSET_INT && sm->submodel[i].grid->smpi->npes>0){
                    for (j=0,jj=sm->submodel[i].grid->my_nnodes; j<sm->submodel[i].grid->smpi->nrecv[sub_pe]; j++,jj++){
                        fmap[i][j+lastpoint[i]] = nnodes_cum_sum + j; // Slightly complicated mapping!
                    }
                    lastpoint[i]   += sm->submodel[i].grid->smpi->nrecv[sub_pe];
                    nnodes_cum_sum += sm->submodel[i].grid->smpi->nrecv[sub_pe];
                }
            }
        }
    }
    
#else
    for (i=0; i<sm->nsubmodels; i++){
        for (j=0; j<sm->submodel[i].grid->nnodes; j++){
            fmap[i][j] = nnodes_cum_sum + j; // Each node is mapped to itslef, for now
        }
        nnodes_cum_sum += sm->submodel[i].grid->nnodes;
        my_nnodes_cum_sum += sm->submodel[i].grid->my_nnodes;
        //sm->total_nnodes_matrix += sm->submodel[i].grid->nnodes_matrix;
        sm->total_macro_nnodes += sm->submodel[i].grid->macro_nnodes;
    }
#endif
    
    
    
    //sm->total_nnodes = nnodes_cum_sum;  /* This is used for (re)allocating the Jacobian and residuals */
    sm->total_my_nnodes = my_nnodes_cum_sum;
    
    for (i=0; i<sm->NumInterfaces; i++){
        ifce = &(sm->interface[i]);
        int master_mod, slave_mod;
        
        // for shallow water to shallow water coupling
        if (sm->submodel[ifce->model_id[0]].flag.SW_FLOW && sm->submodel[ifce->model_id[1]].flag.SW_FLOW) {
            if (sm->submodel[ifce->model_id[0]].flag.SW3_FLOW && sm->submodel[ifce->model_id[1]].flag.SW2_FLOW){
                /* Only if user input was model 0 = 3D, and model 1 = 2D. Then the order is flipped in memory. */
                master_mod = 1;
                slave_mod  = 0;
            }
            else{
                master_mod = 0;
                slave_mod = 1;
            }
        } else {
            // cjt :: not sure order matters for other types of coupling. need to verify this
            master_mod = 0;
            slave_mod = 1;
        }
        
        // Set all Slave model node fmaps to -1 temporarily
        for (j=0; j<ifce->NumNodeColumns; j++){
            ndlist = &(ifce->nodelist[j]);
            if ((ndlist->size[master_mod] != 0)&&(sm->submodel[ifce->model_id[master_mod]].proc_flag>0)){
                for (k=0; k<ndlist->size[slave_mod]; k++){
                    //old_fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][k]] = fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][k]];
                    fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][k]] = -1;
                }
            }/*else if((ndlist->size[slave_mod] != 0)&&(sm->submodel[ifce->model_id[slave_mod]].proc_flag>0)){
              for (k=1; k<ndlist->size[slave_mod]; k++){
              
              fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][k]] = -1;
              }
              }*/
        }
    }
    
    
    char descr[125];
    int reduction_count = 0;
    nnodes_cum_sum = 0;
#ifdef _MESSG
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            for (j=0; j<sm->submodel[i].grid->my_nnodes; j++){
                if (fmap[i][j] != nnodes_cum_sum + j){
                    reduction_count++;
                }
                fmap[i][j]-=reduction_count;
            }
            nnodes_cum_sum += sm->submodel[i].grid->my_nnodes;
        }
    }
    
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            lastpoint[i] = sm->submodel[i].grid->my_nnodes;
        }else{
            lastpoint[i]=0;
        }
    }
    
    
    for (i_pe=0; i_pe<sm->supersmpi->npes; i_pe++){
        for (i=0; i<sm->nsubmodels; i++){
            if(sm->submodel[i].proc_flag==1){
                sub_pe = sm->proc_map[i_pe][i];
                int jj = sm->submodel[i].grid->my_nnodes;
                if((sub_pe != UNSET_INT)&&(sm->submodel[i].grid->smpi->npes>0)){
                    for (j=0,jj=lastpoint[i]; j<sm->submodel[i].grid->smpi->nrecv[sub_pe]; j++,jj++){
                        if (fmap[i][jj] != nnodes_cum_sum + j){
                            reduction_count++;
                        }
                        fmap[i][jj]-=reduction_count;
                    }
                    lastpoint[i]   += sm->submodel[i].grid->smpi->nrecv[sub_pe];
                    nnodes_cum_sum += sm->submodel[i].grid->smpi->nrecv[sub_pe];
                }
            }
        }
    }
    
    nnodes_cum_sum -= reduction_count;
    sm->total_nnodes = nnodes_cum_sum;
    
    lastpoint = (int *) tl_free(sizeof(int), sm->nsubmodels, lastpoint);
#else
    for (i=0; i<sm->nsubmodels; i++){
        for (j=0; j<sm->submodel[i].grid->nnodes; j++){
            if (fmap[i][j] != nnodes_cum_sum + j){
                reduction_count++;
            }
            fmap[i][j]-=reduction_count;
        }
        nnodes_cum_sum += sm->submodel[i].grid->nnodes;
    }
#endif
    
    int rootnode_equation_id;
    // Find the number of master and slave interface nodes
    int n_master_ifc_nodes=0, n_master_ifc_my_nodes=0;  // Number of master interface nodes
    int n_slave_ifc_nodes=0, n_slave_ifc_my_nodes=0;   // Number of slave interface nodes
    for (i=0; i<sm->NumInterfaces; i++){
        ifce = &(sm->interface[i]);
        int master_mod, slave_mod;
        
        
        if (sm->submodel[ifce->model_id[0]].flag.SW_FLOW && sm->submodel[ifce->model_id[1]].flag.SW_FLOW) {
            if (sm->submodel[ifce->model_id[0]].flag.SW3_FLOW && sm->submodel[ifce->model_id[1]].flag.SW2_FLOW){
                /* Only if user input was model 0 = 3D, and model 1 = 2D. Then the order is flipped in memory. */
                /* Caution: This affects the immediate next for loop! */
                master_mod = 1;
                slave_mod  = 0;
            }
            else{
                master_mod = 0;
                slave_mod = 1;
            }
        } else {
            // cjt :: not sure order matters for other types of coupling. need to verify this
            master_mod = 0;
            slave_mod = 1;
        }
        
        for (j=0; j<sm->interface[i].NumNodeColumns; j++){
            if(&(ifce->nodelist[j])!=NULL){ndlist = &(ifce->nodelist[j]);}
            
            if ((ndlist->size[master_mod] != 0)&&(sm->submodel[ifce->model_id[master_mod]].proc_flag>0)){
                n_master_ifc_nodes += ifce->nodelist[j].size[master_mod];
                n_slave_ifc_nodes += ifce->nodelist[j].size[slave_mod];
                if(sm->submodel[sm->interface[i].model_id[master_mod]].grid->node[sm->interface[i].nodelist[j].couplednodes[master_mod][0]].resident_pe == sm->submodel[sm->interface[i].model_id[master_mod]].grid->smpi->myid){
                    n_master_ifc_my_nodes += ifce->nodelist[j].size[master_mod];
                    n_slave_ifc_my_nodes += ifce->nodelist[j].size[slave_mod];
                }
                if (sm->submodel[ifce->model_id[master_mod]].flag.SW3_FLOW){ // ONLY possible if BOTH models are SW3!
                    for (k=0; k<ndlist->size[slave_mod]; k++){
                        fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][k]] =
                        fmap[ifce->model_id[master_mod]][ndlist->couplednodes[master_mod][k]];
                    }
                } else {
                    // printf("\nPE[%i] ndlist[%i].size[%i] = %i", sm->submodel[0].grid->smpi->myid, j, slave_mod, ndlist->size[slave_mod]);
                    rootnode_equation_id = fmap[ifce->model_id[master_mod]][ndlist->couplednodes[master_mod][0]]; // This is the 2D interface node
                    for (k=0; k<ndlist->size[slave_mod]; k++){ // All slave nodes mapped to the single 2D node.
                        fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][k]] = rootnode_equation_id;
                    }
                    //n_master_ifc_nodes += ifce->nodelist[j].size[master_mod];
                    //n_slave_ifc_nodes += ifce->nodelist[j].size[slave_mod];
                    if(sm->submodel[sm->interface[i].model_id[master_mod]].grid->node[sm->interface[i].nodelist[j].couplednodes[master_mod][0]].resident_pe == sm->submodel[sm->interface[i].model_id[master_mod]].grid->smpi->myid){
                        //n_master_ifc_my_nodes += ifce->nodelist[j].size[master_mod];
                        //n_slave_ifc_my_nodes += ifce->nodelist[j].size[slave_mod];
                    }
                }
            }/*else if((ndlist->size[slave_mod] != 0)&&(sm->submodel[ifce->model_id[slave_mod]].proc_flag>0)){
              if (sm->submodel[ifce->model_id[master_mod]].flag.SW3_FLOW){ // ONLY possible if BOTH models are SW3!
              for (k=0; k<ndlist->size[slave_mod]; k++){
              fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][k]] =
              fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][0]];
              }
              n_master_ifc_nodes += ifce->nodelist[j].size[slave_mod];
              n_slave_ifc_nodes += ifce->nodelist[j].size[slave_mod];
              } else {
              rootnode_equation_id = fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][0]]; // This is the 2D interface node
              for (k=0; k<ndlist->size[slave_mod]; k++){ // All slave nodes mapped to the single 2D node.
              fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][k]] = rootnode_equation_id;
              }
              n_master_ifc_nodes += 1;
              n_slave_ifc_nodes += ifce->nodelist[j].size[slave_mod];
              }
              }*/
        }
    }
    
    // Find the number of non-interface nodes - need that for each model
    int *n_non_ifc_nodes, *n_non_ifc_my_nodes;
    n_non_ifc_nodes = (int *) tl_alloc(sizeof(int), sm->nsubmodels); // Number of non-interface nodes
    n_non_ifc_my_nodes = (int *) tl_alloc(sizeof(int), sm->nsubmodels);
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            n_non_ifc_nodes[i] = sm->submodel[i].grid->nnodes;
            n_non_ifc_my_nodes[i] = sm->submodel[i].grid->my_nnodes;
        }else{
            n_non_ifc_nodes[i] = 0;
            n_non_ifc_my_nodes[i] = 0;//if(sm->supersmpi->myid==3)printf("model %d setting n_non nodes %d and %d  \n",i, n_non_ifc_my_nodes[i], n_non_ifc_nodes[i]);
        }
    }
    
    for (i=0; i<sm->NumInterfaces; i++){
        for (j=0; j<sm->interface[i].NumNodeColumns; j++){
            if(sm->interface[i].nodelist[j].size[0]>0){
                if(sm->submodel[sm->interface[i].model_id[0]].grid->node[sm->interface[i].nodelist[j].couplednodes[0][0]].resident_pe ==
                   sm->submodel[sm->interface[i].model_id[0]].grid->smpi->myid){
                    n_non_ifc_my_nodes[sm->interface[i].model_id[0]] -= sm->interface[i].nodelist[j].size[0];
                }
                n_non_ifc_nodes[sm->interface[i].model_id[0]] -= sm->interface[i].nodelist[j].size[0];
            }
            if(sm->interface[i].nodelist[j].size[1]>0){
                if(sm->submodel[sm->interface[i].model_id[1]].grid->node[sm->interface[i].nodelist[j].couplednodes[1][0]].resident_pe ==
                   sm->submodel[sm->interface[i].model_id[1]].grid->smpi->myid){
                    n_non_ifc_my_nodes[sm->interface[i].model_id[1]] -= sm->interface[i].nodelist[j].size[1];        }
                n_non_ifc_nodes[sm->interface[i].model_id[1]] -= sm->interface[i].nodelist[j].size[1];
            }
        }
    }
    
    // Find the number of master and slave interface nodes
    // int n_master_ifc_nodes=0;  // Number of master interface nodes
    // int n_slave_ifc_nodes=0;   // Number of slave interface nodes
    for (i=0; i<sm->NumInterfaces; i++){
        ifce = &(sm->interface[i]);
        int master_mod, slave_mod;
        if (sm->submodel[ifce->model_id[0]].flag.SW_FLOW && sm->submodel[ifce->model_id[1]].flag.SW_FLOW) {
            if (sm->submodel[ifce->model_id[0]].flag.SW3_FLOW && sm->submodel[ifce->model_id[1]].flag.SW2_FLOW){
                /* Only if user input was model 0 = 3D, and model 1 = 2D. Then the order is flipped in memory. */
                master_mod = 1;
                slave_mod  = 0;
            }
            else {
                master_mod = 0;
                slave_mod = 1;
            }
        } else {
            master_mod = 0;
            slave_mod = 1;
        }
        for (j=0; j<ifce->NumNodeColumns; j++){
            if(sm->interface[i].nodelist[j].size[master_mod]>0){
                if(sm->submodel[sm->interface[i].model_id[master_mod]].grid->node[sm->interface[i].nodelist[j].couplednodes[master_mod][0]].resident_pe == sm->submodel[sm->interface[i].model_id[master_mod]].grid->smpi->myid){
                    //      n_master_ifc_nodes += ifce->nodelist[j].size[master_mod];
                    //      n_slave_ifc_nodes += ifce->nodelist[j].size[slave_mod];
                }
            }
        }
    }
    
#ifdef _DEBUG
    if (DEBUG_FMAP){
#ifdef _MESSG
        int proc_count; int myid = sm->supersmpi->myid;
        for(proc_count=0;proc_count<sm->supersmpi->npes;proc_count++){
            if(proc_count==sm->supersmpi->myid){
                for(i=0; i<sm->nsubmodels; i++){
                    printf("[PE%i] Number of HVEL non-interface nodes for submodel %2i : %d\n", myid, i, n_non_ifc_nodes[i]);
                }
                printf("[PE%i] HVEL map: Number of master interface nodes         : %d\n", myid, n_master_ifc_nodes);
                printf("[PE%i] HVEL map: Number of slave interface nodes          : %d\n", myid, n_slave_ifc_nodes);
            }
            //messg_barrier(sm->supersmpi->ADH_COMM);
        }
#else
        for(i=0; i<sm->nsubmodels; i++){
            printf("Number of HVEL non-interface nodes for submodel %i : %d\n", i, n_non_ifc_nodes[i]);
        }
        printf("HVEL map: Number of master interface nodes : %d\n", n_master_ifc_nodes);
        printf("HVEL map: Number of slave interface nodes : %d\n", n_slave_ifc_nodes);
#endif
    }
    //exit(-1);
#endif
    
    sm->nnodes = 0;
    sm->my_nnodes = 0;
    sm->macro_nnodes = sm->total_macro_nnodes;
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            for (j=0; j<sm->submodel[i].grid->nnodes; j++){
                if((sm->my_nnodes<fmap[i][j])&&(j<sm->submodel[i].grid->my_nnodes)) sm->my_nnodes = fmap[i][j];
                if((sm->nnodes<fmap[i][j])) sm->nnodes = fmap[i][j];
            }
        }
    }

#if _DEBUG
    FILE *fp;
    if(DEBUG_FMAP){
        for (i=0; i<sm->nsubmodels; i++){
            if(sm->submodel[i].proc_flag==1){
                sprintf(descr,"FMAP_%d_%d.txt", i, sm->supersmpi->myid);
                fp = fopen(descr, "a");
                for (j=0; j<sm->submodel[i].grid->nnodes; j++){
                    fprintf(fp,"%d    %d \n",j,fmap[i][j]);
                }
                fclose(fp);
            }
        }
    }
#endif
    
    
    sm->my_nnodes++;
    sm->nnodes++;
    //sm->nnodes += n_master_ifc_nodes;    /* Final size of the HVEL matrix, but NOT used for memory allocation! */
    //sm->my_nnodes += n_master_ifc_my_nodes; /* Resident portion of the HVEL matrix */
    
    
#ifdef _MESSG
    sm->macro_nnodes = messg_isum(sm->my_nnodes,sm->supersmpi->ADH_COMM);
#else
    sm->macro_nnodes = sm->my_nnodes;
#endif
    
    sm->macro_nnodes_original = sm->macro_nnodes;
    
    n_non_ifc_nodes = (int *) tl_free(sizeof(int), sm->nsubmodels, n_non_ifc_nodes); // Number of non-interface nodes
    n_non_ifc_my_nodes = (int *) tl_free(sizeof(int), sm->nsubmodels,n_non_ifc_my_nodes);
    
    
    sm->total_nnodes = sm->nnodes;
    
    for (i=0; i<sm->nsubmodels; i++){
        sm->submodel[i].fmap = fmap[i];
    }
    
#ifdef _MESSG
    //tag(MPI_COMM_WORLD); fflush(stdout); MPI_Barrier(MPI_COMM_WORLD);
    comm_set_keys_supermodel(sm);
    //tag(MPI_COMM_WORLD); fflush(stdout); MPI_Barrier(MPI_COMM_WORLD);
    //tl_check_all_pickets(__FILE__,__LINE__); MPI_Barrier(MPI_COMM_WORLD);
    
#endif
    
}

/***********************************************************/
/***********************************************************/

void ssuperModel_forward_map_wvel(SSUPER_MODEL *sm) {
    int i, j, k, sub_pe,ii;
    SINTERFACE *ifce;
    SINTERFACE_NODELIST *ndlist;
    
    // fmap_wvel_node[submodel] :: stores max_nnodes for each submodel
    if (sm->fmap_wvel_nodes == NULL){
        sm->fmap_wvel_nodes = (int *) tl_alloc(sizeof(int ) , sm->nsubmodels);
    }
    
    // fmap_wvel[nsubmodels][submodel[i].grid.nnodes] :: stores the monolithic map for all nnodes in a submodel
    if (sm->fmap_wvel == NULL){
        for (i=0; i<sm->nsubmodels; i++) sm->fmap_wvel_nodes[i]=0;
        sm->fmap_wvel = (int **) tl_alloc(sizeof(int *) , sm->nsubmodels);
    }
    
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            if (sm->submodel[i].grid->max_nnodes > sm->fmap_wvel_nodes[i]){
                sm->submodel[i].fmap_wvel = NULL; /* to prevent memory leak? I think! */
                //printf("nnodes = %i\nnnodes_matrix = %i\n",sm->submodel[i].grid->nnodes, sm->submodel[i].grid->nnodes_matrix);
                sm->fmap_wvel[i] = (int *) tl_realloc(sizeof(int) , sm->submodel[i].grid->max_nnodes, sm->fmap_wvel_nodes[i], sm->fmap_wvel[i]); // realloc for adaptivity
                for(j=0;j<sm->submodel[i].grid->nnodes; j++){
                    sm->fmap_wvel[i][j] = 0;
                }
                sm->fmap_wvel_nodes[i]=sm->submodel[i].grid->max_nnodes;
            }
        }
    }
    
    
    
    int **fmap = sm->fmap_wvel; /* alias */
    int nnodes_cum_sum = 0, my_nnodes_cum_sum = 0;
#ifdef _MESSG
    
    // Map residential nodes 1:1 for now
    for (i=0; i<sm->nsubmodels; i++){
        if (sm->submodel[i].flag.SW3_FLOW && sm->submodel[i].proc_flag==1){
            for (j=0; j<sm->submodel[i].grid->my_nnodes; j++){
                fmap[i][j] = my_nnodes_cum_sum + j;
            }
            my_nnodes_cum_sum += sm->submodel[i].grid->my_nnodes;
        }
    }
    
    int i_pe; /* counter for PE's */
    nnodes_cum_sum = my_nnodes_cum_sum; // will be updated in the following loop
    int *lastpoint = (int *) tl_alloc(sizeof(int), sm->nsubmodels);
    for (i=0; i<sm->nsubmodels; i++){
        if (sm->submodel[i].flag.SW3_FLOW && sm->submodel[i].proc_flag==1){
            lastpoint[i] = sm->submodel[i].grid->my_nnodes;
        }
        else{
            lastpoint[i] = 0;
        }
    }
    for (i_pe=0; i_pe<sm->supersmpi_wvel->npes; i_pe++){
        for (i=0; i<sm->nsubmodels; i++){
            if (sm->submodel[i].flag.SW3_FLOW && sm->submodel[i].proc_flag==1){
                sub_pe =sm->proc_map_wvel[i_pe][i];
                int jj = sm->submodel[i].grid->my_nnodes;
                if(sub_pe != UNSET_INT && sm->submodel[i].grid->smpi->npes>0){
                    for (j=0,jj=sm->submodel[i].grid->my_nnodes; j<sm->submodel[i].grid->smpi->nrecv[sub_pe]; j++,jj++){
                        fmap[i][j+lastpoint[i]] = nnodes_cum_sum + j; // Slightly complicated mapping!
                    }
                    lastpoint[i]   += sm->submodel[i].grid->smpi->nrecv[sub_pe];
                    nnodes_cum_sum += sm->submodel[i].grid->smpi->nrecv[sub_pe];
                }
            }
        }
    }
    
#else
    
    // Map nodes 1:1
    for (i=0; i<sm->nsubmodels; i++){
        if (sm->submodel[i].flag.SW3_FLOW){
            for (j=0; j<sm->submodel[i].grid->nnodes; j++){
                fmap[i][j] = nnodes_cum_sum + j;
            }
            nnodes_cum_sum += sm->submodel[i].grid->nnodes;
            my_nnodes_cum_sum += sm->submodel[i].grid->my_nnodes;
        }
    }
    
#endif
    
    for (i=0; i<sm->NumInterfaces; i++){
        ifce = &(sm->interface[i]);
        if (!(sm->submodel[ifce->model_id[0]].flag.SW3_FLOW && sm->submodel[ifce->model_id[1]].flag.SW3_FLOW)){
            continue;
        }
        
        int master_mod = 0;
        int slave_mod = 1;
        
        // Set all Slave model node fmaps to -1 temporarily
        for (j=0; j<ifce->NumNodeColumns; j++){
            ndlist = &(ifce->nodelist[j]);
            if ((ndlist->size[master_mod] != 0)&&(sm->submodel[ifce->model_id[master_mod]].proc_flag>0)){
                for (k=0; k<ndlist->size[slave_mod]; k++){
                    fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][k]] = -1;
                }
            }/*else if((ndlist->size[slave_mod] != 0)&&(sm->submodel[ifce->model_id[slave_mod]].proc_flag>0)){
              for (k=1; k<ndlist->size[slave_mod]; k++){
              fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][k]] = -1;
              }
              }*/
        }
    }
    
    int reduction_count = 0, total_macro_nnodes = 0;
    nnodes_cum_sum = 0;
#ifdef _MESSG
    for (i=0; i<sm->nsubmodels; i++){
        if (sm->submodel[i].flag.SW3_FLOW && (sm->submodel[i].proc_flag==1)){
            for (j=0; j<sm->submodel[i].grid->my_nnodes; j++){
                if (fmap[i][j] != nnodes_cum_sum + j){
                    reduction_count++;
                }
                fmap[i][j]-=reduction_count;
            }
            nnodes_cum_sum += sm->submodel[i].grid->my_nnodes;
            total_macro_nnodes += sm->submodel[i].grid->macro_nnodes;
        }
    }
    
    for (i=0; i<sm->nsubmodels; i++){
        if (sm->submodel[i].flag.SW3_FLOW && sm->submodel[i].proc_flag==1){
            lastpoint[i] = sm->submodel[i].grid->my_nnodes;
        }
        else{
            lastpoint[i] = 0;
        }
    }
    for (i_pe=0; i_pe<sm->supersmpi_wvel->npes; i_pe++){
        for (i=0; i<sm->nsubmodels; i++){
            if (sm->submodel[i].flag.SW3_FLOW && sm->submodel[i].proc_flag==1){
                sub_pe = sm->proc_map[i_pe][i];
                int jj = sm->submodel[i].grid->my_nnodes;
                if(sub_pe != UNSET_INT&&(sm->submodel[i].grid->smpi->npes>0)){
                    for (j=0,jj=lastpoint[i]; j<sm->submodel[i].grid->smpi->nrecv[sub_pe]; j++,jj++){
                        if (fmap[i][jj] != nnodes_cum_sum + j){
                            reduction_count++;
                        }
                        fmap[i][jj]-=reduction_count;
                    }
                    lastpoint[i]   += sm->submodel[i].grid->smpi->nrecv[sub_pe];
                    nnodes_cum_sum += sm->submodel[i].grid->smpi->nrecv[sub_pe];
                }
            }
        }
    }
    nnodes_cum_sum -= reduction_count;
    lastpoint = (int *) tl_free(sizeof(int), sm->nsubmodels, lastpoint);
#else
    for (i=0; i<sm->nsubmodels; i++){
        if (sm->submodel[i].flag.SW3_FLOW){
            for (j=0; j<sm->submodel[i].grid->nnodes; j++){
                if (fmap[i][j] != nnodes_cum_sum + j){
                    reduction_count++;
                }
                fmap[i][j]-=reduction_count;
            }
            nnodes_cum_sum += sm->submodel[i].grid->nnodes;
            total_macro_nnodes += sm->submodel[i].grid->macro_nnodes;
        }
    }
#endif
    // Find the number of non-interface nodes - need that for each model
    int n_master_ifc_nodes=0, n_master_ifc_my_nodes=0;  // Number of master interface nodes
    int n_slave_ifc_nodes=0, n_slave_ifc_my_nodes=0;   // Number of slave interface nodes
    for (i=0; i<sm->NumInterfaces; i++){
        ifce = &(sm->interface[i]);
        if (!(sm->submodel[ifce->model_id[0]].flag.SW3_FLOW && sm->submodel[ifce->model_id[1]].flag.SW3_FLOW)){
            continue;
        }
        
        int master_mod = 0;
        int slave_mod = 1;
        
        for (j=0; j<sm->interface[i].NumNodeColumns; j++){
            ndlist = &(ifce->nodelist[j]);
            if ((ndlist->size[master_mod] != 0)&&(sm->submodel[ifce->model_id[master_mod]].proc_flag>0)){
                n_master_ifc_nodes += ifce->nodelist[j].size[master_mod];
                n_slave_ifc_nodes += ifce->nodelist[j].size[slave_mod];
                if(sm->submodel[sm->interface[i].model_id[master_mod]].grid->node[sm->interface[i].nodelist[j].couplednodes[master_mod][0]].resident_pe == sm->submodel[sm->interface[i].model_id[master_mod]].grid->smpi->myid){
                    n_master_ifc_my_nodes += ifce->nodelist[j].size[master_mod];
                    n_slave_ifc_my_nodes += ifce->nodelist[j].size[slave_mod];
                }
                for (k=0; k<ndlist->size[slave_mod]; k++){
                    fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][k]] =
                    fmap[ifce->model_id[master_mod]][ndlist->couplednodes[master_mod][k]];
                }
            }/*else if((ndlist->size[slave_mod] != 0)&&(sm->submodel[ifce->model_id[slave_mod]].proc_flag>0)){
              if (sm->submodel[ifce->model_id[master_mod]].flag.SW3_FLOW){ // ONLY possible if BOTH models are SW3!
              for (k=0; k<ndlist->size[slave_mod]; k++){
              fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][k]] =
              fmap[ifce->model_id[slave_mod]][ndlist->couplednodes[slave_mod][0]];
              }
              n_master_ifc_nodes += ifce->nodelist[j].size[slave_mod];
              n_slave_ifc_nodes += ifce->nodelist[j].size[slave_mod];
              }
              }*/
        }
    }
    
    // Find the number of non-interface nodes - need that for each model
    int *n_non_ifc_nodes, *n_non_ifc_my_nodes;
    n_non_ifc_nodes = (int *) tl_alloc(sizeof(int) , sm->nsubmodels); // Number of non-interface nodes
    n_non_ifc_my_nodes = (int *) tl_alloc(sizeof(int) , sm->nsubmodels);
    for (i=0; i<sm->nsubmodels; i++){
        if((sm->submodel[i].proc_flag==1)&&(sm->submodel[i].flag.SW3_FLOW)){
            n_non_ifc_nodes[i] = sm->submodel[i].grid->nnodes;
            n_non_ifc_my_nodes[i] = sm->submodel[i].grid->my_nnodes;
        }else{
            n_non_ifc_nodes[i] = 0;
            n_non_ifc_my_nodes[i] = 0;
        }
    }
    for (i=0; i<sm->NumInterfaces; i++){
        ifce = &(sm->interface[i]);
        if (!(sm->submodel[ifce->model_id[0]].flag.SW3_FLOW && sm->submodel[ifce->model_id[1]].flag.SW3_FLOW)){
            continue;
        }
        for (j=0; j<sm->interface[i].NumNodeColumns; j++){
            if(sm->interface[i].nodelist[j].size[0]>0){
                if(sm->submodel[sm->interface[i].model_id[0]].grid->node[sm->interface[i].nodelist[j].couplednodes[0][0]].resident_pe ==sm->submodel[sm->interface[i].model_id[0]].grid->smpi->myid){
                    n_non_ifc_nodes[sm->interface[i].model_id[0]] -= sm->interface[i].nodelist[j].size[0];
                }
                n_non_ifc_nodes[sm->interface[i].model_id[0]] -= sm->interface[i].nodelist[j].size[0];
            }
            if((sm->interface[i].nodelist[j].size[1]>0)){
                if(sm->submodel[sm->interface[i].model_id[1]].grid->node[sm->interface[i].nodelist[j].couplednodes[1][0]].resident_pe == sm->submodel[sm->interface[i].model_id[1]].grid->smpi->myid){
                    n_non_ifc_my_nodes[sm->interface[i].model_id[1]] -= sm->interface[i].nodelist[j].size[1];
                }
                n_non_ifc_nodes[sm->interface[i].model_id[1]] -= sm->interface[i].nodelist[j].size[1];
            }
        }
    }
    
    // Find the number of master and slave interface nodes
    //int n_master_ifc_nodes=0;  // Number of master interface nodes
    //int n_slave_ifc_nodes=0;   // Number of slave interface nodes
    for (i=0; i<sm->NumInterfaces; i++){
        ifce = &(sm->interface[i]);
        if (!(sm->submodel[ifce->model_id[0]].flag.SW3_FLOW && sm->submodel[ifce->model_id[1]].flag.SW3_FLOW)){ //2D-2D, 2D-3D or 3D-2D only
            if (sm->submodel[ifce->model_id[0]].flag.SW3_FLOW){ // 3D-2D
                for (j=0; j<ifce->NumNodeColumns; j++){
                    if(sm->interface[i].nodelist[j].size[0]>0){
                        if(sm->submodel[sm->interface[i].model_id[0]].grid->node[sm->interface[i].nodelist[j].couplednodes[0][0]].resident_pe == sm->submodel[sm->interface[i].model_id[0]].grid->smpi->myid){
                            //   n_master_ifc_nodes += ifce->nodelist[j].size[0];
                        }
                    }
                }
            }
            else if (sm->submodel[ifce->model_id[1]].flag.SW3_FLOW){ // 2D-3D
                for (j=0; j<ifce->NumNodeColumns; j++){
                    if(sm->interface[i].nodelist[j].size[1]>0){
                        if(sm->submodel[sm->interface[i].model_id[1]].grid->node[sm->interface[i].nodelist[j].couplednodes[1][0]].resident_pe == sm->submodel[sm->interface[i].model_id[1]].grid->smpi->myid){
                            //  n_master_ifc_nodes += ifce->nodelist[j].size[1];
                        }
                    }
                }
            }
        }
        else{ // 3D-3D
            int master_mod = 0;
            int slave_mod = 1;
            for (j=0; j<ifce->NumNodeColumns; j++){
                if(sm->interface[i].nodelist[j].size[master_mod]>0){
                    if(sm->submodel[sm->interface[i].model_id[master_mod]].grid->node[sm->interface[i].nodelist[j].couplednodes[master_mod][0]].resident_pe == sm->submodel[sm->interface[i].model_id[master_mod]].grid->smpi->myid){
                        //  n_master_ifc_nodes += sm->interface[i].nodelist[j].size[master_mod];
                        //  n_slave_ifc_nodes += sm->interface[i].nodelist[j].size[slave_mod];
                    }
                }
            }
        }
    }
    
#ifdef _DEBUG
    if (DEBUG_FMAP){
#ifdef _MESSG
        int proc_count; int myid = sm->supersmpi_wvel->myid;
        for(proc_count=0;proc_count<sm->supersmpi_wvel->npes;proc_count++){
            if(proc_count==sm->supersmpi->myid){
                for(i=0; i<sm->nsubmodels; i++){
                    printf("[PE%i] Number of WVEL non-interface nodes for submodel %2i : %d\n", myid, i, n_non_ifc_nodes[i]);
                }
                printf("[PE%i] WVEL map: Number of master interface nodes         : %d\n", myid, n_master_ifc_nodes);
                printf("[PE%i] WVEL map: Number of slave interface nodes          : %d\n", myid, n_slave_ifc_nodes);
            }
            //messg_barrier(sm->supersmpi->ADH_COMM);
        }
#else
        for(i=0; i<sm->nsubmodels; i++){
            printf("Number of WVEL non-interface nodes for submodel %i : %d\n", i, n_non_ifc_nodes[i]);
        }
        printf("WVEL map: Number of master interface nodes : %d\n", n_master_ifc_nodes);
        printf("WVEL map: Number of slave interface nodes : %d\n", n_slave_ifc_nodes);
#endif
    }
    //exit(-1);
#endif
    sm->wvel_nnodes = 0;
    sm->wvel_my_nnodes = 0;
    for (i=0; i<sm->nsubmodels; i++){
        if (sm->submodel[i].flag.SW3_FLOW && (sm->submodel[i].proc_flag==1)){
            sm->wvel_nnodes += n_non_ifc_nodes[i];
            sm->wvel_my_nnodes += n_non_ifc_my_nodes[i];
        }
    }
    
    sm->wvel_nnodes += n_master_ifc_nodes;  /* Final size of the WVEL matrix, but NOT used for memory allocation! */
    sm->wvel_my_nnodes += n_master_ifc_my_nodes;  /* Resident portion of the WVEL matrix */
    
    n_non_ifc_nodes = (int *) tl_free(sizeof(int), sm->nsubmodels, n_non_ifc_nodes); // Number of non-interface nodes
    n_non_ifc_my_nodes = (int *) tl_free(sizeof(int), sm->nsubmodels,n_non_ifc_my_nodes);
    sm->wvel_nnodes = 0;
    sm->wvel_my_nnodes = 0;
    sm->macro_nnodes = sm->total_macro_nnodes;
    for (i=0; i<sm->nsubmodels; i++){
        if((sm->submodel[i].flag.SW3_FLOW) &&(sm->submodel[i].proc_flag==1)){
            for (j=0; j<sm->submodel[i].grid->nnodes; j++){
                if((sm->wvel_my_nnodes<fmap[i][j])&&(j<sm->submodel[i].grid->my_nnodes)) sm->wvel_my_nnodes = fmap[i][j];
                if((sm->wvel_nnodes<fmap[i][j])) sm->wvel_nnodes = fmap[i][j];
            }
        }
    }
    
    sm->wvel_my_nnodes++;
    sm->wvel_nnodes++;
    for (i=0; i<sm->nsubmodels; i++){
        sm->submodel[i].fmap_wvel = fmap[i];
    }
    
    
#ifdef _MESSG
    sm->wvel_macro_nnodes = messg_isum(sm->wvel_my_nnodes,sm->supersmpi_wvel->ADH_COMM);
#else
    sm->wvel_macro_nnodes = sm->wvel_my_nnodes;
#endif
    
    
    
    if(sm->wvel_nnodes>sm->total_nnodes) sm->total_nnodes = sm->wvel_nnodes;
#if _DEBUG
    if(DEBUG) printf("MYID %d WVEL supermodel my_nnodes %d supermodel nnodes %d total %d macro %d \n",sm->supersmpi->myid,sm->wvel_my_nnodes, sm->wvel_nnodes, sm->total_nnodes, sm->wvel_macro_nnodes);
#endif
#ifdef _MESSG
    comm_set_keys_supermodel_wvel(sm);
#endif
    
}

/***********************************************************/
/***********************************************************/

void ssuperModel_create_interface_data(SSUPER_MODEL *sm) {
    int i, j, k, l;
    ID_LIST_ITEM *ptr;
    SMODEL *mod;
    SGRID *grid;
    SINTERFACE *ifce;
    SINTERFACE_NODELIST *ndlist;
    
    for (i=0; i<sm->NumInterfaces; i++){
        ifce = &(sm->interface[i]);
        for (k=0; k<2; k++){
            mod = &(sm->submodel[ifce->model_id[k]]);
            if(mod->proc_flag==1){
                grid = mod->grid;
                if (mod->flag.SW2_FLOW){
#ifdef _MESSG
                    for (j=0; j<ifce->NumNodeColumns; j++){
                        ndlist = &(ifce->nodelist[j]);
                        if(ndlist->surfnode[k] != UNSET_INT){
                            //ndlist->size[k] = 0;
                            int local_id = 0;
                            for (local_id=0; local_id < grid->nnodes; local_id++){
                                if (grid->node[local_id].gid == ndlist->surfnode[k]){
                                    ndlist->size[k] = 1;
                                    break;
                                }
                            }
                            if(ndlist->size[k]==0){
                                continue;
                            }
                            // printf("\nPE[%i] local_id = %i, gid = %i", sm->submodel[0].grid->smpi->myid, local_id, grid->node[local_id].gid);
                            ndlist->couplednodes[k] = (int *) tl_alloc(sizeof(int), ndlist->size[k]);
                            ndlist->couplednodes[k][0] = local_id;//ndlist->surfnode[k];
                        }
                    }
#else
                    for (j=0; j<ifce->NumNodeColumns; j++){
                        ndlist = &(ifce->nodelist[j]);
                        ndlist->size[k] = 1;
                        ndlist->couplednodes[k] = (int *) tl_alloc(sizeof(int), ndlist->size[k]);
                        ndlist->couplednodes[k][0] = ndlist->surfnode[k];
                    }
#endif
                    
                }
                else if (mod->flag.SW3_FLOW){
                    for (j=0; j<ifce->NumNodeColumns; j++){
                        ndlist = &(ifce->nodelist[j]);
                        if(ndlist->surfnode[k] != UNSET_INT){
#ifdef _MESSG
                            //ndlist->size[k] = 0;
                            int gid = 0;
                            for (gid=0; gid < grid->nnodes; gid++){
                                if (grid->node[gid].gid == ndlist->surfnode[k]){
                                    break;
                                }
                            }
                            
                            int col_id = 0;
                            for (col_id=0; col_id < grid->nnodes_sur; col_id++){
                                if (grid->vertical_list[col_id]->id == gid){
#ifdef _DEBUG
                                    if (DEBUG){
                                        printf("\ngid = %i, grid->vertical_list[%i] = %i", gid, col_id, grid->vertical_list[col_id]->id);
                                    }
#endif
                                    break;
                                }
                            }
                            if (col_id >=grid->nnodes_sur) {
#ifdef _DEBUG
                                if (DEBUG){
                                    printf("[PE%i][Caution] Surface node with ID %i in submodel %i is not present on this processor\n", grid->smpi->myid, ndlist->surfnode[k]+1, ifce->model_id[k]+1);
                                }
#endif
                                continue;
                            }
#else
                            int col_id = 0;
                            for (col_id=0; col_id < grid->nnodes_sur; col_id++){
                                if (grid->vertical_list[col_id]->id == ndlist->surfnode[k]){
                                    break;
                                }
                            }
                            if (col_id >=grid->nnodes_sur) {
                                printf("Surface node with ID %i could not be found in submodel %i\n", ndlist->surfnode[k]+1,ifce->model_id[k]+1);
                                tl_error("Please check interface data");
                            }
#endif
                            /* Get the first node in the vertical list (which is the surface node) */
                            ptr = grid->vertical_list[col_id];
                            //printf("ptr-id = %i\n",ptr->id);
                            int inode= UNSET_INT;
                            int count = 0;
                            
                            while (ptr->next!=NULL){
                                count++;
                                ptr= ptr->next;
                            }
                            ndlist->size[k] = count;
                            ndlist->couplednodes[k] = (int *) tl_alloc(sizeof(int), ndlist->size[k]);
                            ptr = grid->vertical_list[col_id];
                            count = 0;
                            while (ptr->next!=NULL){
                                inode = ptr->id;
                                ndlist->couplednodes[k][count] = inode;
                                count++;
                                ptr = ptr->next; /* skip first node which is surface */
                            }
                        }
                        
                    }
                }
                
                
                //            exit(-1);
            }
            
        }
    }
}

/***********************************************************/
/***********************************************************/

void ssuperModel_free(SSUPER_MODEL **superModel_ptr, SSUPER_INTERFACE **superInterface_ptr, int nsupermodels, int nsuperinterfaces) {
    int i,j,k;
    SSUPER_MODEL *sm = (* superModel_ptr);
    assert(sm != NULL);
    
    for (i=0; i<nsupermodels; i++) {
        int total_nnodes_matrix = sm[i].total_nnodes_matrix;
#ifdef _MESSG
        for(j=0;j<sm[i].supersmpi->npes;j++){
            sm[i].proc_map[j] = tl_free(sizeof(int), sm[i].nsubmodels, sm[i].proc_map[j]);
        }
        sm[i].proc_map = tl_free(sizeof(int *), sm[i].supersmpi->npes, sm[i].proc_map);
#endif
        if(sm[i].proc_map_wvel !=NULL){
            for(j=0;j<sm[i].supersmpi_wvel->npes;j++){;
                sm[i].proc_map_wvel[j] = tl_free(sizeof(int), sm[i].nsubmodels, sm[i].proc_map_wvel[j]);
            }
            sm[i].proc_map_wvel = tl_free(sizeof(int *), sm[i].supersmpi_wvel->npes, sm[i].proc_map_wvel);
        }
#ifdef _PETSC
        int ierr;
        if (sm[i].A != PETSC_NULL) {ierr = MatDestroy(&(sm[i].A));}
        if (sm[i].P != PETSC_NULL) {ierr = MatDestroy(&(sm[i].P));}
        //if (sm[i].bc_mask != NULL) {sm[i].bc_mask = (int *) tl_free(sizeof(int),total_nnodes_matrix * sm[i].max_nsys,sm[i].bc_mask);}
        if (sm[i].bc_mask != NULL) {sm[i].bc_mask = (int *) tl_free(sizeof(int),sm[i].old_bc_mask_size,sm[i].bc_mask);}
        if (sm[i].residual != PETSC_NULL) {ierr = VecDestroy(&(sm[i].residual));}
        if (sm[i].sol != PETSC_NULL) {ierr = VecDestroy(&(sm[i].sol));}
        if (sm[i].ksp != PETSC_NULL) {ierr = KSPDestroy(&(sm[i].ksp));}
#else
        if (sm[i].matrix != NULL) {
            for (j = 0; j < /*sm[i].nnodes_matrix*/ total_nnodes_matrix; j++) {
                sm[i].matrix[j].value = (double *) tl_free(sizeof(double), sm[i].matrix[j].max_size * sm[i].max_nsys_sq, sm[i].matrix[j].value);
                sm[i].matrix[j].index = (int *) tl_free(sizeof(int), sm[i].matrix[j].max_size,sm[i].matrix[j].index);
            }
            sm[i].matrix  = (SPARSE_VECT *) tl_free(sizeof(SPARSE_VECT), total_nnodes_matrix,sm[i].matrix);
        }
        if (sm[i].diagonal != NULL) {sm[i].diagonal = (double *) tl_free(sizeof(double),total_nnodes_matrix * sm[i].max_nsys_sq,sm[i].diagonal);}
        if (sm[i].bc_mask != NULL) {sm[i].bc_mask = (int *) tl_free(sizeof(int),total_nnodes_matrix * sm[i].max_nsys,sm[i].bc_mask);}
        if (sm[i].residual != NULL) {sm[i].residual = (double *) tl_free(sizeof(double),total_nnodes_matrix * sm[i].max_nsys,sm[i].residual);}
        if (sm[i].sol != NULL) {sm[i].sol = (double *) tl_free(sizeof(double),total_nnodes_matrix * sm[i].max_nsys,sm[i].sol);}
        if (sm[i].scale_vect != NULL) {sm[i].scale_vect = (double *) tl_free(sizeof(double),total_nnodes_matrix * sm[i].max_nsys,sm[i].scale_vect);}
#endif

        /* free solver info node block (should be done in Solver_Info (free) ) */
        if (sm[i].solver_info.node_block != NULL){sm[i].solver_info.node_block = (int *) tl_free(sizeof(int),total_nnodes_matrix,sm[i].solver_info.node_block);}
        solv_blk_set_clean(&(sm[i].solver_info.profile));
        // gkc below
        if (sm[i].fmap != NULL){
            for(j=0; j<sm[i].nsubmodels; j++){
                //printf("model %d freeing %d nodes %d \n", j, sm[i].fmap_nodes[j]);
                if ((sm[i].fmap[j] != NULL) && (sm[i].submodel[j].proc_flag==1)) sm[i].fmap[j] = (int *) tl_free(sizeof(int), sm[i].fmap_nodes[j], sm[i].fmap[j]);
            }
            sm[i].fmap = (int **) tl_free(sizeof(int *), sm[i].nsubmodels, sm[i].fmap);
            sm[i].fmap_nodes = (int *) tl_free(sizeof(int ) , sm[i].nsubmodels, sm[i].fmap_nodes);
        }
        if (sm[i].fmap_wvel != NULL){
            for(j=0; j<sm[i].nsubmodels; j++){
                
                if ((sm[i].fmap_wvel[j] != NULL) && (sm[i].submodel[j].proc_flag==1)) sm[i].fmap_wvel[j] = (int *) tl_free(sizeof(int), sm[i].fmap_wvel_nodes[j], sm[i].fmap_wvel[j]);
            }
            sm[i].fmap_wvel = (int **) tl_free(sizeof(int *), sm[i].nsubmodels, sm[i].fmap_wvel);
            sm[i].fmap_wvel_nodes = (int *) tl_free(sizeof(int ) , sm[i].nsubmodels, sm[i].fmap_wvel_nodes);
        }
        
        // gkc above
    }
    for (i=0; i<nsupermodels; i++) {

        // gkc below
        for (j=0; j<sm[i].NumInterfaces; j++) {
            sinterface_free( &(sm[i].interface[j]) );
        }
        if (sm[i].interface != NULL) sm[i].interface = (SINTERFACE *) tl_free(sizeof(SINTERFACE), sm[i].NumInterfaces, sm[i].interface);
        if (sm[i].ntransport > 0) {
            sm[i].con_type = (int *) tl_free(sizeof(int), sm[i].ntransport, sm[i].con_type);
        }
        
        // gkc above
        for (j=0; j<sm[i].nsubmodels; j++) {
            smodel_free( &(sm[i].submodel[j]) );
        }
        sm[i].submodel = (SMODEL *) tl_free(sizeof(SMODEL), sm[i].nsubmodels, sm[i].submodel);
    }
    for (i=0; i<nsupermodels; i++) {
        if (sm[i].supersmpi !=NULL){
            smpi_free(sm[i].supersmpi);
            sm[i].supersmpi = (SMPI *) tl_free(sizeof(SMPI), 1, sm[i].supersmpi);
        }
        if (sm[i].supersmpi_wvel !=NULL){
            smpi_free(sm[i].supersmpi_wvel);
            sm[i].supersmpi_wvel = (SMPI *) tl_free(sizeof(SMPI), 1, sm[i].supersmpi_wvel);
        }
    }
    
    sm = NULL; /* To prevent a memory leak */
    (*superModel_ptr) = (SSUPER_MODEL *) tl_free(sizeof(SSUPER_MODEL), nsupermodels, (*superModel_ptr));
    
    SSUPER_INTERFACE *si = (* superInterface_ptr);
    if (nsuperinterfaces>0){
        assert(si != NULL);
        for (i=0; i<nsuperinterfaces; i++) {
            ssuperinterface_free(&(si[i]));
        }
        si = NULL; /* To prevent a memory leak */
        if (*superInterface_ptr!=NULL) (*superInterface_ptr) = (SSUPER_INTERFACE *) tl_free(sizeof(SSUPER_INTERFACE), nsuperinterfaces, (*superInterface_ptr));
    }
}

/***********************************************************/
/***********************************************************/
// give this routine ONE supermodel pointer
void ssuperModel_printFmap(SSUPER_MODEL *sm) {
    
    int i,isubmodel;
    SMODEL *mod = NULL;
    
    FILE *fp = stdout;
#ifdef _MESSG
    int ierr_code, myid=0;
    ierr_code = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    char fn[30+1];
    snprintf(fn, 30, "model_fmaps_world_pe%d.txt", myid);
    fp = fopen(fn, "w");
#endif
    fflush(fp);
    for (isubmodel=0;isubmodel<sm->nsubmodels;isubmodel++){
        if(sm->submodel[isubmodel].proc_flag==1){
            mod = &(sm->submodel[isubmodel]);
            fprintf(fp,"model: %d \n",isubmodel);
            for (i=0; i<mod->grid->nnodes; i++) {
                fprintf(fp,"fmap[%d] = %d \n",i,mod->fmap[i]);
            }
            fprintf(fp,"\n");
        }
    }
#ifdef _MESSG
    fclose(fp);
#endif
}

void ssuperModel_printScreen(SSUPER_MODEL *sm, SSUPER_INTERFACE *si, int nsupermodels, int nsuperinterfaces) {
    
    FILE *fp = stdout;
#ifdef _MESSG
    int ierr_code, myid=0;
    ierr_code = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    char fn[30+1];
    snprintf(fn, 30, "superModel_info_world_pe%d.dat", myid);
    fp = fopen(fn, "w");
#endif
    
    int i,j,k,l;
    fprintf(fp,"SUPERFILE COUPLING \n");
    fprintf(fp,"----------------------\n");
    fprintf(fp,"Simulation has a total of %d supermodels\n",nsupermodels);
    for (i=0; i<nsupermodels; i++) {
        fprintf(fp,"**********************************************************\n");
        fprintf(fp,"Super Model: %d consists of projects: \n",i);
        for (j=0; j<sm[i].nsubmodels; j++) {
            fprintf(fp,"\tProject %i --> %s\n",j,sm[i].submodel[j].name);
        }
        
        fprintf(fp,"\tmax_nsys: %d \t max_nsys_sq: %d \n",sm[i].max_nsys,sm[i].max_nsys_sq);
        fprintf(fp,"\tmy_nnodes:               %d \n",sm[i].my_nnodes);
        fprintf(fp,"\tnnodes:                  %d \n",sm[i].nnodes);
        fprintf(fp,"\tmacro_nnodes:            %d \n",sm[i].macro_nnodes);
        fprintf(fp,"\twvel_my_nnodes:          %d \n",sm[i].wvel_my_nnodes);
        fprintf(fp,"\twvel_nnodes:             %d \n",sm[i].wvel_nnodes);
        fprintf(fp,"\twvel_macro_nnodes:       %d \n",sm[i].wvel_macro_nnodes);
        
        //        printf("\tmacro_nnodes_original:   %d \n",sm[i].macro_nnodes_original);
        fprintf(fp,"\ttotal_my_nnodes:         %d \n",sm[i].total_my_nnodes);
        fprintf(fp,"\ttotal_nnodes:            %d \n",sm[i].total_nnodes);
        fprintf(fp,"\ttotal_macro_nnodes:      %d \n",sm[i].total_macro_nnodes);
        //        printf("\tnnodes_matrix:           %d \n",sm[i].nnodes_matrix);
        //        printf("\ttotal_nnodes_matrix:     %d \n",sm[i].total_nnodes_matrix);
        
        fprintf(fp,"\tmax_nonlin_it:           %d \n",sm[i].max_nonlin_it);
        fprintf(fp,"\tinc_nonlin:              %.15e \n",sm[i].inc_nonlin);
        fprintf(fp,"\ttol_nonlin:              %.15e \n",sm[i].tol_nonlin);
        fprintf(fp,"\tnblock:                  %d \n",sm[i].nblock);
        
        fprintf(fp,"\tntransport:              %d \n",sm[i].ntransport);
        //        printf("\tnonlinear_it_total: %d \n",sm[i].nonlinear_it_total);
        //        printf("\tnonlinear_it_total_hvel: %d \n",sm[i].nonlinear_it_total_hvel);
        //        printf("\tnonlinear_it_total_wvel: %d \n",sm[i].nonlinear_it_total_wvel);
        
        if (sm[i].NumInterfaces > 0){
            fprintf(fp,"Super Model: %d consists of the following submodel interfaces: \n",i);
            for (j=0; j<sm[i].NumInterfaces; j++) {
                fprintf(fp,"*******************************************\n");
                fprintf(fp,"Interface between submodels %i  and %i contains %i vertical columns of coupled nodes\n",
                        sm[i].interface[j].model_id[0]+1,
                        sm[i].interface[j].model_id[1]+1,
                        sm[i].interface[j].NumNodeColumns);
                for (k=0; k<sm[i].interface[j].NumNodeColumns; k++) {
                    if (sm[i].interface[j].nodelist[k].size[0] > 0){
                        fprintf(fp,"Model %i node(s)[GID]: ", sm[i].interface[j].model_id[0]+1);
                        for(l=0;l<sm[i].interface[j].nodelist[k].size[0];l++){
                            fprintf(fp,"%10i[%i] ",sm[i].interface[j].nodelist[k].couplednodes[0][l]+1,sm[i].submodel[sm[i].interface[j].model_id[0]].grid->node[sm[i].interface[j].nodelist[k].couplednodes[0][l]].gid);
                        }
                        fprintf(fp,"\tcoupled to Model %i node(s)[GID]:", sm[i].interface[j].model_id[1]+1);
                        for(l=0;l<sm[i].interface[j].nodelist[k].size[1];l++){
                            fprintf(fp,"%10i[%i] ",sm[i].interface[j].nodelist[k].couplednodes[1][l]+1,sm[i].submodel[sm[i].interface[j].model_id[1]].grid->node[sm[i].interface[j].nodelist[k].couplednodes[1][l]].gid );
                        }
                        fprintf(fp,"\n");
                    }
                }
            }
            fprintf(fp,"*******************************************\n");
        }
        else{
            fprintf(fp,"There are no submodel-to-submodel interfaces in this project specified by the user\n");
        }
        
#ifdef _DEBUG
        if (DEBUG_FMAP){
            fprintf(fp,"The forward map from submodel nodes to supermodel HVEL system of equations is given below:\n");
            for (j=0; j<sm[i].nsubmodels; j++){
                if(sm[i].submodel[j].proc_flag==1){
                    fprintf(fp,"Map for Submodel[%i]\n\t", j);
                    for (k=0; k<sm[i].submodel[j].grid->nnodes - 6; k=k + 6){
                        for (l=0; l < 6; l++){
                            fprintf(fp,"%12i", sm[i].submodel[j].fmap[k+l]);
                        }
                        fprintf(fp,"\n\t");
                    }
                    for (l=k; l <sm[i].submodel[j].grid->nnodes ; l++){
                        fprintf(fp,"%12i", sm[i].submodel[j].fmap[l]);
                    }
                    fprintf(fp,"\n");
                }
            }
            fprintf(fp,"\n");
            
            fprintf(fp,"The forward map from submodel nodes to supermodel WVEL system of equations is given below:\n");
            for (j=0; j<sm[i].nsubmodels; j++){
                if(sm[i].submodel[j].proc_flag==1){
                    if (sm[i].submodel[j].flag.SW2_FLOW){
                        fprintf(fp,"Submodel[%i] is a 2D model; no map required for it for WVEL.\n", j);
                    }
                    else if (sm[i].submodel[j].flag.SW3_FLOW){
                        fprintf(fp,"Map for Submodel[%i]\n\t", j);
                        for (k=0; k<sm[i].submodel[j].grid->nnodes - 6; k=k + 6){
                            for (l=0; l < 6; l++){
                                fprintf(fp,"%12i", sm[i].submodel[j].fmap_wvel[k+l]);
                            }
                            fprintf(fp,"\n\t");
                        }
                        for (l=k; l <sm[i].submodel[j].grid->nnodes ; l++){
                            fprintf(fp,"%12i", sm[i].submodel[j].fmap_wvel[l]);
                        }
                        fprintf(fp,"\n");
                    }
                }
            }
            fprintf(fp,"\n");
        }
#endif
    }
    fprintf(fp,"**********************************************************\n");
    if (nsuperinterfaces > 0) {
        fprintf(fp,"This simulation has a total of %d superinterfaces\n",nsuperinterfaces);
        for (i=0; i<nsuperinterfaces; i++) {
            int sm1 = si[i].sm_id[0];
            int sm2 = si[i].sm_id[1];
            int mod1 = si[i].model_id[0];
            int mod2 = si[i].model_id[1];
            fprintf(fp,"**********************************************************\n");
            fprintf(fp,"Superinterface %d is the interface between:\n", i+1);
            fprintf(fp,"        supermodel %i, submodel[%i]\n",sm1+1,mod1+1);
            fprintf(fp,"        supermodel %i, submodel[%i]\n",sm2+1,mod2+1);
            fprintf(fp,"***************************************\n");
            fprintf(fp,"This superinterface contains %i vertical columns of coupled edges\n", si[i].NumEdges);
            //for (k=0; k<si[i].NumEdges; k++) {
            //    if (si[i].bdrylist[k].size[0] > 0){
            //        printf("SupMod %i SubMod %i edge/face(s): ", sm1+1, mod1+1);
            //        for(l=0;l<si[i].bdrylist[k].size[0];l++){
            //            printf("%10i ",si[i].bdrylist[k].couplededges[0][l]+1);
            //        }
            //        printf("\tcoupled to SupMod %i SubMod %i edge/face(s):", sm2+1, mod2+1);
            //        for(l=0;l<si[i].bdrylist[k].size[1];l++){
            //            printf("%10i ",si[i].bdrylist[k].couplededges[1][l]+1);
            //        }
            //        printf("\n");
            //    }
            //}
        }
    }
    else {
        fprintf(fp,"There are no supermodel-to-supermodel interfaces specified by the user\n");
    }
    
    fprintf(fp,"**********************************************************\n");
    
    
#ifdef _MESSG
    fclose(fp);
#endif
}

/***********************************************************/
/***********************************************************/

void ssuperModel_readSuperfile(SSUPER_MODEL **superModel_ptr, SSUPER_INTERFACE **superInterface_ptr, int *nsupmods, int *nsupifcs) {
    
    int DEBUG_READ = OFF;
    
    int i,j,id1,id2;
    FILE *fp;
    char name[MAXLINE];
    char line[MAXLINE];
    char *data = NULL, *subdata = NULL, *subsubdata = NULL;
    CARD card, sub_card;
    int dum_int,super_procs,*super_procs_wgt,*sub_procs,**sub_procs_wgt, query;
    int sum1,sum2, val[2], errno, token_len;
    const char *starting;
    char *end;
    SIO dummy;
    double tempdouble;
    SMODEL *model;
    
    int npes = 1, myid = 0, ierr = UNSET_INT;;
#ifdef _MESSG
    ierr = MPI_Comm_size(cstorm_comm, &npes);
    ierr = MPI_Comm_rank(cstorm_comm, &myid);
#endif
    
#ifdef _DEBUG
    if (myid == 0) printf("-----------------------------------------\n");
    if (myid == 0) printf("Reading Superfile -----------------------\n");
    if (myid == 0) printf("-----------------------------------------\n");
#endif
    
    // make sure file exists
    fp = fopen("superfile.in","r");
    if (fp == NULL) {
        tl_error(">> superfile.in does not exist.");
    }
    
    // read the number of super-models
    rewind(fp);
    int nsupermodelcount = 0, nsuperinterfacecount = 0;
    while (fgets(line, MAXLINE, fp) != NULL) {
        if (strip_comments(line) <= 1) {
            continue;
        }
        card = parse_card(line, &data);
        switch (card) {
            case CARD_SMODEL:
                nsupermodelcount++;
                break;
            case CARD_SUPIFC:
                nsuperinterfacecount++;
                break;
        }
    }
    
    (* nsupmods) = nsupermodelcount;
    (* nsupifcs) = nsuperinterfacecount;
    int nsupermodels = nsupermodelcount;
    int nsuperinterfaces = nsuperinterfacecount;
    
    if (myid==0) printf("-- # of superModels: %d\n",nsupermodels);
    if (myid==0) printf("-- # of flux interfaces between superModels: %d\n",nsuperinterfaces);
    
    // CJT :: What are these???
    super_procs_wgt = (int *) tl_alloc(sizeof(int), nsupermodels);
    sub_procs       = (int *) tl_alloc(sizeof(int), nsupermodels);
    sub_procs_wgt   = (int **) tl_alloc(sizeof(int *), nsupermodels);
    
    // subModel interfaces
    int *nsubmodels, *nInterfaces;
    nsubmodels  = (int *) tl_alloc(sizeof(int) , nsupermodels);  // memory freed at the end of this file
    nInterfaces = (int *) tl_alloc(sizeof(int) , nsupermodels);  // memory freed at the end of this file
    sarray_init_int(nsubmodels,nsupermodels);
    sarray_init_int(nInterfaces,nsupermodels);
    
    // superModel interfaces
    SSUPER_INTERFACE *si=NULL;
    if (nsuperinterfaces > 0){
        (*superInterface_ptr) = (SSUPER_INTERFACE *) tl_alloc(sizeof(SSUPER_INTERFACE), nsuperinterfaces);
        si = (*superInterface_ptr); /* alias */
    }
    
    rewind(fp);
    nsupermodelcount = 0;
    nsuperinterfacecount = 0;
    while (fgets(line, MAXLINE, fp) != NULL) {
        if (strip_comments(line) <= 1) {
            continue;
        }
        card = parse_card(line, &data);
        switch (card) {
            case CARD_SMODEL:
                
                //  for this submodel, read the number of submodels
                i=0; j=0;
                val[0]=1; // the # of submodels
                val[1]=1; // this # is MPI weights of the super models (I think)
                dum_int=0;
                for(subdata=data;*subdata!=0;subdata+=j){
                    j=1;
                    if(isdigit(*subdata)){
                        sscanf(subdata, "%d%n", &dum_int, &j);
                        val[i] = dum_int;
                        i++;
                    }
                }
                nsubmodels[nsupermodelcount] = val[0];
                
                // cjt :: not sure what this is for, but not used for example models
                super_procs_wgt[nsupermodelcount] = val[1];
                sub_procs_wgt[nsupermodelcount] = (int *) tl_alloc(sizeof(int ), nsubmodels[nsupermodelcount]);
                for(i=0;i<nsubmodels[nsupermodelcount];i++){
                    sub_procs_wgt[nsupermodelcount][i] = 1;
                    fgets(line,MAXLINE,fp);
                    j=0;
                    dum_int=0;
                    for(subdata=line;*subdata!=0;subdata+=j){
                        j=1;
                        if(isspace(*subdata) && (isdigit(*(subdata+1)))){
                            sscanf(subdata, "%d%n", &dum_int, &j);
                            sub_procs_wgt[nsupermodelcount][i] = dum_int;
                            printf("sub_procs_wgt[nsupermodelcount=%d][%d]: %d\n",nsupermodelcount,i,sub_procs_wgt[nsupermodelcount][i]);
                        }
                    }
                }
                
                
                if (nsubmodels[nsupermodelcount] < 1) {
                    tl_error(">> number of submodels cannot be 0 or negative.");
                }
                nsupermodelcount++;
                break;
                
            case CARD_SUPIFC:
                // read flux interface from one superModel to another
                sscanf(data,"%d%d",&(si[nsuperinterfacecount].sm_id[0]) ,
                       &(si[nsuperinterfacecount].sm_id[1]) );
                si[nsuperinterfacecount].sm_id[0]--;
                si[nsuperinterfacecount].sm_id[1]--;
                printf("\n\nInterface %i: sm_id[0] = %i, sm_id[1] = %i\n\n", nsuperinterfacecount, si[nsuperinterfacecount].sm_id[0],si[nsuperinterfacecount].sm_id[1]);
                if ((si[nsuperinterfacecount].sm_id[0]<0) || (si[nsuperinterfacecount].sm_id[1] < 0)) {
                    tl_error(">> supermodel id's of superinterfaces cannot be 0 or negative.");
                }
                nsuperinterfacecount++;
                break;
                
            case CARD_NUMIFC:
                // read monolithic interface between two submodels of a superModel
                sscanf(data,"%d",&(nInterfaces[nsupermodelcount-1]));
                
                if (nInterfaces[nsupermodelcount-1]<(nsubmodels[nsupermodelcount-1]-1)){
                    tl_error(">> number of interfaces cannot be less than (num_submodels -1)");
                }
                if (nInterfaces[nsupermodelcount-1]>((nsubmodels[nsupermodelcount-1]*(nsubmodels[nsupermodelcount-1]-1))/2)){
                    tl_error(">> number of interfaces cannot be greater than (num_submodels*(num_submodels-1)/2)");
                }
                break;
                // GAJANAN above
        }
    }
    // allocate coupler struct
#ifdef _DEBUG
    if (DEBUG_READ) {
        for(i=0;i<nsupermodels;i++){
            printf("-- pe: %d allocating superModel: %d with %d subModels and %d subModel interfaces\n",myid,i,nsubmodels[i],nInterfaces[i]);
        }
#ifdef _MESSG
        fflush(stdout); messg_barrier(cstorm_comm);
#endif
    }
#endif
    ssuperModel_alloc_init(superModel_ptr, nsupermodels, nsubmodels, nInterfaces);
    SSUPER_MODEL *sm = (* superModel_ptr); /* alias */
    
#ifdef _MESSG
    // CJT NOTES
    // sm[i].proc_flag :: a binary flag determining if a PE works on the ith superModel
    // sm[i].subModel[j].proc_flag :: a binary flag determining if a PE works on the jth subModel of the ith superModel
    
    
    super_procs=0;
    for(i=0;i<nsupermodelcount;i++){
        super_procs=super_procs + super_procs_wgt[i];
        sub_procs[i]=0;
        for(j=0;j<nsubmodels[i];j++){
            sub_procs[i] = sub_procs[i] + sub_procs_wgt[i][j];
        }
    }
    

    int myid2, myid3;
    int color = UNSET_INT;
    int color2 = UNSET_INT;
    MPI_Comm super_comm;
    MPI_Comm sub_comm;
    
    if (npes > 1) {
        sum1=0;
        for(i=0;i<nsupermodelcount;i++) {
            super_procs_wgt[i] = (int)((float)(npes)*((float)(super_procs_wgt[i])/(float)(super_procs)));
            sum1=sum1+super_procs_wgt[i];
            if((i==(nsupermodelcount-1)) && (sum1<npes)){
                super_procs_wgt[i] = super_procs_wgt[i] + (npes - sum1);
                sum1=npes;
            }
            if((myid<sum1)&&(color==UNSET_INT)){
                color=i;
            }
            sum2=0;
            if(myid==0){
                printf("-- supermodel %d is using %d procs of %d \n",i,super_procs_wgt[i],npes);
                if(super_procs_wgt[i] < 1) tl_error("Bad ratio or not enough processors to run all supermodels");
            }
            for (j=0;j<nsubmodels[i];j++) {
                sub_procs_wgt[i][j] = (int)((float)super_procs_wgt[i]*((float)sub_procs_wgt[i][j]/(float)sub_procs[i]));
                sum2 = sum2 + sub_procs_wgt[i][j];
                if((j==(nsubmodels[i]-1)) && (sum2<super_procs_wgt[i])) sub_procs_wgt[i][j] = sub_procs_wgt[i][j] + (super_procs_wgt[i] - sum2);
                if(myid==0){
                    printf("----- submodel %d is using %d procs of %d \n",j,sub_procs_wgt[i][j],super_procs_wgt[i]);
                    if(sub_procs_wgt[i][j] < 1) tl_error("Bad ratio or not enough processors to run all submodels");
                }
            }
        }
        
        //printf("myid: %d color: %d\n",myid,color);
        
        // cjt :: note :: if only one supermodel is used, color is all 0
        MPI_Comm_split(cstorm_comm, color, myid, &(super_comm));
        ierr  = MPI_Comm_rank(super_comm,&myid2);
        
        for(i=0;i<nsupermodelcount;i++) {
            smpi_init(sm[i].supersmpi, super_comm);
            if(color==i){
                sum2=0;
                for (j=0;j<nsubmodels[i];j++) {
                    sum2 = sum2 + sub_procs_wgt[i][j];
                    if((j==(nsubmodels[i]-1)) && (sum2<super_procs_wgt[i])) sum2 = super_procs_wgt[i];
                    if((myid2<sum2)&&(color2==UNSET_INT)) {
                        color2=j;
                    }
                }
                //printf("color2: %d myid2: %d\n",color2,myid2);
                MPI_Comm_split(super_comm, color2, myid2, &sub_comm);
            }
        }
        ierr  = MPI_Comm_rank(sub_comm,&myid3);
        
    } else {
        printf("One processor running all super and submodels\n");
        ierr = MPI_Comm_dup(cstorm_comm,&super_comm);
        ierr = MPI_Comm_dup(cstorm_comm, &sub_comm);
        for(i=0;i<nsupermodelcount;i++) {
            smpi_init(sm[i].supersmpi, super_comm);
            smpi_init(sm[i].supersmpi_wvel, super_comm);
        }
    }
    
    for(i=0;i<nsupermodelcount;i++) {
        for (j=0;j<nsubmodels[i];j++) {
            MPI_Comm_dup(sub_comm, &(sm[i].submodel[j].model_comm));
            
            int myidG,npesG,myidS,npesS,ierr_code;
            ierr_code = MPI_Comm_rank(cstorm_comm, &myidG);
            ierr_code = MPI_Comm_size(cstorm_comm, &npesG);
            ierr_code = MPI_Comm_rank(sm[i].submodel[j].model_comm,    &myidS);
            ierr_code = MPI_Comm_size(sm[i].submodel[j].model_comm,    &npesS);
            //printf("superModel: %d submodel: %d :: world pe %d of %d :: model pe %d of %d\n",i,j,myidG+1,npesG,myidS+1,npesS);
            
        }
    }
    //fflush(stdout); messg_barrier(cstorm_comm); //tl_error("for now");
    
    if (npes > 1) {
        for (i=0; i<nsupermodels; i++) {
            if(color==i){
                sm[i].proc_flag=1;
                for (j=0; j<nsubmodels[i]; j++) {
                    if(color2==j)sm[i].submodel[j].proc_flag=1;
                }
            }else{
                for (j=0; j<nsubmodels[i]; j++) {
                    sm[i].submodel[j].proc_flag=0;
                }
            }
        }
    }else{
        for(i=0;i<nsupermodelcount;i++) {
            sm[i].proc_flag=1;
            for(j=0;j<nsubmodels[i];j++){
                sm[i].submodel[j].proc_flag=1;
            }
        }
    }
    
#else
    /* tells serial proces to run all models */
    for(i=0;i<nsupermodelcount;i++) {
        smpi_init(sm[i].supersmpi);
        smpi_init(sm[i].supersmpi_wvel);
        sm[i].proc_flag=1;
        for(j=0;j<nsubmodels[i];j++){
            sm[i].submodel[j].proc_flag=1;
        }
    }
    
#endif
    
#ifdef _DEBUG
    if (DEBUG_READ) {
        for(i=0;i<nsupermodelcount;i++) {
            printf("-- myid: %d :: superModel[%d] proc flag: %d \n",myid,i,sm[i].proc_flag);
            for(j=0;j<nsubmodels[i];j++){
                printf("---- myid: %d :: subModel[%d] proc flag: %d \n",myid,j,sm[i].submodel[j].proc_flag);
            }
        }
#ifdef _MESSG
        fflush(stdout); messg_barrier(cstorm_comm);
#endif
    }
#endif
    
    // now read main sub-model and initialize
    rewind(fp);
    nsupermodelcount = 0;
    int super,sub;
    while (fgets(line, MAXLINE, fp) != NULL) {
        if (strip_comments(line) <= 1) {
            continue;
        }
        card = parse_card(line, &data);
        switch (card) {
            case CARD_SMODEL:
                sscanf(data,"%d",&id1);
                if(sm[nsupermodelcount].proc_flag==1){
                    super=nsupermodelcount;
                    for (i=0; i<id1; i++) {
                        fgets(line,MAXLINE,fp);
                        sscanf(line,"%s",name);
                        //here the model being worked on is read first to get any interface nodes
                        if(sm[nsupermodelcount].submodel[i].proc_flag==1){
                            printf("\n------------------------------------------------------------------------\n");
                            sub=i;
                            smodel_init(&(sm[nsupermodelcount].submodel[i]),name);
                            printf("------------------------------------------------------------------------\n");
                            printf("------------------------------------------------------------------------\n");
                        }
                    }
                }
                nsupermodelcount++;
                break;
        }
    }
#ifdef _MESSG
    fflush(stdout); messg_barrier(cstorm_comm); messg_barrier(cstorm_comm); //tl_error("for now");
#endif

    // now read other sub-models and initialize
    rewind(fp);
    nsupermodelcount = 0;
    nsuperinterfacecount = 0;
    int model_2d_id, model_3d_id, NumNodeColumns, nInterfacecount;
    int surfnode[2];
    while (fgets(line, MAXLINE, fp) != NULL) {
        if (strip_comments(line) <= 1) {
            continue;
        }
        card = parse_card(line, &data);
        switch (card) {
            case CARD_SMODEL:
                //                sscanf(data,"%d",&id1);
                //                if(sm[nsupermodelcount].proc_flag==1){
                //                    super=nsupermodelcount;
                //                    for (i=0; i<id1; i++) {
                //                        fgets(line,MAXLINE,fp);
                //                        sscanf(line,"%s",name);
                //                        printf("SuperModel Card Found :: SMODEL %s\n",name);
                //                        //here the model being worked on is read first to get any interface nodes
                //                        if(sm[nsupermodelcount].submodel[i].proc_flag==1){
                //                            sub=i;
                //                            smodel_init( &(sm[nsupermodelcount].submodel[i]), name );
                //                        }
                //                    }
                //        }
                
                
                sscanf(data,"%d",&id1);
                for (i=0; i<id1; i++) {
                    fgets(line,MAXLINE,fp);
                    sscanf(line,"%s",name);

                    /*only inittialize models with 0 resident nodes.
                     * smodel_init will only set BC and time series for these models
                     * this way if adaption sends interface or resident nodes for these models,
                     * the parameters are set up */
                    if((sm[nsupermodelcount].proc_flag==1) && (sm[nsupermodelcount].submodel[i].proc_flag==0)){
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d SuperModel Card Found :: SMODEL %s :: SUB_MODEL_PROC_FLAG = 0\n",myid,name);
#endif
                        model = &(sm[nsupermodelcount].submodel[i]);
                        smodel_init(model,name);
                    }
                    if(sm[nsupermodelcount].proc_flag==0){ // cjt :: this is another superModel not associated with my PE
                        strcpy(model->name,name);
                        model->io = (SIO *) tl_alloc(sizeof(SIO), 1);
                        sio_init(model->io, name);
                    }
                    //tl_check_all_pickets(__FILE__,__LINE__);
                }
               
                nsupermodelcount++;
                nInterfacecount = 0;
                
                break;
               
                // now read nonlinear/linear solver information for this supermodel
            case CARD_OP:
                switch (parse_card(data, &subdata)) {
                    case CARD_PRE:
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: OP PRE\n",myid,name);
#endif
                        sscanf(subdata,"%d", &(sm[nsupermodelcount-1].solver_info.prec_value));
                        break;
                        
                    case CARD_BLK:
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: OP BLK\n",myid,name);
#endif
                        sscanf(subdata,"%d", &(sm[nsupermodelcount-1].nblock));
                        if (sm[nsupermodelcount-1].nblock != 1) {
                            tl_error("BLK must be set to 1 to currently run coupled models.\n");
                        }
                        break;
                        
                    case CARD_INC:
                        sscanf(subdata,"%d", &(sm[nsupermodelcount-1].nalloc_inc));
                        break;
                        
                    default:
                        break;
                }
                break;
            case CARD_IP:
                switch (parse_card(data, &subdata)) {
                        
                    case CARD_MIT:
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: IP BLK\n",myid,name);
#endif
                        sscanf(subdata,"%d", &(sm[nsupermodelcount-1].solver_info.max_lin_it));
                        sm[nsupermodelcount-1].solver_info.force_lin_it = NO;
                        break;
                        
                    case CARD_NIT:
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: IP MIT\n",myid,name);
#endif
                        sscanf(subdata,"%d", &(sm[nsupermodelcount-1].max_nonlin_it));
                        sm[nsupermodelcount-1].solver_info.force_nonlin_it = NO;
                        break;
                        
                    case CARD_ITL:
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: IP ITL\n",myid,name);
#endif
                        sscanf(subdata,"%lf", &(sm[nsupermodelcount-1].inc_nonlin));
                        break;
                        
                    case CARD_NTL:
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: IP NTL\n",myid,name);
#endif
                        sscanf(subdata,"%lf", &(sm[nsupermodelcount-1].tol_nonlin));
                        break;
                        
                    default:
                        break;
                }
                break;
                // GAJANAN below
            case CARD_TC:
                switch (parse_card(data, &subdata)) {
                        
                    case CARD_T0:
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: TC T0\n",myid,name);
#endif

                        sscanf(subdata,"%lf", &(tempdouble));
                        for (i=0;i<sm[nsupermodelcount-1].nsubmodels; i++){
                            sm[nsupermodelcount-1].submodel[i].t_prev = tempdouble;
                            sm[nsupermodelcount-1].submodel[i].t_init = tempdouble;
                        }
                        break;
                        
                    case CARD_TF:
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: TC TF\n",myid,name);
#endif
                        sscanf(subdata,"%lf", &(tempdouble));
                        for (i=0;i<sm[nsupermodelcount-1].nsubmodels; i++){
                            sm[nsupermodelcount-1].submodel[i].t_final = tempdouble;
                        }
                        break;
                        
                    case CARD_DT:
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: TC DT\n",myid,name);
#endif
                        sscanf(subdata,"%lf", &(tempdouble));
                        sm[nsupermodelcount-1].dt      = tempdouble;
                        sm[nsupermodelcount-1].dt_prev = tempdouble;
                        printf("\n\nWarning: Overriding DT series values for all submodels with the Supermodel DT supplied in superfile.in. DT Value: %14.6e\n\n", tempdouble);
                        SSERIES *series_dt;
                        for (i=0;i<sm[nsupermodelcount-1].nsubmodels; i++){
                            sm[nsupermodelcount-1].submodel[i].dt = tempdouble;
                            series_dt = sm[nsupermodelcount-1].submodel[i].series_dt;
                            for (j=0; j<series_dt->size; j++){
                                series_dt->entry[j].value[0] = tempdouble;
                            }
                        }
                        break;
                        
                    default:
                        break;
                }
                break;
            case CARD_TEST:
                card = parse_card(data, &subdata);
                switch (card) {
                    case CARD_TIDE23:
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: TEST TIDE23\n",myid,name);
#endif
                        data = subdata;
                        test_case_flag.model_id_2d3d = (int) strtol(data,&subdata,10);
                        test_case_flag.model_id_2d3d--;
                        break;
                    case CARD_SLSH23:
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: TEST SLSH23\n",myid,name);
#endif
    
                        data = subdata;
                        test_case_flag.model_id_2d3d = (int) strtol(data,&subdata,10);
                        test_case_flag.model_id_2d3d--;
                        break;
                    default:
                        tl_error("Only TEST SLSH23 (2D-3D slosh) and TEST TIDE23 (2D-3D tide) are currently supported");
                        break;
                }
                break;
                
            case CARD_CN:
#ifdef _DEBUG
                if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: CN \n",myid,name);
#endif
                sm[nsupermodelcount-1].ntransport = (int) strtol(data,&subdata,10); // number of transport constituents
                sm[nsupermodelcount-1].con_type = (int *) tl_alloc(sizeof(int),sm->ntransport);
                break;
                
            case CARD_SUPIFC:
                fgets(line,MAXLINE,fp);
                int sm_id[2], model_id[2], edgestring[2];
                sm_id[0] = si[nsuperinterfacecount].sm_id[0];
                sm_id[1] = si[nsuperinterfacecount].sm_id[1];
                if (DEBUG) printf("\n\nInterface %i: sm_id[0] = %i, sm_id[1] = %i\n\n", nsuperinterfacecount, si[nsuperinterfacecount].sm_id[0],si[nsuperinterfacecount].sm_id[1]);
                card = parse_card(line,&data);
                switch(card){
                    case CARD_IFCSM:
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: SUPIFC IFCSM\n",myid,name);
#endif
                        sscanf(data,"%d %d %d", &(model_id[0]), &(model_id[1]), &NumNodeColumns);
                        model_id[0]--;
                        model_id[1]--;
                        si[nsuperinterfacecount].model_id[0] = model_id[0];
                        si[nsuperinterfacecount].model_id[1] = model_id[1];
                        si[nsuperinterfacecount].NumEdges = NumNodeColumns;
                        if (DEBUG) printf("\n\nInterface %i: model_id[0] = %i, model_id[1] = %i\n\n", nsuperinterfacecount, si[nsuperinterfacecount].model_id[0],si[nsuperinterfacecount].model_id[1]);
                        if ( model_id[0]>=sm[sm_id[0]].nsubmodels || model_id[1]>=sm[sm_id[1]].nsubmodels ){
                            printf("\nError (in superinterface[%i])", nsuperinterfacecount);
                            tl_error(">> Error: Unacceptable model id's supplied. Check card INTFCE under the SUPIFC card in the superfile.");
                        }
                        break;
                    default:
                        break;
                }
                fgets(line,MAXLINE,fp);
                card = parse_card(line,&data);
                switch(card){
                    case CARD_CPL:
#ifdef _DEBUG
                        if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: SUPIFC CPL\n",myid,name);
#endif
                        sscanf(data,"%d %d", &(edgestring[0]),&(edgestring[1]));
                        edgestring[0]--;
                        edgestring[1]--;
                        if (DEBUG) printf("\n\n\nEdgestrings coupled: %i and %i\n\n", edgestring[0], edgestring[1]);
                        ssuperinterface_alloc_init( &(si[nsuperinterfacecount]), NumNodeColumns, sm_id, model_id, edgestring);
                        //exit(0);
                        break;
                    default:
                        break;
                }
                nsuperinterfacecount++;
                //printf("\n\nsuperinterface %i:",nsuperinterfacecount);
                //printf("\nsupermodel %i flux-coupled to supermodel %i", sm_2d_id+1, sm_3d_id+1);
                //printf("\nsubmodel %i flux-coupled to submodel %i", model_2d_id+1, model_3d_id+1);
                //exit(0);
                break;
                
            case CARD_INTFCE:
#ifdef _DEBUG
                if (DEBUG) printf("GLOBAL MYPE: %d :: SMODEL %s :: Reading in Superfile card :: INTFCE\n",myid,name);
#endif
                if (nInterfacecount > sm[nsupermodelcount-1].NumInterfaces-1){
                    printf("\nError (in supermodel[%i]): More interfaces supplied than specified (%i).\n", nsupermodelcount-1, sm[nsupermodelcount-1].NumInterfaces);
                    tl_error(">> Too many interfaces. Check cards INTFCE and NUMIFC in the input file.");
                }
                
                sscanf(data,"%d %d %d", &model_2d_id, &model_3d_id, &NumNodeColumns);
                model_2d_id--;
                model_3d_id--;
                
                if ( model_2d_id>=sm[nsupermodelcount-1].nsubmodels || model_3d_id>=sm[nsupermodelcount-1].nsubmodels ){
                    int dummyid = (model_2d_id > model_3d_id ? (model_2d_id+1) : (model_3d_id+1)) ;
                    printf("\nError (in supermodel[%i]): Submodel ID %i exceeds the number of submodels (%i)\n", nsupermodelcount-1, dummyid, sm[nsupermodelcount-1].nsubmodels);
                    tl_error(">> Error: Unacceptable model id's supplied. Check card INTFCE in the input file.");
                }
                
                sinterface_alloc_init( &(sm[nsupermodelcount-1].interface[nInterfacecount]), NumNodeColumns, model_2d_id, model_3d_id);
                
                /* The following lines read coupled surface node numbers manually supplied by the user */
                // for (i=0; i<NumNodeColumns; i++) {
                //     fgets(line,MAXLINE,fp);
                //     /* line format is:    <surfnode[0]>  <surfnode[1]> */
                //     data = line;
                //     surfnode[0] = (int) strtol(data,&subdata,10);
                //     surfnode[0]--; // Need --
                //     data = subdata; /* resetting the pointer */
                
                //     surfnode[1]  = (int) strtol(data,&subdata,10);
                //     surfnode[1]--; // Need --
                
                //     if(surfnode[0]>=sm[nsupermodelcount-1].submodel[model_2d_id].grid->macro_nnodes ||
                //        surfnode[1]>=sm[nsupermodelcount-1].submodel[model_3d_id].grid->macro_nnodes ){
                //         tl_error("Error: Node id exceeds the number of nodes in some model. Check interface input data.");
                //     }
                //     sm[nsupermodelcount-1].interface[nInterfacecount].nodelist[i].surfnode[0] = surfnode[0];
                //     sm[nsupermodelcount-1].interface[nInterfacecount].nodelist[i].surfnode[1] = surfnode[1];
                // }
                
                nInterfacecount++;
                break;
            default:
                //                tl_error("Please supply SMODEL, INTFCE, NUMIFC, PRE, MIT, NIT, ITL, and NTL cards in the super file.\n");
                break;
        }
    }
    rewind(fp);
#ifdef _MESSG
    fflush(stdout); messg_barrier(cstorm_comm); //tl_error("for now");
#endif
    int max_nsys=0;
    /* genereate monolithic couplic interfaces */
    for (i=0; i<nsupermodels; i++) {
        for (j=0; j<sm[i].NumInterfaces; j++) {
            int mod_id1 = sm[i].interface[j].model_id[0];
            int mod_id2 = sm[i].interface[j].model_id[1];
            if ( (sm[i].submodel[mod_id1].flag.SW2_FLOW) && (sm[i].submodel[mod_id2].flag.SW2_FLOW) ){ /*2D-2D*/
                generate_2d_2d_interface(&(sm[i]), &(sm[i].interface[j]), mod_id1, mod_id2, nsupermodels);
            }
            else if ( (sm[i].submodel[mod_id1].flag.SW3_FLOW) && (sm[i].submodel[mod_id2].flag.SW3_FLOW) ){ /*3D-3D*/
                generate_3d_3d_interface(&(sm[i]), &(sm[i].interface[j]), mod_id1, mod_id2, nsupermodels);
            }
            else if ( (sm[i].submodel[mod_id1].flag.SW3_FLOW) && (sm[i].submodel[mod_id2].flag.SW2_FLOW) ){ /*3D-2D*/
                sm[i].interface[j].model_id[0] = mod_id2; /* Flipping order from*/
                sm[i].interface[j].model_id[1] = mod_id1;
                mod_id1 = sm[i].interface[j].model_id[0];
                mod_id2 = sm[i].interface[j].model_id[1];
                generate_2d_3d_interface(&(sm[i]), &(sm[i].interface[j]), mod_id1, mod_id2, nsupermodels);
            }
            else if ( (sm[i].submodel[mod_id1].flag.SW2_FLOW) && (sm[i].submodel[mod_id2].flag.SW3_FLOW) ){ /*2D-3D*/
                generate_2d_3d_interface(&(sm[i]), &(sm[i].interface[j]), mod_id1, mod_id2, nsupermodels);
            }
        }
        if(sm[i].submodel[j].max_nsys > max_nsys){
            max_nsys = sm[i].submodel[j].max_nsys;
        }
        sm[i].max_nsys = max_nsys;
        sm[i].max_nsys_sq = max_nsys*max_nsys;
    }
    
    //tag();
    
#ifdef _MESSG
    messg_barrier(cstorm_comm);
    
    int k;
    MPI_Comm part_comm;
    /* fix interface ownership and set super/sub processor map */
    for (i=0; i<nsupermodels; i++) {
        if(sm[i].supersmpi->npes>1){
            fix_interface_ownership(&(sm[i]), nsupermodels);
            sm[i].proc_map = tl_alloc(sizeof(int *), sm[i].supersmpi->npes);
            for(j=0;j<sm[i].supersmpi->npes;j++){
                sm[i].proc_map[j] = tl_alloc(sizeof(int), sm[i].nsubmodels);
            }
            /* initialize the processor map */
            for(j=0;j<sm[i].supersmpi->npes;j++){
                for(k=0;k<sm[i].nsubmodels;k++){
                    sm[i].proc_map[j][k] = UNSET_INT;
                }
            }
            
            /* set myids in the map */
            j=sm[i].supersmpi->myid;
            if(j!=UNSET_INT){
                for(k=0;k<sm[i].nsubmodels;k++){
                    if(sm[i].submodel[k].proc_flag==1){
                        sm[i].proc_map[j][k] = sm[i].submodel[k].grid->smpi->myid;
                    }
                }
            }
            for(j=0;j<sm[i].supersmpi->npes;j++){
                ierr = MPI_Bcast(sm[i].proc_map[j], sm[i].nsubmodels, MPI_INT, j, sm[i].supersmpi->ADH_COMM);
            }
            if(sm[i].supersmpi_wvel !=NULL){
                sm[i].proc_map_wvel = tl_alloc(sizeof(int *), sm[i].supersmpi_wvel->npes);
                for(j=0;j<sm[i].supersmpi_wvel->npes;j++){
                    sm[i].proc_map_wvel[j] = tl_alloc(sizeof(int), sm[i].nsubmodels);
                }
                /* initialize the processor map */
                for(j=0;j<sm[i].supersmpi_wvel->npes;j++){
                    for(k=0;k<sm[i].nsubmodels;k++){
                        sm[i].proc_map_wvel[j][k] = UNSET_INT;
                    }
                }
                
                /* set myids in the map */
                j=sm[i].supersmpi_wvel->myid;
                if(j!=UNSET_INT){
                    for(k=0;k<sm[i].nsubmodels;k++){
                        if(sm[i].submodel[k].proc_flag==1){
                            sm[i].proc_map_wvel[j][k] = sm[i].submodel[k].grid->smpi->myid;
                        }
                    }
                }
                for(j=0;j<sm[i].supersmpi_wvel->npes;j++){
                    ierr = MPI_Bcast(sm[i].proc_map_wvel[j], sm[i].nsubmodels, MPI_INT, j, sm[i].supersmpi_wvel->ADH_COMM);
                }
            }
        }else{
            sm[i].proc_map = tl_alloc(sizeof(int *), sm[i].supersmpi->npes);
            sm[i].proc_map_wvel = tl_alloc(sizeof(int *), sm[i].supersmpi->npes);
            for(j=0;j<sm[i].supersmpi->npes;j++){
                sm[i].proc_map[j] = tl_alloc(sizeof(int), sm[i].nsubmodels);
                sm[i].proc_map_wvel[j] = tl_alloc(sizeof(int), sm[i].nsubmodels);
            }
            for(k=0;k<sm[i].nsubmodels;k++){
                sm[i].proc_map[0][k]=0;
                sm[i].proc_map_wvel[0][k]=0;
            }
        }
    }
#endif
#ifdef _MESSG
    for(i=0;i<nsupermodelcount;i++) {
        for (j=0;j<nsubmodels[i];j++) {
            MPI_Comm_dup(sub_comm, &(sm[i].submodel[j].model_comm));
            
            int myidG,npesG,myidS,npesS,ierr_code;
            ierr_code = MPI_Comm_rank(cstorm_comm, &myidG);
            ierr_code = MPI_Comm_size(cstorm_comm, &npesG);
            ierr_code = MPI_Comm_rank(sm[i].submodel[j].model_comm,    &myidS);
            ierr_code = MPI_Comm_size(sm[i].submodel[j].model_comm,    &npesS);
            //printf("after INTERFACING :: superModel: %d submodel: %d :: world pe %d of %d :: model pe %d of %d\n",i,j,myidG+1,npesG,myidS+1,npesS);
            
        }
    }
    fflush(stdout); messg_barrier(cstorm_comm); //tl_error("for now");
#endif
    
    /* determine flux coupling interfaces */
    int flux_interface=0, ie, istr;
    int **flux_submodels;
    flux_submodels  = (int **)tl_alloc(sizeof(int *), nsupermodels);
    
    for (i=0; i<nsupermodels; i++) {
        flux_submodels[i]  = (int *)tl_alloc(sizeof(int ), sm[i].nsubmodels);
        for(j=0;j<sm[i].nsubmodels;j++){
            flux_submodels[i][j]=0;
            if(sm[i].submodel[j].proc_flag==1){
                if(sm[i].submodel[j].grid->ndim==2){
                    for (ie=0; ie<sm[i].submodel[j].grid->nelems1d; ie++){
                        istr = sm[i].submodel[j].grid->elem1d[ie].string;
                        if (sm[i].submodel[j].str_values[istr].ol_flow.bc_flag == BCT_FLUX_COUPLE){
                            flux_interface=1;
                            flux_submodels[i][j]=1;
                        }
                        if(flux_interface >0) break;
                    }
                }else{
                    for (ie=0; ie<sm[i].submodel[j].grid->nelems2d; ie++){
                        istr = sm[i].submodel[j].grid->elem2d[ie].string;
                        if (sm[i].submodel[j].str_values[istr].flow.bc_flag == BCT_FLUX_COUPLE){
                            flux_interface=1;
                            flux_submodels[i][j]=1;
                        }
                        if(flux_interface >0) break;
                    }
                }
            }
        }
    }
    
    //tag();
    
    /* every processor needs submodel root world ids for flux broadcast */
    if(flux_interface>0){
        root_ids = (int **)tl_alloc(sizeof(int *), nsupermodels);
        int root_recv;
        for (i=0; i<nsupermodels; i++) {
            root_ids[i] = (int *)tl_alloc(sizeof(int), sm[i].nsubmodels);
            for (j=0; j< sm[i].nsubmodels;j++){
#ifdef _MESSG
                if(sm[i].submodel[j].grid->smpi->myid == 0){
                    root_recv = myid;
                }else{
                    root_recv = UNSET_INT;
                }
                ierr = MPI_Allreduce(&root_recv, &root_ids[i][j],1,MPI_INT,MPI_MAX,cstorm_comm);
#else
                root_ids[i][j]=0;
#endif
            }
        }
        
#ifdef _MESSG
        /* insure every procesor knows which submodels are flux coupled */
        for (i=0; i<nsupermodels; i++) {
            int flux_recv;
            for (i=0; i<nsupermodels; i++) {
                for(j=0;j<sm[i].nsubmodels;j++){
                    flux_recv=0;
                    ierr = MPI_Allreduce(&(flux_submodels[i][j]),&flux_recv,1,MPI_INT,MPI_MAX,cstorm_comm);
                    flux_submodels[i][j]=flux_recv;
                }
            }
        }
#endif
        
        /* Find flux coupled submodel elements */
        SMODEL *model = NULL;
        for (i=0; i<nsupermodels; i++) {
            coupled_normal_flux = NULL;
            int jj;
            if(flux_interface==1){
                coupled_normal_flux = (double ***)tl_alloc(sizeof(double **), nsupermodels);
                for (i=0; i<nsupermodels; i++) {
                    coupled_normal_flux[i] = (double **)tl_alloc(sizeof(double *), sm[i].nsubmodels);
                    ie = 0;
                    for(j=0;j<sm[i].nsubmodels;j++){
                        if(flux_submodels[i][j]==1){
                            model=&(sm[i].submodel[j]);
                            for(jj=0;jj<sm[i].nsubmodels;jj++){
                                if(sm[i].submodel[jj].proc_flag==1 &&jj!=j){
                                    read_flux_geo(i,j,model->io ,sm[i].submodel[jj].grid);
                                }
                            }
                            if(model->grid->smpi->myid==0){
                                if(model->grid->ndim == 2 ){
                                    ie=model->grid->macro_nelems1d;
                                }else{
                                    ie=model->grid->macro_nelems2d;
                                }
                            }else{
                                ie=0;
                            }
#ifdef _MESSG
                            MPI_Bcast(&ie, 1, MPI_INT,root_ids[i][j],cstorm_comm);
                            //printf("i %d j %d ie %d \n",i,j,ie);
#endif
                            coupled_normal_flux[i][j] = (double *)tl_alloc(sizeof(double ), ie);
                            sarray_init_dbl(coupled_normal_flux[i][j], ie);
                        }else{
                            coupled_normal_flux[i][j] = NULL;
                        }
                    }
                }
            }
        }
        
        for (i=0; i<nsupermodels; i++) {
            flux_submodels[i]  = (int *)tl_free(sizeof(int ), sm[i].nsubmodels, flux_submodels[i]);
            root_ids[i] = (int *)tl_free(sizeof(int), sm[i].nsubmodels,root_ids[i]);
        }
        flux_submodels  = (int **)tl_free(sizeof(int *), nsupermodels,flux_submodels);
        root_ids = (int **)tl_free(sizeof(int *), nsupermodels,root_ids);
    } else {
        for (i=0; i<nsupermodels; i++) {
            flux_submodels[i]  = (int *)tl_free(sizeof(int ), sm[i].nsubmodels, flux_submodels[i]);
        }
        flux_submodels  = (int **)tl_free(sizeof(int *), nsupermodels,flux_submodels);
    }
    
    //tag();
    
    /* set fmaps, interface data, restrictions for adaption and adaption partition comms */
    
    for (i=0; i<nsupermodels; i++) {
        for(j=0;j<sm[i].nsubmodels;j++){
            if(sm[i].submodel[j].proc_flag==1){
                
                for (ie=0; ie<sm[i].submodel[j].grid->nelems1d; ie++){
                    istr = sm[i].submodel[j].grid->elem1d[ie].string;
                    if (sm[i].submodel[j].str_values[istr].ol_flow.bc_flag == BCT_HYBRID_INTERNAL){
                        sm[i].submodel[j].grid->elem2d[sm[i].submodel[j].grid->elem1d[ie].elem2d].interface=1;
                        if((sm[i].submodel[j].grid->node[sm[i].submodel[j].grid->elem1d[ie].nodes[0]].resident_pe == sm[i].submodel[j].grid->smpi->myid) ||
                           (sm[i].submodel[j].grid->node[sm[i].submodel[j].grid->elem1d[ie].nodes[1]].resident_pe == sm[i].submodel[j].grid->smpi->myid)){
                            sm[i].submodel[j].grid->interface=1;
                        }
                    }
                }
                for (ie=0; ie<sm[i].submodel[j].grid->nelems2d; ie++){
                    istr = sm[i].submodel[j].grid->elem2d[ie].string;
                    if (sm[i].submodel[j].str_values[istr].ol_flow.bc_flag == BCT_HYBRID_INTERNAL && sm[i].submodel[j].grid->nelems3d>0){
                        sm[i].submodel[j].grid->elem3d[sm[i].submodel[j].grid->elem2d[ie].id_3d].interface=1;
                        if((sm[i].submodel[j].grid->node[sm[i].submodel[j].grid->elem2d[ie].nodes[0]].resident_pe == sm[i].submodel[j].grid->smpi->myid) ||
                           (sm[i].submodel[j].grid->node[sm[i].submodel[j].grid->elem2d[ie].nodes[1]].resident_pe == sm[i].submodel[j].grid->smpi->myid) ||
                           (sm[i].submodel[j].grid->node[sm[i].submodel[j].grid->elem2d[ie].nodes[2]].resident_pe == sm[i].submodel[j].grid->smpi->myid)){
                            sm[i].submodel[j].grid->interface=1;
                        }
                    }
                }
                
#ifdef _MESSG
                assert(sm[i].submodel[j].grid->interface>=0);
                assert(sm[i].submodel[j].grid->smpi->myid>=0);
                MPI_Comm_split(sm[i].submodel[j].grid->smpi->ADH_COMM,sm[i].submodel[j].grid->interface,sm[i].submodel[j].grid->smpi->myid,&(part_comm));
                sm[i].submodel[j].grid->part_smpi = (SMPI *) tl_alloc(sizeof(SMPI), 1);
                smpi_init(sm[i].submodel[j].grid->part_smpi, part_comm);
                sm[i].submodel[j].grid->part_map = (int *)tl_alloc(sizeof(int), sm[i].submodel[j].grid->part_smpi->npes);
                if(sm[i].submodel[j].grid->interface==0){
                    ierr = MPI_Allgather(&(sm[i].submodel[j].grid->smpi->myid),1, MPI_INT, sm[i].submodel[j].grid->part_map, 1, MPI_INT, sm[i].submodel[j].grid->part_smpi->ADH_COMM);
                }
#endif
            }
        }
        ssuperModel_create_interface_data(&(sm[i]));
        ssuperModel_forward_map(&(sm[i]));
        if(sm[i].supersmpi_wvel!=NULL)ssuperModel_forward_map_wvel(&(sm[i]));
    }
#ifdef _MESSG    
    messg_barrier(cstorm_comm);
    for(i=0;i<nsupermodelcount;i++) {
        for (j=0;j<nsubmodels[i];j++) {
            MPI_Comm_dup(sub_comm, &(sm[i].submodel[j].model_comm));
            
            int myidG,npesG,myidS,npesS,ierr_code;
            ierr_code = MPI_Comm_rank(cstorm_comm, &myidG);
            ierr_code = MPI_Comm_size(cstorm_comm, &npesG);
            ierr_code = MPI_Comm_rank(sm[i].submodel[j].model_comm,    &myidS);
            ierr_code = MPI_Comm_size(sm[i].submodel[j].model_comm,    &npesS);
            //printf("superModel: %d submodel: %d :: world pe %d of %d :: model pe %d of %d\n",i,j,myidG+1,npesG,myidS+1,npesS);
            
        }
    }
    messg_barrier(cstorm_comm);
#endif

    ssuperModel_printScreen(sm, si, nsupermodels, nsuperinterfaces);
    super_procs_wgt = (int *) tl_free(sizeof(int), nsupermodels, super_procs_wgt);
    sub_procs = (int *) tl_free(sizeof(int), nsupermodels, sub_procs);
    for(i=0;i<nsupermodels;i++)	sub_procs_wgt[i] = (int *) tl_free(sizeof(int ), nsubmodels[i], sub_procs_wgt[i]);
    sub_procs_wgt = (int **) tl_free(sizeof(int *), nsupermodels, sub_procs_wgt);
    fclose(fp);
    nsubmodels = (int *) tl_free(sizeof(int), nsupermodels, nsubmodels);
    nInterfaces = (int *) tl_free(sizeof(int), nsupermodels, nInterfaces);
    //messg_barrier(cstorm_comm);//exit(-1);
}


/***********************************************************/
/***********************************************************/

void generate_2d_2d_interface(SSUPER_MODEL *sm, SINTERFACE *ifce, int mod_id1, int mod_id2, int nsupermodels){
    int i, j, k, l;
    int istr1, istr2;
    int ie1, ie2;
    int check1,check2;
    
    SMODEL *mod1 = &(sm->submodel[mod_id1]);
    SMODEL *mod2 = &(sm->submodel[mod_id2]);
    
    SGRID *grid1 = mod1->grid;
    SGRID *grid2 = mod2->grid;
    
    check1=0;check2=0;
    if((mod1->proc_flag==1)||(mod2->proc_flag==1)){
        for (ie1=0; ie1<grid1->nelems1d; ie1++){
            istr1 = grid1->elem1d[ie1].string;
            int n01 = grid1->elem1d[ie1].nodes[0];
            int n02 = grid1->elem1d[ie1].nodes[1];
            //printf("\nsubmodel[%i].elem1d[%i].string[%i], value = %i",mod_id1,ie1, istr1, mod1->str_values[istr1].ol_flow.bc_flag);
            if (mod1->str_values[istr1].ol_flow.bc_flag == BCT_HYBRID_INTERNAL){
                check1=1;
                for (ie2=0; ie2<grid2->nelems1d; ie2++){
                    istr2 = grid2->elem1d[ie2].string;
                    int n11 = grid2->elem1d[ie2].nodes[0];
                    int n12 = grid2->elem1d[ie2].nodes[1];
                    if (mod2->str_values[istr2].ol_flow.bc_flag == BCT_HYBRID_INTERNAL){
                        check2=1;
                        if   ((projected_node_distance(grid1,n01,grid2,n11)<NOT_QUITE_SMALL)
                              &&  (projected_node_distance(grid1,n02,grid2,n12)<NOT_QUITE_SMALL)){
                            int already_recorded = 0;
                            for (i=0; i<ifce->NumNodeColumns; i++){
                                if (ifce->nodelist[i].surfnode[0] == UNSET_INT){
                                    /* If you are here, then it means n01 and n11 have not been added */
                                    break;
                                }
                                if (ifce->nodelist[i].surfnode[0] == grid1->node[n01].gid){
                                    already_recorded += 1;
                                }
                                if (ifce->nodelist[i].surfnode[0] == grid1->node[n02].gid){
                                    already_recorded += 2;
                                }
                            }
                            switch (already_recorded){
                                case(0):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n01].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n11].gid;
                                    ifce->nodelist[i+1].surfnode[0] = grid1->node[n02].gid;
                                    ifce->nodelist[i+1].surfnode[1] = grid2->node[n12].gid;
                                    break;
                                case(1):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n02].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n12].gid;
                                    break;
                                case(2):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n01].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n11].gid;
                                    break;
                                case (3):
                                    /* Means both nodes have already been recorded */
                                    break;
                                default:
                                    printf("Unexpected behavior in interface creation!");
                                    exit (-1);
                                    break;
                            }
                            break; /* breaks the ie2 for loop */
                        } else if   ((projected_node_distance(grid1,n01,grid2,n12)<NOT_QUITE_SMALL)
                                     &&  (projected_node_distance(grid1,n02,grid2,n11)<NOT_QUITE_SMALL)){
                            int already_recorded = 0;
                            for (i=0; i<ifce->NumNodeColumns; i++){
                                if (ifce->nodelist[i].surfnode[0] == UNSET_INT){
                                    /* If you are here, then it means n01 and n11 have not been added */
                                    break;
                                }
                                if (ifce->nodelist[i].surfnode[0] == grid1->node[n01].gid){
                                    already_recorded += 1;
                                }
                                if (ifce->nodelist[i].surfnode[0] == grid1->node[n02].gid){
                                    already_recorded += 2;
                                }
                            }
                            switch (already_recorded){
                                case(0):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n01].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n12].gid;
                                    ifce->nodelist[i+1].surfnode[0] = grid1->node[n02].gid;
                                    ifce->nodelist[i+1].surfnode[1] = grid2->node[n11].gid;
                                    break;
                                case(1):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n02].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n11].gid;
                                    break;
                                case(2):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n01].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n12].gid;
                                    break;
                                case (3):
                                    /* Means both nodes have already been recorded */
                                    break;
                                default:
                                    printf("Unexpected behavior in interface creation!");
                                    exit (-1);
                                    break;
                            }
                            break; /* breaks the ie2 for loop */
                        }
                    }
                } /* ie2 loop */
            }
        } /* ie1 loop */
    }
    
#ifdef _MESSG
    if((check1==0)||(check2==0)){ /* no interfaces, re-initialize submodel grids */
        for(i=0;i<nsupermodels;i++){
            for(j=0;j<sm[j].nsubmodels;j++){
                if(sm[i].submodel[j].proc_flag==0){
                    /*sgrid_free(sm[i].submodel[j].grid);
                     sgrid_init(sm[i].submodel[j].grid, 0, NULL);*/
                    //smodel_free(&sm[i].submodel[j]);
                }
            }
        }
    }
#endif
    
    
}



/***********************************************************/
/***********************************************************/

void generate_3d_3d_interface(SSUPER_MODEL *sm, SINTERFACE *ifce, int mod_id1, int mod_id2, int nsupermodels){
    int i, j, k, l;
    int istr1, istr2;
    int ie1, ie2;
    int check1, check2;
    
    SMODEL *mod1 = &(sm->submodel[mod_id1]);
    SMODEL *mod2 = &(sm->submodel[mod_id2]);
    
    SGRID *grid1 = mod1->grid;
    SGRID *grid2 = mod2->grid;
    
    check1=0;check2=0;
    
    double node_distance;
    if((mod1->proc_flag==1)||(mod2->proc_flag==1)){
        for (ie1=0; ie1<grid1->nelems2d; ie1++){
            istr1 = grid1->elem2d[ie1].string;
            //printf("\nsubmodel[%i].elem2d[%i].string[%i], value = %i",mod_id1,ie1, istr1, mod1->str_values[istr1].ol_flow.bc_flag);
            if (mod1->str_values[istr1].ol_flow.bc_flag == BCT_HYBRID_INTERNAL){
                check1=1;
                node_distance = projected_node_distance(grid1,grid1->elem2d[ie1].nodes[0], grid1, grid1->elem2d[ie1].nodes[1]);
                int n01, n02;
                if (node_distance < NOT_QUITE_SMALL){
                    n01 = find_vertical_segment(grid1, grid1->elem2d[ie1].nodes[0], grid1->vertical_hash); n01 = grid1->vertical_list[n01]->id;
                    n02 = find_vertical_segment(grid1, grid1->elem2d[ie1].nodes[2], grid1->vertical_hash); n02 = grid1->vertical_list[n02]->id;
                    //printf("\nn01 = %i->%i, n02 = %i->%i", grid1->elem2d[ie1].nodes[0], n01, grid1->elem2d[ie1].nodes[2], n02);
                }
                else{
                    n01 = find_vertical_segment(grid1, grid1->elem2d[ie1].nodes[0], grid1->vertical_hash); n01 = grid1->vertical_list[n01]->id;
                    n02 = find_vertical_segment(grid1, grid1->elem2d[ie1].nodes[1], grid1->vertical_hash); n02 = grid1->vertical_list[n02]->id;
                    //printf("\nn01 = %i->%i, n02 = %i->%i", grid1->elem2d[ie1].nodes[0], n01, grid1->elem2d[ie1].nodes[1], n02);
                }
                
                assert (projected_node_distance(grid1,n01,grid1,n02) > NOT_QUITE_SMALL);
                
                for (ie2=0; ie2<grid2->nelems2d; ie2++){
                    istr2 = grid2->elem2d[ie2].string;
                    if (mod2->str_values[istr2].ol_flow.bc_flag == BCT_HYBRID_INTERNAL){
                        check2=0;
                        node_distance = projected_node_distance(grid2,grid2->elem2d[ie2].nodes[0], grid2, grid2->elem2d[ie2].nodes[1]);
                        int n11, n12;
                        if (node_distance < NOT_QUITE_SMALL){
                            n11 = find_vertical_segment(grid2, grid2->elem2d[ie2].nodes[0], grid2->vertical_hash); n11 = grid2->vertical_list[n11]->id;
                            n12 = find_vertical_segment(grid2, grid2->elem2d[ie2].nodes[2], grid2->vertical_hash); n12 = grid2->vertical_list[n12]->id;
                            //printf("\nn11 = %i->%i, n12 = %i->%i", grid2->elem2d[ie2].nodes[0], n11, grid2->elem2d[ie2].nodes[2], n12);
                        }
                        else{
                            n11 = find_vertical_segment(grid2, grid2->elem2d[ie2].nodes[0], grid2->vertical_hash); n11 = grid2->vertical_list[n11]->id;
                            n12 = find_vertical_segment(grid2, grid2->elem2d[ie2].nodes[1], grid2->vertical_hash); n12 = grid2->vertical_list[n12]->id;
                            //printf("\nn11 = %i->%i, n12 = %i->%i", grid2->elem2d[ie2].nodes[0], n11, grid2->elem2d[ie2].nodes[1], n12);
                        }
                        assert (projected_node_distance(grid2,n11,grid2,n12) > NOT_QUITE_SMALL);
                        
                        if   ((projected_node_distance(grid1,n01,grid2,n11)<NOT_QUITE_SMALL)
                              &&  (projected_node_distance(grid1,n02,grid2,n12)<NOT_QUITE_SMALL)){
                            //printf("\nsubmodel[%i], node[%i] coupled to submodel[%i], node[%i]",
                            //    mod_id1, n01, mod_id2, n11);
                            //printf("\nsubmodel[%i], node[%i] coupled to submodel[%i], node[%i]\n",
                            //    mod_id1, n02, mod_id2, n12);
                            
                            int already_recorded = 0;
                            for (i=0; i<ifce->NumNodeColumns; i++){
                                if (ifce->nodelist[i].surfnode[0] == UNSET_INT){
                                    /* If you are here, then it means n01 and n11 have not been added */
                                    break;
                                }
                                if (ifce->nodelist[i].surfnode[0] == grid1->node[n01].gid){
                                    already_recorded += 1;
                                }
                                if (ifce->nodelist[i].surfnode[0] == grid1->node[n02].gid){
                                    already_recorded += 2;
                                }
                            }
                            switch (already_recorded){
                                case(0):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n01].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n11].gid;
                                    ifce->nodelist[i+1].surfnode[0] = grid1->node[n02].gid;
                                    ifce->nodelist[i+1].surfnode[1] = grid2->node[n12].gid;
                                    break;
                                case(1):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n02].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n12].gid;
                                    break;
                                case(2):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n01].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n11].gid;
                                    break;
                                case (3):
                                    /* Means both nodes have already been recorded */
                                    break;
                                default:
                                    printf("Unexpected behavior in interface creation!");
                                    exit (-1);
                                    break;
                            }
                            break; /* breaks the ie2 for loop */
                        }
                        else if   ((projected_node_distance(grid1,n01,grid2,n12)<NOT_QUITE_SMALL)
                                   &&  (projected_node_distance(grid1,n02,grid2,n11)<NOT_QUITE_SMALL)){
                            //printf("\nsubmodel[%i], node[%i] coupled to submodel[%i], node[%i]",
                            //    mod_id1, n01, mod_id2, n12);
                            //printf("\nsubmodel[%i], node[%i] coupled to submodel[%i], node[%i]\n",
                            //    mod_id1, n02, mod_id2, n11);
                            int already_recorded = 0;
                            for (i=0; i<ifce->NumNodeColumns; i++){
                                if (ifce->nodelist[i].surfnode[0] == UNSET_INT){
                                    /* If you are here, then it means n01 and n11 have not been added */
                                    break;
                                }
                                if (ifce->nodelist[i].surfnode[0] == grid1->node[n01].gid){
                                    already_recorded += 1;
                                }
                                if (ifce->nodelist[i].surfnode[0] == grid1->node[n02].gid){
                                    already_recorded += 2;
                                }
                            }
                            switch (already_recorded){
                                case(0):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n01].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n12].gid;
                                    ifce->nodelist[i+1].surfnode[0] = grid1->node[n02].gid;
                                    ifce->nodelist[i+1].surfnode[1] = grid2->node[n11].gid;
                                    break;
                                case(1):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n02].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n11].gid;
                                    break;
                                case(2):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n01].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n12].gid;
                                    break;
                                case (3):
                                    /* Means both nodes have already been recorded */
                                    break;
                                default:
                                    printf("Unexpected behavior in interface creation!");
                                    exit (-1);
                                    break;
                            }
                            break; /* breaks the ie2 for loop */
                        }
                    }
                } /* ie2 loop */
            }
        } /* ie1 loop */
    }
    
#ifdef _MESSG
    if((check1==0)||(check2==0)){ /* no interfaces, re-initialize submodel grids */
        for(i=0;i<nsupermodels;i++){
            for(j=0;j<sm[j].nsubmodels;j++){
                if(sm[i].submodel[j].proc_flag==0){
                    /*sgrid_free(sm[i].submodel[j].grid);
                     sgrid_init(sm[i].submodel[j].grid, 0, NULL);*/
                    //smodel_free(&sm[i].submodel[j]);
                }
            }
        }
    }
#endif
    
    //    exit(-1);
}


/***********************************************************/
/***********************************************************/

void generate_2d_3d_interface(SSUPER_MODEL *sm, SINTERFACE *ifce, int mod_id1, int mod_id2, int nsupermodels){
    int i, j, k, l;
    int istr1, istr2;
    int ie1, ie2;
    int check1,check2;
    
    SMODEL *mod1 = &(sm->submodel[mod_id1]); /* MUST be 2D */
    SMODEL *mod2 = &(sm->submodel[mod_id2]); /* MUST be 3D */
    assert (mod1->flag.SW2_FLOW && mod2->flag.SW3_FLOW);
    
    SGRID *grid1 = mod1->grid;
    SGRID *grid2 = mod2->grid;
    
    check1=0;check2=0;
    
    double node_distance;
    if((mod1->proc_flag==1)||(mod2->proc_flag==1)){
        for (ie1=0; ie1<grid1->nelems1d; ie1++){
            istr1 = grid1->elem1d[ie1].string;
            int n01 = grid1->elem1d[ie1].nodes[0];
            int n02 = grid1->elem1d[ie1].nodes[1];
            //printf("\nsubmodel[%i].elem1d[%i].string[%i], value = %i",mod_id1,ie1, istr1, mod1->str_values[istr1].ol_flow.bc_flag);
            if (mod1->str_values[istr1].ol_flow.bc_flag == BCT_HYBRID_INTERNAL){
                check1=1;
                for (ie2=0; ie2<grid2->nelems2d; ie2++){
                    istr2 = grid2->elem2d[ie2].string;
                    if (mod2->str_values[istr2].ol_flow.bc_flag == BCT_HYBRID_INTERNAL){
                        check2=1;
                        node_distance = projected_node_distance(grid2,grid2->elem2d[ie2].nodes[0], grid2, grid2->elem2d[ie2].nodes[1]);
                        int n11, n12;
                        if (node_distance < NOT_QUITE_SMALL){
                            n11 = find_vertical_segment(grid2, grid2->elem2d[ie2].nodes[0], grid2->vertical_hash); n11 = grid2->vertical_list[n11]->id;
                            n12 = find_vertical_segment(grid2, grid2->elem2d[ie2].nodes[2], grid2->vertical_hash); n12 = grid2->vertical_list[n12]->id;
                            //printf("\nn11 = %i->%i, n12 = %i->%i", grid2->elem2d[ie2].nodes[0], n11, grid2->elem2d[ie2].nodes[2], n12);
                        }
                        else{
                            n11 = find_vertical_segment(grid2, grid2->elem2d[ie2].nodes[0], grid2->vertical_hash); n11 = grid2->vertical_list[n11]->id;
                            n12 = find_vertical_segment(grid2, grid2->elem2d[ie2].nodes[1], grid2->vertical_hash); n12 = grid2->vertical_list[n12]->id;
                            //printf("\nn11 = %i->%i, n12 = %i->%i", grid2->elem2d[ie2].nodes[0], n11, grid2->elem2d[ie2].nodes[1], n12);
                        }
                        assert (projected_node_distance(grid2,n11,grid2,n12) > NOT_QUITE_SMALL);
                        
                        if   ((projected_node_distance(grid1,n01,grid2,n11)<NOT_QUITE_SMALL)
                              &&  (projected_node_distance(grid1,n02,grid2,n12)<NOT_QUITE_SMALL)){
                            //printf("\nsubmodel[%i], node[%i] coupled to submodel[%i], node[%i]",
                            //    mod_id1, n01, mod_id2, n11);
                            //printf("\nsubmodel[%i], node[%i] coupled to submodel[%i], node[%i]\n",
                            //    mod_id1, n02, mod_id2, n12);
                            
                            int already_recorded = 0;
                            for (i=0; i<ifce->NumNodeColumns; i++){
                                if (ifce->nodelist[i].surfnode[0] == UNSET_INT){
                                    /* If you are here, then it means n01 and n11 have not been added */
                                    break;
                                }
                                if (ifce->nodelist[i].surfnode[0] == grid1->node[n01].gid){
                                    already_recorded += 1;
                                }
                                if (ifce->nodelist[i].surfnode[0] == grid1->node[n02].gid){
                                    already_recorded += 2;
                                }
                            }
                            switch (already_recorded){
                                case(0):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n01].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n11].gid;
                                    ifce->nodelist[i+1].surfnode[0] = grid1->node[n02].gid;
                                    ifce->nodelist[i+1].surfnode[1] = grid2->node[n12].gid;
                                    break;
                                case(1):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n02].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n12].gid;
                                    break;
                                case(2):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n01].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n11].gid;
                                    break;
                                case (3):
                                    /* Means both nodes have already been recorded */
                                    break;
                                default:
                                    printf("Unexpected behavior in interface creation!");
                                    exit (-1);
                                    break;
                            }
                            break; /* breaks the ie2 for loop */
                        }
                        else if   ((projected_node_distance(grid1,n01,grid2,n12)<NOT_QUITE_SMALL)
                                   &&  (projected_node_distance(grid1,n02,grid2,n11)<NOT_QUITE_SMALL)){
                            //printf("\nsubmodel[%i], node[%i] coupled to submodel[%i], node[%i]",
                            //    mod_id1, n01, mod_id2, n12);
                            //printf("\nsubmodel[%i], node[%i] coupled to submodel[%i], node[%i]\n",
                            //    mod_id1, n02, mod_id2, n11);
                            int already_recorded = 0;
                            for (i=0; i<ifce->NumNodeColumns; i++){
                                if (ifce->nodelist[i].surfnode[0] == UNSET_INT){
                                    /* If you are here, then it means n01 and n11 have not been added */
                                    break;
                                }
                                if (ifce->nodelist[i].surfnode[0] == grid1->node[n01].gid){
                                    already_recorded += 1;
                                }
                                if (ifce->nodelist[i].surfnode[0] == grid2->node[n02].gid){
                                    already_recorded += 2;
                                }
                            }
                            switch (already_recorded){
                                case(0):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n01].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n12].gid;
                                    ifce->nodelist[i+1].surfnode[0] = grid1->node[n02].gid;
                                    ifce->nodelist[i+1].surfnode[1] = grid2->node[n11].gid;
                                    break;
                                case(1):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n02].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n11].gid;
                                    break;
                                case(2):
                                    ifce->nodelist[i].surfnode[0] = grid1->node[n01].gid;
                                    ifce->nodelist[i].surfnode[1] = grid2->node[n12].gid;
                                    break;
                                case (3):
                                    /* Means both nodes have already been recorded */
                                    break;
                                default:
                                    printf("Unexpected behavior in interface creation!");
                                    exit (-1);
                                    break;
                            }
                            break; /* breaks the ie2 for loop */
                        }
                    }
                } /* ie2 loop */
            }
        } /* ie1 loop */
    }
#ifdef _MESSG
    if((check1==0)||(check2==0)){ /* no interfaces, re-initialize submodel grids */
        for(i=0;i<nsupermodels;i++){
            for(j=0;j<sm[j].nsubmodels;j++){
                if(sm[i].submodel[j].proc_flag==0){
                    /*sgrid_free(sm[i].submodel[j].grid);
                     sgrid_init(sm[i].submodel[j].grid, 0, NULL);*/
                    
                }
            }
        }
    }
#endif
    
    //    exit(-1);
}

/*************************************************************************************/
/* This function sets the owner of interface nodes to model_id[0] in the interface struct */
/* This arbitrary, but only one processor can own the interface nodes so this is easiest way */
/* After assigning ownership, the processor need to be included in the submodel communicator */
/* Then updates will work for interface nodes */

void fix_interface_ownership(SSUPER_MODEL *sm, int nsupermodels) {
    
#ifdef _MESSG
    int i,j,k,nd,inode,inode2,ierr,my_nnode_max,found,check,kk,position;
    int node1, node0, node;
    int *add_me, *gid_data, *recv_gid_data, color, prime_mod;
    int *rid_data, *recv_rid_data;
    int **old_add_me, **new_add_me, reset_smpi=0;
    int master_mod, slave_mod, old_nnodes;
    SINTERFACE *ifce;
    MPI_Comm new_sub_comm;
    
    ID_LIST_ITEM *ptr;
//    if(DEBUG){
//        for(k=0;k<sm->supersmpi->npes;k++){
//            if(k==sm->supersmpi->myid){
//                for(i=0;i<sm->nsubmodels;i++){
//                    if(sm->submodel[i].proc_flag==1){
//                        printf("MYID %d checking model %d of %d \n", sm->supersmpi->myid,i,sm->nsubmodels);
//                        sgrid_check(sm->submodel[i].grid, __FILE__,__LINE__);
//                    }
//                    if(sm->submodel[i].grid->ndim==2){
//                        if(i==1){
//                            print_grid_to_file(sm->submodel[i].grid, "PRE_FIX_2D_1");
//                        }else{
//                            print_grid_to_file(sm->submodel[i].grid, "PRE_FIX_2D_2");
//                        }
//                    }else{
//                        print_grid_to_file(sm->submodel[i].grid, "PRE_FIX_3D");
//                    }
//                }
//            }
//        }
//    }
    
    /* reset resident id's for interface nodes on higher submodel */
    add_me = (int *) tl_alloc(sizeof(int), sm->nsubmodels);
    old_add_me = (int **) tl_alloc(sizeof(int *), sm->nsubmodels);
    new_add_me = (int **) tl_alloc(sizeof(int *), sm->nsubmodels);
    color = 0;
    for(k=0;k<sm->nsubmodels;k++){
        old_add_me[k] = (int *) tl_alloc(sizeof(int), sm->supersmpi->npes);
        new_add_me[k] = (int *) tl_alloc(sizeof(int), sm->supersmpi->npes);
        add_me[k]=MPI_UNDEFINED;
        if(sm->submodel[k].proc_flag==1) {
            /* make sure all the originally assigned processors stay */
            add_me[k]=1;
            prime_mod=k;
        }else{
            for (inode=0;inode<sm->submodel[k].grid->my_nnodes;inode++){
                sm->submodel[k].grid->node[inode].resident_pe = UNSET_INT;
            }
        }
    }
    for (j=0;j<sm->nsubmodels;j++){
        ierr = MPI_Allgather(&(add_me[j]), 1, MPI_INT, old_add_me[j], 1, MPI_INT, sm->supersmpi->ADH_COMM);
    }
    
    for(j=0;j<sm->NumInterfaces;j++) {
        ifce = &(sm->interface[j]);
        if (sm->submodel[ifce->model_id[0]].flag.SW3_FLOW && sm->submodel[ifce->model_id[1]].flag.SW2_FLOW){
            /* Only if user input was model 0 = 3D, and model 1 = 2D. Then the order is flipped in memory. */
            master_mod = 1;
            slave_mod  = 0;
        }
        else{
            master_mod = 0;
            slave_mod = 1;
        }
        for(k=0;k<sm->interface[j].NumNodeColumns;k++) {
            if((sm->interface[j].nodelist[k].surfnode[0] != UNSET_INT) &&
               (sm->interface[j].nodelist[k].surfnode[1] != UNSET_INT) &&
               (sm->submodel[sm->interface[j].model_id[slave_mod]].proc_flag==1)){
                add_me[sm->interface[j].model_id[slave_mod]]=1;
                /* Interface nodes exist and model_id is "higher" */
                for(inode=0;inode<sm->submodel[sm->interface[j].model_id[slave_mod]].grid->nnodes;inode++){
                    check=0;
                    if(sm->submodel[sm->interface[j].model_id[slave_mod]].grid->node[inode].gid==
                       sm->interface[j].nodelist[k].surfnode[1]){
                        check++;
                        sm->submodel[sm->interface[j].model_id[slave_mod]].grid->node[inode].resident_pe=UNSET_INT;
                        /* if model is 3D, reset all nodes below surface node */
                        if(sm->submodel[sm->interface[j].model_id[slave_mod]].grid->type == COLUMNAR){
                            nd = find_vertical_segment(sm->submodel[sm->interface[j].model_id[slave_mod]].grid, inode, sm->submodel[sm->interface[j].model_id[slave_mod]].grid->vertical_hash);
                            ptr = sm->submodel[sm->interface[j].model_id[slave_mod]].grid->vertical_list[nd];
                            nd = ptr->id;
                            ptr = ptr->next;
                            sm->submodel[sm->interface[j].model_id[slave_mod]].grid->node[nd].resident_pe = UNSET_INT;
                            while (ptr->next != NULL) {
                                nd = ptr->id;
                                sm->submodel[sm->interface[j].model_id[slave_mod]].grid->node[nd].resident_pe = UNSET_INT;
                                ptr = ptr->next;
                                
                            }
                        }
                    }
                    if(check>=2)
                        break; /*if we found both nodes, break from loop */
                }
                /* reset interface surface nodes since they will no longer be owned by this processor */
                //sm->interface[j].nodelist[k].surfnode[0] = UNSET_INT;
                //sm->interface[j].nodelist[k].surfnode[1] = UNSET_INT;
            } else if((sm->interface[j].nodelist[k].surfnode[0] != UNSET_INT) &&
                      (sm->interface[j].nodelist[k].surfnode[1] != UNSET_INT) &&
                      (sm->submodel[sm->interface[j].model_id[slave_mod]].proc_flag==0)){
                /* We have interface nodes and primary model is lower. Need to be added to higher submodel */
                
                add_me[sm->interface[j].model_id[slave_mod]]=1;
                node0=-1;
                node1=-1;
                node=-1;
                for(inode=0;inode<sm->submodel[sm->interface[j].model_id[slave_mod]].grid->my_nnodes;inode++){
                    check=0;
                    if(sm->submodel[sm->interface[j].model_id[slave_mod]].grid->node[inode].gid==
                       sm->interface[j].nodelist[k].surfnode[1]){
                        check++;
                        if(prime_mod==sm->interface[j].model_id[master_mod]){
                            for(inode2=0;inode2<sm->submodel[prime_mod].grid->my_nnodes;inode2++){
                                if(sm->interface[j].nodelist[k].surfnode[0]==sm->submodel[prime_mod].grid->node[inode2].gid){
                                    /* Matched a interface node with one of my nodes in primary mod, set secondary mod nodes to myid */
                                    sm->submodel[sm->interface[j].model_id[slave_mod]].grid->node[inode].resident_pe = sm->submodel[sm->interface[j].model_id[slave_mod]].grid->smpi->myid;
                                    if(sm->submodel[sm->interface[j].model_id[slave_mod]].grid->type == COLUMNAR){
                                        nd = find_vertical_segment(sm->submodel[sm->interface[j].model_id[slave_mod]].grid, inode, sm->submodel[sm->interface[j].model_id[slave_mod]].grid->vertical_hash);
                                        ptr = sm->submodel[sm->interface[j].model_id[slave_mod]].grid->vertical_list[nd];
                                        nd = ptr->id;
                                        
                                        sm->submodel[sm->interface[j].model_id[slave_mod]].grid->node[nd].resident_pe = sm->submodel[sm->interface[j].model_id[slave_mod]].grid->smpi->myid;
                                        while (ptr->next != NULL) {
                                            nd = ptr->id;
                                            sm->submodel[sm->interface[j].model_id[slave_mod]].grid->node[nd].resident_pe = sm->submodel[sm->interface[j].model_id[slave_mod]].grid->smpi->myid;
                                            ptr = ptr->next;
                                            
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                    if(check>=1)
                        break;
                }
            }
        }
    }
    
    for (j=0;j<sm->nsubmodels;j++){
        ierr = MPI_Allgather(&(add_me[j]), 1, MPI_INT, new_add_me[j], 1, MPI_INT, sm->supersmpi->ADH_COMM);
    }
    
    for (j=0;j<sm->nsubmodels;j++) {
        reset_smpi=0;
        for (k=0;k<sm->supersmpi->npes;k++){
            if(old_add_me[j][k]!=new_add_me[j][k]){
                reset_smpi=1;
                break;
            }
        }
        if(reset_smpi==1){
            MPI_Comm_split(sm->supersmpi->ADH_COMM, add_me[j], sm->supersmpi->myid, &new_sub_comm);
            if(add_me[j]==1){
                //MPI_Comm_dup(new_sub_comm, &(sm->submodel[j].model_comm));
                if (sm->submodel[j].grid->start_node_ID !=NULL)
                    sm->submodel[j].grid->start_node_ID = (int *) tl_free(sizeof(int), sm->submodel[j].grid->smpi->npes,sm->submodel[j].grid->start_node_ID);
                if (sm->submodel[j].grid->end_node_ID !=NULL)
                    sm->submodel[j].grid->end_node_ID = (int *) tl_free(sizeof(int), sm->submodel[j].grid->smpi->npes, sm->submodel[j].grid->end_node_ID);
                sm->submodel[j].grid->start_node_ID=NULL;
                sm->submodel[j].grid->end_node_ID=NULL;
                smpi_free(sm->submodel[j].grid->smpi);
                MPI_Comm_dup(new_sub_comm, &(sm->submodel[j].model_comm));
                smpi_init(sm->submodel[j].grid->smpi,sm->submodel[j].model_comm);
                for (i=0; i <  sm->submodel[j].grid->nnodes; i++){
                    sm->submodel[j].grid->smpi->partition_info[i] = sm->submodel[j].grid->node[i].resident_pe;
                }
                if (sm->submodel[j].grid->ndim == 3) {
                    for (i=0; i <  sm->submodel[j].grid->nnodes_sur; i++){
                        sm->submodel[j].grid->smpi->surface_partition_info[i] = sm->submodel[j].grid->node[sm->submodel[j].grid->nodeID_2d_to_3d_sur[i]].resident_pe;
                    }
                }
                sm->submodel[j].proc_flag=1;
                /* reset my nodes to my new id */
                for(inode=0;inode<sm->submodel[j].grid->my_nnodes;inode++){
                    if(sm->submodel[j].grid->node[inode].resident_pe!=UNSET_INT){
                        sm->submodel[j].grid->node[inode].resident_pe=sm->submodel[j].grid->smpi->myid;
                    }
                }
                ierr = MPI_Allreduce(&(sm->submodel[j].grid->my_nnodes), &my_nnode_max, 1, MPI_INT, MPI_MAX, sm->submodel[j].grid->smpi->ADH_COMM);
                /* put all my_nnode gids into arrays in myid */
                recv_gid_data = (int *) tl_alloc(sizeof(int *), sm->submodel[j].grid->smpi->npes*my_nnode_max);
                recv_rid_data = (int *) tl_alloc(sizeof(int *), sm->submodel[j].grid->smpi->npes*my_nnode_max);
                gid_data = (int *) tl_alloc(sizeof(int), my_nnode_max);
                rid_data = (int *) tl_alloc(sizeof(int), my_nnode_max);
                /* load the gid array */
                for (inode=0;inode<my_nnode_max;inode++){
                    if((inode<sm->submodel[j].grid->my_nnodes)
                       &&(sm->submodel[j].grid->node[inode].resident_pe != UNSET_INT)){
                        /* only send gids for nodes I own, interface nodes I gave up
                         * should be set to UNSET_INT above */
                        gid_data[inode]=sm->submodel[j].grid->node[inode].gid;
                        rid_data[inode]=sm->submodel[j].grid->node[inode].resident_id;
                    }else{
                        gid_data[inode]=UNSET_INT;
                        rid_data[inode]=UNSET_INT;
                    }
                    
                }
                
                ierr = MPI_Allgather(gid_data, my_nnode_max, MPI_INT, recv_gid_data, my_nnode_max, MPI_INT, sm->submodel[j].model_comm);
                ierr = MPI_Allgather(rid_data, my_nnode_max, MPI_INT, recv_rid_data, my_nnode_max, MPI_INT, sm->submodel[j].model_comm);
                
                /* recv_gid_data has node gids in blocks of new ID*my_nnode_max Reset resident IDs for cleanup*/
                for(inode=0;inode<sm->submodel[j].grid->nnodes;inode++){
                    if((sm->submodel[j].grid->node[inode].resident_pe == UNSET_INT)||(inode>=sm->submodel[j].grid->my_nnodes)){
                        /* set node resident id to new resident id
                         * my_nnodes should have been set in loop above
                         * Interface nodes will be UNSET_INT and may be inside my_nnodes */
                        found=0;
                        for(k=0;k<sm->submodel[j].grid->smpi->npes;k++){
                            for(kk=0;kk<my_nnode_max;kk++){
                                position = kk+(k*my_nnode_max);
                                if ((recv_gid_data[position]==sm->submodel[j].grid->node[inode].gid)&&(recv_gid_data[position]!=UNSET_INT)) {
                                    /* we have found the node and can reset the resident ID */
                                    sm->submodel[j].grid->node[inode].resident_pe = k;
                                    sm->submodel[j].grid->node[inode].resident_id = recv_rid_data[position];
                                    found++;
                                    break;
                                }
                            }
                            if(found>0) break;
                        }
                    }
                    /* needed for node_cmp */
                    
                    sm->submodel[j].grid->node[inode].myid = sm->submodel[j].grid->smpi->myid;
                }
                
                /* All the resident IDs should be updated and ready to clean up and renumber the mesh
                 * decided to use partition_main with a flag rather than running all the subroutines here */
                
                partition_main(&(sm->submodel[j]), 1);
                recv_gid_data = (int *) tl_free(sizeof(int *), sm->submodel[j].grid->smpi->npes*my_nnode_max, recv_gid_data);
                gid_data = (int *) tl_free(sizeof(int), my_nnode_max, gid_data);
                recv_rid_data = (int *) tl_free(sizeof(int *), sm->submodel[j].grid->smpi->npes*my_nnode_max, recv_rid_data);
                rid_data = (int *) tl_free(sizeof(int), my_nnode_max, rid_data);
            }else{
                sm->submodel[j].proc_flag=0;
                if(sm->submodel[j].grid->ndim >0){
                    smodel_free(&(sm->submodel[j]));
                    smodel_defaults(&(sm->submodel[j]));
                    sm->submodel[j].grid = (SGRID *) tl_alloc(sizeof(SGRID), 1);
                    sgrid_init(sm->submodel[j].grid, 0, NULL);
                }
            }
        }
    }
    
    
    color=0;
    for(i=0;i<sm->nsubmodels;i++){
        if((sm->submodel[i].proc_flag==1) && sm->submodel[i].flag.SW3_FLOW){
            color=1;
        }
    }
    int ierr_code;
    ierr_code = MPI_Comm_split(sm->supersmpi->ADH_COMM, color, sm->supersmpi->myid, &new_sub_comm);
    if (ierr_code != MPI_SUCCESS){
        messg_err(ierr_code);
    }
    if(color>0){
        smpi_init(sm->supersmpi_wvel, new_sub_comm);
    }else{
        sm->supersmpi_wvel=(SMPI *) tl_free(sizeof(SMPI), 1, sm->supersmpi_wvel);
        sm->supersmpi_wvel=NULL;
    }
    /* update proc map */
    
    add_me = (int *) tl_free(sizeof(int), sm->nsubmodels, add_me);
    for(k=0;k<sm->nsubmodels;k++){
        old_add_me[k] = (int *) tl_free(sizeof(int), sm->supersmpi->npes, old_add_me[k]);
        new_add_me[k] = (int *) tl_free(sizeof(int), sm->supersmpi->npes, new_add_me[k]);
    }
    old_add_me = (int **) tl_free(sizeof(int *), sm->nsubmodels, old_add_me);
    new_add_me = (int **) tl_free(sizeof(int *), sm->nsubmodels, new_add_me);
    for(i=0;i<sm->nsubmodels;i++){
        if(sm->submodel[i].proc_flag==1){
            for(j=0;j<sm->submodel[i].grid->my_nnodes;j++){
                sm->submodel[i].grid->node[j].block=sm->supersmpi->myid;
                
            }
            messg_barrier(sm->submodel[i].grid->smpi->ADH_COMM);
            comm_update_snode(sm->submodel[i].grid);
            if(sm->submodel[i].flag.SW2_FLOW) comm_update_sw2 (sm->submodel[i].sw->d2, sm->submodel[i].grid);
            if(sm->submodel[i].flag.SW3_FLOW) comm_update_sw3 (sm->submodel[i].sw->d3, sm->submodel[i].grid);
#ifdef _ADH_GROUNDWATER
            if(sm->submodel[i].flag.GW_FLOW) comm_update_gw (sm->submodel[i].sgw, sm->submodel[i].grid);
#endif
            
//            if(DEBUG){
//#ifdef _DEBUG
//                tl_check_all_pickets(__FILE__,__LINE__);
//#endif
//                sgrid_check(sm->submodel[i].grid, __FILE__,__LINE__);
//                comm_check(sm->submodel[i].grid);
//                if(sm->submodel[i].grid->ndim==2){
//                    if(i==1){
//                        print_grid_to_file(sm->submodel[i].grid, "POST_FIX_2D_1");
//                    }else{
//                        print_grid_to_file(sm->submodel[i].grid, "POST_FIX_2D_2");
//                    }
//                }else{
//                    print_grid_to_file(sm->submodel[i].grid, "POST_FIX_3D");
//                    print_sw3_to_file("FIX_IDS", sm->submodel[i].sw,sm->submodel[i].grid,sm->submodel[i].ntransport,sm->submodel[i].con);
//                }
//            }
        }
    }
#endif
    return;
}

