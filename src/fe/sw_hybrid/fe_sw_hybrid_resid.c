/* This routine loops over submodels to calculate the shallow water residuals */
#include "global_header.h"
/*****************************************************************************************/
/*****************************************************************************************/
static int DEBUG           = OFF;
static int DEBUG_INTERFACE = OFF; /* Only for coupled models */
/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/

void fe_sw_hybrid_resid(SSUPER_MODEL *sm, int idum) {
#ifndef _PETSC
    int i,j,idof,ndof;
    SMODEL *mod;
    
    for (idof=0, ndof = sm->nnodes * sm->max_nsys; idof<ndof; idof++){
        sm->residual[idof] = 0.0;
    }
    
    // CJT :: DO ALL 2D DOMAINS FIRST!
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            mod = &(sm->submodel[i]);
            if(mod->flag.SW2_FLOW){
                fe_sw2_resid(sm,i);
                
                // CJT :: now save for 3D distribution
                for (j=0; j<mod->grid->nnodes; j++) {
                    mod->sw->d2->dacontResid[j] = sm->residual[sm->nsys*mod->fmap[j]+2];
                }
#ifdef _DEBUG
                if (DEBUG){
                    printf("\n2D Model %i partial HVEL residuals (Caution: may contain contributions from other submodels!):", i+1);
                    int j;
                    for (j=0; j<mod->grid->nnodes; j++) {
                        printf("\n    Node % 5i: R_mx = % 9.6e\t R_my = % 9.6e\t R_c = % 9.6e", j+1
                               , sm->residual[sm->nsys*mod->fmap[j]]
                               , sm->residual[sm->nsys*mod->fmap[j]+1]
                               , sm->residual[sm->nsys*mod->fmap[j]+2]);
                    }
                    // printScreen_resid("HVEL residual", sm->residual, sm->nnodes, sm->nsys, __LINE__, __FILE__);
                }
#endif
                
            }
        }
    }
    
    // NOW 3D DOMAINS
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            mod = &(sm->submodel[i]);
            if(mod->flag.SW3_FLOW){
                fe_hvel_resid(sm,i);
#ifdef _DEBUG
                if (DEBUG){
                    printf("\n3D Model %i partial HVEL residuals (Caution: may contain contributions from other submodels!):", i+1);
                    int j;
                    for (j=0; j<mod->grid->nnodes; j++) {
                        printf("\n    Node % 5i: R_mx = % 9.6e\t R_my = % 9.6e\t R_c = % 9.6e", j+1
                               , sm->residual[sm->nsys*mod->fmap[j]]
                               , sm->residual[sm->nsys*mod->fmap[j]+1]
                               , sm->residual[sm->nsys*mod->fmap[j]+2]);
                    }
                    // printScreen_resid("HVEL residual", sm->residual, sm->nnodes, sm->nsys, __LINE__, __FILE__);
                }
#endif
            }
        }
    }
    
    /*****************************************************************************************/
    /*****************************************************************************************/
    
#ifdef _DEBUG
    if (DEBUG_INTERFACE){
        int j, k, l;
        SMODEL *mod1;
        SMODEL *mod2;
        SINTERFACE *ifce;
        SINTERFACE_NODELIST *ndlist;
        
        printf("\n************************************************************\nPrinting Interface Residuals\n");
        /* Print to check equality of interface variables */
        for (j=0; j<sm->NumInterfaces; j++){
            ifce = &(sm->interface[j]);
            
            if (DEBUG_INTERFACE){
                printf("\nPrinting HVEL residuals for interface %i\n",j+1);
            }
            
            mod1 = &(sm->submodel[ifce->model_id[0]]);
            mod2 = &(sm->submodel[ifce->model_id[1]]);
            
            if (sm->submodel[ifce->model_id[0]].flag.SW2_FLOW){
                if(sm->submodel[ifce->model_id[1]].flag.SW2_FLOW){   /* 2D 2D coupling */
                    for (k=0; k<ifce->NumNodeColumns; k++){
                        ndlist = &(ifce->nodelist[k]);
                        printf("\nColumn %i :",k+1);
                        printf("\n  Model 1 (2D) :");
                        printf("\n    Node % 5i : R_mx = % 9.6e\t R_my = % 9.6e\t R_c = % 9.6e",
                               ndlist->couplednodes[0][0]+1,
                               sm->residual[sm->nsys*mod1->fmap[ndlist->couplednodes[0][0]]  ],
                               sm->residual[sm->nsys*mod1->fmap[ndlist->couplednodes[0][0]]+1],
                               sm->residual[sm->nsys*mod1->fmap[ndlist->couplednodes[0][0]]+2]
                               );
                        printf("\n  Model 2 (2D) :");
                        printf("\n    Node % 5i : R_mx = % 9.6e\t R_my = % 9.6e\t R_c = % 9.6e",
                               ndlist->couplednodes[1][0]+1,
                               sm->residual[sm->nsys*mod2->fmap[ndlist->couplednodes[1][0]]  ],
                               sm->residual[sm->nsys*mod2->fmap[ndlist->couplednodes[1][0]]+1],
                               sm->residual[sm->nsys*mod2->fmap[ndlist->couplednodes[1][0]]+2]
                               );
                        printf("\n");
                    }
                    printf("\n");
                }
                else{                                                /* 2D 3D coupling */
                    for (k=0; k<ifce->NumNodeColumns; k++){
                        ndlist = &(ifce->nodelist[k]);
                        printf("\nColumn %i :",k+1);
                        printf("\n  Model 1 (2D) :");
                        printf("\n    Node % 5i : R_mx = % 9.6e\t R_my = % 9.6e\t R_c = % 9.6e",
                               ndlist->couplednodes[0][0]+1,
                               sm->residual[sm->nsys*mod1->fmap[ndlist->couplednodes[0][0]]  ],
                               sm->residual[sm->nsys*mod1->fmap[ndlist->couplednodes[0][0]]+1],
                               sm->residual[sm->nsys*mod1->fmap[ndlist->couplednodes[0][0]]+2]
                               );
                        printf("\n  Model 2 (3D) :");
                        for (l=0; l<ndlist->size[1]; l++){
                            printf("\n    Node % 5i : R_mx = % 9.6e\t R_my = % 9.6e\t R_c = % 9.6e",
                                   ndlist->couplednodes[1][l]+1,
                                   sm->residual[sm->nsys*mod2->fmap[ndlist->couplednodes[1][l]]  ],
                                   sm->residual[sm->nsys*mod2->fmap[ndlist->couplednodes[1][l]]+1],
                                   sm->residual[sm->nsys*mod2->fmap[ndlist->couplednodes[1][l]]+2]
                                   );
                        }
                        printf("\n");
                    }
                    printf("\n");
                }
            }
            
            else if (sm->submodel[ifce->model_id[0]].flag.SW3_FLOW){
                if(sm->submodel[ifce->model_id[1]].flag.SW3_FLOW){   /* 3D 3D coupling */
                    for (k=0; k<ifce->NumNodeColumns; k++){
                        ndlist = &(ifce->nodelist[k]);
                        printf("\nColumn %i :",k+1);
                        printf("\n  Model 1 (3D) :");
                        for (l=0; l<ndlist->size[0]; l++){
                            printf("\n    Node % 5i : R_mx = % 9.6e\t R_my = % 9.6e\t R_c = % 9.6e",
                                   ndlist->couplednodes[0][l]+1,
                                   sm->residual[sm->nsys*mod1->fmap[ndlist->couplednodes[0][l]]  ],
                                   sm->residual[sm->nsys*mod1->fmap[ndlist->couplednodes[0][l]]+1],
                                   sm->residual[sm->nsys*mod1->fmap[ndlist->couplednodes[0][l]]+2]
                                   );
                        }
                        printf("\n  Model 2 (3D) :");
                        for (l=0; l<ndlist->size[1]; l++){
                            printf("\n    Node % 5i : R_mx = % 9.6e\t R_my = % 9.6e\t R_c = % 9.6e",
                                   ndlist->couplednodes[1][l]+1,
                                   sm->residual[sm->nsys*mod2->fmap[ndlist->couplednodes[1][l]]  ],
                                   sm->residual[sm->nsys*mod2->fmap[ndlist->couplednodes[1][l]]+1],
                                   sm->residual[sm->nsys*mod2->fmap[ndlist->couplednodes[1][l]]+2]
                                   );
                        }
                        printf("\n");
                    }
                    printf("\n");
                }
                else{                                                /* 3D 2D coupling */
                    for (k=0; k<ifce->NumNodeColumns; k++){
                        ndlist = &(ifce->nodelist[k]);
                        printf("\nColumn %i :",k+1);
                        printf("\n  Model 2 (2D) :");
                        printf("\n    Node % 5i : R_mx = % 9.6e\t R_my = % 9.6e\t R_c = % 9.6e",
                               ndlist->couplednodes[1][0]+1,
                               sm->residual[sm->nsys*mod2->fmap[ndlist->couplednodes[1][0]]  ],
                               sm->residual[sm->nsys*mod2->fmap[ndlist->couplednodes[1][0]]+1],
                               sm->residual[sm->nsys*mod2->fmap[ndlist->couplednodes[1][0]]+2]
                               );
                        printf("\n  Model 1 (3D) :");
                        for (l=0; l<ndlist->size[0]; l++){
                            printf("\n    Node % 5i : R_mx = % 9.6e\t R_my = % 9.6e\t R_c = % 9.6e",
                                   ndlist->couplednodes[0][l]+1,
                                   sm->residual[sm->nsys*mod1->fmap[ndlist->couplednodes[0][l]]  ],
                                   sm->residual[sm->nsys*mod1->fmap[ndlist->couplednodes[0][l]]+1],
                                   sm->residual[sm->nsys*mod1->fmap[ndlist->couplednodes[0][l]]+2]
                                   );
                        }
                        printf("\n");
                    }
                    printf("\n");
                }
            }
            printf("\n");
        } /* end for (j=0; j<sm->NumInterfaces; j++){} */
    } /* end if (DEBUG) */
    
    if (DEBUG) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<sm->supersmpi->npes;proc_count++){
            if(proc_count==sm->supersmpi->myid){
                printf("***********myid %d",sm->supersmpi->myid);
#endif
                printScreen_resid("Final HVEL residual", sm->residual, sm->nnodes, sm->nsys, __LINE__, __FILE__);
#ifdef _MESSG
            }
            messg_barrier(MPI_COMM_WORLD);
        }
        //exit(-1);
#endif
    }
#endif
    
#endif
}
