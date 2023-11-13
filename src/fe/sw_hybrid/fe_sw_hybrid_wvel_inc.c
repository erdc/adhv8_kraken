/* This routine loops over submodels to increment the shallow water solution */
#include "global_header.h"

/***********************************************************/
/***********************************************************/
static int DEBUG           = OFF;
static int DEBUG_INTERFACE = OFF; /* Only for coupled models */
/***********************************************************/
/***********************************************************/

void fe_sw_hybrid_wvel_inc(SSUPER_MODEL *sm, int idum) {
#ifndef _PETSC
    
    int i, node_count = 0;
    SMODEL *mod;
    
    for (i=0; i<sm->nsubmodels; i++){
        mod = &(sm->submodel[i]);
        if (mod->flag.SW3_FLOW) {
            fe_wvel_inc(sm,i);
#ifdef _DEBUG
            if (DEBUG){
                printf("\n3D Model %i solution (WVEL stage):", i+1);
                int j;
                for (j = 0; j < mod->grid->nnodes; j++) {
                    printf("\n    Node % 5i: vel.x = % 9.6e\t vel.y = % 9.6e\t vel.z = % 9.6e\t displ = % 9.6e", j+1
                           , mod->sw->d3->vel[j].x
                           , mod->sw->d3->vel[j].y
                           , mod->sw->d3->vel[j].z
                           , mod->sw->d3->displacement[j]);
                }
                printf("\n");
            }
#endif
            node_count = node_count + mod->grid->nnodes;
        }
    }
    
#ifdef _DEBUG
    if (DEBUG_INTERFACE){
        int j, k, l;
        SMODEL *mod1;
        SMODEL *mod2;
        SINTERFACE *ifce;
        SINTERFACE_NODELIST *ndlist;
        
        printf("\n******************************\nPrinting Interface Solution\n");
        /* Print to check equality of interface variables */
        for (j=0; j<sm->NumInterfaces; j++){
            ifce = &(sm->interface[j]);
            ifce = &(sm->interface[j]);
            mod1 = &(sm->submodel[ifce->model_id[0]]);
            mod2 = &(sm->submodel[ifce->model_id[1]]);
            
            if (sm->submodel[ifce->model_id[0]].flag.SW2_FLOW){
                if(sm->submodel[ifce->model_id[1]].flag.SW2_FLOW){   /* 2D 2D coupling */
                    /* Do nothing */
                }
                else{                                                /* 2D 3D coupling */
                    for (k=0; k<ifce->NumNodeColumns; k++){
                        ndlist = &(ifce->nodelist[k]);
                        printf("\nColumn %i :",k+1);
                        printf("\n  Model 1 (2D) :");
                        printf("\n    Node % 5i : vel.x = % 9.6e\t vel.y = % 9.6e\t vel.z = % 9.6e\t depth = % 9.6e",
                               ndlist->couplednodes[0][0]+1,
                               mod1->sw->d2->vel[ndlist->couplednodes[0][0]].x,
                               mod1->sw->d2->vel[ndlist->couplednodes[0][0]].y,
                               0.0,
                               mod1->sw->d2->head[ndlist->couplednodes[0][0]]
                               );
                        printf("\n  Model 2 (3D) :");
                        for (l=0; l<ndlist->size[1]; l++){
                            printf("\n    Node % 5i : vel.x = % 9.6e\t vel.y = % 9.6e\t vel.z = % 9.6e\t depth = % 9.6e",
                                   ndlist->couplednodes[1][l]+1,
                                   mod2->sw->d3->vel[ndlist->couplednodes[1][l]].x,
                                   mod2->sw->d3->vel[ndlist->couplednodes[1][l]].y,
                                   mod2->sw->d3->vel[ndlist->couplednodes[1][l]].z,
                                   mod2->sw->d3->displacement[ndlist->couplednodes[1][l]]
                                   + (mod2->grid->node[ndlist->couplednodes[1][l]].z
                                      - mod2->grid->node[ndlist->couplednodes[1][ndlist->size[1]-1]].z)
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
                            printf("\n    Node % 5i : vel.x = % 9.6e\t vel.y = % 9.6e\t vel.z = % 9.6e\t depth = % 9.6e",
                                   ndlist->couplednodes[0][l]+1,
                                   mod1->sw->d3->vel[ndlist->couplednodes[0][l]].x,
                                   mod1->sw->d3->vel[ndlist->couplednodes[0][l]].y,
                                   mod1->sw->d3->vel[ndlist->couplednodes[0][l]].z,
                                   mod1->sw->d3->displacement[ndlist->couplednodes[0][l]]
                                   + (mod1->grid->node[ndlist->couplednodes[0][l]].z
                                      - mod1->grid->node[ndlist->couplednodes[0][ndlist->size[0]-1]].z)
                                   );
                        }
                        printf("\n  Model 2 (3D) :");
                        for (l=0; l<ndlist->size[1]; l++){
                            printf("\n    Node % 5i : vel.x = % 9.6e\t vel.y = % 9.6e\t vel.z = % 9.6e\t depth = % 9.6e",
                                   ndlist->couplednodes[1][l]+1,
                                   mod2->sw->d3->vel[ndlist->couplednodes[1][l]].x,
                                   mod2->sw->d3->vel[ndlist->couplednodes[1][l]].y,
                                   mod2->sw->d3->vel[ndlist->couplednodes[1][l]].z,
                                   mod2->sw->d3->displacement[ndlist->couplednodes[1][l]]
                                   + (mod2->grid->node[ndlist->couplednodes[1][l]].z
                                      - mod2->grid->node[ndlist->couplednodes[1][ndlist->size[1]-1]].z)
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
                        printf("\n    Node % 5i : vel.x = % 9.6e\t vel.y = % 9.6e\t vel.z = % 9.6e\t depth = % 9.6e",
                               ndlist->couplednodes[1][0]+1,
                               mod2->sw->d2->vel[ndlist->couplednodes[1][0]].x,
                               mod2->sw->d2->vel[ndlist->couplednodes[1][0]].y,
                               0.0,
                               mod2->sw->d2->head[ndlist->couplednodes[1][0]]
                               );
                        printf("\n  Model 1 (3D) :");
                        for (l=0; l<ndlist->size[0]; l++){
                            printf("\n    Node % 5i : vel.x = % 9.6e\t vel.y = % 9.6e\t vel.z = % 9.6e\t depth = % 9.6e",
                                   ndlist->couplednodes[0][l]+1,
                                   mod1->sw->d3->vel[ndlist->couplednodes[0][l]].x,
                                   mod1->sw->d3->vel[ndlist->couplednodes[0][l]].y,
                                   mod1->sw->d3->vel[ndlist->couplednodes[0][l]].z,
                                   mod1->sw->d3->displacement[ndlist->couplednodes[0][l]]
                                   + (mod1->grid->node[ndlist->couplednodes[0][l]].z
                                      - mod1->grid->node[ndlist->couplednodes[0][ndlist->size[0]-1]].z)
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
#endif
#endif
}
