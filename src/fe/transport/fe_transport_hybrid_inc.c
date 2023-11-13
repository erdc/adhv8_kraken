/* This routine loops over submodels to increment the transport solution */
#include "global_header.h"

/***********************************************************/
/***********************************************************/
static int DEBUG = OFF;

/***********************************************************/
/***********************************************************/

void fe_transport_hybrid_inc(SSUPER_MODEL *sm, int dum) {
#ifndef _PETSC
    int i;	    		/* loop counter */
    // aliases
    SMODEL *mod;
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            fe_transport_inc(sm, i);
        }
    }

#ifdef _DEBUG

    if (DEBUG){
        int j, k, l;
        SMODEL *mod1;
        SMODEL *mod2;
        SINTERFACE *ifce;
        SINTERFACE_NODELIST *ndlist;

        /* Print to check equality of interface variables */
        for (j=0; j<sm->NumInterfaces; j++){
            ifce = &(sm->interface[j]);
            ifce = &(sm->interface[j]);
            mod1 = &(sm->submodel[ifce->model_id[0]]);
            mod2 = &(sm->submodel[ifce->model_id[1]]);

            if (sm->submodel[ifce->model_id[0]].flag.SW2_FLOW){
                if(sm->submodel[ifce->model_id[1]].flag.SW2_FLOW){   /* 2D 2D coupling */
                    for (k=0; k<ifce->NumNodeColumns; k++){
                        ndlist = &(ifce->nodelist[k]);
            	    printf("\nColumn %i :",k+1);
//                        printf("\nMod1 (2D) node %5i : vel.x = %9.6e\t vel.y = %9.6e\t depth = %9.6e",
//                               ndlist->couplednodes[0][0]+1,
//                               mod1->sw->d2->vel[ndlist->couplednodes[0][0]].x,
//                               mod1->sw->d2->vel[ndlist->couplednodes[0][0]].y,
//                               mod1->sw->d2->head[ndlist->couplednodes[0][0]]
//                              );
//                        printf("\nMod2 (2D) node %5i : vel.x = %9.6e\t vel.y = %9.6e\t depth = %9.6e",
//                               ndlist->couplednodes[1][0]+1,
//                               mod2->sw->d2->vel[ndlist->couplednodes[1][0]].x,
//                               mod2->sw->d2->vel[ndlist->couplednodes[1][0]].y,
//                               mod2->sw->d2->head[ndlist->couplednodes[1][0]]
//                              );
                        printf("\n");
                    }
                    printf("\n");
                }
                else{                                                /* 2D 3D coupling */
                    for (k=0; k<ifce->NumNodeColumns; k++){
                        ndlist = &(ifce->nodelist[k]);
            	    printf("\nColumn %i :",k+1);
//                        printf("\nMod1 (2D) node %5i : vel.x = %9.6e\t vel.y = %9.6e\t depth = %9.6e",
//                               ndlist->couplednodes[0][0]+1,
//                               mod1->sw->d2->vel[ndlist->couplednodes[0][0]].x,
//                               mod1->sw->d2->vel[ndlist->couplednodes[0][0]].y,
//                               mod1->sw->d2->head[ndlist->couplednodes[0][0]]
//                              );
//                        for (l=0; l<ndlist->size[1]; l++){
//                            printf("\nMod2 (3D) node %5i : vel.x = %9.6e\t vel.y = %9.6e\t displ = %9.6e",
//                                   ndlist->couplednodes[1][l]+1,
//                                   mod2->sw->d3->vel[ndlist->couplednodes[1][l]].x,
//                                   mod2->sw->d3->vel[ndlist->couplednodes[1][l]].y,
//                                   mod2->sw->d3->displacement[ndlist->couplednodes[1][l]]
//                                  );
//                        }
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
//                        for (l=0; l<ndlist->size[0]; l++){
//                            printf("\nMod1 (3D) node %5i : vel.x = %9.6e\t vel.y = %9.6e\t displ = %9.6e",
//                                   ndlist->couplednodes[0][l]+1,
//                                   mod1->sw->d3->vel[ndlist->couplednodes[0][l]].x,
//                                   mod1->sw->d3->vel[ndlist->couplednodes[0][l]].y,
//                                   mod1->sw->d3->displacement[ndlist->couplednodes[0][l]]
//                                  );
//                        }
//                        for (l=0; l<ndlist->size[1]; l++){
//                            printf("\nMod2 (3D) node %5i : vel.x = %9.6e\t vel.y = %9.6e\t displ = %9.6e",
//                                   ndlist->couplednodes[1][l]+1,
//                                   mod2->sw->d3->vel[ndlist->couplednodes[1][l]].x,
//                                   mod2->sw->d3->vel[ndlist->couplednodes[1][l]].y,
//                                   mod2->sw->d3->displacement[ndlist->couplednodes[1][l]]
//                                  );
//                        }
                        printf("\n");
                    }
                    printf("\n");
                }
                else{                                                /* 3D 2D coupling */
                    for (k=0; k<ifce->NumNodeColumns; k++){
                        ndlist = &(ifce->nodelist[k]);
            	    printf("\nColumn %i :",k+1);
//                        printf("\nMod2 (2D) node %5i : vel.x = %9.6e\t vel.y = %9.6e\t depth = %9.6e",
//                               ndlist->couplednodes[1][0]+1,
//                               mod2->sw->d2->vel[ndlist->couplednodes[1][0]].x,
//                               mod2->sw->d2->vel[ndlist->couplednodes[1][0]].y,
//                               mod2->sw->d2->head[ndlist->couplednodes[1][0]]
//                              );
//                        for (l=0; l<ndlist->size[1]; l++){
//                            printf("\nMod1 (3D) node %5i : vel.x = %9.6e\t vel.y = %9.6e\t displ = %9.6e",
//                                   ndlist->couplednodes[0][l]+1,
//                                   mod1->sw->d3->vel[ndlist->couplednodes[0][l]].x,
//                                   mod1->sw->d3->vel[ndlist->couplednodes[0][l]].y,
//                                   mod1->sw->d3->displacement[ndlist->couplednodes[0][l]]
//                                  );
//                        }
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
