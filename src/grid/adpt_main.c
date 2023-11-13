/* adapts the mesh */
/*  A. Liu and B. Joe, Siam J. Sci. Comput. Vol. 16, No. 6, pp 1269-1291, Nov 1995
 We don't follow the cite exactly.  We do borrow extensively from their ideas.
 We use longest edge with a non standard metric:
 1) compare the level of the edge - the higher of the two node levels of the edge.
 2) then compare lengths.
 To avoid possible conflicts from roundoff in the calculation of edge lengths, the edges are shorted
 so that the length comparison becomes an integer comparison.
 */

#include "global_header.h"

static int DEBUG = ON;
static int DEBUG_REFINE = ON;
static int DEBUG_UNREFINE = ON;

void adpt_main(SMODEL *mod) {
    
    char str[25];
    
    //print_parents_to_file(mod->grid, str);
    //print_sw2_to_file(mod->sw->d2,mod->grid,str);
    //print_grid_to_file(mod->grid,str);
    
    /* static variables */
    static int cycle = 0;         /* counter for the refinement cycles */
    
    /* store old nnodes and element counts */
    mod->grid->nnodes_old = mod->grid->nnodes;
    mod->grid->nelems2d_old = mod->grid->nelems2d;
    mod->grid->nelems3d_old = mod->grid->nelems3d;
    
    /* local variables */
    int change_nnode=0;
    int REF_FLAG_RECV = 0;
    int UNREF_FLAG_RECV = 0;
    
#ifdef _MESSG
    int i;
    int min_nnodes = 0;
    int max_nnodes = 0;
    double ratio = 0.0;
    double tval=0, tsum=0, tmax=0;      /* times */
#endif
    
    int init_nnode = 0;     // the initial number of my_nnodes over my pe
    int init_tot_nnode = 0; // the initial total my_nnodes over all pe's
    int tot_nnode = 0;      // the total my_nnodes over all pe's after adaption
    FILE *fp;
    SGRID *g=mod->grid;
    SCON *con=mod->con;
    SSW_3D *sw;
    SNS_3D *ns;
    
    if (mod->flag.SW_FLOW
#ifdef _ADH_GROUNDWATER
#ifdef _DWGW_COUPLING
            && (mod->flag.GW_FLOW==OFF)
#endif
#endif
            ) {
        sw = mod->sw->d3;
    } else if (mod->flag.NS_FLOW) {
        ns = mod->ns->d3;
    }
    
    int itrns;
    char filename[MAXLINE];
#ifdef _MESSG
    int myid = g->smpi->myid;
#else
    int myid = 0;
#endif
    
    /* get flags */
    error_main(mod);
    
    /* check on the node information */
    mod->flag.ADAPTED_THE_GRID = NO;
    init_nnode = mod->grid->my_nnodes;
    tot_nnode = init_nnode;
#ifdef _MESSG
    if(DEBUG) tot_nnode = messg_isum(init_nnode, mod->grid->smpi->ADH_COMM);
#endif
    init_tot_nnode = tot_nnode;
    
    
#ifdef _MESSG
    int REF_FLAG = mod->flag.GRID_REFINED;
    int UNREF_FLAG = mod->flag.GRID_UNREFINED;
    MPI_Allreduce(&REF_FLAG, &REF_FLAG_RECV, 1, MPI_INT, MPI_MAX, mod->grid->smpi->ADH_COMM);
    MPI_Allreduce(&UNREF_FLAG, &UNREF_FLAG_RECV, 1, MPI_INT, MPI_MAX, mod->grid->smpi->ADH_COMM);
#else
    REF_FLAG_RECV = mod->flag.GRID_REFINED;
    UNREF_FLAG_RECV = mod->flag.GRID_UNREFINED;
#endif
    
    int ierr;
    
    
    /********************************************************************/
    /* unrefinement pass ************************************************/
    if(UNREF_FLAG_RECV == YES) {
        
        if (mod->grid->type == COLUMNAR) {
            column_adpt_unref(mod);
        } else {
            adpt_unref(mod);
        }
        
#ifdef _MESSG
        if(DEBUG) tot_nnode = messg_isum(mod->grid->my_nnodes, mod->grid->smpi->ADH_COMM);
        i = abs(init_nnode - mod->grid->my_nnodes);
        change_nnode = messg_imax(i, mod->grid->smpi->ADH_COMM);
#else
        tot_nnode = mod->grid->my_nnodes;
        change_nnode = abs(init_nnode - mod->grid->my_nnodes);
#endif
        
        if (change_nnode != 0) {
            mod->flag.ADAPTED_THE_GRID = YES;
        }
        
        if (DEBUG) {
            if (DEBUG_UNREFINE){
                if (mod->flag.SW2_FLOW
#ifdef _ADH_GROUNDWATER
#ifdef _DWGW_COUPLING
                        && (mod->flag.GW_FLOW==OFF)
#endif
#endif
                        ){
                    printf("\nUNREFINEMENT RESULTS\n");
                    
                    for (ierr=0; ierr<mod->grid->nnodes; ierr++) {
                        printf("node: %d string: %d parent1: %d parent2: %d head: %20.10f vel: %20.10f %20.10f\n",mod->grid->node[ierr].id, mod->grid->node[ierr].string,mod->grid->node[ierr].parent[0],mod->grid->node[ierr].parent[1],mod->sw->d2->head[ierr],mod->sw->d2->vel[ierr].x,mod->sw->d2->vel[ierr].y);
                    }
                    
                    for (ierr=0; ierr<mod->grid->nelems1d; ierr++) {
                        printf("elem1d: %d nodes: %d %d string: %d mat: %d djac: %20.10f\n",mod->grid->elem1d[ierr].id, mod->grid->elem1d[ierr].nodes[0], mod->grid->elem1d[ierr].nodes[1], mod->grid->elem1d[ierr].string,mod->grid->elem1d[ierr].mat,mod->grid->elem1d[ierr].djac);
                    }
                    
                    for (ierr=0; ierr<mod->grid->nelems2d; ierr++) {
                        if (mod->grid->elem2d[ierr].id != ierr) {printf("WHOA id and loop are != for 2d element ierr = %d id = %d\n",ierr,mod->grid->elem2d[ierr].id); exit(-1);}
                        printf("elem2d: %d nodes: %d %d %d string: %d mat: %d djac: %20.10f error: %20.10f\n",mod->grid->elem2d[ierr].id, mod->grid->elem2d[ierr].nodes[0], mod->grid->elem2d[ierr].nodes[1], mod->grid->elem2d[ierr].nodes[2], mod->grid->elem2d[ierr].string,mod->grid->elem2d[ierr].mat,mod->grid->elem2d[ierr].djac,mod->grid->elem_error[ierr]);
                    }
                    
                    //exit(-1);
                }
                else if (mod->flag.SW3_FLOW || mod->flag.NS3_FLOW){
                    printf("\nUNREFINEMENT RESULTS\n");
                    
                    if (mod->grid->havePrisms==TRUE){
                        for (ierr=0; ierr<mod->grid->nnodes; ierr++) {
                            printf("node: %d string: %d parent1: %d parent2: %d dpl: %20.10f  vel: %20.10f %20.10f %20.10f\n",mod->grid->node[ierr].id, mod->grid->node[ierr].string,mod->grid->node[ierr].parent[0],mod->grid->node[ierr].parent[1],mod->sw->d3->displacement[ierr],mod->sw->d3->vel[ierr].x,mod->sw->d3->vel[ierr].y, mod->sw->d3->vel[ierr].z);
                        }
                        
                        for (ierr=0; ierr<mod->grid->nelems2d; ierr++) {
                            /* Gajanan gkc - May not work for triangular faces. Segmentation fault possible. */
                            printf("elem2d: %d nodes: %d %d %d %d string: %d mat: %d djac: %20.10f\n",mod->grid->elem2d[ierr].id, mod->grid->elem2d[ierr].nodes[0], mod->grid->elem2d[ierr].nodes[1], mod->grid->elem2d[ierr].nodes[2], mod->grid->elem2d[ierr].nodes[3], mod->grid->elem2d[ierr].string,mod->grid->elem2d[ierr].mat,mod->grid->elem2d[ierr].djac);
                        }
                        
                        for (ierr=0; ierr<mod->grid->nelems3d; ierr++) {
                            if (mod->grid->elem3d[ierr].id != ierr) {printf("WHOA id and loop are != for 3d element ierr = %d id = %d\n",ierr,mod->grid->elem3d[ierr].id); exit(-1);}
                            printf("elem3d: %d nodes: %d %d %d %d %d %d string: %d mat: %d djac: %20.10f error: %20.10f\n",mod->grid->elem3d[ierr].id, mod->grid->elem3d[ierr].nodes[0], mod->grid->elem3d[ierr].nodes[1], mod->grid->elem3d[ierr].nodes[2], mod->grid->elem3d[ierr].nodes[3], mod->grid->elem3d[ierr].nodes[4], mod->grid->elem3d[ierr].nodes[5], mod->grid->elem3d[ierr].string,mod->grid->elem3d[ierr].mat,mod->grid->elem3d[ierr].djac,mod->grid->elem_error[ierr]);
                        }
                    }
                    //exit(-1);
                }
#ifdef _ADH_GROUNDWATER
                else if (mod->flag.GW_FLOW){
                    //printf("\nUNREFINEMENT RESULTS\n");
                    //
                    //for (ierr=0; ierr<mod->grid->nnodes; ierr++) {
                    //    printf("node: %d string: %d parent1: %d parent2: %d phead: %20.10f \n",mod->grid->node[ierr].id, mod->grid->node[ierr].string,mod->grid->node[ierr].parent[0],mod->grid->node[ierr].parent[1],mod->sgw->gw_phead[ierr]);
                    //}
                    //
                    //for (ierr=0; ierr<mod->grid->nelems2d; ierr++) {
                    //    printf("elem2d: %d nodes: %d %d %d string: %d mat: %d djac: %20.10f\n",mod->grid->elem2d[ierr].id, mod->grid->elem2d[ierr].nodes[0], mod->grid->elem2d[ierr].nodes[1], mod->grid->elem2d[ierr].nodes[2], mod->grid->elem2d[ierr].string,mod->grid->elem2d[ierr].mat,mod->grid->elem2d[ierr].djac);
                    //}
                    //
                    //for (ierr=0; ierr<mod->grid->nelems3d; ierr++) {
                    //    if (mod->grid->elem3d[ierr].id != ierr) {printf("WHOA id and loop are != for 3d element ierr = %d id = %d\n",ierr,mod->grid->elem3d[ierr].id); exit(-1);}
                    //    printf("elem3d: %d nodes: %d %d %d %d string: %d mat: %d djac: %20.10f error: %20.10f\n",mod->grid->elem3d[ierr].id, mod->grid->elem3d[ierr].nodes[0], mod->grid->elem3d[ierr].nodes[1], mod->grid->elem3d[ierr].nodes[2], mod->grid->elem3d[ierr].nodes[3], mod->grid->elem3d[ierr].string,mod->grid->elem3d[ierr].mat,mod->grid->elem3d[ierr].djac,mod->grid->elem_error[ierr]);
                    //}
                }
#endif
            }
#ifdef _MESSG
            if (mod->grid->smpi->myid <= 0) {
#endif
                printf("\nUNRefinement statistics\n");
                printf("max_nnodes: %d nnodes: %d total nodes:  %8d (before)  %8d (after)  %8d (maximum change on a pe)  %6d (cycle)\n", mod->grid->max_nnodes,mod->grid->nnodes,init_tot_nnode, tot_nnode, change_nnode, cycle);
                printf("max_nelems3d: %d nelems3d: %d nelems3d_old: %d\n", mod->grid->max_nelems3d,mod->grid->nelems3d, mod->grid->nelems3d_old);
                printf("max_nelems2d: %d nelems2d: %d nelems2d_old: %d\n", mod->grid->max_nelems2d,mod->grid->nelems2d, mod->grid->nelems2d_old);
                printf("max_nelems1d: %d nelems1d: %d\n", mod->grid->max_nelems1d,mod->grid->nelems1d);
                
                
#ifdef _MESSG
            }
            printf("id: %4d  nodes:  before: %8d  after: %8d  cycle: %8d\n", mod->grid->smpi->myid, init_nnode, mod->grid->my_nnodes, cycle);
            messg_barrier(mod->grid->smpi->ADH_COMM);
#endif
        }
        
    }
    
    
    init_nnode = mod->grid->my_nnodes;
    init_tot_nnode = tot_nnode;
    
    
    /********************************************************************/
    /* refinement pass **************************************************/
    if(REF_FLAG_RECV == YES) {
        
        if (mod->grid->type == COLUMNAR) {
            cycle = column_adpt_ref(mod, cycle);
        } else {
            cycle = adpt_ref(mod, cycle);
        }
        
        int ie,i;
        for(ie = 0; ie < mod->grid->nelems3d; ie++) {
            if (mod->grid->node[mod->grid->elem3d[ie].nodes[0]].x < UNSET_FLT + 1e-6 ||
                mod->grid->node[mod->grid->elem3d[ie].nodes[0]].y < UNSET_FLT + 1e-6 ||
                mod->grid->node[mod->grid->elem3d[ie].nodes[0]].z < UNSET_FLT + 1e-6) {
                selem3d_printScreen(&mod->grid->elem3d[ie]);
                for (i=0; i<mod->grid->elem3d[ie].nnodes; i++) printf("nd: %d \t x: %20.10f y: %20.10f z: %20.10f\n",i,
                                                                      mod->grid->node[mod->grid->elem3d[ie].nodes[i]].x,
                                                                      mod->grid->node[mod->grid->elem3d[ie].nodes[i]].y,
                                                                      mod->grid->node[mod->grid->elem3d[ie].nodes[i]].z);
                tl_error("One of the nodes has UNSET locations after REfinement.");
            }
        }
        
#ifdef _MESSG
        if(DEBUG) tot_nnode = messg_isum(mod->grid->my_nnodes, mod->grid->smpi->ADH_COMM);
        i = abs(mod->grid->my_nnodes - init_nnode);
        change_nnode = messg_imax(i, mod->grid->smpi->ADH_COMM);
#else
        tot_nnode = mod->grid->my_nnodes;
        change_nnode = abs(mod->grid->my_nnodes - init_nnode);
#endif
        
        if (change_nnode != 0) {
            mod->flag.ADAPTED_THE_GRID = YES;
            
        }
        
        if (DEBUG) {
            if (DEBUG_REFINE) {
                if (mod->flag.SW2_FLOW
#ifdef _ADH_GROUNDWATER
#ifdef _DWGW_COUPLING
                        && (mod->flag.GW_FLOW==OFF)
#endif
#endif
                        ){
                    printf("\nREFINEMENT RESULTS\n");
                    
                    for (ierr=0; ierr<mod->grid->nnodes; ierr++) {
                        printf("node: %d string: %d parent1: %d parent2: %d head: %20.10f  vel: %20.10f %20.10f\n",mod->grid->node[ierr].id, mod->grid->node[ierr].string,mod->grid->node[ierr].parent[0],mod->grid->node[ierr].parent[1],mod->sw->d2->head[ierr],mod->sw->d2->vel[ierr].x,mod->sw->d2->vel[ierr].y);
                    }
                    
                    for (ierr=0; ierr<mod->grid->nelems1d; ierr++) {
                        printf("elem1d: %d nodes: %d %d string: %d mat: %d djac: %20.10f\n",mod->grid->elem1d[ierr].id, mod->grid->elem1d[ierr].nodes[0], mod->grid->elem1d[ierr].nodes[1], mod->grid->elem1d[ierr].string,mod->grid->elem1d[ierr].mat,mod->grid->elem1d[ierr].djac);
                    }
                    
                    for (ierr=0; ierr<mod->grid->nelems2d; ierr++) {
                        printf("elem2d: %d nodes: %d %d %d string: %d mat: %d djac: %20.10f error: %20.10f\n",mod->grid->elem2d[ierr].id, mod->grid->elem2d[ierr].nodes[0], mod->grid->elem2d[ierr].nodes[1], mod->grid->elem2d[ierr].nodes[2], mod->grid->elem2d[ierr].string,mod->grid->elem2d[ierr].mat,mod->grid->elem2d[ierr].djac,mod->grid->elem_error[ierr]);
                    }
                }
                else if (mod->flag.SW3_FLOW || mod->flag.NS3_FLOW){
                    printf("\nREFINEMENT RESULTS\n");
                    
                    if (mod->grid->havePrisms==TRUE){
                        for (ierr=0; ierr<mod->grid->nnodes; ierr++) {
                            printf("node: %d string: %d parent1: %d parent2: %d dpl: %20.10f  vel: %20.10f %20.10f %20.10f\n",mod->grid->node[ierr].id, mod->grid->node[ierr].string,mod->grid->node[ierr].parent[0],mod->grid->node[ierr].parent[1],mod->sw->d3->displacement[ierr],mod->sw->d3->vel[ierr].x,mod->sw->d3->vel[ierr].y, mod->sw->d3->vel[ierr].z);
                        }
                        
                        for (ierr=0; ierr<mod->grid->nelems2d; ierr++) {
                            /* Gajanan gkc - May not work for triangular faces. Segmentation fault possible. */
                            printf("elem2d: %d nodes: %d %d %d %d string: %d mat: %d djac: %20.10f\n",mod->grid->elem2d[ierr].id, mod->grid->elem2d[ierr].nodes[0], mod->grid->elem2d[ierr].nodes[1], mod->grid->elem2d[ierr].nodes[2], mod->grid->elem2d[ierr].nodes[3], mod->grid->elem2d[ierr].string,mod->grid->elem2d[ierr].mat,mod->grid->elem2d[ierr].djac);
                        }
                        
                        for (ierr=0; ierr<mod->grid->nelems3d; ierr++) {
                            printf("elem3d: %d nodes: %d %d %d %d %d %d string: %d mat: %d djac: %20.10f error: %20.10f\n",mod->grid->elem3d[ierr].id, mod->grid->elem3d[ierr].nodes[0], mod->grid->elem3d[ierr].nodes[1], mod->grid->elem3d[ierr].nodes[2], mod->grid->elem3d[ierr].nodes[3], mod->grid->elem3d[ierr].nodes[4], mod->grid->elem3d[ierr].nodes[5], mod->grid->elem3d[ierr].string,mod->grid->elem3d[ierr].mat,mod->grid->elem3d[ierr].djac,mod->grid->elem_error[ierr]);
                        }
                    }
                    //exit(-1);
                }
#ifdef _ADH_GROUNDWATER
                else if (mod->flag.GW_FLOW){
                    //printf("\nREFINEMENT RESULTS\n");
                    //
                    //for (ierr=0; ierr<mod->grid->nnodes; ierr++) {
                    //    printf("node: %d string: %d parent1: %d parent2: %d phead: %20.10f \n",mod->grid->node[ierr].id, mod->grid->node[ierr].string,mod->grid->node[ierr].parent[0],mod->grid->node[ierr].parent[1],mod->sgw->gw_phead[ierr]);
                    //}
                    //
                    //for (ierr=0; ierr<mod->grid->nelems2d; ierr++) {
                    //    printf("elem2d: %d nodes: %d %d %d string: %d mat: %d djac: %20.10f\n",mod->grid->elem2d[ierr].id, mod->grid->elem2d[ierr].nodes[0], mod->grid->elem2d[ierr].nodes[1], mod->grid->elem2d[ierr].nodes[2], mod->grid->elem2d[ierr].string,mod->grid->elem2d[ierr].mat,mod->grid->elem2d[ierr].djac);
                    //}
                    //
                    //for (ierr=0; ierr<mod->grid->nelems3d; ierr++) {
                    //    if (mod->grid->elem3d[ierr].id != ierr) {printf("WHOA id and loop are != for 3d element ierr = %d id = %d\n",ierr,mod->grid->elem3d[ierr].id); exit(-1);}
                    //    printf("elem3d: %d nodes: %d %d %d %d string: %d mat: %d djac: %20.10f error: %20.10f\n",mod->grid->elem3d[ierr].id, mod->grid->elem3d[ierr].nodes[0], mod->grid->elem3d[ierr].nodes[1], mod->grid->elem3d[ierr].nodes[2], mod->grid->elem3d[ierr].nodes[3], mod->grid->elem3d[ierr].string,mod->grid->elem3d[ierr].mat,mod->grid->elem3d[ierr].djac,mod->grid->elem_error[ierr]);
                    //}
                }
#endif
            }
#ifdef _MESSG
            if (mod->grid->smpi->myid <= 0) {
#endif
                printf("\nRefinement statistics\n");
                printf("max_nnodes: %d nnodes: %d, total nodes:  %8d (before)  %8d (after)  %8d (maximum change on a pe)  %6d (cycle)\n", mod->grid->max_nnodes, mod->grid->nnodes, init_tot_nnode, tot_nnode, change_nnode, cycle);
                printf("max_nelems3d: %d nelems3d: %d nelems3d_old: %d\n", mod->grid->max_nelems3d,mod->grid->nelems3d, mod->grid->nelems3d_old);
                printf("max_nelems2d: %d nelems2d: %d nelems2d_old: %d\n", mod->grid->max_nelems2d,mod->grid->nelems2d, mod->grid->nelems2d_old);
                printf("max_nelems1d: %d nelems1d: %d\n", mod->grid->max_nelems1d,mod->grid->nelems1d);
                
                
#ifdef _MESSG
            }
            printf("REFINE id: %4d  nodes:  before: %8d  after: %8d  cycle: %8d\n", mod->grid->smpi->myid, init_nnode, mod->grid->my_nnodes, cycle);
            messg_barrier(mod->grid->smpi->ADH_COMM);
#endif
        }
    }
    
    
    /*   Determine if repartitioning should be performed */
    /*   If no adaption occurred, do not repartition */
    /*   If the node ratio over the processors is greater than 1/10, repartition */
    
    if (mod->flag.ADAPTED_THE_GRID == YES) {
        if (mod->grid->ndim==2) {
            if (mod->flag.SW_FLOW) {
                mod->sw->d2->elem_rhs_realloc = 0;
            } else if (mod->flag.NS_FLOW) {
                mod->ns->d2->elem_rhs_realloc = 0;
            }
        } else if (mod->grid->ndim==3) {
            if (mod->flag.SW_FLOW
#ifdef _ADH_GROUNDWATER
                    && (mod->flag.GW_FLOW==OFF)
#endif
                    ) {
                mod->sw->d3->elem_rhs_realloc = 0;
#ifdef _ADH_GROUNDWATER
#ifdef _DWGW_COUPLING
            } else if (mod->flag.SW_FLOW && (mod->flag.GW_FLOW==ON) ) {
                mod->sw->d2->elem_rhs_realloc = 0;
#endif
#endif
            } else if (mod->flag.NS_FLOW) {
                mod->ns->d3->elem_rhs_realloc = 0;
            }
        }
        /*might only need this for outputting adapted grid but for now be safe */
        adpt_fix_global_ids(mod->grid);
        
        
        
#ifdef _MESSG
        if(mod->grid->interface==0){
          if(mod->grid->part_smpi != NULL){
            max_nnodes = messg_imax(mod->grid->my_nnodes, mod->grid->part_smpi->ADH_COMM);
            min_nnodes = messg_imin(mod->grid->my_nnodes, mod->grid->part_smpi->ADH_COMM);
          }else{
            max_nnodes = messg_imax(mod->grid->my_nnodes, mod->grid->smpi->ADH_COMM);
            min_nnodes = messg_imin(mod->grid->my_nnodes, mod->grid->smpi->ADH_COMM);
          }
          ratio = (1.0 * min_nnodes) / max_nnodes;
        } else {
          ratio=1.0;
        }
        ratio = messg_dmin(ratio, mod->grid->smpi->ADH_COMM);
        
        
        if (mod->grid->smpi->myid <= 0 && DEBUG == ON) {
            printf("\nAdaptive Repartitioning Analysis:\n");
            printf("number of nodes:  %8d (max)  %8d (min)  %8.2f (min/max)\n", max_nnodes, min_nnodes, ratio);
        }
        if (ratio < 0.80) {
            
            if (DEBUG) {
                if(mod->grid->smpi->myid <= 0)printf("\n\nUsing Parmetis Adaptive Repartitioning.\n");
                tval = MPI_Wtime();
            }
            
            
            partition_main(mod,0);
            
            
            if (DEBUG) {
                tval = MPI_Wtime() - tval;
                MPI_Allreduce(&tval, &tsum, 1, MPI_DOUBLE, MPI_SUM, mod->grid->smpi->ADH_COMM);
                MPI_Allreduce(&tval, &tmax, 1, MPI_DOUBLE, MPI_MAX, mod->grid->smpi->ADH_COMM);
                if (mod->grid->smpi->myid <= 0) printf("timings:  average: %8.4f  max: %8.4f  max/avg: %6.2f\n", tsum / mod->grid->smpi->npes, tmax, tmax * mod->grid->smpi->npes / tsum);
            }
        }
        
        comm_set_keys(mod->grid);
        
        //comm_check(mod->grid);
        comm_update_snode(mod->grid);
#endif
        int ie3d = UNSET_INT, icol=UNSET_INT;
        ID_LIST_ITEM *ptr;
        if(mod->grid->ndim == 3) {
            if (mod->grid->type == COLUMNAR) {
                build_columns(mod->grid, NO);
                for(icol=0; icol<mod->grid->ncolumns; icol++)  {
                    ptr = mod->grid->column_list[icol];
                    while(ptr->next != NULL) {
                        ie3d = ptr->id;
                        mod->grid->elem3d[ie3d].icol = icol;
                        mod->grid->elem3d[ie3d].elem2d_sur = mod->grid->elem2d_sur[icol];
                        mod->grid->elem3d[ie3d].elem2d_bed = mod->grid->elem2d_bed[icol];
                        ptr = ptr->next;
                    }
                }
            } else {
                classify_2d_elements(&mod->grid);
            }
        }
#ifdef _MESSG
        if (mod->flag.SW2_FLOW) comm_update_sw2(mod->sw->d2, mod->grid);
        if(mod->flag.SW3_FLOW) {
            comm_update_sw3(mod->sw->d3, mod->grid);
            comm_update_sw3_surface(mod->sw->d3, mod->grid);
        }
        
#ifdef _SEDIMENT
        if(mod->sed != NULL) comm_update_sed(mod->sed, mod->grid, mod->nconti, mod->nsed, mod->nlayers);
#endif
        if(mod->ntransport > 0) comm_update_con(mod->con, mod->grid, mod->ntransport);
#endif
    }
    // cjt :: make sure flag is on when # of nodes/elements change
    if (mod->grid->nnodes != mod->grid->nnodes_old) assert(mod->flag.ADAPTED_THE_GRID == YES);
    if (mod->grid->ndim == 2 && mod->grid->nelems2d != mod->grid->nelems2d_old) assert(mod->flag.ADAPTED_THE_GRID == YES);
    if (mod->grid->ndim == 3 && mod->grid->nelems3d != mod->grid->nelems3d_old) assert(mod->flag.ADAPTED_THE_GRID == YES);
    
#ifdef _DEBUG
    //tl_check_all_pickets(__FILE__,__LINE__);
#endif
}
