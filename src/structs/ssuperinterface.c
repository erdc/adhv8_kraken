#include "global_header.h"
//void couple_edge_to_edge(SSUPER_MODEL *sm);

/***********************************************************/
/***********************************************************/
static int DEBUG = OFF;

/***********************************************************/
/***********************************************************/

void ssuperinterface_alloc_init(SSUPER_INTERFACE *si, int NumEdges, int *sm_id, int *mod_id, int *edgestring) {
    int i,j;
    si->NumEdges = NumEdges;
    for (i=0;i<2;i++){
        si->sm_id[i] = sm_id[i];
        si->model_id[i] = mod_id[i];
        si->bdrystrings[i] = edgestring[i];
    }
    si->bdrylist = (SSUP_IFCE_BOUNDARY_LIST *) tl_alloc(sizeof(SSUP_IFCE_BOUNDARY_LIST), NumEdges);
    for (i=0; i<NumEdges; i++){
        for (j=0;j<2;j++){
            si->bdrylist[i].surfnode[j][0] = UNSET_INT;
            si->bdrylist[i].surfnode[j][1] = UNSET_INT;
            si->bdrylist[i].size[j] = 0;
        }
    }
}

/***********************************************************/
/***********************************************************/

void ssuperinterface_boundarylist_alloc_init(SSUP_IFCE_BOUNDARY_LIST *bdrylist, int *size) {
    int i,j;
    for (j=0;j<2;j++){
        bdrylist->size[j] = size[j];
        bdrylist->couplededges[j] = (int *) tl_alloc(sizeof(int), bdrylist->size[j]);
        bdrylist->coupled_normal_flux[j] = (double *) tl_alloc(sizeof(double), bdrylist->size[j]);
        for (i=0;i<size[j];i++){
            bdrylist->couplededges[j][i] = UNSET_INT;
            bdrylist->coupled_normal_flux[j][i] = 0.0;
            //if (DEBUG) printf("\ncouplededges: %i",bdrylist->couplededges[j][i]);
        }
    }
}

/***********************************************************/
/***********************************************************/

void ssuperinterface_free(SSUPER_INTERFACE *si) {
    int i,j;
    if (si != NULL){
        for (i=0; i<si->NumEdges; i++) {
            for (j=0;j<2;j++){
                if (si->bdrylist[i].size[j] > 0) {
                    //printf("\n\nNumEdges = %i, i=%i, size = %i\n", si->NumEdges, i, si->bdrylist[i].size[j]);
                    //printf("\n\n j = %i, cnf = %20.5f \n", j, si->bdrylist[i].coupled_normal_flux[j][0]);
                    si->bdrylist[i].couplededges[j] = (int *) tl_free(sizeof(int), si->bdrylist[i].size[j], si->bdrylist[i].couplededges[j]);
                    si->bdrylist[i].coupled_normal_flux[j] = (double *) tl_free(sizeof(double), si->bdrylist[i].size[j], si->bdrylist[i].coupled_normal_flux[j]);
                }
            }
        }
        si->bdrylist = (SSUP_IFCE_BOUNDARY_LIST *) tl_free(sizeof(SSUP_IFCE_BOUNDARY_LIST), si->NumEdges, si->bdrylist);
    }
}

/***********************************************************/
/***********************************************************/

// void ssuperinterface_create_data(SSUPER_INTERFACE *si, SSUPER_MODEL *sm) {
//     int i, j, k, l;
//     ID_LIST_ITEM *ptr;
//     SMODEL *mod;
//     SGRID *grid;
//     SSUP_IFCE_BOUNDARY_LIST *ndlist;
// 
//     for (k=0; k<2; k++){
//         mod = submodel[k];
//         grid = mod->grid;
//         if (mod->flag.SW2_FLOW){
// #ifdef _MESSG
//             for (j=0; j<si->NumEdges; j++){
//                 ndlist = &(si->bdrylist[j]);
//                 //ndlist->size[k] = 0;
//                 int local_id = 0;
//                 for (local_id=0; local_id < grid->nnodes; local_id++){
//                     if (grid->node[local_id].gid == ndlist->surfnode[k]){
//                         ndlist->size[k] = 1;
//                         break;
//                     }
//                 }
//                 if(ndlist->size[k]==0){
//                     continue;
//                 }
//                 // printf("\nPE[%i] local_id = %i, gid = %i", sm->submodel[0].grid->smpi->myid, local_id, grid->node[local_id].gid);
//                 ndlist->couplededges[k] = (int *) tl_alloc(sizeof(int), ndlist->size[k]);
//                 ndlist->couplededges[k][0] = local_id;//ndlist->surfnode[k];
//             }
// #else
//             for (j=0; j<si->NumEdges; j++){
//                 ndlist = &(si->bdrylist[j]);
//                 ndlist->size[k] = 1;
//                 ndlist->couplededges[k] = (int *) tl_alloc(sizeof(int), ndlist->size[k]);
//                 ndlist->couplededges[k][0] = ndlist->surfnode[k];
//             }
// #endif
//         }
//         else if (mod->flag.SW3_FLOW){
//             for (j=0; j<si->NumEdges; j++){
//                 ndlist = &(si->bdrylist[j]);
// #ifdef _MESSG
//                 //ndlist->size[k] = 0;
//                 int gid = 0;
//                 for (gid=0; gid < grid->nnodes; gid++){
//                     if (grid->node[gid].gid == ndlist->surfnode[k]){
//                         break;
//                     }
//                 }
// 
//                 int col_id = 0;
//                 for (col_id=0; col_id < grid->nnodes_sur; col_id++){
//                     if (grid->vertical_list[col_id]->id == gid){
// #ifdef _DEBUG
//                         if (DEBUG){
//                             printf("\ngid = %i, grid->vertical_list[%i] = %i", gid, col_id, grid->vertical_list[col_id]->id);
//                         }
// #endif
//                         break;
//                     }
//                 }
//                 if (col_id >=grid->nnodes_sur) {
// #ifdef _DEBUG
//                     if (DEBUG){
//                         printf("[PE%i][Caution] Surface node with ID %i in submodel %i is not present on this processor\n", grid->smpi->myid, ndlist->surfnode[k]+1, si->model_id[k]+1);
//                     }
// #endif
//                     continue;
//                 }
// #else
//                 int col_id = 0;
//                 for (col_id=0; col_id < grid->nnodes_sur; col_id++){
//                     if (grid->vertical_list[col_id]->id == ndlist->surfnode[k]){
//                         break;
//                     }
//                 }
//                 if (col_id >=grid->nnodes_sur) {
//                     printf("Surface node with ID %i could not be found in submodel %i\n", ndlist->surfnode[k]+1,si->model_id[k]+1);
//                     tl_error("Please check interface data");
//                 }
// #endif
//                 /* Get the first node in the vertical list (which is the surface node) */
//                 ptr = grid->vertical_list[col_id];
//                 //printf("ptr-id = %i\n",ptr->id);
//                 int inode= UNSET_INT;
//                 int count = 0;
//                 while (ptr->next!=NULL){
//                     count++;
//                     ptr= ptr->next;
//                 }
//                 ndlist->size[k] = count; 
//                 ndlist->couplededges[k] = (int *) tl_alloc(sizeof(int), ndlist->size[k]);
//                 ptr = grid->vertical_list[col_id];
//                 count = 0;
//                 while (ptr->next!=NULL){
//                     inode = ptr->id;
//                     ndlist->couplededges[k][count] = inode;
//                     count++;
//                     ptr = ptr->next; /* skip first node which is surface */
//                 }
//                 //printf("\n");
//             }
//         }
// //        exit(-1);
//     }
// }

/***********************************************************/
/***********************************************************/

//double projected_node_distance(SGRID *grid1, int n1, SGRID *grid2, int n2){
//    return ( fabs(grid2->node[n2].x - grid1->node[n1].x) + fabs(grid2->node[n2].y - grid1->node[n1].y) );
//}


/***********************************************************/
/***********************************************************/

void generate_2d2d_superinterface(SSUPER_INTERFACE *si, SSUPER_MODEL *sm, int mastersupermodel, int slavesupermodel){
    int i, j, k, l;
    int istr1, istr2;
    int ie1, ie2;
    SMODEL *mod1, *mod2;
    SSUPER_MODEL *mastersm = &(sm[mastersupermodel]);
    SSUPER_MODEL *slavesm  = &(sm[slavesupermodel]);

    int mastermodel, slavemodel;
    if (si->sm_id[0] == mastersupermodel){
        mastermodel = 0;
        slavemodel = 1;
    }
    else /* if (si->sm_id[1] == mastersupermodel) */ {
        mastermodel = 1;
        slavemodel = 0;
    }
    mod1 = &(mastersm->submodel[si->model_id[mastermodel]]);
    mod2 = &(slavesm->submodel[si->model_id[slavemodel]]);
#ifdef _DEBUG
    if (DEBUG) {
        printf("\nmastersupermodel : %i, slavesupermodel : %i", mastersupermodel, slavesupermodel);
        printf("\nmastermodel : %i, slavemodel : %i", mastermodel, slavemodel);
        printf("\nsm[%i].submodel[%i] and sm[%i].submodel[%i]\n",
                     si->sm_id[mastersupermodel], si->model_id[mastermodel],
                     si->sm_id[slavesupermodel], si->model_id[slavemodel]);
    }
#endif

    SGRID *grid1 = mod1->grid;
    SGRID *grid2 = mod2->grid;

    SSUP_IFCE_BOUNDARY_LIST *edgelist;
    int size[2]={1, 1}; // Easy since this is a 2D-2D interface

    int edgecount = 0;
    for (ie1=0; ie1<grid1->nelems1d; ie1++){
        istr1 = grid1->elem1d[ie1].string;
        if (istr1 == si->bdrystrings[mastermodel]) {
            int n01 = grid1->elem1d[ie1].nodes[0];
            int n02 = grid1->elem1d[ie1].nodes[1];
            ssuperinterface_boundarylist_alloc_init(&(si->bdrylist[edgecount]), size);
            edgelist = &(si->bdrylist[edgecount]);
            edgelist->couplededges[mastermodel][0] = ie1;
            edgelist->surfnode[mastermodel][0] = n01;
            edgelist->surfnode[mastermodel][1] = n02;
#ifdef _DEBUG
            if (DEBUG){
                printf("\nsupermodel[%i]: elem1d[%i].string[%i]",mastersupermodel,ie1, istr1);
                printf("\nbdrystrings: %i and %i", si->bdrystrings[mastermodel], si->bdrystrings[slavemodel]);
                printf("\nsizes: %i and %i", edgelist->size[mastermodel], edgelist->size[slavemodel]);
                printf("\nsurfnodes[0]: %i and %i", edgelist->surfnode[mastermodel][0], edgelist->surfnode[mastermodel][1]);
                printf("\ncouplededges: %i and %i\n", edgelist->couplededges[mastermodel][0], edgelist->couplededges[slavemodel][0]);
            }
#endif
            /* Finding coupled edge below: */
            for (ie2=0; ie2<grid2->nelems1d; ie2++){
                istr2 = grid2->elem1d[ie2].string;
                if (istr2 == si->bdrystrings[slavemodel]) {
                    int n11 = grid2->elem1d[ie2].nodes[0];
                    int n12 = grid2->elem1d[ie2].nodes[1];
                    if   ((projected_node_distance(grid1,n01,grid2,n11)<NOT_QUITE_SMALL)
                      &&  (projected_node_distance(grid1,n02,grid2,n12)<NOT_QUITE_SMALL)){
                        mod2->str_values[istr2].ol_flow.bc_flag = BCT_HYBRID_EXTERNAL;
                        edgelist->couplededges[slavemodel][0] = ie2;
                        edgelist->surfnode[slavemodel][0] = n11;
                        edgelist->surfnode[slavemodel][1] = n12;
                        break; /* breaks the ie2 for loop */
                    }
                    else if   ((projected_node_distance(grid1,n01,grid2,n12)<NOT_QUITE_SMALL)
                           &&  (projected_node_distance(grid1,n02,grid2,n11)<NOT_QUITE_SMALL)){
                        mod2->str_values[istr2].ol_flow.bc_flag = BCT_HYBRID_EXTERNAL;
                        edgelist->couplededges[slavemodel][0] = ie2;
                        edgelist->surfnode[slavemodel][0] = n12;
                        edgelist->surfnode[slavemodel][1] = n11;
                        break; /* breaks the ie2 for loop */
                    }
                    if (DEBUG) printf("\nIE1: %3i, n1 %3i, n2 %3i\nIE2: %3i, n1 %3i, n2 %3i\n", ie1, n01, n02, ie2, n11, n12);
                }
            } /* ie2 loop */
#ifdef _DEBUG
            if (DEBUG) {
                printf("\nbdrystrings: %i and %i", si->bdrystrings[mastermodel], si->bdrystrings[slavemodel]);
                printf("\nsizes: %i and %i", edgelist->size[mastermodel], edgelist->size[slavemodel]);
                printf("\nsurfnodes[0]: %i and %i", edgelist->surfnode[slavemodel][0], edgelist->surfnode[slavemodel][1]);
                printf("\ncouplededges: %i and %i\n", edgelist->couplededges[mastermodel][0], edgelist->couplededges[slavemodel][0]);
            }
#endif
            if (edgelist->couplededges[slavemodel][0] == UNSET_INT){
                tl_error("Problem in superinterface! Pairing edge for coupling not found!");
            }
            edgecount++;
        }
    } /* ie1 loop */
    assert (edgecount == si->NumEdges);
}
