#include "global_header.h"

/*
   \brief Update "global_id" Info in Global_Node entries
 */
void comm_update_GN(int flag,    /* 1 -> updating node_map, 0 for node_ladj or node_radj */
                    SGRID *g
  )
{
#ifdef _MPI
  int i, j, k;
  int *resident_pe;
  int *resident_id;
  int *old_resident_pe;
  int *old_resident_id;
  int *new_resident_pe;
  int *new_resident_id;

  resident_pe = (int *)tl_alloc(sizeof(int), g->nnodes);
  resident_id = (int *)tl_alloc(sizeof(int), g->nnodes);
  if (flag==1){
    old_resident_pe = (int *)tl_alloc(sizeof(int), g->nnodes);
    old_resident_id = (int *)tl_alloc(sizeof(int), g->nnodes);
    new_resident_pe = (int *)tl_alloc(sizeof(int), g->nnodes);
    new_resident_id = (int *)tl_alloc(sizeof(int), g->nnodes);
  }

  for (i=0;i<g->nnodes;i++){
    resident_pe[i]=g->node[i].resident_pe;
    resident_id[i]=g->node[i].resident_id;
  }
  comm_update_int(resident_pe, 1, g->smpi);
  comm_update_int(resident_id, 1, g->smpi);

  k=0;
  for (i=g->my_nnodes;i<g->nnodes;i++){
    if(flag==1){
      if((g->node[i].resident_id != UNSET_INT) && 
          ((g->node[i].resident_pe != resident_pe[i]) || (g->node[i].resident_id != resident_id[i]))){
        old_resident_id[k] = g->node[i].resident_id;
        old_resident_pe[k] = g->node[i].resident_pe;
        new_resident_id[k] = resident_id[i];
        new_resident_pe[k] = resident_pe[i];
        k++;
      }
    }
    g->node[i].resident_pe=resident_pe[i];
    g->node[i].resident_id=resident_id[i];
  }
  if((flag==1) && (k!=0)){
    for(i=0;i<k;i++){
      for(j=0;j<g->my_nnodes;j++){
        if((old_resident_id[i] == g->node[j].parent_res_id[0]) && (old_resident_pe[i] == g->node[j].parent_res_pe[0])){
          g->node[j].parent_res_id[0]=new_resident_id[i];
          g->node[j].parent_res_pe[0]=new_resident_pe[i];
        }else if ((old_resident_id[i] == g->node[j].parent_res_id[1]) && (old_resident_pe[i] == g->node[j].parent_res_pe[1])){
          g->node[j].parent_res_id[1]=new_resident_id[i];
          g->node[j].parent_res_pe[1]=new_resident_pe[i];
        }
      }
    }
  }
	if (flag==1){
    old_resident_pe = (int *)tl_free(sizeof(int), g->nnodes, old_resident_pe);
    old_resident_id = (int *)tl_free(sizeof(int), g->nnodes, old_resident_id);
    new_resident_pe = (int *)tl_free(sizeof(int), g->nnodes, new_resident_pe);
    new_resident_id = (int *)tl_free(sizeof(int), g->nnodes, new_resident_id);
  }
  resident_pe = (int *)tl_free(sizeof(int), g->nnodes, resident_pe);
  resident_id = (int *)tl_free(sizeof(int), g->nnodes, resident_id);
#endif
  return;
}
