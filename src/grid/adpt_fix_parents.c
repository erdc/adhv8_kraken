/* These routine switch parent values between global IDs and local IDs. */
/* It is required when repartiioning the grid, otherwise, local parnet information is sent/recevied */


#include "global_header.h"

void set_parents_global(SGRID *g) {
  int i, parent0, parent1;

  for(i=0;i<g->nnodes;i++){
    if(g->node[i].original_id == UNSET_INT){//if(g->smpi->myid==2)printf("myid %d i %d of %d nodes parnets %d %d \n", g->smpi->myid, i, g->nnodes, g->node[i].parent[0], g->node[i].parent[1]);
      parent0 = g->node[i].parent[0];
      parent1 = g->node[i].parent[1];
      if((parent0 != UNSET_INT) && (parent0 < g->nnodes) && (parent0 >= 0)) g->node[i].parent[0] = g->node[parent0].gid;
      else g->node[i].parent[0] = UNSET_INT;
      if((parent1 != UNSET_INT) && (parent1 < g->nnodes) && (parent1 >= 0)) g->node[i].parent[1] = g->node[parent1].gid;
      else g->node[i].parent[1] = UNSET_INT;
     
    }
  }

  return;
}

void set_parents_local(SGRID *g) {
  int i, j, parent0, parent1, break1, break0;
  int ierr=0;
  for(i=0;i<g->nnodes;i++){
    if(g->node[i].original_id == UNSET_INT){
      parent0 = g->node[i].parent[0];
      parent1 = g->node[i].parent[1];
      break0=0;
      break1=0;
      if(parent0 == UNSET_INT) break0=1;
      if(parent1 == UNSET_INT) break1=1;
        
      for(j=0;j<g->nnodes;j++){
        if(break0 == 0){
          if(g->node[j].gid==parent0){
            g->node[i].parent[0] = j;
            break0=1;
          }
        }
        if(break1 == 0){
          if(g->node[j].gid==parent1){
            g->node[i].parent[1] = j;
            break1=1;
          }
        }
        if((break0>0)&&(break1>0))
          break;
      }
    }
  }
  
  return;

}

