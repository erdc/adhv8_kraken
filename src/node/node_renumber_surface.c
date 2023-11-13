/* this routine renumbers the nodes */
/* NOTE:  THE GLOBAL NODE NUMBERS NEED TO BE UPDATED AFTER CALLING THIS ROUTINE !! */

#include "global_header.h"

void node_renumber_surface(SMODEL *mod) {
     
    int i, j, k;
    int num_alt_nodes;              /* the total number of nodes that have been changed */
    
    SGRID *grid = mod->grid; // alias
    int max_nnode = mod->grid->max_nnodes_sur;
    int nnode = mod->grid->nnodes_sur;
 
    // temporary storage arrays
    SVECT *vtmp = (SVECT *) tl_alloc(sizeof(SVECT), max_nnode);
    SVECT2D *v2dtmp = (SVECT2D *) tl_alloc(sizeof(SVECT2D), max_nnode);
    int *itmp = (int *) tl_alloc(sizeof(int), max_nnode);
    double *dtmp = (double *) tl_alloc(sizeof(double), max_nnode);
    int *order_tmp = (int *) tl_alloc(sizeof(int), max_nnode);
    int *new_node_number = (int *) tl_alloc(sizeof(int), max_nnode);

    /* Initialize the order array so we don't do it every sub call */
    j=0;
    for (i = 0; i < max_nnode; i++) {
      if(grid->node[grid->nodeID_2d_to_3d_sur[j]].global_surf_id == UNSET_INT) j++;
        order_tmp[i]=0;
        new_node_number[i]=UNSET_INT;
    }
    

    /* check to see if any nodes have changed numbers */
    k=0;
    num_alt_nodes = 0;
    for (i = 0; i < max_nnode; i++) {
      if(i<nnode){
        if(grid->old_global_surf[i]==grid->node[grid->nodeID_2d_to_3d_sur[i]].global_surf_id) {
          new_node_number[i]=i; 
        }
        else {
          for(j=0;j<nnode;j++){
            if(grid->old_global_surf[i]==grid->node[grid->nodeID_2d_to_3d_sur[j]].global_surf_id) {
              new_node_number[i]=j;
              num_alt_nodes++;
              if (new_node_number[i] > i) order_tmp[new_node_number[i]] = 1;
              break;
            }
          }
          if(new_node_number[i]==UNSET_INT){
            new_node_number[i]=nnode + k;
            if (new_node_number[i] > i) order_tmp[new_node_number[i]] = 1;
            k++;
          }
        }
      }
      else{
        if(grid->old_global_surf[i] != UNSET_INT){
          for(j=0;j<nnode;j++){
            if(grid->old_global_surf[i]==grid->node[grid->nodeID_2d_to_3d_sur[j]].global_surf_id) {
              new_node_number[i]=j;
              num_alt_nodes++;
              if (new_node_number[i] > i) order_tmp[new_node_number[i]] = 1;
              break;
            }
          }
        }
       if(new_node_number[i]==UNSET_INT){
         new_node_number[i]=nnode + k;
         k++;
       } 
      }
     }
    
    

    if (num_alt_nodes > 0) {
      
        /* renumber physics variables*/ 
        if (mod->flag.SW3_FLOW)  {ssw_3d_renumber_surface(mod->sw->d3, max_nnode, new_node_number, order_tmp, dtmp, v2dtmp);}
        //if (mod->flag.NS3_FLOW)  {sns_3d_renumber_surface(mod->ns->d3, max_nnode, new_node_number, order_tmp, dtmp, v2dtmp);}
#ifdef _SEDIMENT
        if (mod->flag.SW3_FLOW && mod->flag.SEDIMENT)  {ssediment_renumber_3d_bed(mod->sed, max_nnode, new_node_number, order_tmp, itmp, dtmp, v2dtmp);}
#endif
       
#ifdef _ADH_ICM
        if (mod->flag.ICM) {}
#endif

    }
    // free temporary storage arrays
    vtmp = (SVECT *) tl_free(sizeof(SVECT), max_nnode, vtmp);
    v2dtmp = (SVECT2D *) tl_free(sizeof(SVECT2D), max_nnode, v2dtmp);
    itmp = (int *) tl_free(sizeof(int), max_nnode, itmp);
    dtmp = (double *) tl_free(sizeof(double), max_nnode, dtmp);
    order_tmp = (int *) tl_free(sizeof(int), max_nnode, order_tmp);
    new_node_number = (int *) tl_free(sizeof(int), max_nnode, new_node_number);

}
