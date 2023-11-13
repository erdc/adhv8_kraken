/* returns the 2D gradient of the scalar field for at a node */

#include "global_header.h"

//***************************************************************************************************//
//***************************************************************************************************//
//***************************************************************************************************//

void tl_nodal_grad_avg_2d(SGRID *grid, STR_VALUE *string, double *sfield, SVECT2D *nodal_grad_avg) {

    double elem_sfield[NDONTRI];

    SVECT2D nodal_grad;
    svect2d_init_array(nodal_grad_avg, grid->nnodes);

    int inode = 0, gnode = UNSET_INT, ie = 0;
    int n0=0, n1=0, n2=0;
    double coef = 0;
    for(ie = 0; ie < grid->nelems2d; ie++) {
        if (string[grid->elem2d[ie].string].phys_flag != OFF) {
	
            /* elemental scalar field */ 
            ELEM2D_GET_LOCAL(sfield, elem_sfield, grid->elem2d[ie].nodes);
            
            /* find the gradient of the solution */
            GRAD2D_FUNC(grid->elem2d[ie].grad_shp, elem_sfield, nodal_grad);
            
            /* weight */
            coef = 1. / 3. * grid->elem2d[ie].djac;

            for (inode=0; inode<NDONTRI; inode++) {
                gnode = grid->elem2d[ie].nodes[inode]; // global node id
                grid->node[gnode].area += coef;
                nodal_grad_avg[gnode].x += nodal_grad.x * coef;
                nodal_grad_avg[gnode].y += nodal_grad.y * coef;
            }
        }
    }
  
    /* divides the nodal bedslope by the areas to get a local value for bedslope */ 
    for(inode=0; inode < grid->my_nnodes; inode++)  {
      if(grid->node[inode].area > SMALL) {
          nodal_grad_avg[inode].x /= grid->node[inode].area;
          nodal_grad_avg[inode].y /= grid->node[inode].area;
      } else if(grid->node[inode].area <= SMALL) {
          nodal_grad_avg[inode].x = 0.0;
          nodal_grad_avg[inode].y = 0.0;
      }
    }

}

//***************************************************************************************************//
//***************************************************************************************************//
//***************************************************************************************************//

void tl_nodal_grad_avg_3d(SGRID *grid, STR_VALUE *string, double *sfield, SVECT2D *nodal_grad_avg) {
    
    double elem_sfield[NDONTRI];
    
    SVECT2D nodal_grad;
    svect2d_init_array(nodal_grad_avg, grid->nnodes_bed);
    
    double *total_area = (double *) tl_alloc(sizeof(double), grid->nnodes_bed);
    sarray_init_dbl(total_area, grid->nnodes_bed);
    
    int inode = 0, gnode = UNSET_INT, ie = 0, id2d = 0;
    int n0=0, n1=0, n2=0;
    double coef = 0;
    for(ie = 0; ie < grid->nelems2d_bed; ie++) {
        
        id2d = grid->elem2d_bed[ie]; // the 2d element id (all 2d elements, sidewall, surface, etc.
        
        if (string[grid->elem2d[id2d].string].phys_flag != OFF) {
            
            /* elemental scalar field */
            for (inode=0; inode<NDONTRI; inode++) {
                elem_sfield[inode] = sfield[ grid->nodeID_3d_to_2d_bed[ grid->elem2d[id2d].nodes[inode] ] ];
            }
            
            /* find the gradient of the solution */
            GRAD2D_FUNC(grid->elem2d[id2d].grad_shp, elem_sfield, nodal_grad);
            
            /* weight */
            coef = 1. / 3. * grid->elem2d[id2d].djac;
            
            for (inode=0; inode<NDONTRI; inode++) {
                gnode = grid->nodeID_3d_to_2d_bed[ grid->elem2d[id2d].nodes[inode] ]; // global 2d node id
                total_area[gnode] += coef;
                nodal_grad_avg[gnode].x += nodal_grad.x * coef;
                nodal_grad_avg[gnode].y += nodal_grad.y * coef;
            }
        }
    }
    
    /* divides the nodal bedslope by the areas to get a local value for bedslope */
    for(inode=0; inode < grid->nnodes_bed; inode++)  {
        if(total_area[inode] > SMALL) {
            nodal_grad_avg[inode].x /= total_area[inode];
            nodal_grad_avg[inode].y /= total_area[inode];
        } else if(total_area[inode] <= SMALL) {
            nodal_grad_avg[inode].x = 0.0;
            nodal_grad_avg[inode].y = 0.0;
        }
    }
    
    total_area = (double *) tl_free(sizeof(double), grid->nnodes_bed, total_area);
    
}











