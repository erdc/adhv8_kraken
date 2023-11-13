/* checks comm_update by comparing the ghost values to the actual values in the 
   overlap region */

#include "global_header.h"

void comm_check(SGRID *grid)
{
#ifdef _MESSG
  int ierror = 0;               /* error indicator */
  char msg[MAXLINE];
  int ii = 0;                   /* loop counter */
  double *vx = NULL;            /* the coordinates to be communicated */
  double *vy = NULL;            /* the coordinates to be communicated */
  double *vz = NULL;            /* the coordinates to be communicated */
  double *vd_surf = NULL;       /* sum of coordinates for surface comm test*/
  SVECT2D *v2 = NULL;
  SVECT *v3 = NULL;
  int *v = NULL;
  int *vi_surf = NULL;
  double temp = DZERO;          /* temporary variable */
  double diff = DZERO;          /* difference between actual grid and communicated grid 
                                   (should be zero) */

  /* allocate space */
  vx = (double *) tl_alloc(sizeof(double), grid->nnodes);
  vy = (double *) tl_alloc(sizeof(double), grid->nnodes);
  vz = (double *) tl_alloc(sizeof(double), grid->nnodes);
  vd_surf = (double *) tl_alloc(sizeof(double), grid->nnodes_sur);
  v = (int *) tl_alloc(sizeof(int), grid->nnodes);
  vi_surf = (int *) tl_alloc(sizeof(int), grid->nnodes_sur);
  v2 = (SVECT2D *) tl_alloc(sizeof(SVECT2D), grid->nnodes);
  if (grid->nelems3d > 0) v3 = (SVECT *) tl_alloc(sizeof(SVECT), grid->nnodes);

  /* load the coordinate data for communication */
  int i,jj;
  for (ii = 0, jj = 0; ii < grid->my_nnodes; ii++)
    {
      vx[ii] = grid->node[ii].x;
      vy[ii] = grid->node[ii].y;
      vz[ii] = grid->node[ii].z;
      v[ii] = 1;
      if (grid->node[ii].global_surf_id != UNSET_INT) {
        vi_surf[jj]=grid->node[ii].global_surf_id;
        vd_surf[jj]=grid->node[ii].x+grid->node[ii].y+grid->node[ii].z;
        jj++;
      }
      v2[ii].x = grid->node[ii].x;
      v2[ii].y = grid->node[ii].y;
      if (grid->nelems3d > 0) {
        v3[ii].x = grid->node[ii].x;
        v3[ii].y = grid->node[ii].y;
        v3[ii].z = grid->node[ii].z;
      } 
    }
  for (ii = grid->my_nnodes; ii < grid->nnodes; ii++)
    {
      vx[ii] = UNSET_FLT;
      vy[ii] = UNSET_FLT;
      vz[ii] = UNSET_FLT;
      v[ii] = UNSET_INT;
      if (grid->node[ii].global_surf_id != UNSET_INT) {
        vi_surf[jj] = UNSET_INT;
        vd_surf[jj] = UNSET_FLT;
        jj++;
      }
      v2[ii].x = UNSET_FLT;
      v2[ii].y = UNSET_FLT;
      if (grid->nelems3d > 0) {
        v3[ii].x = UNSET_FLT;
        v3[ii].y = UNSET_FLT;
        v3[ii].z = UNSET_FLT;
      }
      //printf("myid %d ii %d vx %f grid->node[ii].x %f vy %f grid->node[ii].y %f \n", grid->smpi->myid, ii, v2[ii].x, grid->node[ii].x,v2[ii].y,grid->node[ii].y);
    }

  /* communicate the coordinates */
  comm_update_double(vx, 1, grid->smpi);
  comm_update_double(vy, 1, grid->smpi);
  comm_update_double(vz, 1, grid->smpi);
  comm_update_int(v, 1, grid->smpi);
  //comm_update_double_surf(vd_surf, 1, grid->smpi);
  //comm_update_int_surf(vi_surf, 1, grid->smpi);
  comm_update_VECT2D(v2, grid->smpi);
  //comm_update_double(v2, 2, grid->smpi);
  if (grid->nelems3d > 0) comm_update_VECT(v3, grid->smpi);
  messg_barrier(grid->smpi->ADH_COMM);
  /* calculate the difference between my actual grid (node[ii].x) and 
     the grid my neighbor thinks I have (vx), etc. */
  int first_error =0;
  for (ii = 0; ii < grid->nnodes; ii++)
    {
      
      temp = fabs(vx[ii] - grid->node[ii].x);
      diff += temp;
      temp += fabs(vy[ii] - grid->node[ii].y);
      diff += temp;
      temp += fabs(vz[ii] - grid->node[ii].z);
      diff += temp;
      if (diff >SMALL && first_error==0){printf("myid %d ii %d vx[ii] %f grid->node[ii].x %f vy[ii] %f grid->node[ii].y %f vz[ii] %f grid->node[ii].z %f diff %7.5e DZERO %7.5e SMALL %7.5e\n",grid->smpi->myid, ii, vx[ii], grid->node[ii].x, vy[ii], grid->node[ii].y, vz[ii], grid->node[ii].z, diff, DZERO,SMALL );first_error++;}
    }
  if (diff > SMALL)
    {
      printf("ghost border MISMATCH found in DOUBLE.\n");ierror++;
    }
  else
    {
      printf("===============myid=%d, comm_check DOUBLE was SUCCESSFUL=============\n", grid->smpi->myid);
    }
  messg_barrier(grid->smpi->ADH_COMM);
  diff = DZERO;first_error=0;
  for (ii = 0; ii < grid->nnodes; ii++)
    {
      temp = fabs(v[ii] - 1);
      diff += temp;
      if (diff >SMALL&&first_error==0){printf("myid %d ii %d v[ii] %d \n",grid->smpi->myid, ii, v[ii]);first_error++;}
    }
  if (diff > SMALL)
    {
      printf("ghost border MISMATCH found in INT.\n");ierror++;
    }
  else
    {
      printf("===============myid=%d, comm_check was INT SUCCESSFUL=============\n", grid->smpi->myid);
    }
  messg_barrier(grid->smpi->ADH_COMM);
  diff = DZERO;
  for (ii = 0; ii < grid->nnodes_sur; ii++)
    {
      
      if(grid->ndim==3) {jj = grid->nodeID_2d_to_3d_sur[ii];}
      else{jj=ii;}
      temp = fabs(vd_surf[ii] - (grid->node[jj].x + grid->node[jj].y + grid->node[jj].z));
      diff += temp;
    //  if (diff >SMALL)printf("myid %d ii %d jj %d vd_surf[ii] %f coor_diff %f \n",grid->smpi->myid, ii, jj, vd_surf[ii], grid->node[jj].x - grid->node[jj].y - grid->node[jj].z);
    }
  if (diff > SMALL)
    {
     // printf("ghost border MISMATCH found in SURFACE DOUBLE.\n");ierror++;
    }
  else
    {
      printf("===============myid=%d, comm_check SURFACE DOUBLE was SUCCESSFUL=============\n", grid->smpi->myid);
    }
  messg_barrier(grid->smpi->ADH_COMM);
  diff = DZERO;
  for (ii = 0; ii < grid->nnodes_sur; ii++)
    {
      if(grid->ndim==3) {jj = grid->nodeID_2d_to_3d_sur[ii];}
      else{jj = ii;}
      temp = fabs(vi_surf[ii] - grid->node[jj].global_surf_id);
      diff += temp;
     // if (diff >SMALL)printf("myid %d ii %d vi_surf[ii] %d surf_id %d \n",grid->smpi->myid, ii, vi_surf[ii], grid->node[jj].global_surf_id);
    }
  if (diff > SMALL)
    {
     // printf("ghost border MISMATCH found in INT.\n");ierror++;
    }
  else
    {
      printf("===============myid=%d, comm_check was SURFACE INT SUCCESSFUL=============\n", grid->smpi->myid);
    }
  messg_barrier(grid->smpi->ADH_COMM);
  for (ii = 0; ii < grid->nnodes; ii++)
    {
      temp = fabs(v2[ii].x - grid->node[ii].x);
      diff += temp;
      temp += fabs(v2[ii].y - grid->node[ii].y);
      diff += temp;
      if (diff >SMALL)printf("myid %d ii %d v2[ii].x %f grid->node[ii].x %f v2[ii].y %f grid->node[ii].y %f\n",grid->smpi->myid, ii, v2[ii].x, grid->node[ii].x, v2[ii].y, grid->node[ii].y);diff = SMALL;
    }
  if (diff > SMALL)
    {
      printf("ghost border MISMATCH found in VECT2D.\n");ierror++;
    }
  else
    {
      printf("===============myid=%d, comm_check VECT2D was SUCCESSFUL=============\n", grid->smpi->myid);
    }
  messg_barrier(grid->smpi->ADH_COMM);
  if (grid->nelems3d > 0) {
  diff = DZERO;
  for (ii = 0; ii < grid->nnodes; ii++)
    {
      temp = fabs(v3[ii].x - grid->node[ii].x);
      diff += temp;
      temp += fabs(v3[ii].y - grid->node[ii].y);
      diff += temp;
      temp += fabs(v3[ii].z - grid->node[ii].z);
      diff += temp;
      if (diff >SMALL)printf("myid %d ii %d v3[ii].x %f grid->node[ii].x %f v3[ii].y %f grid->node[ii].y %f v3[ii].z %f grid->node[ii].z %f \n",grid->smpi->myid, ii, v3[ii].x, grid->node[ii].x, v3[ii].y, grid->node[ii].y, v3[ii].z, grid->node[ii].z );
    }
  if (diff > SMALL)
    {
      printf("ghost border MISMATCH found in VECT.\n");ierror++;
    }
  else
    {
      printf("===============myid=%d, comm_check VECT was SUCCESSFUL=============\n", grid->smpi->myid);
    }
  }
  messg_barrier(grid->smpi->ADH_COMM);
  if (ierror >0) {
    sprintf(msg, "there were %d communication errors!", ierror);
    tl_error(msg);
  }
  
  vx = (double *) tl_free(sizeof(double), grid->nnodes, vx);
  vy = (double *) tl_free(sizeof(double), grid->nnodes, vy);
  vz = (double *) tl_free(sizeof(double), grid->nnodes, vz);
  v = (int *) tl_free(sizeof(int), grid->nnodes, v);
  v2 = (SVECT2D *) tl_free(sizeof(SVECT2D), grid->nnodes, v2);
  if (grid->nelems3d > 0) v3 = (SVECT *) tl_free(sizeof(SVECT), grid->nnodes, v3);
#endif
  return;
}
