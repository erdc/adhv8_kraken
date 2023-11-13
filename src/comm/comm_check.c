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
  double temp_surf;

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
  int jj;
  for (ii = 0, jj = 0; ii < grid->my_nnodes; ii++)
    {
      vx[ii] = grid->node[ii].x;
      vy[ii] = grid->node[ii].y;
      vz[ii] = grid->node[ii].z;
      v[ii] = grid->node[ii].gid;
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
    }

  /* communicate the coordinates */
  comm_update_double(vx, 1, grid->smpi);
  comm_update_double(vy, 1, grid->smpi);
  comm_update_double(vz, 1, grid->smpi);
  comm_update_int(v, 1, grid->smpi);
  comm_update_double_surf(vd_surf, 1, grid->smpi);
  comm_update_int_surf(vi_surf, 1, grid->smpi);
  comm_update_VECT2D(v2, grid->smpi);messg_barrier(grid->smpi->ADH_COMM);
  comm_update_double(v2, 2, grid->smpi);
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
      if (diff >SMALL && first_error==0){first_error++;}
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
      temp = abs(v[ii] - grid->node[ii].gid);
      diff += temp;
      if (diff >SMALL&&first_error==0){first_error++;}
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
  diff = DZERO;first_error=0;
  for (ii = 0; ii < grid->nnodes_sur; ii++)
    {
      
      if(grid->type == COLUMNAR) {jj = grid->nodeID_2d_to_3d_sur[ii];}
      else{jj=ii;}
      temp = fabs(vd_surf[ii] - (grid->node[jj].x + grid->node[jj].y + grid->node[jj].z));
      if(isnan(temp)){printf("MYID %d temp anan for node %d surf %d\n",grid->smpi->myid, jj, ii);}
      else{diff += temp;}
      if ((diff > SMALL) && (first_error==0)){printf("myid %d ii %d jj %d vd_surf[ii] %14.10f coor_diff %14.10f temp %14.10e diff %10.6e \n",grid->smpi->myid, ii, jj, vd_surf[ii], grid->node[jj].x + grid->node[jj].y + grid->node[jj].z, (vd_surf[ii] - (grid->node[jj].x + grid->node[jj].y + grid->node[jj].z)), diff);first_error++;}
    }
  
  if (diff > SMALL)
    {
      printf("ghost border MISMATCH found in SURFACE DOUBLE.\n");ierror++;
    }
  else
    {
      printf("===============myid=%d, comm_check SURFACE DOUBLE was SUCCESSFUL=============\n", grid->smpi->myid);
    }
  messg_barrier(grid->smpi->ADH_COMM);
  diff = DZERO;first_error=0;
  for (ii = 0; ii < grid->nnodes_sur; ii++)
    {
      if(grid->type == COLUMNAR) {jj = grid->nodeID_2d_to_3d_sur[ii];}
      else{jj = ii;}
      temp = abs(vi_surf[ii] - grid->node[jj].global_surf_id);
      diff += temp;
      if (diff >SMALL&& first_error==0){printf("myid %d ii %d vi_surf[ii] %d surf_id %d \n",grid->smpi->myid, ii, vi_surf[ii], grid->node[jj].global_surf_id);first_error++;}
    }
  if (diff > SMALL)
    {
      printf("ghost border MISMATCH found in SURFACE INT.\n");ierror++;
    }
  else
    {
      printf("===============myid=%d, comm_check was SURFACE INT SUCCESSFUL=============\n", grid->smpi->myid);
    }
  messg_barrier(grid->smpi->ADH_COMM);first_error=0;
  diff = DZERO;first_error=0;
  for (ii = 0; ii < grid->nnodes; ii++)
    {
      temp = fabs(v2[ii].x - grid->node[ii].x);
      diff += temp;
      temp += fabs(v2[ii].y - grid->node[ii].y);
      diff += temp;
      if (diff >SMALL&& first_error==0){printf("myid %d ii %d v2[ii].x %f grid->node[ii].x %f v2[ii].y %f grid->node[ii].y %f diff %7.5e SMALL %7.5e\n",grid->smpi->myid, ii, v2[ii].x, grid->node[ii].x, v2[ii].y, grid->node[ii].y, diff, SMALL);first_error++;}
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
  diff = DZERO;first_error=0;
  for (ii = 0; ii < grid->nnodes; ii++)
    {
      temp = fabs(v3[ii].x - grid->node[ii].x);
      diff += temp;
      temp += fabs(v3[ii].y - grid->node[ii].y);
      diff += temp;
      temp += fabs(v3[ii].z - grid->node[ii].z);
      diff += temp;
      if (diff >SMALL && first_error==0){printf("myid %d ii %d v3[ii].x %f grid->node[ii].x %f v3[ii].y %f grid->node[ii].y %f v3[ii].z %f grid->node[ii].z %f \n",grid->smpi->myid, ii, v3[ii].x, grid->node[ii].x, v3[ii].y, grid->node[ii].y, v3[ii].z, grid->node[ii].z );first_error++;}
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
  vi_surf = (int *) tl_free(sizeof(int), grid->nnodes_sur, vi_surf);
  vd_surf = (double *) tl_free(sizeof(double), grid->nnodes_sur, vd_surf);
  if (grid->nelems3d > 0) v3 = (SVECT *) tl_free(sizeof(SVECT), grid->nnodes, v3);
#endif
  return;
}
