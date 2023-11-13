/* this routine packs the double information for the given node in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void print_sw3_to_file(char *description, SSW *sw,
                   SGRID *grid,
                   int ntransport,
                   SCON *con )
{
  int itrn;                     /* loop counter over the transported quantities */
  int i,ind, k;                     /* counters */
  int nnode = grid->my_nnodes;
  FILE *fp;
  char filename[MAXLINE];
  
  fp = io_fopen(build_filename2(filename, MAXLINE, description, "_sw3_", grid->smpi->myid, ".txt", UNSET_INT), "w", TRUE);
  fprintf(fp,"%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n", "node_x","node_y","node_z","vel_x","vel_y","vel_z","old_vel_x","old_vel_y","old_vel_z","older_vel_x","older_vel_y","older_vel_z", "prs","prs_plus","prs_minus", "disp", "old_disp", "older_disp","grid_speed","old_grid");
  for(i=0;i<nnode;i++){
    
    fprintf(fp,"%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n",grid->node[i].x, grid->node[i].y,grid->node[i].z,sw->d3->vel[i].x,sw->d3->vel[i].y,sw->d3->vel[i].z,sw->d3->old_vel[i].x,sw->d3->old_vel[i].y,sw->d3->old_vel[i].z,sw->d3->older_vel[i].x,sw->d3->older_vel[i].y,sw->d3->older_vel[i].z, sw->d3->prs[i],sw->d3->prs_plus[i],sw->d3->prs_minus[i], sw->d3->displacement[i], sw->d3->old_displacement[i], sw->d3->older_displacement[i],sw->d3->grid_speed[i],sw->d3->old_grid_speed[i]);
  }
  fclose(fp);

fp = io_fopen(build_filename2(filename, MAXLINE, description, "_con_", grid->smpi->myid, ".txt", UNSET_INT), "w", TRUE);
    /* cjt */
    if (con != NULL) {
      
      for (itrn = 0; itrn < ntransport; itrn++) {
        fprintf(fp,"%12s%12s","nod_source", "decay_coef");
        fprintf(fp,"%11s1%11s1%11s1","con","old_con","older_con");
      }
      fprintf(fp,"\n");
      
         for (itrn = 0; itrn < ntransport; itrn++) {
           for(i=0;i<nnode;i++){
            fprintf(fp,"%12.4e%12.4e",con[itrn].source[i], con[itrn].nodal_decay_coefficient[i]);
            fprintf(fp,"%12.4e%12.4e%12.4e",con[itrn].concentration[i]/con[itrn].property[0],con[itrn].old_concentration[i]/con[itrn].property[0],con[itrn].older_concentration[i]/con[itrn].property[0]);
            fprintf(fp,"\n");
          }
       }
    }
    fclose(fp);

    if (sw->d3->waves != NULL ) printf("NEED WAVE PRINT\n");
      

    if (sw->d3->winds != NULL) {
      printf("NEED WIND PRINT\n");
    }
 }
#endif
