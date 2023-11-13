#include "global_header.h"

/* Updates the SCON struct after partitioning */

void comm_update_con(SCON *con, SGRID *grid, int ntransport){

  int i,j,itrn; //counters
  int num_doubles_con = 10; //from scon.h not moving properties
  double *packed_con;
 
  packed_con = (double *)tl_alloc(sizeof(double), grid->nnodes * num_doubles_con * ntransport);

  /* pack the array */
  for(j = 0, itrn = 0; itrn < ntransport; itrn++){
    for(i = 0; i < grid->nnodes;i++){
      packed_con[j++] = con[itrn].concentration[i];
      packed_con[j++] = con[itrn].old_concentration[i];
      packed_con[j++] = con[itrn].older_concentration[i];
      packed_con[j++] = con[itrn].sink[i];
      packed_con[j++] = con[itrn].source[i];
      packed_con[j++] = con[itrn].nodal_decay_coefficient[i];
      packed_con[j++] = con[itrn].error[i];
      packed_con[j++] = con[itrn].mfcf[i];
      packed_con[j++] = con[itrn].vcf[i].x;
      packed_con[j++] = con[itrn].vcf[i].y;
    }
  }

  /*update the packed array */
  comm_update_double(packed_con, num_doubles_con*ntransport, grid->smpi);

  /*update the con struct (only ghost nodes to save time and be safe */
  for(itrn = 0; itrn < ntransport; itrn++){
    for(j = ((grid->my_nnodes * num_doubles_con) + (grid->nnodes * num_doubles_con * itrn)), i = grid->my_nnodes; i < grid->nnodes;i++){
      con[itrn].concentration[i] = packed_con[j++];
      con[itrn].old_concentration[i] = packed_con[j++];
      con[itrn].older_concentration[i] = packed_con[j++];
      con[itrn].sink[i] = packed_con[j++];
      con[itrn].source[i] = packed_con[j++];
      con[itrn].nodal_decay_coefficient[i] = packed_con[j++];
      con[itrn].error[i] = packed_con[j++];
      con[itrn].mfcf[i] = packed_con[j++];
      con[itrn].vcf[i].x = packed_con[j++];
      con[itrn].vcf[i].y = packed_con[j++];
    }
  } 

  packed_con = (double *)tl_free(sizeof(double), grid->nnodes * num_doubles_con * ntransport, packed_con);
}
