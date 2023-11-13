#include "global_header.h"

void print_parents_to_file(SGRID *g, char *description)
{
  FILE *fp;
  int i;
  char filename[MAXLINE];
#ifdef _MESSG
  int myid = g->smpi->myid;
#else
  int myid = 0;
#endif
   
  fp = io_fopen(build_filename2(filename, MAXLINE, description, "parents", myid, ".txt", UNSET_INT), "w", TRUE);
  
  fprintf(fp, "my_nnodes %d nnodes %d nelems3d %d nelems2d %d nelems1d %d\n", g->my_nnodes,g->nnodes,g->nelems3d,g->nelems2d,g->nelems1d);

  fprintf(fp, "ND       i    r_pe    r_id   0r_pe   0r_id   1r_pe   1r_id\n");
  for (i=0;i<g->nnodes;i++){
#ifdef _MESSG
    if(g->node[i].original_id == UNSET_INT)
  fprintf(fp,"ND%8d%8d%8d%8d%8d%8d%8d \n", i, g->node[i].resident_pe, g->node[i].resident_id,  g->node[i].parent_res_pe[0], g->node[i].parent_res_id[0],  g->node[i].parent_res_pe[1], g->node[i].parent_res_id[1]);
#else
  fprintf(fp,"ND %d     %d          %d           %d                  %d \n",i, g->node[i].id, g->node[i].id, myid, g->node[i].id);
#endif
  }
  
  fclose(fp);
  return;
}
