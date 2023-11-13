#include "global_header.h"

void print_grid_to_file(SGRID *g, char *description)
{
  FILE *fp;
  int i;
  char filename[MAXLINE];
#ifdef _MESSG
  int myid = g->smpi->myid;
#else
  int myid = 0;
#endif
   
  fp = io_fopen(build_filename2(filename, MAXLINE, description, "", myid, ".txt", UNSET_INT), "w", TRUE);
  
  fprintf(fp, "my_nnodes %d nnodes %d nelems3d %d nelems2d %d nelems1d %d\n", g->my_nnodes,g->nnodes,g->nelems3d,g->nelems2d,g->nelems1d);

  if(g->nelems3d > 0){
    fprintf(fp, "E4T              i   elem3d[i].id elem3d[i].gid          node0         node1         node2         node3          mat          string\n");
    for (i=0;i<g->nelems3d;i++){
#ifdef _MESSG
    fprintf(fp, "E4T %14d%14d%14d%14d%14d%14d%14d%14d%14d\n", i, g->elem3d[i].id, g->elem3d[i].gid, g->elem3d[i].nodes[0], g->elem3d[i].nodes[1], g->elem3d[i].nodes[2], g->elem3d[i].nodes[3], g->elem3d[i].mat, g->elem3d[i].string); 
#else
    fprintf(fp, "E4T %d     %d         %d          %d    %d     %d    %d\n", i, g->elem3d[i].id, g->elem3d[i].id, g->elem3d[i].nodes[0], g->elem3d[i].nodes[1], g->elem3d[i].nodes[2], g->elem3d[i].nodes[3]);
#endif
    }
  }

  if(g->nelems2d > 0){
    fprintf(fp, "E3T              i  elem2d[i].id   node0        node1         node2         djac\n");
    for (i=0;i<g->nelems2d;i++){
#ifdef _MESSG
    fprintf(fp, "E3T %14d%14d%14d%14d%14d   %10.3e \n", i, g->elem2d[i].id, g->elem2d[i].nodes[0], g->elem2d[i].nodes[1], g->elem2d[i].nodes[2], g->elem2d[i].djac);
#else
    fprintf(fp, "E3T %d     %d           %d            %d     %d    %d\n", i, g->elem2d[i].id, g->elem2d[i].id, g->elem2d[i].nodes[0], g->elem2d[i].nodes[1], g->elem2d[i].nodes[2]);
#endif
    }
  }

  if(g->nelems1d > 0){
    fprintf(fp, "E2T i elem1d[i].id node0 node1 \n");
    for (i=0;i<g->nelems1d;i++){
    fprintf(fp, "E2T %d      %d      %d    %d \n", i, g->elem1d[i].id, g->elem1d[i].nodes[0], g->elem1d[i].nodes[1]);
    }
  }
  
  fprintf(fp, "ND      i    id   gid  gsd  r_pe  r_id   blk eg_str level el_flg\n");
  for (i=0;i<g->nnodes;i++){
#ifdef _MESSG
  fprintf(fp,"ND %6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",i, g->node[i].id,  g->node[i].gid,g->node[i].global_surf_id, g->node[i].resident_pe, g->node[i].resident_id, g->node[i].      block, g->node[i].edge_string, g->node[i].level, g->node[i].els_flag);
#else
  fprintf(fp,"ND %d     %d          %d           %d                  %d \n",i, g->node[i].id, g->node[i].id, myid, g->node[i].id);
#endif
  }
  if(g->ndim ==2){
    fprintf(fp, "\n\nMESH2D\n");
    for (i=0;i<g->nelems2d;i++) fprintf(fp, "E3T  %8d %8d %8d %8d %4d\n", g->elem2d[i].id, g->elem2d[i].nodes[0], g->elem2d[i].nodes[1], g->elem2d[i].nodes[2], g->elem2d[i].mat);
    for (i=0;i<g->nnodes;i++) fprintf(fp, "ND %d %16.8e %16.8e %16.8e\n", g->node[i].id,g->node[i].x,g->node[i].y,g->node[i].z);
  }
  else{
    fprintf(fp, "\n\nMESH3D\n");
    for (i=0;i<g->nelems3d;i++) fprintf(fp, "E4T  %8d %8d %8d %8d %8d %4d\n", g->elem3d[i].id, g->elem3d[i].nodes[0], g->elem3d[i].nodes[1], g->elem3d[i].nodes[2], g->elem3d[i].nodes[3], g->elem3d[i].mat);
    for (i=0;i<g->nnodes;i++) fprintf(fp, "ND %d %16.8e %16.8e %16.8e\n", g->node[i].id,g->node[i].x,g->node[i].y,g->node[i].z);
  }
  fclose(fp);
  return;
}
