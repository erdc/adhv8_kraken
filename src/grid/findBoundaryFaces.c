#include "global_header.h"

void findBoundaryFaces(SGRID *grid)
{

  int *faceNormal;
  int nfaces = grid->nelems3d * NDONTER;  /* number of faces on a tetrahedron */
  SNODE *node;

  faceNormal = (int *) tl_alloc(sizeof(int), nfaces);
  node = (SNODE *) tl_alloc(sizeof(int), NDONTRI);

  for (ielem = 0; ielem < grid->nelems3d; ielem++) {

  for (ielem = 0; ielem < grid->nelems3d; ielem++) {
    // FACE 1
    node[0] = grid->elem3d[ielem]->node[0];
    node[1] = grid->elem3d[ielem]->node[1];
    node[2] = grid->elem3d[ielem]->node[2];
    nrml[] = elem2d_normal[node];

    // FACE 2
    node[0] = grid->elem3d[ielem]->node[0];
    node[3] = grid->elem3d[ielem]->node[1];
    node[1] = grid->elem3d[ielem]->node[2];
    nrml[] = elem2d_normal[node];

    // FACE 3
    node[0] = grid->elem3d[ielem]->node[0];
    node[3] = grid->elem3d[ielem]->node[1];
    node[2] = grid->elem3d[ielem]->node[2];
    nrml[] = elem2d_normal[node];

    // FACE 4
    node[2] = grid->elem3d[ielem]->node[0];
    node[3] = grid->elem3d[ielem]->node[1];
    node[1] = grid->elem3d[ielem]->node[2];
    nrml[] = elem2d_normal[node];

  }
