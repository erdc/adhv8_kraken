#ifndef H_EXTRUSION_
#define H_EXTRUSION_

/***********************************************************/
/***********************************************************/
/***********************************************************/

#include "global_header.h"
#include "xdmf_generator.h"

#define STANDARD_MESH 0
#define HYBRID_MESH 1
#define PRISM_MESH 2
#define TET_MESH 3
#define MIXED_ELEMENT_MESH 4
#define ALLOC_INC 5
#define VALID 1
#define INVALID 0

// data structures
typedef struct{
    int n, flag;
    double *ztop;
    double *dz;
} zbins;

typedef struct {
    int node1;
    int valid;         /* 1 means the node is valid; 0 means ignore */
} NODE;

typedef struct {
    int count;         /* the number of nodes in the string */
    int allocated;     /* the number of spots allocated (count must be <= allocated) */
    NODE *nodes;       /* the nodes in the string */
    int id_2d;         /* the new id of the string */
    int id_3d;
    int valid;         /* 1 means that this node string exits, 0 means ignore */
} NODE_STRING;

typedef struct {
    int node1;
    int node2;
    int valid;         /* 1 means the edge is valid; 0 means ignore */
} EDGE;

typedef struct {
    int count;         /* the number of edges in the string */
    int allocated;     /* the number of spots allocated in the string (count must be <= allocated) */
    EDGE *edges;       /* the edges in the string */
    int id_2d;         /* the new id of the string */
    int valid;         /* 1 means that this edge string exists, 0 means ignore */
} EDGE_STRING;

typedef struct {
    int elem_id;
    int face_num;
    int node1;
    int node2;       /* node1 and node2 are the original 2D surface nodes that comprise the edge that was extruded to created this face */
    int nnodes;
    int node[4];
} FACE;

typedef struct {
    int string_id;   /* the original id of the edge string in the 2D surface mesh */
    int id_3d;
    int count;
    int allocated;
    FACE *faces;
    int valid;        /* 1 means the edge is valid; 0 means ignore */
} FACE_STRING;

// prototypes
int extrudeAdH(int flag_bins,int flag_min_layer, int minimum_layers, char *adh_root, int mesh_type);
void create_3d_geometry_adh(SGRID *geo2d, SGRID *geo3d, int *node_count, int *node_start, int *node_bin, double dz, double *head, zbins *bin, int bin_flag, int mesh_type);
void create_3d_bc_adh(SMODEL *mod, SGRID *geo3d, int *node_count, int *node_start, int mesh_type, int *node_flags2d, int *node_flags3d, int *new_node_number2d, int *new_node_number3d);
void create_sw3_hot(SMODEL *mod, SGRID *grid3d, int *node_start, int *node_count);
double find_dz_bin(zbins *bin, double z, int mat);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif

