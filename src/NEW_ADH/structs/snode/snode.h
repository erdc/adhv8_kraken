#ifndef H_SNODE_
#define H_SNODE_

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++//

typedef struct {

    int id;                 /* local id */
    int string;             /* string attached to node if there is one */
    int edge_string;        /* edge string associated with node if there is one */
	int node_string;        /* node string associated with node if there is one */
    int original_id;        /* original node number - UNSET_INT if not an original      node -- different than global_node?*/
    int parent[2];          /* (node_ladj/radj) when a node is adapted, it is averaged   based on these two original nodes. */
    int parent_res_pe[2];   /* (node_ladj/radj) when a node is adapted, it is averaged   based on these two original nodes. */
    int parent_res_id[2];   /* (node_ladj/radj) when a node is adapted, it is averaged   based on these two original nodes. */
    int level;              /* the level of the node */
    int block;              /* which block a node is assigned */
    int els_flag;           /* ?? */
    int bflag;              /* boundary flag :: 0 = surface, 1 = bottom, 2 = sidewall (3d) */
    double x, y, z;         /* node cartesian postion */
    double area;            /* the total area of elements surrounded the node */

    int nelems_connected;   /* the total number of elements connected to this node */
    int *elemID;            /* an array of element IDs connected to this node */
    int gid;                /* (global_num) global id */
    int global_surf_id;     /* for winds/waves I/O global surface id = gid in 2D or UNSET_INT in 3D if not a surface node */
    int global_bed_id;      /* for sediment I/O global surface id = gid in 2D or UNSET_INT in 3D if not a bed node */

    int myid;                /* myid needed for node_order */
    int resident_pe;        /* (sd) owning processor */
    int resident_id;        /* (rnode) node number on OWNING precessor - if this is     ghost, this will be different than calling node */

} SNODE;


/*********************************************************/
/* struct methods -------------------------------------- */
void snode_alloc_array(SNODE **node, int nnodes);
void snode_init(SNODE *);
void snode_init_array(SNODE *, int);
void snode_init_alloc_array(SNODE **node, int nnodes);
void snode_copy(SNODE *, SNODE);
void snode_copy_array(SNODE *, SNODE *, int);
void snode_free(SNODE *);
void snode_free_array(SNODE *node, int nnodes);
void snode_printScreen(SNODE);
void snode2svect(SNODE *nd, SVECT *v, int nnodes);
void snode2svect2d(SNODE *nd, SVECT2D *v, int nnodes);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
