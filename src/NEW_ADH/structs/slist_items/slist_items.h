#ifndef H_SLISTITEM_
#define H_SLISTITEM_

// dependencies ::
//  svect
//  svect2d

struct FACE_LIST_ITEM
{
  int nd1;          /* the first node on the face */
  int nd2;          /* the second node on the face */
  int nd3;          /* the third node on the face */
  int ie1;          /* the first element on the face */
  int ie2;          /* the second element on the face */
  struct FACE_LIST_ITEM *next;  /* the next item in the linked list */
};
typedef struct FACE_LIST_ITEM FACE_LIST_ITEM;   /* special line for linked list structure - refers to itself */

struct ELEM3D_LIST_ITEM
{
  int nd1;          /* the first node */
  int nd2;          /* the second node */
  int nd3;          /* the third node */
  int nd4;          /* the fourth node */
  int ielem;            /* the element number */
  int mat;          /* material id*/
  struct ELEM3D_LIST_ITEM *next;    /* the next item in the linked list */
};              /* a linked list item for 3D elements */
typedef struct ELEM3D_LIST_ITEM ELEM3D_LIST_ITEM;   /* special line for linked list structure - refers to itself */

struct ELEM2D_LIST_ITEM
{
  int nd1;          /* the first node */
  int nd2;          /* the second node */
  int nd3;          /* the third node */
  int mat;          /* material id*/
  int ielem;            /* the element number */
  struct ELEM2D_LIST_ITEM *next;    /* the next item in the linked list */
};              /* a linked list item for 2D elements */
typedef struct ELEM2D_LIST_ITEM ELEM2D_LIST_ITEM;   /* special line for linked list structure - refers to itself */

struct ELEM1D_LIST_ITEM
{
  int nd1;          /* the first node */
  int nd2;          /* the second node */
  int ielem;            /* the element number */
  struct ELEM1D_LIST_ITEM *next;    /* the next item in the linked list */
};              /* a linked list item for 1D elements */
typedef struct ELEM1D_LIST_ITEM ELEM1D_LIST_ITEM;   /* special line for linked list structure - refers to itself */

struct EDGE_LIST_ITEM
{
  int nd1;          /* the first node on the edge */
  int nd2;          /* the second node on the edge */
  int new_node;         /* the new node */
  int rank;         /* the rank of the edge */
  struct EDGE_LIST_ITEM *next;  /* the next item in the linked list */
};              /* a linked list item for edges */
typedef struct EDGE_LIST_ITEM EDGE_LIST_ITEM;   /* special line for linked list structure - refers to itself */


struct NODE_LIST_ITEM
{
  struct NODE_LIST_ITEM *next;  /* the next item in the linked list */
  int rnode;            /* relative node number on owning processor */
  int sd;           /* owning processor */
  int local;            /* local node number on myid */
//#ifdef _MESSG
  int global;
  int global_surf;
  SVECT coords;
//#endif
};              /* a linked list item for global node numbers */
typedef struct NODE_LIST_ITEM NODE_LIST_ITEM;   /* special line for linked list structure - refers to itself */


struct ELEM_REF_LIST_ITEM
{
      struct ELEM_REF_LIST_ITEM *next;  /* the next item in the linked list */
        int ielem;            /* the element number */
};              /* a linked list item for element numbers */
typedef struct ELEM_REF_LIST_ITEM ELEM_REF_LIST_ITEM;   /* special line for linked list structure - refers to itself */


/* a linked list item for centroids (points) */
struct CENT_LIST_ITEM
{
  SVECT2D vect;                  /* the coordinates of the point */
  int index;                    /* an index for this column (e.g., column id) */
  struct CENT_LIST_ITEM *next;  /* the next item in the linked list */
};
typedef struct CENT_LIST_ITEM CENT_LIST_ITEM;   /* special line for linked list structure - refers to itself */


struct ID_LIST_ITEM
{
  int id;           /* the id of the member */
  struct ID_LIST_ITEM *next;    /* the next item in the linked list */
  struct ID_LIST_ITEM *prev;    /* the next item in the linked list */
};
typedef struct ID_LIST_ITEM ID_LIST_ITEM;   /* special line for linked list structure - refers to itself */


// this is an interesting list:
// -- an array of num_midpts of these are allocated in the code
// -- each one points to the top midpoint of the column
// -- this midpoint may not be on the surface, 
// ---- for example, if the first midpoint falls between the top two vertical nodes.
struct MIDPT_LIST_ITEM
{
  SVECT vect;              /* the coordinates of the point */
  int index;               /* an index for this column (e.g., column id) */
  int node1;               /* this midpoint lies on the edge formed by node1 and node2 */
  int node2;
  int surf_node1;          /* the two surface nodes that lie above node1 and node2 */
  int surf_node2;
  int vertical;            /* 0 = not a vertical edge; 1 = this is a vertical edge */
  int column1;             /* the indices of the 2 columns that share this edge */
  int column2;

  int elem_upper[2];          /* the element that lies above this edge/midpt in each of the two columns */
  int elem_lower[2];          /* the element that lies below this edge/midpt in each of the two columns */
  double value[5];            /* data associated with this midpoint */
                              /* value[0] = pressure at midpt;
                               * value[1] = pressure from positive perturbation at the "left" surface node
                               * value[2] = pressure from positive perturbation at the "right" surface node
                               * value[3] = pressure from negative perturbation at the "left" surface node
                               * value[4] = pressure from negative perturbation at the "right" surface node
                               */
  struct MIDPT_LIST_ITEM *next; /* the next item in the linked list */
  struct MIDPT_LIST_ITEM *prev; /* the next item in the linked list */
};              /* a linked list item for centroids (points) */
typedef struct MIDPT_LIST_ITEM MIDPT_LIST_ITEM;


/* structure used to rank the edges */
typedef struct
{
  int number;           /* which edge */
  double length;        /* the length of the edge */
} EDGE_RANK;            /* structure to sort the edges */

/*********************************************************/
/* struct methods -------------------------------------- */
int add_node(NODE_LIST_ITEM ** head, int local_index, int global_index, SVECT coords, int sd, int rnode, int global_surf);
int add_node_tmp(NODE_LIST_ITEM ** head, int local_index, int global_index);
int add_ghost_node(NODE_LIST_ITEM **head, int index, int my_nnode, int *num_ghosts);
int add_surf_node(NODE_LIST_ITEM **, int);
int search_ghost_node(NODE_LIST_ITEM * head, int index);
int add_elem2d_list(ELEM2D_LIST_ITEM ** head, int index, int *node, int iel_mat);
int add_elem3d_list(ELEM3D_LIST_ITEM ** head, int global_index, int *node, int iel_mat);
int free_elem3d_list(ELEM3D_LIST_ITEM * elem3d_list);
int free_elem2d_list(ELEM2D_LIST_ITEM * elem2d_list);
int free_node_list(NODE_LIST_ITEM * head);
int free_ghost_nodes(NODE_LIST_ITEM *** ghost_nodes, int npes);

#endif
