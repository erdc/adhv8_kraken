
// Constants.h

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__ 


/* the nodes on an edge of a 2d element */
//static const int nd_on_2dedge[NEDGEPRELM][NDPREDG] = { {0, 1}, {0, 2}, {1, 2} };

/* the nodes on an edge of a 3d element */
//static const int nd_on_3dedge[NEDGEPRELM][NDPREDG] = { {0, 1}, {0, 2}, {0, 3}, {1, 2}, {3, 2}, {1, 3} };

/* the nodes on a face */
//static const int nd_on_fc[NFCPRELM][NDPRFC] = { {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1} };

/* face opposite a node - redundant with this numbering but maintained for readability */
//static const int fc_opp_nd[NDPRELM] = { 0, 1, 2, 3 };

/* node opposite a face - also redundant */
//static const int nd_opp_fc[NFCPRELM] = { 0, 1, 2, 3 };

/* edge opposite a node (2d element) - redundant with this numbering but maintained for readability */
//static const int edge_opp_nd[NDPRFC] = { 0, 1, 2 };

/* node opposite an edge (2d element) - also redundant */
//static const int nd_opp_edge[NEDGEPRFC] = { 0, 1, 2 };

/* finds an edge given the local */
//static const int find_edge[NDPRELM][NDPRELM] = { {UNSET_INT, 0, 1, 2}, {0, UNSET_INT, 3, 5}, {1, 3, UNSET_INT, 4}, {2, 5, 4, UNSET_INT} };

#endif

