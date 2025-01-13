#ifndef H_DOFMAPS_
#define H_DOFMAPS_


//cg_maps

//uses array look up
//void get_cell_dofs(int *local_dofs, int *fmaplocal, int nnodes, int *local_node_ids ,int elem_nvars, int *elem_vars, SMAT_PHYSICS *node_physics_mat, int *nodal_physics_mat_id);
void get_cell_dofs(int *local_dofs, int *fmaplocal, int nnodes, int *local_node_ids ,int elem_nvars, int *elem_vars, SMAT_PHYSICS **node_physics_mat);
//int get_cg_dof(int var, int NodeID, int *fmaplocal, SMAT_PHYSICS *node_physics_mat, int *nodal_physics_mat_id);
int get_cg_dof(int var, int NodeID, int *fmaplocal, SMAT_PHYSICS **node_physics_mat);
//void global_to_local_dbl_cg(double *local, double *global, int *nodes, int nnodes, int var, int *fmaplocal, SMAT_PHYSICS *node_physics_mat, int *nodal_physics_mat_id);
void global_to_local_dbl_cg(double *local, double *global, int *nodes, int nnodes, int var, int *fmaplocal, SMAT_PHYSICS **node_physics_mat);
//fully implicit, uses redundant calculations
void get_cell_dofs_2(int *local_dofs, int nnodes, int *local_node_ids ,int elem_nvars, int *elem_vars, SMAT_PHYSICS *node_physics_mat, int *nodal_physics_mat_id);
int get_cg_dof_2(int var, int NodeID, SMAT_PHYSICS *node_physics_mat, int *nodal_physics_mat_id);
void global_to_local_dbl_cg_2(double *local, double *global, int *nodes, int nnodes, int var, SMAT_PHYSICS *node_physics_mat, int *nodal_physics_mat_id);
void global_to_local_SVECT2D_cg(SVECT2D *local, double *global, int *nodes, int nnodes, int varx, int vary, int *fmaplocal, SMAT_PHYSICS **node_physics_mat);
void global_to_local_dbl_cg_arr_map(double *local_vals, int *nodeIDs, int nnodes, int *map_array, double *global_vals);
//general_maps
void local_dofs_to_global_dofs(int *global_dofs,int ndofs_on_ele,int *dofs,int *local_range,int *ghosts);
int global_to_local(int global_col_no,int local_size,int *ghosts, int nghost);


#endif
