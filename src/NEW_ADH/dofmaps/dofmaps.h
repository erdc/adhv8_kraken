#ifndef H_DOFMAPS_
#define H_DOFMAPS_


//cg_maps
void get_cell_dofs(int *local_dofs, int *fmaplocal, int nnodes, int *local_node_ids ,int elem_nvars, int *elem_vars, int *nodal_nvars, int **nodal_vars);
void get_cell_dofs_2(int *local_dofs, int nnodes, int *local_node_ids ,int elem_nvars, int *elem_vars, SMAT_PHYSICS *node_physics_mat, int *nodal_physics_mat_id);
int get_cg_dof(int var, int NodeID, SMAT_PHYSICS *node_physics_mat, int *nodal_physics_mat_id);

//general_maps
void local_dofs_to_global_dofs(int *global_dofs,int ndofs_on_ele,int *dofs,int *local_range,int *ghosts);
int global_to_local(int global_col_no,int local_size,int *ghosts, int nghost);


#endif