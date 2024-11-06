#ifndef _H_JACOBIAN_
#define _H_JACOBIAN_


void assemble_jacobian(SMODEL_SUPER *sm, SGRID *grid);
void load_global_mat_split_CSR(double *vals, int *indptr, int *indices, double *off_diag_vals, int *off_diag_indptr, int *off_diag_indices, double **elem_mat, int ndofs_ele, int *dofs, int *global_dofs, int *local_range);
void load_global_mat_CSR(double *vals, int *indptr, int *indices, double **elem_mat, int ndofs_ele, int *global_dofs, int *local_range);
int binary_search_part(int *arr, int start, int end, int target);
void perturb_var(double **elem_mat, SMODEL_SUPER *sm, SELEM_PHYSICS *elem_physics, 
    int ie, int nodes_on_element, int nvar_ele, int *elem_vars ,int perturb_var_code, int nsubModels, int ele_var_no, int *NodeIDs, int DEBUG);
void elem_matrix_deriv(int node_no, int dof_no, int nnodes, int elem_nvars, double *local1, double *local2, double **mat, double diff_ep);

#endif