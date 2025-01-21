#ifndef H_SDVAR_
#define H_SDVAR_
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
typedef struct {

    // node-based surface water variables
    int n_dvar; //number of active dependent variables
    int nnode_dvar; //number of active nodes
    double **nodal_dvar; //[n_dvar][nnode_dvar]

    // element-based surface water variables
    int n_dvar_elem_dbl;
    int n_dvar_elem_int;
    int nelem_dvar;//number of active elements
    double **elem_dvar;//[n_dvar_elem_dbl][nelem_dvar]
    int **elem_flags;//[n_dvar_elem_int][nelem_dvar]

    //maps
    int *dvar_node_map; //[nnode] array goes NodeID->index within nodal_dvar
    int *dvar_active; //n_dvar array that has the NodeID, convenient for printing out stuff
    int *dvar_elem_map; //nelem array takes elem # -> index within elem_var

} SDVAR;
void sdvar_alloc_init(SDVAR *sdvar, int nnode, int n_dvar, int n_dvar_elem_dbl, int n_dvar_elem_int, int nnode_dvar, int nelem_dvar);
#endif

