#include "adh.h"
void sdvar_alloc_init(SDVAR *sdvar, int nnode, int n_dvar, int n_dvar_elem_dbl, int n_dvar_elem_int, int nnode_dvar, int nelem_dvar){
	int j;
    printf("In sdvar alloc init\n");
	//allocate and set up the sdvar object
	// node-based surface water variables
    sdvar->n_dvar = n_dvar; //number of active dependent variables
    sdvar->nnode_dvar = nnode_dvar; //number of active nodes
    printf("set 2 vals\n");
    //allocate the double array
    //double **nodal_dvar //[n_dvar][nnode_dvar]
    sdvar->nodal_dvar = NULL;
    if (sdvar->n_dvar > 0){
    	sdvar->nodal_dvar = (double**) tl_alloc(sizeof(double*), n_dvar);
    	for(j=0;j<n_dvar;j++){
        	sdvar->nodal_dvar[j] = (double*) tl_alloc(sizeof(double), nnode_dvar);
    	}
    	//intialize all data to 0
    	sarray_init_double_2d(sdvar->nodal_dvar, n_dvar, nnode_dvar);
    }
    printf("In sdvar alloc init 2\n");
    // element-based surface water variables
    sdvar->n_dvar_elem_dbl = n_dvar_elem_dbl;
    sdvar->nelem_dvar = nelem_dvar;//number of active elements

    sdvar->elem_dvar = NULL;
    if(sdvar->n_dvar_elem_dbl > 0){
    	sdvar->elem_dvar = (double**) tl_alloc(sizeof(double*), n_dvar_elem_dbl);
    	for(j=0;j<n_dvar_elem_dbl;j++){
        	sdvar->elem_dvar[j] = (double*) tl_alloc(sizeof(double), nelem_dvar);
    	}
    	//intialize all data to 0
    	sarray_init_double_2d(sdvar->elem_dvar, n_dvar_elem_dbl, nelem_dvar);
    }

    sdvar->n_dvar_elem_int = n_dvar_elem_int;
    sdvar->elem_flags = NULL;
    if(sdvar->n_dvar_elem_int > 0){
    	sdvar->elem_flags = (int**) tl_alloc(sizeof(int*), n_dvar_elem_int);
    	for(j=0;j<n_dvar_elem_int;j++){
        	sdvar->elem_flags[j] = (int*) tl_alloc(sizeof(int), nelem_dvar);
    	}
    	//intialize all data to 0
    	sarray_init_int_2d(sdvar->elem_flags, n_dvar_elem_int, nelem_dvar);
    }


    //maps, only need one if nnodes != nnode_dvar
    sdvar->dvar_node_map = NULL; //[nnode] array goes NodeID->index within nodal_dvar
    //if (nnode!=nnode_dvar){
    	//create map here?, will also need the nodal physics to get this actually
    //}
    //same for dvar_active
    sdvar->dvar_active = NULL;
    //int *dvar_active; //n_dvar array that has the NodeID, convenient for printing out stuff
    //and for elem_map
    sdvar->dvar_elem_map = NULL;
    //int *dvar_elem_map; //nelem array takes elem # -> index within elem_var


}