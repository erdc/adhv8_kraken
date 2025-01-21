void global_to_local_dbl_cg_arr_map(double *local_vals, int *nodeIDs, int nnodes, int *map_array, double *global_vals){
        for(int i =0;i<nnodes;i++){
                local_vals[i] = global_vals[map_array[nodeIDs[i]]];
        }
}
