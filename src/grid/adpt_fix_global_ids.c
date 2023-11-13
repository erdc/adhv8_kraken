#include "global_header.h"

void adpt_fix_global_ids(SGRID *grid) {
    
    int n_old, n_new, n_strt, n_end, n_orig;
    int n_old_bed, n_new_bed, n_strt_bed, n_end_bed, n_orig_bed;
    int n_old_surf, n_new_surf, n_strt_surf, n_end_surf, n_orig_surf;
    int i, ierr, inode, *gid, *surf, *bed;
    
    gid = (int *) tl_alloc(sizeof(int), grid->nnodes);
    surf = (int *) tl_alloc(sizeof(int), grid->nnodes);
    bed = (int *) tl_alloc(sizeof(int), grid->nnodes);
    n_orig = 0;
    n_orig_bed = 0;
    n_orig_surf = 0;
    /*count the number of new nnodes that I own */
    for (i=0;i<grid->my_nnodes;i++){
        if(grid->node[i].original_id != UNSET_INT){
            n_orig++;
            if(grid->ndim == 3){
                if(grid->nodeID_3d_to_2d_sur[i] != UNSET_INT){
                    n_orig_surf++;
                }
                if(grid->nodeID_3d_to_2d_bed[i] != UNSET_INT){
                    n_orig_bed++;
                }
            }
        }
    }
    
    /* number of new nodes that need gids fixed to a contigous set */
    n_new = grid->my_nnodes - n_orig;
    n_new_surf = grid->my_nnodes_sur - n_orig_surf;
    n_new_bed = grid->my_nnodes_bed - n_orig_bed;
    
#ifdef _MESSG
    ierr = MPI_Scan(&n_new, &n_strt, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM);
    if(ierr != MPI_SUCCESS) messg_err(ierr);
    n_strt = n_strt - n_new;
    if(grid->ndim==3){
        ierr = MPI_Scan(&n_new_surf, &n_strt_surf, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM);
        if(ierr != MPI_SUCCESS) messg_err(ierr);
        n_strt_surf = n_strt_surf - n_new_surf;
        ierr = MPI_Scan(&n_new_bed, &n_strt_bed, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM);
        if(ierr != MPI_SUCCESS) messg_err(ierr);
        n_strt_bed = n_strt_bed - n_new_bed;
    }
#else
    n_strt = grid->initial_nnodes;
    if(grid->ndim == 3){
        n_strt_surf = grid->initial_nnodes_bed; /* bed and surface node counts shold be equal and there is no initial surface nodes */
        n_strt_bed = grid->initial_nnodes_bed;
    }
#endif
    
    
    
    /* might not need to start at 0, but for safety we will */
    for(i=0; i<grid->my_nnodes; i++){
        if(grid->node[i].original_id == UNSET_INT){
            gid[i] = grid->orig_macro_nnodes + n_strt++;
            if(grid->ndim == 2){
                surf[i] = gid[i];
                bed[i] = gid[i];
            }else{
                if(grid->nodeID_3d_to_2d_bed[i] != UNSET_INT){
                    bed[i] = grid->orig_macro_nnodes_bed + n_strt_bed++;
                }else{
                    bed[i] = UNSET_INT;
                }
                if(grid->nodeID_3d_to_2d_sur[i] != UNSET_INT){
                    surf[i] = grid->orig_macro_nnodes_sur + n_strt_surf++;
                }else{
                    surf[i] = UNSET_INT;
                }
            }
        }else{
            gid[i] = grid->node[i].gid;
            surf[i] = grid->node[i].global_surf_id;
            bed[i] = grid->node[i].global_bed_id;
        }
    }
    
    /* set new total nodes, elements on mesh and owning elements for printing, etc. */
    int my_elem = 0;
    int elem_tot = 0;
    for (i = 0; i < grid->nelems2d; i++) {
#ifdef _MESSG
        grid->elem2d[i].my_pe = 0;
        if ((grid->node[grid->elem2d[i].nodes[0]].resident_pe != grid->node[grid->elem2d[i].nodes[1]].resident_pe) &&
            (grid->node[grid->elem2d[i].nodes[0]].resident_pe != grid->node[grid->elem2d[i].nodes[2]].resident_pe) &&
            (grid->node[grid->elem2d[i].nodes[1]].resident_pe != grid->node[grid->elem2d[i].nodes[2]].resident_pe)){
            for (inode = 0; inode < 3; inode++){
                if (grid->node[grid->elem2d[i].nodes[inode]].resident_pe > grid->smpi->myid) grid->elem2d[i].my_pe++;
            }
        }
        else{
            for (inode = 0; inode < 3; inode++){
                if (grid->node[grid->elem2d[i].nodes[inode]].resident_pe == grid->smpi->myid) grid->elem2d[i].my_pe++;
            }
        }
        if (grid->elem2d[i].my_pe > 1) my_elem++;
    }
    
    comm_set_keys(grid);
    
    comm_update_int(gid,1,grid->smpi);
    comm_update_int(surf,1,grid->smpi);
    comm_update_int(bed,1,grid->smpi);
    for(i=0;i<grid->nnodes;i++){
        grid->node[i].gid = gid[i] ;
        grid->node[i].global_surf_id = surf[i];
        grid->node[i].global_bed_id = bed[i];
    }
    
    
    n_old = grid->my_nnodes;
    n_old_surf = grid->my_nnodes_sur;
    n_old_bed = grid->my_nnodes_bed;
    
    ierr = MPI_Allreduce(&n_old, &n_end, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM);
    if(ierr != MPI_SUCCESS) messg_err(ierr);
    ierr = MPI_Allreduce(&n_old_surf, &n_end_surf, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM);
    if(ierr != MPI_SUCCESS) messg_err(ierr);
    ierr = MPI_Allreduce(&n_old_bed, &n_end_bed, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM);
    if(ierr != MPI_SUCCESS) messg_err(ierr);
    ierr = MPI_Allreduce(&my_elem, &elem_tot, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM);
    if(ierr != MPI_SUCCESS) messg_err(ierr);
    
#else
        grid->elem2d[i].my_pe = 3;
    }

    for(i=grid->orig_macro_nnodes;i<grid->my_nnodes;i++){
        grid->node[i].gid = i ;
        grid->node[i].global_surf_id = surf[i]; /* gkc modifying */
        grid->node[i].global_bed_id = bed[i]; /* gkc adding */
    }
    elem_tot=grid->nelems2d;
    n_end = grid->my_nnodes;
#endif

grid->macro_nelems2d = elem_tot;
grid->macro_nnodes = n_end; /* sets for printing */
grid->macro_nnodes_sur = n_end_surf; /* sets for printing */
grid->macro_nnodes_bed = n_end_bed; /* sets for printing */


gid = (int *) tl_free(sizeof(int), grid->nnodes, gid);
surf = (int *) tl_free(sizeof(int), grid->nnodes, surf);
bed = (int *) tl_free(sizeof(int), grid->nnodes, bed);
return;
}
