#include "global_header.h"

/*! 
 \brief This routine initializes the communications data
 */
void comm_set_keys(SGRID *grid)
{
    int ii = 0, jj = 0;           /* loop counters */
    int ie = 0;                   /* loop counter over the elements */
    int i_processor = 0;          /* loop counter over the processors */
    int *is_border = NULL;        /* flags for border nodes */
    int local_owner = 0;          /* check */
    
    /*! Allocate memory for the message set up */
    is_border = (int *) tl_alloc(sizeof(int), grid->nnodes);
    
    /*! Clear the Communication Buffers, etc. */
    for (i_processor = 0; i_processor < grid->smpi->npes; i_processor++)
    {
        
        messg_buffer_free(grid->smpi->send_msg + i_processor);
        messg_buffer_init(grid->smpi->send_msg + i_processor, i_processor);
        messg_buffer_free(grid->smpi->recv_msg + i_processor);
        messg_buffer_init(grid->smpi->recv_msg + i_processor, i_processor);
        
        if (grid->smpi->send_key[i_processor].key != NULL){
            grid->smpi->send_key[i_processor].key = (int *) tl_free(sizeof(int), grid->smpi->send_key[i_processor].size, grid->smpi->send_key[i_processor].key);
            
        }
        grid->smpi->send_key[i_processor].size = 0;
        grid->smpi->recv_init[i_processor] = UNSET_INT;
        grid->smpi->nsend[i_processor] = 0;
        grid->smpi->nrecv[i_processor] = 0;
        
        if (grid->smpi->send_key_surf[i_processor].key != NULL){
            grid->smpi->send_key_surf[i_processor].key = (int *) tl_free(sizeof(int), grid->smpi->send_key_surf[i_processor].size, grid->smpi->send_key_surf[i_processor].key);
            
        }
        grid->smpi->send_key_surf[i_processor].size = 0;
        grid->smpi->recv_init_surf[i_processor] = UNSET_INT;
        grid->smpi->nsend_surf[i_processor] = 0;
        grid->smpi->nrecv_surf[i_processor] = 0;
        
    }
    
    /*! Set up the Number of Receive Messages */
    for (ii = grid->my_nnodes; ii < grid->nnodes; ii++)
    {
        local_owner = grid->node[ii].resident_pe;
        if(local_owner<0 || local_owner>=grid->smpi->npes) printf("myid %d problem with node %d local_owner %d gid %d\n", grid->smpi->myid,ii,local_owner, grid->node[ii].gid);
        if(grid->smpi->recv_init[local_owner] == UNSET_INT) grid->smpi->recv_init[local_owner] = ii;
        if(grid->node[ii].global_surf_id != UNSET_INT) {
            if (grid->smpi->recv_init_surf[local_owner] == UNSET_INT) {
                if(grid->type == COLUMNAR) {grid->smpi->recv_init_surf[local_owner] = grid->nodeID_3d_to_2d_sur[ii];}
                else{grid->smpi->recv_init_surf[local_owner] = ii;}
            }
            grid->smpi->nrecv_surf[local_owner]++;
        }
        grid->smpi->nrecv[local_owner]++;
    }
    
    /*! Set up the Number of Receive Messages */
    /* for (ii = grid->my_nnodes; ii < grid->nnodes; ii++)
     {
     local_owner = grid->node[ii].resident_pe;
     if(grid->smpi->recv_init[local_owner] == 0){
     grid->smpi->recv_init[local_owner] = ii;
     grid->smpi->nrecv[local_owner] += ii;
     }
     
     }
     */
    
    /*! We need to set up the send_key. */
    /*! The nodes that I own can be ghost nodes on other processors. */
    /*! They might be scattered all over the place within my memory. */
    /*! But we know that they will be in ascending order on the other PE. */
    /*! So we need to make a quick list, or a key, so that during a comm_update_xxxx, */
    /*! we can pluck out the right values */
    for (i_processor = 0; i_processor < grid->smpi->npes; i_processor++)
    {
        if (grid->smpi->nrecv[i_processor] > 0)
        {
            /* Re-initialize the is_border array for each processor */
            for (ii = 0; ii < grid->nnodes; ii++)
            {
                is_border[ii] = NO;
            }
            /* loop over the 3D elements identifying the border nodes */
            for (ie = 0; ie < grid->nelems3d; ie++)
            {
                if (grid->node[grid->elem3d[ie].nodes[0]].resident_pe == i_processor ||
                    grid->node[grid->elem3d[ie].nodes[1]].resident_pe == i_processor ||
                    grid->node[grid->elem3d[ie].nodes[2]].resident_pe == i_processor ||
                    grid->node[grid->elem3d[ie].nodes[3]].resident_pe == i_processor)
                {
                    is_border[grid->elem3d[ie].nodes[0]] = YES;
                    is_border[grid->elem3d[ie].nodes[1]] = YES;
                    is_border[grid->elem3d[ie].nodes[2]] = YES;
                    is_border[grid->elem3d[ie].nodes[3]] = YES;
                }
            }
            /* loop over the 2D elements identifying the border nodes */
            for (ie = 0; ie < grid->nelems2d; ie++)
            {
                if (grid->node[grid->elem2d[ie].nodes[0]].resident_pe == i_processor ||
                    grid->node[grid->elem2d[ie].nodes[1]].resident_pe == i_processor ||
                    grid->node[grid->elem2d[ie].nodes[2]].resident_pe == i_processor)
                {
                    is_border[grid->elem2d[ie].nodes[0]] = YES;
                    is_border[grid->elem2d[ie].nodes[1]] = YES;
                    is_border[grid->elem2d[ie].nodes[2]] = YES;
                }
            }
            /* loop over the 1D elements identifying the border nodes */
            for (ie = 0; ie < grid->nelems1d; ie++)
            {
                if (grid->node[grid->elem1d[ie].nodes[0]].resident_pe == i_processor ||
                    grid->node[grid->elem1d[ie].nodes[1]].resident_pe == i_processor)
                {
                    is_border[grid->elem1d[ie].nodes[0]] = YES;
                    is_border[grid->elem1d[ie].nodes[1]] = YES;
                }
            }
            /*! Count the border nodes. Note that I only look in the nodes */
            /*! that this processor owns. Ghost nodes will get flagged as well */
            /*! but since the other processor owns it, this PE will not be */
            /*! responsible for it. */
            
            for (ii = 0; ii < grid->my_nnodes; ii++)
            {
                if (is_border[ii] == YES)
                {
                    grid->smpi->nsend[i_processor]++;
                    if (grid->node[ii].global_surf_id != UNSET_INT) {
                        grid->smpi->nsend_surf[i_processor]++;
                    }
                }
            }
            //if(grid->smpi->nsend[i_processor] != grid->smpi->nsend_surf[i_processor]) printf("MYID %d nsend error %d %d \n", grid->smpi->nsend[i_processor],grid->smpi->nsend_surf[i_processor]);
            grid->smpi->send_key[i_processor].size=grid->smpi->nsend[i_processor];
            grid->smpi->send_key_surf[i_processor].size=grid->smpi->nsend_surf[i_processor];
            
            /*! Allocate the send_key Buffer Array */
            if (grid->smpi->send_key[i_processor].size != 0)
            {
                grid->smpi->send_key[i_processor].key = (int *) tl_alloc(sizeof(int), grid->smpi->send_key[i_processor].size);
            }
            if (grid->smpi->send_key_surf[i_processor].size != 0)
            {
                grid->smpi->send_key_surf[i_processor].key = (int *) tl_alloc(sizeof(int), grid->smpi->send_key_surf[i_processor].size);
            }
            
            /*! Now go through the nodes and copy down the IDs */
            /*! so that we have our key, in ascending order */
            for (ii = 0, jj = 0, ie=0; ii < grid->my_nnodes; ii++)
            {
                if (is_border[ii] == YES)
                {
                    grid->smpi->send_key[i_processor].key[jj++] = ii;
                    if (grid->node[ii].global_surf_id != UNSET_INT) {
                        if(grid->type == COLUMNAR){grid->smpi->send_key_surf[i_processor].key[ie++] = grid->nodeID_3d_to_2d_sur[ii];}
                        else{grid->smpi->send_key_surf[i_processor].key[ie++] = ii;}
                    }
                }
            }
            assert(jj == grid->smpi->send_key[i_processor].size);
            assert(ie == grid->smpi->send_key_surf[i_processor].size);
        }  /* Skip All Processors with which I am not commmunicating */
    }   /* Loop over Processors */
    
    /* JLH
     * I'm putting this in to get adaption working again;
     * other pieces of the code need to have nrecv and nsend
     */
    /* count the ghost nodes and set nrecv */
    //    for(ii = my_nnode; ii < nnode; ii++)
    //      {
    //        nrecv[node_map[ii].owner_pe]++;
    //      }
    //    for(i_processor = 0; i_processor < npes; i_processor++)
    //      {
    //        nsend[i_processor] = send_key[i_processor].size;
    //      }
    
    /*! Free the local memory */
    is_border = (int *) tl_free(sizeof(int), grid->nnodes, is_border);
    return;
}
