#include "adh.h"

void smpi_init(SMPI *smpi
#ifdef _MESSG
               , MPI_Comm input_comm
#endif
) {
    
#ifdef _MESSG
    int i_processor = 0;          /* Loop Counter over processors */
    int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */
    
    ierr_code = MPI_Comm_dup(input_comm, &(smpi->ADH_COMM));
    if (ierr_code != MPI_SUCCESS) {
        messg_err(ierr_code);
    }
    
    smpi->msg_status = (MPI_Status *) NULL;
    
    smpi->msg_request = (MPI_Request *) NULL;
    
    ierr_code = MPI_Comm_size(smpi->ADH_COMM, &(smpi->npes));
    
    ierr_code = MPI_Comm_rank(smpi->ADH_COMM, &(smpi->myid));
    
    smpi->msg_starttime = MPI_Wtime();
    
    /* allocates the message arrays */
    
    smpi->recv_init = (int *) tl_alloc(sizeof(int), smpi->npes);
    
    smpi->nsend = (int *) tl_alloc(sizeof(int), smpi->npes);
    
    smpi->nrecv = (int *) tl_alloc(sizeof(int), smpi->npes);
    
    smpi->recv_init_surf = (int *) tl_alloc(sizeof(int), smpi->npes);
    
    smpi->nsend_surf = (int *) tl_alloc(sizeof(int), smpi->npes);
    
    smpi->nrecv_surf = (int *) tl_alloc(sizeof(int), smpi->npes);
    
    smpi->send_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), smpi->npes);
    
    smpi->recv_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), smpi->npes);
    
    smpi->send_edge_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), smpi->npes);
    
    smpi->recv_edge_msg = (MESSG_BUFFER *) tl_alloc(sizeof(MESSG_BUFFER), smpi->npes);
    
    smpi->send_key = (MESSG_KEY *) tl_alloc(sizeof(MESSG_KEY), smpi->npes);
    
    smpi->send_key_surf = (MESSG_KEY *) tl_alloc(sizeof(MESSG_KEY), smpi->npes);
    
    smpi->nsend_edge = (int *) tl_alloc(sizeof(int), smpi->npes);
    
    smpi->nrecv_edge = (int *) tl_alloc(sizeof(int), smpi->npes);
    
    smpi->send_edge_key = (EDGE_LIST_ITEM ***) tl_alloc(sizeof(EDGE_LIST_ITEM **), smpi->npes);
    
    smpi->recv_edge_key = (EDGE_LIST_ITEM ***) tl_alloc(sizeof(EDGE_LIST_ITEM **), smpi->npes);
    
    for (i_processor = 0; i_processor < smpi->npes; i_processor++) {
        messg_buffer_init(smpi->send_msg + i_processor, i_processor);
        messg_buffer_init(smpi->recv_msg + i_processor, i_processor);
        messg_buffer_init(smpi->send_edge_msg + i_processor, i_processor);
        messg_buffer_init(smpi->recv_edge_msg + i_processor, i_processor);
        smpi->send_key[i_processor].size = 0;
        smpi->send_key[i_processor].key = NULL;
        smpi->recv_init[i_processor] = UNSET_INT;
        smpi->nsend[i_processor] = 0;
        smpi->nrecv[i_processor] = 0;
        smpi->send_key_surf[i_processor].size = 0;
        smpi->send_key_surf[i_processor].key = NULL;
        smpi->recv_init_surf[i_processor] = UNSET_INT;
        smpi->nsend_surf[i_processor] = 0;
        smpi->nrecv_surf[i_processor] = 0;
    }
    
    smpi->nmsg_counter = 0;
    smpi->nmsg_status = 0;
    smpi->nmsg_request = 0;
    
    smpi->max_nmsg_status = 1;
    smpi->msg_status = (MPI_Status *) tl_alloc(sizeof(MPI_Status), smpi->max_nmsg_status);
    smpi->max_nmsg_request = 0;
    smpi->partition_flag = 0;
    //smpi->partition_info = NULL; // CJT :: Was commented out in Lucas/Gajanan coupling branch
    //smpi->surface_partition_info = NULL; // CJT :: Was commented out in Lucas/Gajanan coupling branch
#else
    smpi->npes = 1;
    smpi->myid = 0;
#endif
    return;
}

/***********************************************************************************/
void smpi_free(SMPI *smpi) {
#ifdef _MESSG
    int i_processor = 0;          /* Loop Counter over processors */
    /* free the message arrays */
    
    if(smpi->recv_init != NULL) smpi->recv_init = (int *) tl_free(sizeof(int), smpi->npes, smpi->recv_init);
    
    if(smpi->nsend != NULL) smpi->nsend = (int *) tl_free(sizeof(int), smpi->npes, smpi->nsend);
    
    if(smpi->nrecv != NULL) smpi->nrecv = (int *) tl_free(sizeof(int), smpi->npes, smpi->nrecv);
    
    if(smpi->recv_init_surf != NULL) smpi->recv_init_surf = (int *) tl_free(sizeof(int), smpi->npes, smpi->recv_init_surf);
    
    if(smpi->nsend_surf != NULL) smpi->nsend_surf = (int *) tl_free(sizeof(int), smpi->npes, smpi->nsend_surf);
    
    if(smpi->nrecv_surf != NULL) smpi->nrecv_surf = (int *) tl_free(sizeof(int), smpi->npes, smpi->nrecv_surf);
    
    //smpi->send_msg = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), smpi->npes, smpi->send_msg);
    for (i_processor = 0; i_processor < smpi->npes;i_processor++){
        if(smpi->send_msg != NULL) if (smpi->send_key[i_processor].size > 0) messg_buffer_free(smpi->send_msg + i_processor);
        if(smpi->recv_edge_msg != NULL) if (smpi->recv_edge_msg[i_processor].size > 0) messg_buffer_free(smpi->recv_edge_msg + i_processor);
        if(smpi->send_edge_msg != NULL) if (smpi->send_edge_msg[i_processor].size > 0) messg_buffer_free(smpi->send_edge_msg + i_processor);
    }
    
    if(smpi->recv_msg != NULL) smpi->send_msg = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), smpi->npes, smpi->send_msg);
    
    if(smpi->recv_msg != NULL) smpi->recv_msg = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), smpi->npes, smpi->recv_msg);
    
    if(smpi->send_edge_msg != NULL) smpi->send_edge_msg = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), smpi->npes, smpi->send_edge_msg);
    
    if(smpi->recv_edge_msg != NULL) smpi->recv_edge_msg = (MESSG_BUFFER *) tl_free(sizeof(MESSG_BUFFER), smpi->npes, smpi->recv_edge_msg);
    
    for (i_processor = 0; i_processor < smpi->npes; i_processor++)
    {
        if(smpi->send_key != NULL){
            if (smpi->send_key[i_processor].key != NULL){
                smpi->send_key[i_processor].key = (int *) tl_free(sizeof(int), smpi->send_key[i_processor].size, smpi->send_key[i_processor].key);
            }
        }
        if(smpi->send_key_surf != NULL){
            if (smpi->send_key_surf[i_processor].key != NULL){
                smpi->send_key_surf[i_processor].key = (int *) tl_free(sizeof(int), smpi->send_key_surf[i_processor].size, smpi->send_key_surf[i_processor].key);
            }
        }
    }
    
    if(smpi->msg_request != NULL) smpi->msg_request = (MPI_Request *) tl_free(sizeof(MPI_Request), smpi->max_nmsg_request, smpi->msg_request);
    
    if(smpi->send_key != NULL) smpi->send_key = (MESSG_KEY *) tl_free(sizeof(MESSG_KEY), smpi->npes, smpi->send_key);
    
    if(smpi->send_key_surf != NULL) smpi->send_key_surf = (MESSG_KEY *) tl_free(sizeof(MESSG_KEY), smpi->npes, smpi->send_key_surf);
    
    if(smpi->nsend_edge != NULL) smpi->nsend_edge = (int *) tl_free(sizeof(int), smpi->npes, smpi->nsend_edge);
    
    if(smpi->nrecv_edge != NULL) smpi->nrecv_edge = (int *) tl_free(sizeof(int), smpi->npes, smpi->nrecv_edge);
    
    if(smpi->send_edge_key != NULL) smpi->send_edge_key = (EDGE_LIST_ITEM ***) tl_free(sizeof(EDGE_LIST_ITEM **), smpi->npes, smpi->send_edge_key);
    
    if(smpi->recv_edge_key != NULL) smpi->recv_edge_key = (EDGE_LIST_ITEM ***) tl_free(sizeof(EDGE_LIST_ITEM **), smpi->npes, smpi->recv_edge_key);
    
    if(smpi->msg_status != NULL) smpi->msg_status = (MPI_Status *) tl_free(sizeof(MPI_Status), smpi->max_nmsg_status, smpi->msg_status);
#endif
    return;
}

/*******************************************************************************************************************/

void smpi_defaults(SMPI *smpi){
    smpi->npes = 0;        /* number of processors */
    smpi->myid = -1;        /* local processor identification # */
#ifdef _MESSG
    smpi->send_msg = NULL;  /* the send message buffers */
    smpi->recv_msg = NULL;  /* the receive message buffers */
    smpi->send_edge_msg = NULL; /* send buffer for the edge messages in adaption */
    smpi->recv_edge_msg = NULL; /* receive buffer for the edge messages in adaption */
    smpi->send_key = NULL; /* the key to load the send message */
    smpi->send_key_surf = NULL; /* the key to load the send message */
    smpi->recv_init = NULL;      /* the beginning position to unload the message */
    smpi->recv_init_surf = NULL;      /* the beginning position to unload the message */
    smpi->nsend = NULL;      /* the number of items to send */
    smpi->nrecv = NULL;      /* the number of items to receive */
    smpi->nsend_surf = NULL;      /* the number of items to send */
    smpi->nrecv_surf = NULL;      /* the number of items to receive */
    smpi->nsend_edge = NULL;     /* the number of edges to send to the given pe */
    smpi->nrecv_edge = NULL;     /* the number of edges to receive from the given pe */
    smpi->send_edge_key = NULL; /* the key for the edges being sent to the given pe
                                 NOTE:  the key is actually a pointer into the hash table
                                 for the given edge */
    smpi->recv_edge_key = NULL; /* the key for the edges being received from the given pe
                                 NOTE:  the key is actually a pointer into the hash table
                                 for the given edge */
    smpi->msg_status = NULL;  /* return flag for MPI receives */
    smpi->msg_request = NULL;    /* request handle for asynchronous communication */
    
    smpi->partition_info = NULL; /* partition number for each node */
    smpi->surface_partition_info = NULL; /* 3d surface partition info */
#endif
    
    return;
}

/*****************************************************************************************************/

void smpi_realloc(SMPI *smpi
#ifdef _MESSG
                  , MPI_Comm input_comm
#endif
) {
    
#ifdef _MESSG
    
    int i_processor = 0;          /* Loop Counter over processors */
    int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */
    int old_npes=0;
    
    if(smpi->npes >=0) old_npes = smpi->npes;
    
    ierr_code = MPI_Comm_dup(input_comm, &(smpi->ADH_COMM));
    if (ierr_code != MPI_SUCCESS)
    {
        messg_err(ierr_code);
    }
    
    ierr_code = MPI_Comm_size(smpi->ADH_COMM, &(smpi->npes));
    
    ierr_code = MPI_Comm_rank(smpi->ADH_COMM, &(smpi->myid));
    
    smpi->msg_starttime = MPI_Wtime();
    
    /* allocates the message arrays */
    
    smpi->recv_init = (int *) tl_realloc(sizeof(int), smpi->npes, old_npes, smpi->recv_init);
    
    smpi->nsend = (int *) tl_realloc(sizeof(int), smpi->npes, old_npes, smpi->nsend);
    
    smpi->nrecv = (int *) tl_realloc(sizeof(int), smpi->npes, old_npes, smpi->nrecv);
    
    smpi->recv_init_surf = (int *) tl_realloc(sizeof(int), smpi->npes, old_npes, smpi->recv_init_surf);
    
    smpi->nsend_surf = (int *) tl_realloc(sizeof(int), smpi->npes, old_npes, smpi->nsend_surf);
    
    smpi->nrecv_surf = (int *) tl_realloc(sizeof(int), smpi->npes, old_npes, smpi->nrecv_surf);
    
    smpi->send_msg = (MESSG_BUFFER *) tl_realloc(sizeof(MESSG_BUFFER), smpi->npes, old_npes, smpi->send_msg);
    
    smpi->recv_msg = (MESSG_BUFFER *) tl_realloc(sizeof(MESSG_BUFFER), smpi->npes, old_npes, smpi->recv_msg);
    
    smpi->send_edge_msg = (MESSG_BUFFER *) tl_realloc(sizeof(MESSG_BUFFER), smpi->npes, old_npes, smpi->send_edge_msg);
    
    smpi->recv_edge_msg = (MESSG_BUFFER *) tl_realloc(sizeof(MESSG_BUFFER), smpi->npes, old_npes, smpi->recv_edge_msg);
    
    smpi->send_key = (MESSG_KEY *) tl_realloc(sizeof(MESSG_KEY), smpi->npes, old_npes, smpi->send_key);
    
    smpi->send_key_surf = (MESSG_KEY *) tl_realloc(sizeof(MESSG_KEY), smpi->npes, old_npes, smpi->send_key_surf);
    
    smpi->nsend_edge = (int *) tl_realloc(sizeof(int), smpi->npes, old_npes, smpi->nsend_edge);
    
    smpi->nrecv_edge = (int *) tl_realloc(sizeof(int), smpi->npes, old_npes, smpi->nrecv_edge);
    
    smpi->send_edge_key = (EDGE_LIST_ITEM ***) tl_realloc(sizeof(EDGE_LIST_ITEM **), smpi->npes, old_npes, smpi->send_edge_key);
    
    smpi->recv_edge_key = (EDGE_LIST_ITEM ***) tl_realloc(sizeof(EDGE_LIST_ITEM **), smpi->npes, old_npes, smpi->recv_edge_key);
    
    for (i_processor = 0; i_processor < smpi->npes; i_processor++) {
        messg_buffer_init(smpi->send_msg + i_processor, i_processor);
        messg_buffer_init(smpi->recv_msg + i_processor, i_processor);
        messg_buffer_init(smpi->send_edge_msg + i_processor, i_processor);
        messg_buffer_init(smpi->recv_edge_msg + i_processor, i_processor);
        smpi->send_key[i_processor].size = 0;
        smpi->send_key[i_processor].key = NULL;
        smpi->recv_init[i_processor] = UNSET_INT;
        smpi->nsend[i_processor] = 0;
        smpi->nrecv[i_processor] = 0;
        smpi->send_key_surf[i_processor].size = 0;
        smpi->send_key_surf[i_processor].key = NULL;
        smpi->recv_init_surf[i_processor] = UNSET_INT;
        smpi->nsend_surf[i_processor] = 0;
        smpi->nrecv_surf[i_processor] = 0;
    }
    
    smpi->nmsg_counter = 0;
    smpi->nmsg_status = 0;
    smpi->nmsg_request = 0;
    if(smpi->msg_status != NULL ){
        smpi->msg_status = (MPI_Status *) tl_realloc(sizeof(MPI_Status), 1, smpi->max_nmsg_status, smpi->msg_status);
    } else {
        smpi->msg_status = (MPI_Status *) tl_alloc(sizeof(MPI_Status), 1);
    }
    smpi->max_nmsg_status = 1;
    smpi->max_nmsg_request = 0;
    smpi->partition_flag = 0;
#else
    smpi->npes = 1;
    smpi->myid = 0;
#endif
    return;
}

