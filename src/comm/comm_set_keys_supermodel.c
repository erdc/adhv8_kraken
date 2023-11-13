#include "global_header.h"

/*! 
  \brief This routine initializes the communications data for the supermodel
  */
int sort_key(const void *p, const void *q){
  int key1 = *(const int *)p;
  int key2 = *(const int *)q;
  if(key1<key2){
    return(-1);
  }else if(key2<key1){
    return(1);
  }else{
    return(0);
  }
}
static int DEBUG = OFF;

void comm_set_keys_supermodel(SSUPER_MODEL *sm)
{
  SMPI *smpi = sm->supersmpi;   /* alias */
  int **fmap = sm->fmap;

  int i, j, k, sub_pe, check, inode;
  int ii = 0, jj = 0;           /* loop counters */
  int ie = 0;                   /* loop counter over the elements */
  int i_processor = 0;          /* loop counter over the processors */
  int *is_border = NULL;        /* flags for border nodes */
  int local_owner = 0;          /* check */
  SINTERFACE *ifce;
  SINTERFACE_NODELIST *ndlist;
  int master_mod, slave_mod, mod_id;
  int first_non_ifce, recv_size;
  int first_non_ifce_surf, recv_size_surf;
  first_non_ifce=UNSET_INT;
  first_non_ifce_surf=UNSET_INT;
  int *already_recorded;
  int *temp_key;

  already_recorded = (int *) tl_alloc(sizeof(int), sm->nnodes);



  /*! Clear the Communication Buffers, etc. */
  for (i_processor = 0; i_processor < smpi->npes; i_processor++){

    messg_buffer_free(smpi->send_msg + i_processor); 
    messg_buffer_init(smpi->send_msg + i_processor, i_processor);
    messg_buffer_free(smpi->recv_msg + i_processor);      
    messg_buffer_init(smpi->recv_msg + i_processor, i_processor);

    if (smpi->send_key[i_processor].key != NULL){
      smpi->send_key[i_processor].key = (int *) tl_free(sizeof(int), smpi->send_key[i_processor].size, smpi->send_key[i_processor].key);

    }
    smpi->send_key[i_processor].size = 0;
    smpi->recv_init[i_processor] = UNSET_INT;
    smpi->nsend[i_processor] = 0;
    smpi->nrecv[i_processor] = 0;

    if (smpi->send_key_surf[i_processor].key != NULL){
      smpi->send_key_surf[i_processor].key = (int *) tl_free(sizeof(int), smpi->send_key_surf[i_processor].size, smpi->send_key_surf[i_processor].key);

    }
    smpi->send_key_surf[i_processor].size = 0;
    smpi->recv_init_surf[i_processor] = UNSET_INT;
    smpi->nsend_surf[i_processor] = 0;
    smpi->nrecv_surf[i_processor] = 0;

  }

  /*Determine if I own any interface nodes */
  check=0;
  for (jj=0; jj<sm->NumInterfaces; jj++){
    ifce = &(sm->interface[jj]);
    if((sm->submodel[ifce->model_id[0]].proc_flag>0)||(sm->submodel[ifce->model_id[1]].proc_flag>0)){
      for (j=0; j<ifce->NumNodeColumns; j++){
        ndlist = &(ifce->nodelist[j]);
        for(i=0;i<2;i++){
          if(sm->submodel[ifce->model_id[i]].proc_flag>0){
            for(inode=0;inode<sm->submodel[ifce->model_id[i]].grid->my_nnodes;inode++){
              if(sm->submodel[ifce->model_id[i]].grid->node[inode].gid==ndlist->surfnode[i]){ 
                check=1;
                break;
              }
            }
          }if(check>0)break;
        }if(check>0)break;
      }if(check>0)break;
    }
    if(check>0)break;
  }
  int interface_check[smpi->npes]; /* array to determine if processors own interface nodes */

  /* see which procs own interface nodes
   * if interface nodes are owned, we go through the complicated key setup
   * if not, then key set up is just model keys fmapped */
  for(i=0;i<smpi->npes;i++) interface_check[i]=0;

  MPI_Allgather(&(check),1,MPI_INT,interface_check,1,MPI_INT,smpi->ADH_COMM);


  check=0;
  sarray_init_value_int(already_recorded, sm->nnodes, UNSET_INT);
  
  /* Setting up variables for recieving data */
  for (i_processor=0; i_processor<smpi->npes; i_processor++){ 
    for(i=0;i<sm->nsubmodels;i++){
      sub_pe = sm->proc_map[i_processor][i];
      if((sm->submodel[i].proc_flag>0)&&(sub_pe!=UNSET_INT)){
        if(sm->submodel[i].grid->smpi->nrecv[sub_pe]>0){
          for (ii = sm->submodel[i].grid->smpi->recv_init[sub_pe]; ii < (sm->submodel[i].grid->smpi->recv_init[sub_pe] + sm->submodel[i].grid->smpi->nrecv[sub_pe]); ii++){
            j = fmap[i][ii]; 
            if(already_recorded[j] == UNSET_INT&&(interface_check[smpi->myid]==1)){
              smpi->nrecv[i_processor]++;
              already_recorded[j]=1;
            }else{
              smpi->nrecv[i_processor]++;
              already_recorded[j]=1;
            }
            if(((smpi->recv_init[i_processor]==UNSET_INT)||(j<smpi->recv_init[i_processor]))&&(already_recorded[j]==1)){
              smpi->recv_init[i_processor]=j;
            }
          }
        }
      }
    }
  }

  sarray_init_value_int(already_recorded, sm->nnodes, UNSET_INT);
    

  for (i_processor=0; i_processor<smpi->npes; i_processor++){
    for(i=0;i<sm->nsubmodels;i++){
      sub_pe = sm->proc_map[i_processor][i];
      if((sm->submodel[i].proc_flag>0)&&(sub_pe!=UNSET_INT)){
        if(sm->submodel[i].grid->smpi->nrecv_surf[sub_pe]>0){
          for (ii = sm->submodel[i].grid->smpi->recv_init_surf[sub_pe]; ii < (sm->submodel[i].grid->smpi->recv_init_surf[sub_pe] + sm->submodel[i].grid->smpi->nrecv_surf[sub_pe]); ii++){
            j = fmap[i][ii];
            if(already_recorded[j] == UNSET_INT&&(interface_check[smpi->myid]==1)){
              smpi->nrecv_surf[i_processor]++;
              already_recorded[j]=1;
            }else{
              smpi->nrecv_surf[i_processor]++;
              already_recorded[j]=1;
            }
            if(((smpi->recv_init_surf[i_processor]==UNSET_INT)||(j<smpi->recv_init_surf[i_processor]))&&(already_recorded[j]==1)){
              smpi->recv_init_surf[i_processor]=j;
            }
          }
        }
      }
    }
  }
    

  sarray_init_value_int(already_recorded, sm->nnodes, UNSET_INT);

  /* Setting up variables for sending data. */
  int new_size, new_size_surf;
  int cum_keysize, cum_keysurfsize;
  int iii;

  for (i_processor=0; i_processor<smpi->npes; i_processor++){
    sarray_init_value_int(already_recorded, sm->nnodes, UNSET_INT);
    smpi->nsend[i_processor] = 0;
    smpi->send_key[i_processor].size = 0;
    for(i=0;i<sm->nsubmodels;i++){
      sub_pe = sm->proc_map[i_processor][i];
      if((sm->submodel[i].proc_flag>0)&&(sub_pe!=UNSET_INT)){
        if(sm->submodel[i].grid->smpi->nsend[sub_pe]>0){
          for (ii = 0; ii < sm->submodel[i].grid->smpi->send_key[sub_pe].size; ii++){
            j = fmap[i][sm->submodel[i].grid->smpi->send_key[sub_pe].key[ii]];     
            if((already_recorded[j] == UNSET_INT)&&(interface_check[i_processor]==1)){
              smpi->nsend[i_processor]++;
              smpi->send_key[i_processor].size++;
              already_recorded[j]=1;
            }else if(interface_check[i_processor]==0){
              smpi->nsend[i_processor]++;
              smpi->send_key[i_processor].size++;
            }
          }
        }
      }
    }
    if (smpi->send_key[i_processor].size > 0){
      smpi->send_key[i_processor].key = (int *) tl_alloc(sizeof(int), smpi->send_key[i_processor].size);
    }
  }

  sarray_init_value_int(already_recorded, sm->nnodes, UNSET_INT);

  for (i_processor=0; i_processor<smpi->npes; i_processor++){
    sarray_init_value_int(already_recorded, sm->nnodes, UNSET_INT);
    smpi->nsend_surf[i_processor] = 0;
    for(i=0;i<sm->nsubmodels;i++){
      sub_pe = sm->proc_map[i_processor][i];
      if((sm->submodel[i].proc_flag>0)&&(sub_pe!=UNSET_INT)){
        if(sm->submodel[i].grid->smpi->nsend_surf[sub_pe]>0){
          for (ii = 0; ii < sm->submodel[i].grid->smpi->send_key_surf[sub_pe].size; ii++){
            j = fmap[i][sm->submodel[i].grid->smpi->send_key_surf[sub_pe].key[ii]];
            if((already_recorded[j] == UNSET_INT)&&(interface_check[i_processor]==1)){
              smpi->nsend_surf[i_processor]++;
              smpi->send_key_surf[i_processor].size++;
              already_recorded[j]=1;
            }else if(interface_check[i_processor]==0){
              smpi->nsend_surf[i_processor]++;
              smpi->send_key_surf[i_processor].size++;
            }
          }
        }
      }
    }
    if (smpi->send_key_surf[i_processor].size > 0){
      smpi->send_key_surf[i_processor].key = (int *) tl_alloc(sizeof(int), smpi->send_key_surf[i_processor].size);
    }
  }

  sarray_init_value_int(already_recorded, sm->nnodes, UNSET_INT);
  
  for (i_processor=0; i_processor<smpi->npes; i_processor++){
    sarray_init_value_int(already_recorded, sm->nnodes, UNSET_INT);
    new_size = 0;
    if(smpi->send_key[i_processor].size >0) temp_key = (int *) tl_alloc(sizeof(int),smpi->send_key[i_processor].size);
    for(i=0;i<sm->nsubmodels;i++){
      sub_pe = sm->proc_map[i_processor][i];
      if((sm->submodel[i].proc_flag>0)&&(sub_pe!=UNSET_INT)){
        if(sm->submodel[i].grid->smpi->nsend[sub_pe]>0){
          for (ii = 0; ii < sm->submodel[i].grid->smpi->send_key[sub_pe].size; ii++){
            j = fmap[i][sm->submodel[i].grid->smpi->send_key[sub_pe].key[ii]];
            if((already_recorded[j] == UNSET_INT)&&(interface_check[i_processor]==1)){
              temp_key[new_size++] = j;
              already_recorded[j]=1;
            }else if(interface_check[i_processor]==0){
              temp_key[new_size++] = j;
            }
          }
        }
      }
    }
    qsort((void *) temp_key, smpi->send_key[i_processor].size, sizeof(int), sort_key);
    for(i=0;i<smpi->send_key[i_processor].size;i++){
      smpi->send_key[i_processor].key[i] = temp_key[i];
    }
    if(smpi->send_key[i_processor].size >0) temp_key = (int *) tl_free(sizeof(int),smpi->send_key[i_processor].size, temp_key);
    
  }
  
  
  sarray_init_value_int(already_recorded, sm->nnodes, UNSET_INT);
  
  for (i_processor=0; i_processor<smpi->npes; i_processor++){
    sarray_init_value_int(already_recorded, sm->nnodes, UNSET_INT);
    new_size = 0;
    if(smpi->send_key_surf[i_processor].size >0)temp_key = (int *) tl_alloc(sizeof(int),smpi->send_key_surf[i_processor].size);
    for(i=0;i<sm->nsubmodels;i++){
      sub_pe = sm->proc_map[i_processor][i];
      if((sm->submodel[i].proc_flag>0)&&(sub_pe!=UNSET_INT)){
        if(sm->submodel[i].grid->smpi->nsend_surf[sub_pe]>0){
          for (ii = 0; ii < sm->submodel[i].grid->smpi->send_key_surf[sub_pe].size; ii++){
            j = fmap[i][sm->submodel[i].grid->smpi->send_key_surf[sub_pe].key[ii]];
            if((already_recorded[j] == UNSET_INT)&&(interface_check[i_processor]==1)){
              temp_key[new_size++] = j;
              already_recorded[j]=1;
            }else if(interface_check[i_processor]==0){
              temp_key[new_size++] = j;
            }
          }
        }
      }
    }
    qsort((void *) temp_key, smpi->send_key_surf[i_processor].size, sizeof(int), sort_key);
    for(i=0;i<smpi->send_key_surf[i_processor].size;i++){
      smpi->send_key_surf[i_processor].key[i] = temp_key[i];
    }
    if(smpi->send_key_surf[i_processor].size >0)temp_key = (int *) tl_free(sizeof(int),smpi->send_key_surf[i_processor].size, temp_key);
    
  }
  already_recorded = (int *) tl_free(sizeof(int), sm->nnodes, already_recorded); 

  return;
}
