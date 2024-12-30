/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     hard coded for comm update double for the standalone test, will delete later
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 * \note The general routine should update ghost values on each process from the owners
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void comm_update_double(double *vec,  /* the vector to be updated */
                        int size_v,
                        int npe,
                        int rank //number of processors
  )
{
#ifdef _MESSG

  //double *bpntr = NULL;         /* pointer to allow for de-referencing */
  int i_processor = 0, j_processor = 0, ii=0;  /* loop counter */

  //just send whole vector i guess?
  double temp_send[7],temp_recv[7];
  for (ii = 0;ii<7;ii++){
    if(ii<size_v){
      temp_send[ii] = vec[ii];
      temp_recv[ii] = 0;
    }else{temp_send[ii] = 0;
          temp_recv[ii] = 0;}
  }


  /* Post the Receives */
  for (i_processor = 0; i_processor < npe; i_processor++){
    if (rank ==i_processor){
      for (j_processor = 0; j_processor < npe; j_processor++){
        if (rank != j_processor){
          MPI_Send(temp_send, 7, MPI_DOUBLE, j_processor, j_processor, MPI_COMM_WORLD);
        }
      }
    }else{
       MPI_Recv(temp_recv, 7, MPI_DOUBLE, i_processor, rank, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
       //now move temp entries to appropriate vec!!!
       if (rank==0){
          if(i_processor == 1){
            vec[3] = temp_recv[0];
            vec[4] = temp_recv[1];
          }else if(i_processor == 2){
            vec[5] = temp_recv[0];
            vec[6] = temp_recv[1];

          }
       }else if(rank==1){
          if(i_processor == 0){
            vec[3] = temp_recv[0];
            vec[4] = temp_recv[1];
            vec[5] = temp_recv[2];
          }else if(i_processor == 2){
            vec[6] = temp_recv[0];

          }
       }else if(rank==2){
        if(i_processor == 0){
            vec[2] = temp_recv[0];
            vec[3] = temp_recv[1];
            vec[4] = temp_recv[2];
          }else if(i_processor == 1){
            vec[5] = temp_recv[0];
            vec[6] = temp_recv[1];
            vec[7] = temp_recv[2];

          }

       }
    }
  }
#endif
  return;
}
