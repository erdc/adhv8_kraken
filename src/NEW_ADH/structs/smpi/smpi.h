#ifdef _MESSG
#include "mpi.h"
#include <time.h>
#endif
#ifndef H_SMPI_
#define H_SMPI_

/*!
   Every Send/Receive uses a "message buffer" structure that stores info about
   the send/receive.

   size - Allocated Size in Bytes of *buffer.
   It should only be non-zero if *buffer is specifically allocated with its
   "own" memory. For instance, in a send, it may be necessary
   to gather pieces of non-contiguous information. In that case the *buffer
   is allocated its own memory, pieces are packed into it, in some order, and then
   a contiguous chunk of information is ready to send. If *buffer is acting as
   a simple pointer to (other) contiguous memory, then size needs to be zero. Any
   previously allocated memory to *buffer should be freed first.

   nitem - Number of items being sent/received. Needed by the MPI call.

   type - Type of data being sent (MPI_INT, MPI_DOUBLE, etc.)

   sd - The "other" processor with which processor myid is communicating with

   pos - An offset within *buffer, used for packing/unpacking mixed-type messages
   where (MPI_PACK) is called to pack the buffer.

   *buffer - Is the buffer that will be sent/received. It is a (void *) because
   many different types of information may be sent (int, double, packed/mixed, etc.)
   It may be its own allocated space, or it may be a pointer to some other
   (separately allocated memory).
 */


typedef struct {
  int size;         /* the size of the buffer - in bytes */
  int nitem;            /* the number of items in the buffer */
  int type;         /* the type of buffer */
  int sd;           /* the processor the message is being send to */
  int pos;          /* current position in the buffer */
  void *buffer;         /* the buffer */
} MESSG_BUFFER;         /* message buffer */

//************************************************************//
//************************************************************//

typedef struct {
  int size;         /* the size of the key */
  int *key;         /* the key */
} MESSG_KEY;            /* message key */

//************************************************************//
//************************************************************//

typedef struct {
  int rnode;            /* the relative node number - local node number on owning processor */
  int sd;           /* the owning processor */
  int global_num;       /* global node number */
} GLOBAL_NODE;          /* the global node number */

/***********************************************************/
/* struct methods ---------------------------------------- */

/* Message Passing Data */
typedef struct {
  int npes;        /* number of processors */
  int myid;        /* local processor identification # */
#ifdef _MESSG
  double msg_starttime;    /* the start time of the model */
  time_t msg_timestart;    /* the start time of the model */
  MESSG_BUFFER *send_msg;  /* the send message buffers */
  MESSG_BUFFER *recv_msg;  /* the receive message buffers */
  MESSG_BUFFER *send_edge_msg; /* send buffer for the edge messages in adaption */
  MESSG_BUFFER *recv_edge_msg; /* receive buffer for the edge messages in adaption */
  MESSG_KEY *send_key; /* the key to load the send message */
  MESSG_KEY *send_key_surf; /* the key to load the send message */
  int *recv_init;      /* the beginning position to unload the message */
  int *recv_init_surf;      /* the beginning position to unload the message */
  int *nsend;      /* the number of items to send */
  int *nrecv;      /* the number of items to receive */
  int *nsend_surf;      /* the number of items to send */
  int *nrecv_surf;      /* the number of items to receive */
  int *nsend_edge;     /* the number of edges to send to the given pe */
  int *nrecv_edge;     /* the number of edges to receive from the given pe */
  EDGE_LIST_ITEM ***send_edge_key; /* the key for the edges being sent to the given pe
                       NOTE:  the key is actually a pointer into the hash table
                       for the given edge */
  EDGE_LIST_ITEM ***recv_edge_key; /* the key for the edges being received from the given pe
                       NOTE:  the key is actually a pointer into the hash table
                       for the given edge */
  int nmsg_counter;    /* the number of outstanding messages */
  int nmsg_status;     /* the number of status structures allocated */
  int nmsg_request;    /* the number of request structures allocated */
  int max_nmsg_status; /* the number of status structures allocated */
  int max_nmsg_request;    /* the number of request structures allocated */

  MPI_Comm ADH_COMM;   /* Alias for actual communicator used */
  MPI_Status *msg_status;  /* return flag for MPI receives */
  MPI_Request *msg_request;    /* request handle for asynchronous communication */

  int *partition_info; /* partition number for each node */
  int *surface_partition_info; /* 3d surface partition info */
  int partition_flag; /* flags if the grid has been previously partitioned */
#endif
} SMPI;

#endif

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void smpi_init(SMPI *smpi
#ifdef _MESSG
               , MPI_Comm input_comm
#endif
);
void smpi_free(SMPI *smpi);
void smpi_defaults(SMPI *smpi);
void smpi_realloc(SMPI *smpi
#ifdef _MESSG
                  , MPI_Comm input_comm
#endif
);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
