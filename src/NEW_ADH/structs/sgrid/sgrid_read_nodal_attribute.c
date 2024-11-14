#include "adh.h"
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;
#define MATRIX_RANK 2
#define VECTOR_RANK 1
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Read a nodal attribute
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] g           (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void sgrid_read_nodal_attribute(SGRID *g){



#ifdef _ADH_HDF5

    int i,j=0;
    int nsurf=0;
    int lnode3d;



    for (i=0;i<g->nnodes;i++){
        if (g->node[i].bflag==SURFACE || g->node[i].bflag==BODY2D){
            //includes ghost nodes
            //printf("%d\n",g->node[i].global_surf_id);
            nsurf+=1;
    }
    }
    //Note, need these to match to avoid above loop
    printf("PE %d, nnodes_sur = %d, garbage %d\n",g->smpi->myid,nsurf,g->nnodes_sur_old);
    int mapping[nsurf];
    int nodeID_2d_to_3d_sur[nsurf];
    //set up map to go from local to global
    for (i=0;i<g->nnodes;i++){
        if (g->node[i].bflag==SURFACE || g->node[i].bflag==BODY2D){
            nodeID_2d_to_3d_sur[j] = i;
            mapping[j]=g->node[i].global_surf_id;
            j+=1;
        }
    }
    //key attribute here,
    //g->node[local_index].global_surf_id
    //localSurf2GlobalMap[nnodes_surf]; initialize to UNSET_INT
    //count = 0
    //for (i=0; i<nnodes_sur; i++) {
     //lnode3d = nodeID_2d_to_3d_sur[i]
     //gnode3d = g->node[lnode3d].gid
    //localSurf2GlobalMap[count] = gnode3d;
     //count ++;
    //}


    //attempt to read in H5 file in parallel
    hid_t     file_id, plist_id, filespace, dset_id, memspace;
    hid_t     grp1,grp2;
    char      fname[50];

    hsize_t   dims[MATRIX_RANK];
    herr_t    status;
    hsize_t   count[MATRIX_RANK],offset[MATRIX_RANK];
    int nt=0;
    // for adaptive grid writing
    char number[50] = "";
    sprintf(number, "%d", nt);

    // create object to store h5 config
    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    #ifdef _MESSG
        MPI_Info info  = MPI_INFO_NULL;
        H5Pset_fapl_mpio(plist_id, g->smpi->ADH_COMM, info);
    #endif


    // open file

    strcpy(fname,g->filename);
    strcat(fname, ".h5");
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);



    //++++++++++++++++++++++++++++++++++++++
    // Writing Nodes
    //++++++++++++++++++++++++++++++++++++++

    // declare global sizes of data set nodes first. This is local data.
    //Nodally based scalar data, maybe outsource this to other routine later

    int *scalardata;
    scalardata = (int *)malloc(sizeof(int) * g->my_nnodes);



    //create a parallel dataset object and close filespace
    dset_id = H5Dopen(file_id, "/Data/NodalScalar/0", H5P_DEFAULT);


    //create a local dataspace where each process has a few of the rows
    count[0]  = g->my_nnodes;
    count[1]  = dims[1];


    // give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(VECTOR_RANK, count, NULL);

    // determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);

    // the 2 is dimension of our data structure
    hsize_t   coord1[g->my_nnodes];

    // the values in coord1 should be the global node numbers, this time it is trivial
    int l=0;
    for (i=0; i<g->nnodes; i++){
        //temporarily store PE
        if(g->node[i].resident_pe == g->smpi->myid){
            //scalardata[l] = g->node[i].resident_pe;
            coord1[l] = g->node[i].gid;
            l+=1;

    }
    }

    assert(l==g->my_nnodes);
    H5Sselect_elements(filespace, H5S_SELECT_SET,g->my_nnodes,(const hsize_t *)&coord1);

    // Create property list for collective dataset write.
    // declare collextive data file writing

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #ifdef _MESSG
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    //collective write
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, scalardata);



    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);

    // +++++++++++++++++++++++++++++++++++++++
    // Nodal write complete
    // +++++++++++++++++++++++++++++++++++++++

    //close mesh file
    H5Fclose(file_id);


    for (i=0; i<g->nnodes; i++){
        if(g->node[i].resident_pe == g->smpi->myid){
            printf("PE %d, gid %d, scalardat = %d\n",g->smpi->myid,g->node[i].gid,scalardata[i]);
        }
    }
    free(scalardata);

#endif

}

