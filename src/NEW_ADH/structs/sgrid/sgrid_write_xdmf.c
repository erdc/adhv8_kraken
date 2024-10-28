#include "adh.h"
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;
#define SPATIAL_DIM 3
#define MATRIX_RANK 2
#define VECTOR_RANK 1
void sgrid_write_xdmf(SGRID *g){

#ifdef _MESSG
    if(g->smpi->myid==0)
#endif
    {
        FILE *xmf = 0;
        char fname[50];
        float t=0;
        int mesh_no=0;
        strcpy(fname, g->filename);
        strcat(fname, ".xmf");
        xmf = fopen(fname, "w");

        //call from xdmf utils
        write_xdmf_header(xmf);

        fprintf(xmf, "\t\t\t<Grid Name=\"2D Unstructured Mesh\">\n");
        fprintf(xmf, "\t\t\t<Time Value=\"%f\" />\n",t);
        //Geometry start, this is nodes
        fprintf(xmf, "\t\t\t\t<Geometry GeometryType=\"XYZ\">\n");
        fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d %d\" Format=\"HDF\">\n", g->macro_nnodes,SPATIAL_DIM);
        fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Mesh/XY/%d\n",g->filename,mesh_no);
        fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
        fprintf(xmf, "\t\t\t\t</Geometry>\n");
        //Geometry finish, this is nodes

        //Topology start
        fprintf(xmf, "\t\t\t\t<Topology TopologyType=\"Mixed\" NumberOfElements=\"%d\">\n", g->macro_nelems1d+ g->macro_nelems2d+g->macro_nelems3d);
        //Data Item is the mixed connectivity
        fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", g->macro_nQuads * 5 + g->macro_nTris * 4 + g->macro_nTets * 5 + g->macro_nPrisms * 7);
        fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Mesh/Elements/%d\n",g->filename,mesh_no);
        fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
        fprintf(xmf, "\t\t\t\t</Topology>\n");
        //Topology finish, this contains connectivities

        //call from xdmf utils
        write_xdmf_tail(xmf);
        //XDMF File finished

        fclose(xmf);

    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Print an AdH grid to XDMF
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


void sgrid_write_xdmf_nodal_pe(SGRID *g){

if(g->smpi->myid==0){
    FILE *xmf = 0;
    char fname[50];
    float t=0;
    int mesh_no=0;
    strcpy(fname, g->filename);
    strcat(fname, ".xmf");
    xmf = fopen(fname, "r+");


  int lines = 0;
  int nchar = 0,nline=0;
  char line[256],c;

  //get line count first
  // Extract characters from file and store in character c
  for (c = getc(xmf); c != EOF; c = getc(xmf))
        if (c == '\n') // Increment count if this character is newline
            nline+= 1;
  printf("Nline,%d\n",nline);
  rewind(xmf);
  // Count lines and move to the second last line (considering 0-based indexing)
  while (fgets(line, sizeof(line), xmf) != NULL) {
    lines++;
    if(lines==nline-4){break;}
    }
    //get the pointer
    fseek(xmf, 0, SEEK_CUR);

  // Insert text before the current position
  //Any nodal data we can also write
            fprintf(xmf, "\t\t\t\t<Attribute Name=\"NodalData\" AttributeType=\"Scalar\" Center=\"Node\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", g->macro_nnodes);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/NodalScalar/%d\n",g->filename,mesh_no);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Attribute>\n");

  // append file tail
  write_xdmf_tail(xmf);



    //XDMF File finished

    fclose(xmf);
    printf("Text inserted successfully!\n");

}


}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Print an AdH grid to XDMF
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


void sgrid_write_xdmf_elemental_pe(SGRID *g){

if(g->smpi->myid==0){
    FILE *xmf = 0;
    char fname[50];
    float t=0;
    int mesh_no=0;
    strcpy(fname, g->filename);
    strcat(fname, ".xmf");
    xmf = fopen(fname, "r+");


  int lines = 0;
  int nchar = 0,nline=0;
  char line[256],c;

  //get line count first
  // Extract characters from file and store in character c
  for (c = getc(xmf); c != EOF; c = getc(xmf))
        if (c == '\n') // Increment count if this character is newline
            nline+= 1;
  printf("Nline,%d\n",nline);
  rewind(xmf);
  // Count lines and move to the second last line (considering 0-based indexing)
  while (fgets(line, sizeof(line), xmf) != NULL) {
    lines++;
    if(lines==nline-4){break;}
    }
    //get the pointer
    fseek(xmf, 0, SEEK_CUR);


  // Insert text before the current position
  //Any nodal data we can also write
  //Should also be way to write elemental data here too just change Center=\"Node" to Center=\"Cell"
            fprintf(xmf, "\t\t\t\t<Attribute Name=\"ElementalData\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", g->macro_nelems1d+ g->macro_nelems2d+g->macro_nelems3d);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/ElementalScalar/%d\n",g->filename,mesh_no);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Attribute>\n");

  // append file tail
  write_xdmf_tail(xmf);

    //XDMF File finished

    fclose(xmf);


}


}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Print an AdH grid to XDMF
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


void sgrid_write_nodal_pe(SGRID *g){



#ifdef _ADH_HDF5

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

    // mesh properties
    int *connectivity;

    // create object to store h5 config
    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    #ifdef _MESSG
        MPI_Info info  = MPI_INFO_NULL;
        H5Pset_fapl_mpio(plist_id, g->smpi->ADH_COMM, info);
    #endif


    // open file
    strcpy(fname,g->filename);
    strcat(fname, ".h5");
    file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    // group 1 is anything with mesh nodes
    grp1 = H5Gopen(file_id, "/Data/NodalScalar", H5P_DEFAULT); // a group is like a folder

    //++++++++++++++++++++++++++++++++++++++
    // Writing Nodes
    //++++++++++++++++++++++++++++++++++++++

    // declare global sizes of data set nodes first. This is local data.
    //Nodally based scalar data, maybe outsource this to other routine later

    int *scalardata;
    scalardata = (int *)malloc(sizeof(int) * g->my_nnodes);




    // create a dataspace.  This is a global quantity
    dims[0]  = g->macro_nnodes;
    dims[1]  = 0;
    filespace = H5Screate_simple(VECTOR_RANK, dims, NULL);

    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp1, number, H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

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
    int i,j,k=0,l=0;
    for (i=0; i<g->nnodes; i++){
        //temporarily store PE
        if(g->node[i].resident_pe == g->smpi->myid){
            scalardata[l] = g->node[i].resident_pe;
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
    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, scalardata);
    free(scalardata);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp1);

    // +++++++++++++++++++++++++++++++++++++++
    // Nodal write complete
    // +++++++++++++++++++++++++++++++++++++++

    //close mesh file
    H5Fclose(file_id);

#endif

}




/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Print an AdH grid to XDMF
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


void sgrid_write_elemental_pe(SGRID *g){



#ifdef _ADH_HDF5

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

    // mesh properties
    int *connectivity;

    // create object to store h5 config
    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    #ifdef _MESSG
        MPI_Info info  = MPI_INFO_NULL;
        H5Pset_fapl_mpio(plist_id, g->smpi->ADH_COMM, info);
    #endif


    // open file
    strcpy(fname,g->filename);
    strcat(fname, ".h5");
    file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    // group 1 is anything with mesh nodes
    grp1 = H5Gopen(file_id, "/Data/ElementalScalar", H5P_DEFAULT); // a group is like a folder

    //++++++++++++++++++++++++++++++++++++++
    // Writing Nodes
    //++++++++++++++++++++++++++++++++++++++

    // declare global sizes of data set nodes first. This is local data.
    //Nodally based scalar data, maybe outsource this to other routine later

    int *scalardata;
    scalardata = (int *)malloc(sizeof(int) * (g->my_nelems3d + g->my_nelems2d + g->my_nelems1d));




    // create a dataspace.  This is a global quantity
    dims[0]  = g->macro_nelems3d + g->macro_nelems2d + g->macro_nelems1d;
    dims[1]  = 0;
    filespace = H5Screate_simple(VECTOR_RANK, dims, NULL);

    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp1, number, H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    //create a local dataspace where each process has a few of the rows
    count[0]  = g->my_nelems3d + g->my_nelems2d + g->my_nelems1d;
    count[1]  = dims[1];


    // give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(VECTOR_RANK, count, NULL);

    // determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);

    // the 2 is dimension of our data structure
    hsize_t   coord1[g->my_nelems3d + g->my_nelems2d + g->my_nelems1d];

    // the values in coord1 should be the global node numbers, this time it is trivial
    int i,j,k=0,l=0;
    for (i=0; i<g->nelems3d; i++){
        //temporarily store PE
        if(g->elem3d[i].resident_pe == g->smpi->myid){
            scalardata[l] = g->elem3d[i].resident_pe;
            coord1[l] = g->elem3d[i].gid;
            l+=1;
    }
    }

    for (i=0; i<g->nelems2d; i++){
        //temporarily store PE
        if(g->elem2d[i].resident_pe == g->smpi->myid){
            scalardata[l] = g->elem2d[i].resident_pe;
            coord1[l] = g->elem2d[i].gid + g->macro_nelems3d;
            l+=1;
    }
    }

    //currently 1d not available to have PE plotted

    assert(l==g->my_nelems3d + g->my_nelems2d + g->my_nelems1d);
    H5Sselect_elements(filespace, H5S_SELECT_SET,  g->my_nelems3d + g->my_nelems2d + g->my_nelems1d  ,(const hsize_t *)&coord1);

    // Create property list for collective dataset write.
    // declare collextive data file writing

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #ifdef _MESSG
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, scalardata);
    free(scalardata);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp1);

    // +++++++++++++++++++++++++++++++++++++++
    // Nodal write complete
    // +++++++++++++++++++++++++++++++++++++++

    //close mesh file
    H5Fclose(file_id);

#endif

}

