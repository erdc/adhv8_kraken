#include "adh.h"
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;
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

#define MATRIX_RANK 2 // Tensor dimension :: Always be 2
#define VECTOR_RANK 1 // Always be 1
#define SPATIAL_DIM 3 // Always be 3



void sgrid_write_hdf5(SGRID *g){

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
    grp1 = H5Gopen(file_id, "/Mesh/XY", H5P_DEFAULT); // a group is like a folder

    //++++++++++++++++++++++++++++++++++++++
    // Writing Nodes
    //++++++++++++++++++++++++++++++++++++++

    // declare global sizes of data set nodes first. This is local data.
    double (*xyz)[SPATIAL_DIM] = malloc(sizeof(double[g->my_nnodes][SPATIAL_DIM]));




    // create a dataspace.  This is a global quantity
    dims[0]  = g->macro_nnodes;
    dims[1]  = SPATIAL_DIM;
    filespace = H5Screate_simple(MATRIX_RANK, dims, NULL);

    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp1, number, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    //create a local dataspace where each process has a few of the rows
    count[0]  = g->my_nnodes;
    count[1]  = dims[1];
    //offset[0] = g->node[0].gid;
    //offset[1] = 0;

    // give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(MATRIX_RANK, count, NULL);

    // determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    
    // the 2 is dimension of our data structure, RANK is spatial dimension of mesh
    hsize_t   coord1[g->my_nnodes * SPATIAL_DIM * MATRIX_RANK];
    
    // the values in coord1 should be the global node numbers, this time it is trivial
    int i,j,k=0,l=0;
    for (i=0; i<g->nnodes; i++){
        //temporarily store xyz if it is a residential node
#ifdef _MESSG
        if(g->node[i].resident_pe == g->smpi->myid)
#endif
        {
            xyz[l][0] = g->node[i].x; 
            xyz[l][1] = g->node[i].y;
            xyz[l][2] = g->node[i].z;
            l+=1;

            for(j=0;j<SPATIAL_DIM;j++) {
                coord1[k] = g->node[i].gid;
                coord1[k+1] = j;
                k += MATRIX_RANK;
                
            }

        }
    }
    

    assert(l==g->my_nnodes);
    H5Sselect_elements(filespace, H5S_SELECT_SET,g->my_nnodes * SPATIAL_DIM,(const hsize_t *)&coord1);
    
    // Create property list for collective dataset write.
    // declare collextive data file writing

    plist_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef _MESSG
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif
    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xyz);
    free(xyz);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp1);
    // +++++++++++++++++++++++++++++++++++++++
    // Nodal write complete
    // +++++++++++++++++++++++++++++++++++++++


    // +++++++++++++++++++++++++++++++++++++++
    // Now write elemental data
    // +++++++++++++++++++++++++++++++++++++++

    // Connectivity data
    //size of array should be numQuad*(4+1) + numTri*(3+1) + numTet*(4+1) + numPrism*(6+1)
    //+1 is for element code
    grp2 = H5Gopen(file_id, "/Mesh/Elements", H5P_DEFAULT);

    //global attributes
    int nentry = g->macro_nQuads * 5 + g->macro_nTris * 4 + g->macro_nTets * 5 + g->macro_nPrisms * 7;

    //local data, different on each PE
    int nentry_local = g->my_nQuads * 5 + g->my_nTris * 4 + g->my_nTets * 5 + g->my_nPrisms * 7;
    connectivity = (int *)malloc(sizeof(int) * nentry_local);
    count[0] = nentry_local; // local dim[0]
    count[1] = 0; // ignores for elemental since just dim 1 array
    

    int ie;
    k = 0;
    for (ie=0; ie<g->nelems3d; ie++) {
#ifdef _MESSG
        if (g->elem3d[ie].resident_pe == g->smpi->myid) 
#endif
        {

            if      (g->elem3d[ie].nnodes == 4) {connectivity[k] = 6;} // XDMF CODE
            else if (g->elem3d[ie].nnodes == 6) {connectivity[k] = 7;} // XDMF CODE
            k += 1;
            for (i=0; i<g->elem3d[ie].nnodes; i++) {
                connectivity[k] = g->node[g->elem3d[ie].nodes[i]].gid;
                k++;
            }
        }
    }



    for (ie=0; ie<g->nelems2d; ie++) {
        //if (g->elem2d[ie].bflag == BODY) { // only body elements??  // MAKE SURE THEY ARE RESIDENTIAL!!!
#ifdef _MESSG
        if (g->elem2d[ie].resident_pe == g->smpi->myid)
#endif
         {
            if      (g->elem2d[ie].nnodes == 3) {connectivity[k] = 4;} // XDMF CODE
            else if (g->elem2d[ie].nnodes == 4) {connectivity[k] = 5;} // XDMF CODE
            k += 1;
            for (i=0; i<g->elem2d[ie].nnodes; i++) {
                connectivity[k] = g->node[g->elem2d[ie].nodes[i]].gid;
                k++;
            }
        }
    }



    //add 1d element too
     for (ie=0; ie<g->nelems1d; ie++) {
        //if (g->elem2d[ie].bflag == BODY) { // only body elements??  // MAKE SURE THEY ARE RESIDENTIAL!!!
        //should be g->elem1d[ie].resident_pe == g->smpi->myid
        if (true) {
            if      (g->elem1d[ie].nnodes == 2) {connectivity[k] = 2;} // XDMF CODE
            k += 1;
            for (i=0; i<g->elem1d[ie].nnodes; i++) {
                connectivity[k] = g->node[g->elem1d[ie].nodes[i]].gid;
                k++;
            }
        }
    }


    assert(k == nentry_local);

    // create dataspace and add
    //give global size
    dims[0] = nentry;
    filespace = H5Screate_simple(VECTOR_RANK, dims, NULL);
    
    // create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp2, number, H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    
    // give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(VECTOR_RANK, count, NULL);

    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL); // only contiguous arrays
    
    // select elements
    // to do this we need to loop through each element and find out what is before it
    // the 2 is dimension of our data structure, RANK is spatial dimension of mesh
    hsize_t   coord2[nentry_local];
    
    // global element id number always will follow this convention:
    // 3d -> 2d -> 1d (descending order in terms of # nodes per element)
    // and prisms -> Tets
    // Quads -> Tris
    k=0;
    int gid,nnodes;
    for (ie=0; ie<g->nelems3d; ie++){

        if (g->elem3d[ie].resident_pe == g->smpi->myid) {
        gid = g->elem3d[ie].gid;
        nnodes = g->elem3d[ie].nnodes;
        //use gid to calculate global coord numbers

        //Prisms first
        if(nnodes == 6) {
            for(j=0;j<nnodes+1;j++){
                coord2[k] = gid*(nnodes+1) + j;
                k+=1;
            }
        //Tets second
        }else if(nnodes == 4){
            for(j=0;j<nnodes+1;j++){
                coord2[k] = g->macro_nPrisms*(7) + (gid-g->macro_nPrisms)*(nnodes+1) + j;
                k+=1;
            }

        }
        } 

    }



    for (ie=0; ie<g->nelems2d; ie++){
        gid = g->elem2d[ie].gid;
        nnodes = g->elem2d[ie].nnodes;
        //use gid to calculate global coord numbers
#ifdef _MESSG
        if (g->elem2d[ie].resident_pe == g->smpi->myid)
#endif
        {
        //Quads first
        if(nnodes == 4) {
            for(j=0;j<nnodes+1;j++){
                coord2[k] = g->macro_nPrisms*7 + g->macro_nTets*5 + (gid)*(nnodes+1) + j;
                k+=1;
            }
        //Tris second
        }else if(nnodes == 3){
            for(j=0;j<nnodes+1;j++){
                coord2[k] = g->macro_nPrisms*7 + g->macro_nTets*5 + g->macro_nQuads*5 + (gid-g->macro_nQuads)*(nnodes+1) + j;
                k+=1;
            }

        } 
        }

    }

    //add 1d here
    for (ie=0; ie<g->nelems1d; ie++){
        gid = g->elem1d[ie].gid;
        nnodes = g->elem1d[ie].nnodes;
        //use gid to calculate global coord numbers
        //RN no mpi infor for 1d elements, need to fix
        //should be g->elem1d[ie].resident_pe == g->smpi->myid
        if (true) {
        //line seg only
        if(nnodes == 2) {
            for(j=0;j<nnodes+1;j++){
                coord2[k] = g->macro_nPrisms*7 + g->macro_nTets*5 + g->macro_nQuads*5 + g->macro_nTris*4 + (gid)*(nnodes+1) + j;
                k+=1;
            }
        //Tris second
        }
        }

    }

    //put into appropriate part of global array
    H5Sselect_elements(filespace, H5S_SELECT_SET,nentry_local,(const hsize_t *)&coord2);


    //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #ifdef _MESSG
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif    
    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, connectivity);
    free(connectivity);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp2);
    ////////////////////////////////////
    ///Elemental connections complete///

    //close mesh file
    H5Fclose(file_id);
    
#endif

}
