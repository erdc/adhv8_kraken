
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <hdf5.h>
#include <mpi.h>
#include <assert.h>
#define NumElements 4
#define NumNodes 8
#define NEntry 18
#define RANK 2
//RANK is dimesnion of data set, 2 is a matrix, 1 is a vector

void
write_hdf5_data(MPI_Comm comm, MPI_Info info, int mpi_rank, int mpi_size)
{

    //for this simple test case assert the size is only 4
    assert(mpi_size==4);


    hid_t     file_id, dset_id, memspace, filespace, plist_id;
    hsize_t   dims[RANK];
    //float *xyz;
    int i;
    hsize_t   offset[RANK];
    herr_t    status; 
    hsize_t   count[RANK];
    int nodes_per_pe,indx;


    //mesh properties
    int numQuad;
    int numTri;
    int *connectivity;
    numQuad = 2;
    numTri = 2;

    //tutorial at
    //https://docs.hdfgroup.org/hdf5/v1_14/v1_14_4/_ex_a_p_i.html


    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //setting options for parallel
    H5Pset_fapl_mpio(plist_id, comm, info);
    //H5Pset_all_coll_metadata_ops(plist_id, true);
    //H5Pset_coll_metadata_write(plist_id, true);


    //create file
    file_id = H5Fcreate("2DUnstructuredMesh_parallel.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);


    //////////////////////////////
    //Writing Nodes//////////////

    //declare global sizes of data set
    //Nodes first
    //this is local data
    nodes_per_pe=2;
    //xyz = (float *)malloc(sizeof(float) * nodes_per_pe * 2);
    //xyz = (float *)malloc(sizeof(float[nodes_per_pe][2]));
    float (*xyz)[2] = malloc(sizeof(float[nodes_per_pe][2]));

    //Hard code partition this time
    // 2 nodes per partition
    if(mpi_rank==0){
        xyz[0][0] = 0.0;xyz[0][1] = 0.0;
        xyz[1][0] = 1.0;xyz[1][1] = 0.0; }
    else if(mpi_rank==1){   
        xyz[0][0] = 1.0;xyz[0][1] = 1.0;
        xyz[1][0] = 0.0;xyz[1][1] = 1.0;} 
    else if(mpi_rank==2){   
        xyz[0][0] = 2.0;xyz[0][1] = 0.0;
        xyz[1][0] = 2.0;xyz[1][1] = 1.0; }
    else if(mpi_rank==3){
        xyz[0][0] = 2.0;xyz[0][1] = 2.0; 
        xyz[1][0] = 0.0;xyz[1][1] = 2.0;} 



    


    //create a dataspace
    // this is a global quantity
    dims[0]  = NumNodes;
    dims[1]  = 2;
    filespace = H5Screate_simple(2, dims, NULL);

    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(file_id, "/XY", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    //create a simple dataspace where each process has a few of the rows
    count[0]  = nodes_per_pe;
    count[1]  = dims[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;

    //give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(RANK, count, NULL);


    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);


    //create data pointer
    //data = (int *)malloc(sizeof(int) * count[0] * count[1]);
    //for (i = 0; i < count[0] * count[1]; i++) {
    //    data[i] = mpi_rank + 10;
    //}

    //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, xyz);
    free(xyz);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);

    //////////////////////////
    ////Nodal write complete//


    //////////////////////////
    //Element data set///////
    // Connectivity data 
    //size of array should be numQuad*(4+1) + numTri*(3+1) + numTet*(4+1) + numPrism*(6+1)
    //+1 is for element code

    //global attributes
    int nentry;
    nentry = numQuad*(4+1) + numTri*(3+1);
    
    //local data, different on each PE
    //first two PEs hold one quad. second PEs hold one triangle
    if(mpi_rank == 0 || mpi_rank ==1){ 
        connectivity = (int *)malloc(sizeof(int) * 5);}
    else if (mpi_rank ==2 || mpi_rank ==3){
        connectivity = (int *)malloc(sizeof(int) * 4);
    }
    
    //load in connectivity values
    if (mpi_rank==0){
        count[0] = 5;
        offset[0]=0;
        connectivity[0] = 5;
        connectivity[1] = 0;connectivity[2] = 1;connectivity[3] = 2;connectivity[4] = 3;   
    }
    else if(mpi_rank==1){
        count[0]=5;
        offset[0]=5;
        connectivity[0] = 5;
        connectivity[1] = 1;connectivity[2] = 4;connectivity[3] = 5;connectivity[4] = 2;
    }
    else if(mpi_rank==2){
        count[0]=4;
        offset[0]=10;
        connectivity[0] = 4;
        connectivity[1] = 2;connectivity[2] = 5;connectivity[3] = 6;   
        
    }
    else if(mpi_rank==3){
        count[0]=4;
        offset[0]=14;
        connectivity[0] = 4;
        connectivity[1] = 3;connectivity[2] = 2;connectivity[3] = 7;
    }

    //create dataspace and add
    //give global size
    dims[0] = nentry;
    filespace = H5Screate_simple(1, dims, NULL);
    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(file_id, "/Elements", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    //give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(1, count, NULL);

    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

     //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, connectivity);
    free(connectivity);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);

    ////////////////////////////////////
    ///Elemental connections complete///


    //////Scalar data write////////////////
    //Nodally based data first

    float *scalardata; 
    scalardata = (float *)malloc(sizeof(float) * nodes_per_pe);
        
    for(indx = 0; indx < nodes_per_pe; indx++)
    {   
        scalardata[indx] = (mpi_rank*nodes_per_pe+indx)*100;
    }

    dims[0] = NumNodes;
    count[0] = nodes_per_pe;
    offset[0] = mpi_rank * count[0];
    filespace = H5Screate_simple(1, dims, NULL);
    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(file_id, "/Scalar1", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    //give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(1, count, NULL);

    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

     //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, scalardata);
    free(scalardata);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);


    //End of Nodal based scalar


    //Write some elementally based data ////////////
    //Elementally based data

     
    float *scalardata2; 
    scalardata2 = (float *)malloc(sizeof(float) * 1);
        
    for(indx = 0; indx < 1; indx++)
    {   
        scalardata2[indx] = (mpi_rank)*100;
    }

    dims[0] = NumElements;
    count[0] = 1;
    offset[0] = mpi_rank * count[0];
    filespace = H5Screate_simple(1, dims, NULL);
    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(file_id, "/Scalar2", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    //give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(1, count, NULL);

    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

     //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, scalardata2);
    free(scalardata2);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);


    H5Fclose(file_id);

    if (mpi_rank == 0)
        printf("PHDF5 example finished with no errors\n");

    /*
    hid_t     file_id;
    hid_t     dataset_id, dataspace_id;
    hsize_t   dims[2];
    herr_t    status;   
    int       indx;
    int numQuad;
    int numTri;
    int numNodes;
    int j,temp;

    numQuad = 2;
    numTri = 2;
    numNodes = 8;

    // Create the HDF5 file 

    file_id = H5Fcreate("2DUnstructuredMesh.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


    // Node coordinates data 
       
    float xyz[numNodes][2];          
    
    xyz[0][0] = 0.0;xyz[0][1] = 0.0;
    xyz[1][0] = 1.0;xyz[1][1] = 0.0;    
    xyz[2][0] = 1.0;xyz[2][1] = 1.0;
    xyz[3][0] = 0.0;xyz[3][1] = 1.0;    
    xyz[4][0] = 2.0;xyz[4][1] = 0.0;
    xyz[5][0] = 2.0;xyz[5][1] = 1.0; 
    xyz[6][0] = 2.0;xyz[6][1] = 2.0; 
    xyz[7][0] = 0.0;xyz[7][1] = 2.0;    
      
    // Write Node coordinates data to file 

    dims[0] = NumNodes;
    dims[1] = 2;    
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/XY", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz); 
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id);     




    // Connectivity data 
    //size of array should be numQuad*(4+1) + numTri*(3+1) + numTet*(4+1) + numPrism*(6+1)
    //+1 is for element code
    int nentry;
    nentry = numQuad*(4+1) + numTri*(3+1);
    int connectivity[nentry];




    int quadConnectivity[numQuad][4];
    quadConnectivity[0][0] = 0;quadConnectivity[0][1] = 1;quadConnectivity[0][2] = 2;quadConnectivity[0][3] = 3;    
    quadConnectivity[1][0] = 1;quadConnectivity[1][1] = 4;quadConnectivity[1][2] = 5;quadConnectivity[1][3] = 2;

           

    // Write Connectivity data to file 
    //want to use Mixed topology type
    //https://xdmf.org/index.php/XDMF_Model_and_Format#Arbitrary
    //connectivity will be one long array with elem code then connnections for each eleemnt
    //element code for Triangles: 4
    //element code for Quads: 5
    //element code for Tets: 6
    //element code for Prisms: 8

    // loop through quads first
    for(indx =0;indx < numQuad;indx++){
        connectivity[indx*5] = 5;
        for(j=1;j<5;j++){
            connectivity[indx*5+j]=quadConnectivity[indx][j-1];
        }
    }

    temp = (indx-1)*5+j;


    //loop through triangles next
    int triConnectivity[numTri][3];

    triConnectivity[0][0] = 2;triConnectivity[0][1] = 5;triConnectivity[0][2] = 6;   
    triConnectivity[1][0] = 3;triConnectivity[1][1] = 2;triConnectivity[1][2] = 7;


    for(indx =0;indx < numTri;indx++){
        connectivity[indx*4+temp] = 4;
        for(j=1;j<4;j++){
            connectivity[indx*4+j+temp]=triConnectivity[indx][j-1];
        }
    }

    //Debugging
    
    //for (indx=0;indx<nentry;indx++){
    //    printf("connectivity entry %d = %d\n",indx,connectivity[indx]);
    //}
    

    dims[0] = nentry;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/Elements", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, connectivity);
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id);



    // Scalar data 

    //Nodally based dat

    float scalardata[NumNodes]; 
        
    for(indx = 0; indx < numNodes; indx++)
    {   
        scalardata[indx] = indx*100;
    }

    //Write Scalar data to file

    dims[0] = NumNodes; 
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/Scalar1", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, scalardata); 
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id);



    //Elementally based data

    float scalardata2[NumElements]; 
        
    for(indx = 0; indx < NumElements; indx++)
    {   
        scalardata2[indx] = indx*100;
    }

    //Write Scalar data to file

    dims[0] = NumElements; 
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/Scalar2", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, scalardata2); 
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id);



    //Vector data ?




    // Close HDF entities 
      
    status = H5Fclose(file_id);
    */
}
 
void
write_xdmf_xml()
{
    FILE *xmf = 0;
 
    /*
     * Open the file and write the XML description of the mesh.
     */
    xmf = fopen("2DUnstructuredMesh_parallel.xmf", "w");

    //breaks stuff?
    //fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    //fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"3.0\">\n");
    //Domain information
    fprintf(xmf, "<Domain>\n");
        //Grid start
        fprintf(xmf, "<Grid Name=\"2D Unstructured Mesh\">\n");
            //Geometry start, this is nodes
            fprintf(xmf, "<Geometry GeometryType=\"XY\">\n");
            fprintf(xmf, "<DataItem Dimensions=\"%d %d\" Format=\"HDF\">\n", NumNodes,2);
            fprintf(xmf, "2DUnstructuredMesh_parallel.h5:/XY\n");
            fprintf(xmf, "</DataItem>\n");
            fprintf(xmf, "</Geometry>\n");
            //Geometry finish, this is nodes

            //Topology start
            fprintf(xmf, "<Topology TopologyType=\"Mixed\" NumberOfElements=\"%d\">\n", NumElements);
            //Data Item is the mixed connectivity
            fprintf(xmf, "<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", NEntry);
            fprintf(xmf, "2DUnstructuredMesh_parallel.h5:/Elements\n");    
            fprintf(xmf, "</DataItem>\n");
            fprintf(xmf, "</Topology>\n");
            //Topology finish, this contains connectivities

            //Any nodal data we can also write
            fprintf(xmf, "<Attribute Name=\"NodalData\" AttributeType=\"Scalar\" Center=\"Node\">\n");
            fprintf(xmf, "<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", NumNodes);
            fprintf(xmf, "2DUnstructuredMesh_parallel.h5:/Scalar1\n");
            fprintf(xmf, "</DataItem>\n");
            fprintf(xmf, "</Attribute>\n");
            //Nodal data finish

            //Should also be way to write elemental data here too just change Center=\"Node" to Center=\"Cell"
            fprintf(xmf, "<Attribute Name=\"ElementalData\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
            fprintf(xmf, "<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", NumElements);
            fprintf(xmf, "2DUnstructuredMesh_parallel.h5:/Scalar2\n");
            fprintf(xmf, "</DataItem>\n");
            fprintf(xmf, "</Attribute>\n");

        
        fprintf(xmf, "</Grid>\n");
        //Grid Finish
    fprintf(xmf, "</Domain>\n");
    //Domain finished
    fprintf(xmf, "</Xdmf>\n");
    //XDMF File finished
    
    fclose(xmf);
}
 
int
main(int argc, char *argv[])
{

    MPI_Comm comm  = MPI_COMM_WORLD;
    MPI_Info info  = MPI_INFO_NULL;

    int mpi_size,mpi_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
    write_hdf5_data(comm,info,mpi_rank,mpi_size);
    if (mpi_rank ==0){
        write_xdmf_xml();
    }
    MPI_Finalize();
 
    return 0;
}