
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include <hdf5.h>
#include <mpi.h>
#include <assert.h>

void
init_hdf5_file(MPI_Comm comm, MPI_Info info, int mpi_rank, int mpi_size, char *fbase)
{
    assert(mpi_size==4);
    //this function should do the following
    // 1. Create HDf5 File for the simulation
    // 2. Create all necessary data groups for the run

    hid_t     file_id, dset_id, memspace, filespace, plist_id;
    hid_t     grp1,grp2,grp3,grp4,grp5,grp6,grp7,grp8,grp9;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //setting options for parallel
    H5Pset_fapl_mpio(plist_id, comm, info);
    //H5Pset_all_coll_metadata_ops(plist_id, true);
    //H5Pset_coll_metadata_write(plist_id, true);


    //create file
    char fname[50];
    strcpy(fname,fbase);
    strcat(fname, ".h5");
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);


    //1st 2 groups always active even if data is empty

    //group 1 is anything with mesh
    grp1 = H5Gcreate(file_id, "/Mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


    //group 2 contains any info with data
    grp2 = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //Eventually want a routine to handle this based on active models but for
    

    //now just hard code some exmamples
    grp3 = H5Gcreate(grp2, "Time", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    grp4 = H5Gcreate(grp2, "NodalScalar", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    grp5 = H5Gcreate(grp2, "ElementalScalar", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    grp6 = H5Gcreate(grp2, "NodalVector", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    grp7 = H5Gcreate(grp2, "ElementalVector", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //these groups on the other hand will always be there so this can be hard coded
    grp8 = H5Gcreate(grp1, "XY", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    grp9 = H5Gcreate(grp1, "Elements", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);



    H5Gclose(grp1);
    H5Gclose(grp2);
    H5Gclose(grp3);
    H5Gclose(grp4);
    H5Gclose(grp5);
    H5Gclose(grp6);
    H5Gclose(grp7);
    H5Gclose(grp8);
    H5Gclose(grp9);
    H5Fclose(file_id);


}

void write_hdf5_mesh(MPI_Comm comm, MPI_Info info, int mpi_rank, int mpi_size, char *fbase, int nt){
    //given a time step number, rights a mesh to the hdf5 file

    hid_t     file_id, plist_id, filespace, dset_id, memspace;
    hid_t     grp1,grp2;
    char      fname[50];
    hsize_t   dims[RANK];
    hsize_t   offset[RANK];
    herr_t    status; 
    hsize_t   count[RANK];
    int nodes_per_pe,indx;
    int elem_per_pe=1;


    // for adaptive grid writing
    char number[50] = "";
    sprintf(number, "%d", nt);

     //mesh properties
    int numQuad;
    int numTri;
    int *connectivity;
    numQuad = 2;
    numTri = 2;




    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //setting options for parallel
    H5Pset_fapl_mpio(plist_id, comm, info);
    //H5Pset_all_coll_metadata_ops(plist_id, true);
    //H5Pset_coll_metadata_write(plist_id, true);

    //open file
    strcpy(fname,fbase);
    strcat(fname, ".h5");
    file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);


    //group 1 is anything with mesh nodes
    grp1 = H5Gopen(file_id, "/Mesh/XY", H5P_DEFAULT);

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
        xyz[0][0] = 2.0;xyz[0][1] = 2.0+nt; 
        xyz[1][0] = 0.0;xyz[1][1] = 2.0+nt;} 



    


    //create a dataspace
    // this is a global quantity
    dims[0]  = NumNodes;
    dims[1]  = 2;
    filespace = H5Screate_simple(2, dims, NULL);

    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp1, number, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
    H5Gclose(grp1);
    //////////////////////////
    ////Nodal write complete//



    //Now write elemental data
    //////////////////////////
    //Element data set///////
    // Connectivity data 
    //size of array should be numQuad*(4+1) + numTri*(3+1) + numTet*(4+1) + numPrism*(6+1)
    //+1 is for element code
    grp2 = H5Gopen(file_id, "/Mesh/Elements", H5P_DEFAULT);




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
    dset_id = H5Dcreate(grp2, number, H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
    H5Gclose(grp2);

    ////////////////////////////////////
    ///Elemental connections complete///

    //close mesh file
    H5Fclose(file_id);

}

void
write_hdf5_data(MPI_Comm comm, MPI_Info info, int mpi_rank, int mpi_size, char *fbase, int nt, float t)
{
    //writes hdf5 data at particular time step nt
    //for this simple test case assert the size is only 4
    assert(mpi_size==4);


    hid_t     file_id, dset_id, memspace, filespace, plist_id, grp;
    hsize_t   dims[RANK];
    //float *xyz;
    int i;
    hsize_t   offset[RANK];
    herr_t    status; 
    hsize_t   count[RANK];
    int nodes_per_pe,indx;
    nodes_per_pe=2;
    int elem_per_pe=1;


    //mesh properties
    int numQuad;
    int numTri;
    int *connectivity;
    numQuad = 2;
    numTri = 2;


    char number[50] = "";
    sprintf(number, "%d", nt);

    //tutorial at
    //https://docs.hdfgroup.org/hdf5/v1_14/v1_14_4/_ex_a_p_i.html



    //Open file for parallel i/o
    ////////////////////////////////////////////////
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //setting options for parallel
    H5Pset_fapl_mpio(plist_id, comm, info);
    //H5Pset_all_coll_metadata_ops(plist_id, true);
    //H5Pset_coll_metadata_write(plist_id, true);


    //open file
    char fname[50];
    strcpy(fname,fbase);
    strcat(fname, ".h5");
    file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);
    ////////////////////////////////////////////////



    
    //Open up data groups
    //group 2 contains any info with data
    grp = H5Gopen(file_id, "/Data/NodalScalar", H5P_DEFAULT);




    
    //////Scalar data write////////////////
    //Nodally based data first

    float *scalardata; 
    scalardata = (float *)malloc(sizeof(float) * nodes_per_pe);
        
    for(indx = 0; indx < nodes_per_pe; indx++)
    {   
        scalardata[indx] = (mpi_rank*nodes_per_pe+indx)*100 + t*50;
    }

    dims[0] = NumNodes;
    count[0] = nodes_per_pe;
    offset[0] = mpi_rank * count[0];
    filespace = H5Screate_simple(1, dims, NULL);
    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp, number, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
    H5Gclose(grp);



    //End of Nodal based scalar

    


    //Write some elementally based data ////////////
    //Elementally based data
    grp = H5Gopen(file_id, "/Data/ElementalScalar", H5P_DEFAULT);
     
    float *scalardata2; 
    scalardata2 = (float *)malloc(sizeof(float) * elem_per_pe);
        
    for(indx = 0; indx < 1; indx++)
    {   
        scalardata2[indx] = (mpi_rank)*100 + t*50;
    }

    dims[0] = NumElements;
    count[0] = 1;
    offset[0] = mpi_rank * count[0];
    filespace = H5Screate_simple(1, dims, NULL);
    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp, number, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
    H5Gclose(grp);

    //End of elementally based data
    /////////////////////////////////////////////

    
    ////////////////////////////////////////////
    
    //vector based nodal data
    grp = H5Gopen(file_id, "/Data/NodalVector", H5P_DEFAULT);
    float (*vectordata)[2] = malloc(sizeof(float[nodes_per_pe][2]));
        
    for(indx = 0; indx < nodes_per_pe; indx++)
    {   
        vectordata[indx][0] = 1+t;
        vectordata[indx][1] = 1;
    }

    dims[0] = NumNodes;
    dims[1] = 2;
    count[0]  = nodes_per_pe;
    count[1]  = dims[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    filespace = H5Screate_simple(2, dims, NULL);
    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp, number, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    //give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(RANK, count, NULL);

    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

     //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, vectordata);
    free(vectordata);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp);
    //end of nodal vectors
    /////////////////////////////////////////////////////


    //Elementally based vectors
    //////////////////////////////////////////////////////
    grp =  H5Gopen(file_id, "/Data/ElementalVector", H5P_DEFAULT);
    float (*vectordata2)[2] = malloc(sizeof(float[elem_per_pe][2]));
        
    for(indx = 0; indx < elem_per_pe; indx++)
    {   
        vectordata2[indx][0] = 1-t;
        vectordata2[indx][1] = 1;
    }

    dims[0] = NumElements;
    dims[1] = 2;
    count[0]  = elem_per_pe;
    count[1]  = dims[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    filespace = H5Screate_simple(2, dims, NULL);
    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp, number, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    //give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(RANK, count, NULL);

    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

     //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, vectordata2);
    free(vectordata2);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp);

    



    /////////////////////////////////////////////////////
    H5Fclose(file_id);

    if (mpi_rank == 0)
        printf("PHDF5 example finished with no errors\n");
  
}
 
void xdmf_init(FILE *xmf, char *fbase){
    char fname[50];
    strcpy(fname, fbase);
    strcat(fname, ".xmf");
    xmf = fopen(fname, "w");
    //breaks stuff?
    //fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    //fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"3.0\">\n");
    //Domain information
    fprintf(xmf, "\t <Domain Name=\"Adh Sim\">\n");
    //Time based grid start
    fprintf(xmf, "\t\t<Grid Name=\"MeshTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
    fclose(xmf);
}


void xdmf_write_dataset(FILE *xmf, char *fbase, int nt, int mesh_no, float t){
    char fname[50];
    strcpy(fname, fbase);
    strcat(fname, ".xmf");
    xmf = fopen(fname, "a");




    //Grid start
        fprintf(xmf, "\t\t\t<Grid Name=\"2D Unstructured Mesh\">\n");
            fprintf(xmf, "\t\t\t<Time Value=\"%f\" />\n",t);
            //Geometry start, this is nodes
            fprintf(xmf, "\t\t\t\t<Geometry GeometryType=\"XY\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d %d\" Format=\"HDF\">\n", NumNodes,2);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Mesh/XY/%d\n",fbase,mesh_no);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Geometry>\n");
            //Geometry finish, this is nodes

            //Topology start
            fprintf(xmf, "\t\t\t\t<Topology TopologyType=\"Mixed\" NumberOfElements=\"%d\">\n", NumElements);
            //Data Item is the mixed connectivity
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", NEntry);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Mesh/Elements/%d\n",fbase,mesh_no);    
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Topology>\n");
            //Topology finish, this contains connectivities
            
            
            //Any nodal data we can also write
            fprintf(xmf, "\t\t\t\t<Attribute Name=\"NodalData\" AttributeType=\"Scalar\" Center=\"Node\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", NumNodes);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/NodalScalar/%d\n",fbase,nt);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Attribute>\n");
            //Nodal data finish
            

            //Should also be way to write elemental data here too just change Center=\"Node" to Center=\"Cell"
            fprintf(xmf, "\t\t\t\t<Attribute Name=\"ElementalData\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", NumElements);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/ElementalScalar/%d\n",fbase,nt);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Attribute>\n");
            
            
            //Should also be way to write elemental data here too just change Center=\"Node" to Center=\"Cell"
            fprintf(xmf, "\t\t\t\t<Attribute Name=\"NodalVector\" AttributeType=\"Vector\" Center=\"Node\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d %d\" Format=\"HDF\">\n", NumNodes,RANK);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/NodalVector/%d\n",fbase,nt);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Attribute>\n");

            //Should also be way to write elemental data here too just change Center=\"Node" to Center=\"Cell"
            fprintf(xmf, "\t\t\t\t<Attribute Name=\"ElementalVector\" AttributeType=\"Vector\" Center=\"Cell\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d %d\" Format=\"HDF\">\n", NumElements,RANK);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/ElementalVector/%d\n",fbase,nt);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Attribute>\n");
            
            
            
        
        fprintf(xmf, "\t\t\t</Grid>\n");
        fclose(xmf);
}

void xdmf_finalize(FILE *xmf, char *fbase){
    char fname[50];
    strcpy(fname,fbase);
    strcat(fname, ".xmf");
    xmf = fopen(fname, "a");


     //Grid Finish
    fprintf(xmf, "\t\t</Grid>\n");
    fprintf(xmf, "\t</Domain>\n");
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
    char *fbase;
    int idx;
    fbase="Timseries_parallel2";
    //printf("%s \n",fbase);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    // loop through time and write a vector as well as scalar function
    int nt = 5;
    float t=0.0;
    float dt = 2.0;




    //split up hdf5 into basic parts
    //always initialize file and write initial mesh, and initial conditions
    init_hdf5_file(comm,info,mpi_rank,mpi_size,fbase);
    write_hdf5_mesh(comm,info,mpi_rank,mpi_size,fbase,0);
    write_hdf5_data(comm,info,mpi_rank,mpi_size,fbase,0,t);

    if (mpi_rank==0){printf("Initializing xmf file \n");
                    FILE *xmf = 0;
                    xdmf_init(xmf,fbase);
                    xdmf_write_dataset(xmf,fbase,idx,idx,t);}

    
    for (idx=1;idx<nt+1;idx++){
        
        t+=dt;
        //adaptive meshing
        write_hdf5_mesh(comm,info,mpi_rank,mpi_size,fbase,idx);
        write_hdf5_data(comm,info,mpi_rank,mpi_size,fbase,idx,t);

        if (mpi_rank ==0){
            printf("Time step %d\n",idx);
            //
            // Open the file and write the XML description of the mesh.
            //
            FILE *xmf = 0;
            xdmf_write_dataset(xmf, fbase,idx,idx,t);

        }
    }
    
    if (mpi_rank==0){printf("closing xmf file\n");
                    FILE *xmf = 0;
                    xdmf_finalize(xmf, fbase);}
    MPI_Finalize();
 
    return 0;
}

// CJT for writing mesh from memory
void sgrid_write_hdf5(SGRID *g){
    //given a time step number, rights a mesh to the hdf5 file
    
    int RANK = g->smpi->myid;

    hid_t     file_id, plist_id, filespace, dset_id, memspace;
    hid_t     grp1,grp2;
    char      fname[50];
    hsize_t   dims[RANK];
    hsize_t   offset[RANK];
    herr_t    status;
    hsize_t   count[RANK];
    
    // for adaptive grid writing
    char number[50] = "";
    sprintf(number, "%d", nt);

    //mesh properties
    int numQuad = g->nQuads;
    int numTri = g->nTris;
    int numTet = g->nTets;
    int numPrism; = g->nPrism;
    int *connectivity;

    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //setting options for parallel
    H5Pset_fapl_mpio(plist_id, comm, info);
    //H5Pset_all_coll_metadata_ops(plist_id, true);
    //H5Pset_coll_metadata_write(plist_id, true);

    //open file
    strcpy(fname,g->filename);
    strcat(fname, ".h5");
    file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);


    //group 1 is anything with mesh nodes
    grp1 = H5Gopen(file_id, "/Mesh/XY", H5P_DEFAULT);

    //////////////////////////////
    //Writing Nodes//////////////

    //declare global sizes of data set nodes first. This is local data.
    float (*xyz)[2] = malloc(sizeof(float[g->my_nnodes][3]));

    //create a dataspace.  This is a global quantity
    dims[0]  = g->macro_nnodes;
    dims[1]  = 3;
    filespace = H5Screate_simple(2, dims, NULL);

    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp1, number, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    //create a simple dataspace where each process has a few of the rows
    count[0]  = g->my_nnodes;
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
    H5Gclose(grp1);
    //////////////////////////
    ////Nodal write complete//



    //Now write elemental data
    //////////////////////////
    //Element data set///////
    // Connectivity data
    //size of array should be numQuad*(4+1) + numTri*(3+1) + numTet*(4+1) + numPrism*(6+1)
    //+1 is for element code
    grp2 = H5Gopen(file_id, "/Mesh/Elements", H5P_DEFAULT);




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
    dset_id = H5Dcreate(grp2, number, H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
    H5Gclose(grp2);

    ////////////////////////////////////
    ///Elemental connections complete///

    //close mesh file
    H5Fclose(file_id);

}
