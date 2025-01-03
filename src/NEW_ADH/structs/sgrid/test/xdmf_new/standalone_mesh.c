
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <hdf5.h>
#define NumElements 4
#define NumNodes 8
#define NEntry 18

void
write_hdf5_data()
{
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

    /* Create the HDF5 file */

    file_id = H5Fcreate("2DUnstructuredMesh.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


    /* Node coordinates data */
       
    float xyz[numNodes][2];          
    
    xyz[0][0] = 0.0;xyz[0][1] = 0.0;
    xyz[1][0] = 1.0;xyz[1][1] = 0.0;    
    xyz[2][0] = 1.0;xyz[2][1] = 1.0;
    xyz[3][0] = 0.0;xyz[3][1] = 1.0;    
    xyz[4][0] = 2.0;xyz[4][1] = 0.0;
    xyz[5][0] = 2.0;xyz[5][1] = 1.0; 
    xyz[6][0] = 2.0;xyz[6][1] = 2.0; 
    xyz[7][0] = 0.0;xyz[7][1] = 2.0;    
      
    /* Write Node coordinates data to file */

    dims[0] = NumNodes;
    dims[1] = 2;    
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/XY", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz); 
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id);     




    /* Connectivity data */
    /*size of array should be numQuad*(4+1) + numTri*(3+1) + numTet*(4+1) + numPrism*(6+1)*/
    /*+1 is for element code*/
    int nentry;
    nentry = numQuad*(4+1) + numTri*(3+1);
    int connectivity[nentry];




    int quadConnectivity[numQuad][4];
    quadConnectivity[0][0] = 0;quadConnectivity[0][1] = 1;quadConnectivity[0][2] = 2;quadConnectivity[0][3] = 3;    
    quadConnectivity[1][0] = 1;quadConnectivity[1][1] = 4;quadConnectivity[1][2] = 5;quadConnectivity[1][3] = 2;

           

    /* Write Connectivity data to file */
    /*want to use Mixed topology type*/
    /* https://xdmf.org/index.php/XDMF_Model_and_Format#Arbitrary */
    /* connectivity will be one long array with elem code then connnections for each eleemnt*/
    /* element code for Triangles: 4*/
    /* element code for Quads: 5*/
    /* element code for Tets: 6*/
    /* element code for Prisms: 8*/

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
    /*
    for (indx=0;indx<nentry;indx++){
        printf("connectivity entry %d = %d\n",indx,connectivity[indx]);
    }
    */

    dims[0] = nentry;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/Elements", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, connectivity);
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id);



    /* Scalar data */

    /*Nodally based dat*/

    float scalardata[NumNodes]; 
        
    for(indx = 0; indx < numNodes; indx++)
    {   
        scalardata[indx] = indx*100;
    }

    /* Write Scalar data to file */

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

    /* Write Scalar data to file */

    dims[0] = NumElements; 
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "/Scalar2", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, scalardata2); 
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id);



    /*Vector data ?*/




    /* Close HDF entities */
      
    status = H5Fclose(file_id);

}
 
void
write_xdmf_xml()
{
    FILE *xmf = 0;
 
    /*
     * Open the file and write the XML description of the mesh.
     */
    xmf = fopen("2DUnstructuredMesh.xmf", "w");

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
            fprintf(xmf, "2DUnstructuredMesh.h5:/XY\n");
            fprintf(xmf, "</DataItem>\n");
            fprintf(xmf, "</Geometry>\n");
            //Geometry finish, this is nodes

            //Topology start
            fprintf(xmf, "<Topology TopologyType=\"Mixed\" NumberOfElements=\"%d\">\n", NumElements);
            //Data Item is the mixed connectivity
            fprintf(xmf, "<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", NEntry);
            fprintf(xmf, "2DUnstructuredMesh.h5:/Elements\n");    
            fprintf(xmf, "</DataItem>\n");
            fprintf(xmf, "</Topology>\n");
            //Topology finish, this contains connectivities

            //Any nodal data we can also write
            fprintf(xmf, "<Attribute Name=\"NodalData\" AttributeType=\"Scalar\" Center=\"Node\">\n");
            fprintf(xmf, "<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", NumNodes);
            fprintf(xmf, "2DUnstructuredMesh.h5:/Scalar1\n");
            fprintf(xmf, "</DataItem>\n");
            fprintf(xmf, "</Attribute>\n");
            //Nodal data finish

            //Should also be way to write elemental data here too just change Center=\"Node" to Center=\"Cell"
            fprintf(xmf, "<Attribute Name=\"ElementalData\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
            fprintf(xmf, "<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", NumElements);
            fprintf(xmf, "2DUnstructuredMesh.h5:/Scalar2\n");
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
    write_hdf5_data();
    write_xdmf_xml();
 
    return 0;
}