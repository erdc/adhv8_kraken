#include "adh.h"
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;
void init_hdf5_file(SGRID *g)
{
    #ifdef _ADH_HDF5
    //this function should do the following
    // 1. Create HDf5 File for the simulation
    // 2. Create all necessary data groups for the run



    hid_t     file_id, dset_id, memspace, filespace, plist_id;
    hid_t     grp1,grp2,grp3,grp4,grp5,grp6,grp7,grp8,grp9;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //setting options for parallel

    #ifdef _MESSG
    MPI_Info info  = MPI_INFO_NULL;
    H5Pset_fapl_mpio(plist_id, g->smpi->ADH_COMM, info);
    #endif
    //H5Pset_all_coll_metadata_ops(plist_id, true);
    //H5Pset_coll_metadata_write(plist_id, true);


    //create file
    char fname[50];
    strcpy(fname,g->filename);
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

    #endif

}

