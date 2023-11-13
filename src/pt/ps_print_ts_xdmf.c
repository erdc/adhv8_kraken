/* writes a time step to a HDF5 file */

#include "global_header.h"
#ifdef _ADH_HDF5
#include "adh_hdf5.h"
#endif

void ps_print_ts_xdmf(
#ifdef _ADH_HDF5
                      HDF5 *hdf5,
                      char *data_name,  /* String to be written as the data set identifier */
                      void *data,   /* the data */
                      //int nnodes, /* gkc commenting for elem centered data */
                      int nsize,    /* gkc adding for elem centered data */
                      double time,  /* the time of the print */
                      int ntime,    /* index of time step */
                      int center_type,  /* 0=nodal based data, 1=3d-element based data, 
                                           2=2d-element based data */
                      int data_type /* 0=scalar, 1=vector */
#endif
  )
{
#ifdef _ADH_HDF5
  int i;                        /* loop counter */
  hid_t group_id;               /* HDF5 group for data */
  hid_t dataspace, dataset, aid, attr, dataspace_times, filespace;
  hid_t cparms;
  herr_t status;
  hsize_t ndims = 0;
  hsize_t dims[2], maxdims[2], chunkdims[2], size[2];
  hsize_t offset[2], slabdims[2];
  double *scalar_ptr;
  double *vector_data;
  //SVECT2D *vector_ptr;
  SVECT *vector3d_ptr;
  double *data_ptr = NULL;
  int *write_node_index;
  char number[MAXLINE] = "";
  char buff1[MAXLINE] = "";
  xmlNodePtr node1 = NULL, node2 = NULL;
  int compression;
  //int nsize;            /* gkc commenting for elem centered data */
  char type_string[MAXLINE] = "";
  

  // All Data is currently node centered. Leaving this here as a basis if needed later
  // Will only work with PRN_ADPT = TRUE if adaption is turned on, otherwise lots of work for element based data
  if (center_type == NODE_CENTERED)
   {
     //nsize = nnodes;
     sprintf(type_string, "Node");
   }
  else if (center_type == ELEM_3D_CENTERED)
   {
     //nsize = nelems3d;
     sprintf(type_string, "Cell");
     fprintf(stderr, "Cell Output type not supported yet for Particle Tracker");
     exit(-1);
   }
  else if (center_type == ELEM_2D_CENTERED)
   {
     sprintf(type_string, "Cell");
     //tl_error("Not yet supported: 3D element centered data.\n");
     //nsize = nelem2d;
     fprintf(stderr, "Cell Output type not supported yet for Particle Tracker");
     exit(-1);
   }
  else
   {
     //tl_error("Incorrect tag for data centering.\n");
     fprintf(stderr,"Incorrect tag for data centering.; %s:%d\n",__FILE__, __LINE__);
   }
  write_node_index = (int *)malloc(sizeof(int) * nsize);
  if (write_node_index==NULL){fprintf(stderr,"Memory Error: %s:%d",__FILE__,__LINE__);}

  for (i = 0; i < nsize; i++){
    write_node_index[i]=i; /* gkc: May have to be changed later. */
  }

  /* for now, no compression, but can be added as a feature later*/
  compression = 0;

  /* We use the time step index in the name of the HDF group */
  sprintf(number, "%d", ntime);
  /* If this is the first time step, we need to create the HDF group */
  if (ntime == 0)
    {
      group_id = H5Gcreate(hdf5->group_data, data_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      /* create "Times" data set which contains the time stamps */
      dims[0] = 1;
      maxdims[0] = H5S_UNLIMITED;
      dataspace_times = H5Screate_simple(1, dims, maxdims);
      cparms = H5Pcreate(H5P_DATASET_CREATE);
      chunkdims[0] = 20;
      status = H5Pset_chunk(cparms, 1, chunkdims);
      dataset = H5Dcreate(group_id, "Times", H5T_NATIVE_DOUBLE, dataspace_times, H5P_DEFAULT, cparms, H5P_DEFAULT);
      status = H5Pclose(cparms);
      status = H5Sclose(dataspace_times);
      status = H5Dclose(dataset);
    }
  else
    {
      group_id = H5Gopen(hdf5->group_data, data_name, 0);
    }

  /* Create the HDF dataspace */
  if (data_type == SCALAR_DATA)
    {
      ndims = 1;
      dims[0] = nsize;
      scalar_ptr = (double *) data; 
      vector_data = (double *)malloc(sizeof(double)* nsize);
      if (vector_data==NULL){fprintf(stderr,"Memory Error: %s:%d",__FILE__,__LINE__);}
      for (i = 0; i < nsize; i++)
        {
          vector_data[i] = scalar_ptr[write_node_index[i]];
        }        
      data_ptr = (double *) vector_data;
    }
  //else if (data_type == VECTOR2D_DATA)
  //  {
  //    ndims = 2;
  //    dims[0] = nsize;
  //    dims[1] = 3;

  //    /* need to bundle up the data */
  //    vector_ptr = (SVECT2D *) data;
  //    vector_data = (double *)malloc(sizeof(double)* 3*nsize);
  //    if (vector_data==NULL){fprintf(stderr,"Memory Error: %s:%d",__FILE__,__LINE__);}
  //    for (i = 0; i < nsize; i++)
  //      {
  //        vector_data[3 * i] = vector_ptr[write_node_index[i]].x;
  //        vector_data[3 * i + 1] = vector_ptr[write_node_index[i]].y;
  //        vector_data[3 * i + 2] = 0.0;
  //      }
  //    data_ptr = vector_data;
  //  }
  else if (data_type == VECTOR3D_DATA)
    {
      ndims = 2;
      dims[0] = nsize;
      dims[1] = 3;

      /* need to bundle up the data */
      vector3d_ptr = (SVECT *) data;
      vector_data = (double *)malloc(sizeof(double)* 3*nsize);
      if (vector_data==NULL){fprintf(stderr,"Memory Error: %s:%d",__FILE__,__LINE__);}
      for (i = 0; i < nsize; i++)
        {
          vector_data[3 * i] = vector3d_ptr[write_node_index[i]].x;
          vector_data[3 * i + 1] = vector3d_ptr[write_node_index[i]].y;
          vector_data[3 * i + 2] = vector3d_ptr[write_node_index[i]].z;
        }
      data_ptr = vector_data;
    }

  /* Write the data to the HDF file */
  xdmf_write_dataset(group_id, number, H5T_NATIVE_DOUBLE, compression, ndims, dims, data_ptr);

  /* Add attribute to data set which contains the corresponding mesh id */
  xdmf_write_attribute(group_id, number, "Mesh_ID", H5T_NATIVE_INT, &(hdf5->current_mesh));

  /* Add attribute to data set which contains the time stamp */
  xdmf_write_attribute(group_id, number, "time", H5T_NATIVE_DOUBLE, &time);

  /* write the time stamp to the Times dataset */
  dataset = H5Dopen(group_id, "Times", H5P_DEFAULT);
  if (ntime != 0)
    {
      size[0] = ntime + 1;
      status = H5Dset_extent(dataset, size);
    }
  filespace = H5Dget_space(dataset);
  offset[0] = ntime;
  slabdims[0] = 1;
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, slabdims, NULL);

  dims[0] = 1;
  dataspace_times = H5Screate_simple(1, dims, NULL);
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace_times, filespace, H5P_DEFAULT, &time);

  status = H5Sclose(filespace);
  status = H5Dclose(dataset);
  status = H5Sclose(dataspace_times);
  status = H5Gclose(group_id);

  /* Add data entries to XML tree */
  node1 = xmlNewNode(NULL, BAD_CAST "Attribute");
  if (data_type == SCALAR_DATA)
    {
      sprintf(buff1, "Scalar");
    }
  else
    {
      sprintf(buff1, "Vector");
    }
  xmlNewProp(node1, BAD_CAST "AttributeType", BAD_CAST buff1);
  xmlNewProp(node1, BAD_CAST "Center", BAD_CAST type_string);
  xmlNewProp(node1, BAD_CAST "Name", BAD_CAST data_name);
  sprintf(buff1, "%s:Data/%s/%d", hdf5->h5_filename, data_name, ntime);
  node2 = xmlNewChild(node1, NULL, BAD_CAST "DataItem", BAD_CAST buff1);
  xmlNewProp(node2, BAD_CAST "DataType", BAD_CAST "Float");

  if (data_type == SCALAR_DATA)
    {
      sprintf(buff1, "%d", nsize);
    }
  else
    {
      sprintf(buff1, "%d 3", nsize);
    }

  xmlNewProp(node2, BAD_CAST "Dimensions", BAD_CAST buff1);
  xmlNewProp(node2, BAD_CAST "Format", BAD_CAST "HDF");
  xmlNewProp(node2, BAD_CAST "Precision", BAD_CAST "8");
  xmlAddChild(hdf5->current_grid_node, node1);
    if (data_type == SCALAR_DATA)
    {
      //vector_data = (double *)tl_free(sizeof(double), nsize, vector_data);
      free(vector_data);
    }
  else if (data_type == VECTOR2D_DATA)
    {
      //vector_data = (double *)tl_free(sizeof(double), 3*nsize, vector_data);
      free(vector_data);
    }
  else if (data_type == VECTOR3D_DATA)
    {
      //vector_data = (double *)tl_free(sizeof(double), 3*nsize, vector_data);
      free(vector_data);
    }
  //write_node_index = (int *)tl_free(sizeof(int), nnodes, write_node_index);
  //write_node_index = (int *)tl_free(sizeof(int), nsize, write_node_index);
  free(write_node_index);

#endif
  return;
}
