/* writes a 3D time step to a HDF5 file */

#include "global_header.h"
#ifdef _ADH_HDF5
#include "adh_hdf5.h"
#endif

void ps_print_surf_ts_xdmf(
#ifdef _ADH_HDF5
                           HDF5 *hdf5,   /* The HDF5 file */
                           char *data_name, /* String written as the data set identifier */
                           void *data,  /* the data */
                           int nnodes_sur,
                           double time, /* the time of the print */
                           int ntime,   /* index of time step */
                           int center_type, /* 0=nodal based data, 1=3d-element based data, 
                                               2=2d-element based data */
                           int data_type    /* 0=scalar, 1=vector */
#endif
  )
{
#ifdef _ADH_HDF5
  int i;                        /* loop counter */
  hid_t group_id;               /* HDF5 group for data */
  hid_t dataset, dataspace_times, filespace;
  hid_t cparms;
  herr_t status;
  hsize_t ndims;
  hsize_t dims[2], maxdims[2], chunkdims[2], size[2];
  hsize_t offset[2], slabdims[2];
  double *vector_data = NULL;
  SVECT2D *vector_ptr = NULL;
  SVECT *vector3d_ptr = NULL;
  double *data_ptr = NULL;
  char number[MAXLINE] = "";
  char buff1[MAXLINE] = "";
  xmlNodePtr node1 = NULL, node2 = NULL;
  int compression;
  int nsize;
  int index;
  char type_string[MAXLINE] = "";
  
  /* No compression yet LP */
  compression = 0;

  /* We use the time step index in the name of the HDF group */
  sprintf(number, "%d", ntime);

  /* If this is the first time step, we need to create the HDF group */
  if (ntime == 0)
    {
      group_id = H5Gcreate(hdf5->group_surf_data, data_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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
      group_id = H5Gopen(hdf5->group_surf_data, data_name, 0);
    }

  // /* All data is currently node centered. */
  // /* Set the data size and type */
  // if (center_type == NODE_CENTERED)
  //   {
  //     nsize = number_surface_nodes;
  //     sprintf(type_string, "Node");
  //   }
  // else if (center_type == ELEM_3D_CENTERED)
  //   {
  //     tl_error("Surface print should be 2D or 1D.\n");
  //   }
  // else if (center_type == ELEM_2D_CENTERED)
  //   {
  //     nsize = nelem2d;
  //     sprintf(type_string, "Cell");
  //   }
  // else
  //   {
  //     tl_error("Incorrect tag for data centering.\n");
  //   }
  if (center_type != NODE_CENTERED)
    {
      tl_error("Incorrect tag for data centering. Only NODE_CENTER is supported. No data will be written.\n");
      return;
    }
  nsize = nnodes_sur;
  sprintf(type_string, "Node");

  if (nsize > 0)
    {
      /* Create the HDF dataspace */
      if (data_type == SCALAR_DATA) {
          ndims = 1;
          dims[0] = nsize;
          data_ptr = (double *) data;
      }
      else if (data_type == VECTOR2D_DATA) {
          ndims = 2;
          dims[0] = nsize;
          dims[1] = 3;

          /* need to bundle up the data */
          vector_ptr = (SVECT2D *) data;
          vector_data = (double *) tl_alloc(sizeof(double), 3 * nsize);
          index = 0;
          for (i = 0; i < nnodes_sur; i++) {
              vector_data[index] = vector_ptr[i].x;
              vector_data[index + 1] = vector_ptr[i].y;
              vector_data[index + 2] = 0.0;
              index += 3;
          }
          data_ptr = vector_data;
      }
      else if (data_type == VECTOR3D_DATA) {
          ndims = 2;
          dims[0] = nsize;
          dims[1] = 3;

          /* need to bundle up the data */
          vector3d_ptr = (SVECT *) data;
          vector_data = (double *) tl_alloc(sizeof(double), 3 * nsize);
          index = 0;
          for (i = 0; i < nnodes_sur; i++) {
              vector_data[index] = vector3d_ptr[i].x;
              vector_data[index + 1] = vector3d_ptr[i].y;
              vector_data[index + 2] = vector3d_ptr[i].z;
              index += 3;
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
    }

  status = H5Gclose(group_id);

  /* Add data entries to XML tree */
  node1 = xmlNewNode(NULL, BAD_CAST "Attribute");

  if (data_type == SCALAR_DATA) sprintf(buff1, "Scalar");
  else sprintf(buff1, "Vector");

  xmlNewProp(node1, BAD_CAST "AttributeType", BAD_CAST buff1);
  xmlNewProp(node1, BAD_CAST "Center", BAD_CAST type_string);
  xmlNewProp(node1, BAD_CAST "Name", BAD_CAST data_name);
  sprintf(buff1, "%s:Data/%s/%d", hdf5->h5_surf_filename, data_name, ntime);
  node2 = xmlNewChild(node1, NULL, BAD_CAST "DataItem", BAD_CAST buff1);
  xmlNewProp(node2, BAD_CAST "DataType", BAD_CAST "Float");

  if (data_type == SCALAR_DATA) sprintf(buff1, "%d", nsize);
  else sprintf(buff1, "%d 3", nsize);

  xmlNewProp(node2, BAD_CAST "Dimensions", BAD_CAST buff1);
  xmlNewProp(node2, BAD_CAST "Format", BAD_CAST "HDF");
  xmlNewProp(node2, BAD_CAST "Precision", BAD_CAST "8");
  xmlAddChild(hdf5->current_grid_node_surf, node1);

  if (vector_data != NULL) vector_data = (double *) tl_free(sizeof(double), 3 * nsize, vector_data);
#endif

  return;
}
