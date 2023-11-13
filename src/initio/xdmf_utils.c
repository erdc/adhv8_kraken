/*!
   \file xdmf_utils.c
   \brief Functions to interface to HDF5 calls
 */
#include "global_header.h"
#ifdef _ADH_HDF5
#include "adh_hdf5.h"
#endif

#ifndef ADH_INIT
static int myid;
#endif

/*!
   /brief Create and write an HDF dataset.

   \param group_id containing group for the dataset
   \param name the name of the dataset to be created
   \param dtype_id is the datatype (e.g., H5T_NATIVE_INT)
   \param compression is a flag --- 0 no  compression, 1 turn on compression
   \param ndims is the number of dimensions to the data
   \param dims is a vector with the extent of data in each dimension
   \param data is a pointer to the data to be written

   Returns the dataset identifier.
 */
#ifdef _ADH_HDF5
herr_t xdmf_write_dataset(hid_t group_id, char *name, hid_t dtype_id, int compression, int ndims, hsize_t * dims, void *data)
{

  hid_t dataspace, dataset, dcpl;
  hsize_t chunk[2];
  htri_t avail;
  herr_t status;
  hsize_t chunk_size = 400;

  dataspace = H5Screate_simple(ndims, dims, NULL);
  HDF5_ID(dataspace);

  /* check to see if gzip is available */
  avail = H5Zfilter_avail(H5Z_FILTER_DEFLATE);

  /* turn off compression is calling function has flagged no compression */
  if (compression == 0)
    {
      avail = 0;
    }

  /* set chunk dimensions. If data dimensions are too small, then turn off compression */
  if (ndims == 2)
    {
      chunk[0] = chunk_size / dims[1];
      chunk[1] = dims[1];   /* this should be 3, since all vector data is 3-dimensional */
    }
  else
    {
      chunk[0] = chunk_size;
    }

  if (chunk[0] > dims[0])
    {
      avail = 0;
    }

  if (avail)
    {
      dcpl = H5Pcreate(H5P_DATASET_CREATE);
      HDF5_ID(dcpl);
      status = H5Pset_shuffle(dcpl);    /* use shuffle to try to improve compression */
      HDF5_ERR(status);
      status = H5Pset_deflate(dcpl, 9); /* set compression lvl to max of 9 */
      HDF5_ERR(status);
      status = H5Pset_chunk(dcpl, ndims, chunk);
      HDF5_ERR(status);
    }
  else
    {
      /* printf("Compression not available\n"); */
      dcpl = H5P_DEFAULT;
    }

  dataset = H5Dcreate(group_id, name, dtype_id, dataspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  HDF5_ID(dataset);

  status = H5Dwrite(dataset, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  HDF5_ERR(status);

  status = H5Sclose(dataspace);
  HDF5_ERR(status);
  status = H5Pclose(dcpl);
  HDF5_ERR(status);
  status = H5Dclose(dataset);
  HDF5_ERR(status);
  return status;
}

/*!
   /brief Read an HDF dataset.

   \param group_id containing group for the dataset
   \param name the name of the dataset to be read
   \param dtype_id is the datatype (e.g., H5T_NATIVE_INT)
   \param data is a pointer to the data to be read
 */
herr_t xdmf_read_dataset(hid_t group_id, char *name, hid_t dtype_id, void *data)
{
  hid_t dataset;
  herr_t status;

  dataset = H5Dopen(group_id, name, H5P_DEFAULT);
  HDF5_ID(dataset);
  status = H5Dread(dataset, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  HDF5_ERR(status);
  status = H5Dclose(dataset);
  HDF5_ERR(status);

  return status;
}

/*!
   /brief Read the number of points in a data set.

   The number of points is the number of total number of entries in
   the array of data. For example, if data is a 100x20 array, then
   number of points will be 2000.

   \param group_id containing group for the dataset
   \param name the name of the dataset to be created
 */
hsize_t xdmf_read_npoints(hid_t group_id, char *name)
{
  hid_t dataset, filespace;
  hsize_t npoints;

  dataset = H5Dopen(group_id, name, H5P_DEFAULT);
  HDF5_ID(dataset);
  filespace = H5Dget_space(dataset);
  HDF5_ID(filespace);

  npoints = H5Sget_simple_extent_npoints(filespace);

  HDF5_ERR(H5Sclose(filespace));
  HDF5_ERR(H5Dclose(dataset));
  return npoints;
}

/*!
   /brief Read the extents of the data array

   This determines the dimensionality of the data. For example, if the
   data is a 100x20 array, then dims = {100, 20}. The input array dims
   must be at least as long as the number of dimensions in the data.

   \param group_id containing group for the dataset
   \param name the name of the dataset to be created
   \param dims on return contains a vector with the extent in each dimension
 */
herr_t xdmf_read_extents(hid_t group_id, char *name, hsize_t * dims)
{
  hid_t dataset, filespace;
  herr_t status;

  dataset = H5Dopen(group_id, name, H5P_DEFAULT);
  HDF5_ID(dataset);
  filespace = H5Dget_space(dataset);
  HDF5_ID(filespace);
  status = H5Sget_simple_extent_dims(filespace, dims, NULL);
  HDF5_ERR(status);

  status = H5Sclose(filespace);
  HDF5_ERR(status);
  status = H5Dclose(dataset);
  HDF5_ERR(status);

  return status;
}

/*!
   /brief Determine the number of dimensions in the dataset

   This returns the number of dimensions in the dataset. For example, if the
   data is a 100x20x5 array, then ndims = 3. 

   \param group_id containing group for the dataset
   \param name the name of the dataset to be created
 */
hsize_t xdmf_read_ndims(hid_t group_id, char *name)
{
  hid_t dataset, filespace;
  hsize_t ndims;

  dataset = H5Dopen(group_id, name, H5P_DEFAULT);
  HDF5_ID(dataset);
  filespace = H5Dget_space(dataset);
  HDF5_ID(filespace);
  ndims = H5Sget_simple_extent_ndims(filespace);

  HDF5_ERR(H5Sclose(filespace));
  HDF5_ERR(H5Dclose(dataset));

  return ndims;
}

/*!
   /brief Ascribe an attribute to an existing dataset

   This creates an attribute and writes it to a dataset.

   \param group_id containing group for the dataset
   \param name the name of the dataset to be created
   \param attribute_name the name of the attribute to be created
   \param dtype_id is the datatype (e.g., H5T_NATIVE_INT)
   \param data the value for the attribute
 */
herr_t xdmf_write_attribute(hid_t group_id, char *name, char *attribute_name, hid_t dtype_id, void *data)
{
  hid_t dataset, aid, attr;
  herr_t status;
  hsize_t dims = 1;

  aid = H5Screate_simple(1, &dims, NULL);
  HDF5_ID(aid);
  dataset = H5Dopen(group_id, name, H5P_DEFAULT);
  HDF5_ID(dataset);
  attr = H5Acreate(dataset, attribute_name, dtype_id, aid, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_ID(attr);
  status = H5Awrite(attr, dtype_id, data);
  HDF5_ERR(status);

  status = H5Aclose(attr);
  HDF5_ERR(status);
  status = H5Dclose(dataset);
  HDF5_ERR(status);
  status = H5Sclose(aid);
  HDF5_ERR(status);

  return status;
}

/*!
   /brief Reads an attribute

   \param group_id containing group for the dataset
   \param name the name of the dataset to be created
   \param attribute_name the name of the attribute to be read
   \param dtype_id is the datatype (e.g., H5T_NATIVE_INT)
   \param data the value of the attribute
 */
herr_t xdmf_read_attribute(hid_t group_id, char *name, char *attribute_name, hid_t dtype_id, void *data)
{
  hid_t dataset, filespace, attribute_id;
  herr_t status;

  dataset = H5Dopen(group_id, name, H5P_DEFAULT);
  HDF5_ID(dataset);
  filespace = H5Dget_space(dataset);
  HDF5_ID(filespace);
  attribute_id = H5Aopen(dataset, attribute_name, H5P_DEFAULT);
  HDF5_ID(attribute_id);
  status = H5Aread(attribute_id, dtype_id, data);
  HDF5_ERR(status);

  status = H5Aclose(attribute_id);
  HDF5_ERR(status);
  status = H5Sclose(filespace);
  HDF5_ERR(status);
  status = H5Dclose(dataset);
  HDF5_ERR(status);

  return status;
}
#else
void xdmf_write_dataset()
  {
    return;
  }
void xdmf_read_dataset()
  {
    return;
  }
void xdmf_read_npoints()
  {
    return;
  }
void xdmf_read_extents()
  {
    return;
  }
void xdmf_read_ndims()
  {
    return;
  }
void xdmf_write_attribute()
  {
    return;
  }
void xdmf_read_attribute()
  {
    return;
  }
#endif
