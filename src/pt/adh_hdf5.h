/************************************************************************/
/************************************************************************/
#ifndef __ADH_HDF5_H__
#define __ADH_HDF5_H__
#ifdef HDF5_INIT
#define EXTERN_HDF5
#else
#define EXTERN_HDF5 extern
#endif
#ifdef _ADH_HDF5
#include "hdf5.h"
#include "libxml/tree.h"

#define HDF5_NAME_LEN 180 /* the maximum number of characters for filenames */


/* error reporting macros */
#define HDF5_ERR(res) \
{ \
  int _error_code = (int)(res); \
  if (_error_code < 0) \
    { \
      fprintf(stderr, "\nHDF5 call returned error code %d (%s:%d, processor = %d).\n\n", \
        _error_code, __FILE__, __LINE__, myid); \
    } \
}

#define HDF5_ID(id) \
  if ((id) < 0) \
    { \
    fprintf(stderr, "\nHDF5 call returned invalid ID (%s:%d, processor = %d).\n\n", \
      __FILE__, __LINE__, myid); \
    }

typedef struct {
    int current_mesh;   /* counter for the mesh; incremented each time mesh is adapted */
    int current_ts;     /* counter for the time-steps */
    int iprint_step;
    char h5_filename[HDF5_NAME_LEN];
    char xmf_filename[HDF5_NAME_LEN];
    
    hid_t h5_file_ptr;
    hid_t group_mesh;       /* top-level container group for the meshes */
    hid_t group_data;       /* top-level container group for the data */
    hid_t group_nodes;      /* the group containing the nodes */
    hid_t group_elements;   /* the group containing the elements */
    
    xmlDocPtr doc;
    xmlNodePtr root_node;
    xmlNodePtr collection_node;
    xmlNodePtr current_grid_node;
} HDF5;
void ps_print_ts_xdmf(HDF5 *, char *, void *, int, double, int, int, int);
herr_t xdmf_write_dataset(hid_t, char *, hid_t, int, int, hsize_t *, void *);
herr_t xdmf_write_attribute(hid_t, char *, char *, hid_t, void *);

#ifdef _ADH_DEFINE_COPY
/*LP XDMF print control flags */
#define NODE_CENTERED 0
#define ELEM_3D_CENTERED 1
#define ELEM_2D_CENTERED 2

#define SCALAR_DATA 0
#define VECTOR2D_DATA 1
#define VECTOR3D_DATA 2

#define MAXLINE 300
#endif

#endif
#endif
/************************************************************************/
/************************************************************************/
