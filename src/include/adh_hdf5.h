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


//EXTERN_HDF5 char h5_filename[HDF5_NAME_LEN];
//EXTERN_HDF5 char xmf_filename[HDF5_NAME_LEN];
//EXTERN_HDF5 hid_t h5_file_ptr;
//EXTERN_HDF5 hid_t group_mesh;   /* top-level container group for the meshes */
//EXTERN_HDF5 hid_t group_data;   /* top-level container group for the data */
//EXTERN_HDF5 hid_t group_nodes;  /* the group containing the nodes */
//EXTERN_HDF5 hid_t group_ghost_nodes;    /* the group containing the ghost nodes */
//EXTERN_HDF5 hid_t group_elements;   /* the group containing the elements */
//EXTERN_HDF5 hid_t group_orig_nd;    /* the group containing the original nodes numbers */
//EXTERN_HDF5 hid_t group_mat_id; /* the group containing the material id */
//EXTERN_HDF5 hid_t group_node_map;   /* the group containing info from node_pair */
//
//EXTERN_HDF5 int current_mesh;   /* counter for the mesh; incremented each time mesh is adapted */
//EXTERN_HDF5 int current_ts;     /* counter for the time-steps */
//EXTERN_HDF5 xmlDocPtr doc;
//EXTERN_HDF5 xmlNodePtr root_node;
//EXTERN_HDF5 xmlNodePtr collection_node;
//EXTERN_HDF5 xmlNodePtr current_grid_node;
//
//#if defined (_ADH_HEAT) || defined (_ADH_OVERLAND)
//EXTERN_HDF5 char h5_surf_filename[HDF5_NAME_LEN];
//EXTERN_HDF5 char xmf_surf_filename[HDF5_NAME_LEN];
//EXTERN_HDF5 hid_t h5_surf_file_ptr;
//EXTERN_HDF5 hid_t group_surf_mesh;  /* top-level container group for the meshes */
//EXTERN_HDF5 hid_t group_surf_data;  /* top-level container group for the data */
//EXTERN_HDF5 hid_t group_surf_nodes; /* the group containing the nodes */
//EXTERN_HDF5 hid_t group_surf_elements;  /* the group containing the elements */
//
//EXTERN_HDF5 xmlDocPtr doc_surf;
//EXTERN_HDF5 xmlNodePtr root_node_surf;
//EXTERN_HDF5 xmlNodePtr collection_node_surf;
//EXTERN_HDF5 xmlNodePtr current_grid_node_surf;
//#endif

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

/* Gajanan gkc adding HDF5 parts. May 2017 */

typedef struct {
    int current_mesh;   /* counter for the mesh; incremented each time mesh is adapted */
    int current_ts;     /* counter for the time-steps */
    int iprint_step;
    char h5_filename[HDF5_NAME_LEN];
    char xmf_filename[HDF5_NAME_LEN];
    hid_t h5_file_ptr;
    hid_t group_mesh;   /* top-level container group for the meshes */
    hid_t group_data;   /* top-level container group for the data */
    hid_t group_nodes;  /* the group containing the nodes */
    hid_t group_ghost_nodes;    /* the group containing the ghost nodes */
    hid_t group_elements;   /* the group containing the elements */
    hid_t group_orig_nd;    /* the group containing the original nodes numbers */
    hid_t group_mat_id; /* the group containing the material id */
    hid_t group_part_id; /* the group containing the partition id */
    hid_t group_node_map;   /* the group containing info from node_map */
    
    xmlDocPtr doc;
    xmlNodePtr root_node;
    xmlNodePtr collection_node;
    xmlNodePtr current_grid_node;
    
    char h5_surf_filename[HDF5_NAME_LEN];
    char xmf_surf_filename[HDF5_NAME_LEN];
    hid_t h5_surf_file_ptr;
    hid_t group_surf_mesh;  /* top-level container group for the meshes */
    hid_t group_surf_data;  /* top-level container group for the data */
    hid_t group_surf_nodes; /* the group containing the nodes */
    hid_t group_surf_ghost_nodes;    /* the group containing the ghost nodes */
    hid_t group_surf_elements;   /* the group containing the elements */
    hid_t group_surf_orig_nd;    /* the group containing the original nodes numbers */
    hid_t group_surf_mat_id; /* the group containing the material id */
    hid_t group_surf_node_map;   /* the group containing info from node_map */
    
    xmlDocPtr doc_surf;
    xmlNodePtr root_node_surf;
    xmlNodePtr collection_node_surf;
    xmlNodePtr current_grid_node_surf;
} HDF5;
void ps_print_ts_xdmf(HDF5 *, char *, void *, int, double, int, int, int);
void ps_print_surf_ts_xdmf(HDF5 *, char *, void *, int, double, int, int, int);
herr_t xdmf_write_dataset(hid_t, char *, hid_t, int, int, hsize_t *, void *);
herr_t xdmf_read_dataset(hid_t group_id, char *, hid_t, void *);
hsize_t xdmf_read_npoints(hid_t, char *);
herr_t xdmf_read_extents(hid_t, char *, hsize_t *);
hsize_t xdmf_read_ndims(hid_t, char *);
herr_t xdmf_write_attribute(hid_t, char *, char *, hid_t, void *);
herr_t xdmf_read_attribute(hid_t group_id, char *, char *, hid_t, void *);


#endif
#endif
/************************************************************************/
/************************************************************************/
