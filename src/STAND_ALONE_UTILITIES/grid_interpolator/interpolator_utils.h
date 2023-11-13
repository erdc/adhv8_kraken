#ifndef READ_WRITE_INPUT_FILES_H
#define READ_WRITE_INPUT_FILES_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

struct grib_record{

    // Lat/Lon grid data
    double lat1;     // Latitude of first point
    double lat2;     // Latitude of last point
    double lon1;     // Longitude of first point
    double lon2;     // Longitude of last point
    int ni;          // Number of points along parallel latitude line
    int nj;          // Number of points along a meridian line
    int num_pts;     // Total number of points in the grid (ni*nj)
    int scan_order;  // Order in which values are listed in the grid file
    double lat_res;  // Constant spacing between latitudes of points
    double lon_res;  // Constant spacing between latitudes of points
    
    // XY grid data
    double x1,x2,y1,y2;
    double dx, dy;

    // Time
    int date;        // Format is:YYYYMMDDHH

    // Variable type
    int var;

    // Data origin
    int center_id;
    int table_version;

    // Type of grid
    int data_rep_type;

    // Grid values
    double *vals;

    // PDS and GDS numbers
    int *PDS, *GDS;
    int vals_length, PDS_length, GDS_length;
    int max_PDS_length, max_GDS_length;

};

// Entry struct functions
void init_entry(struct grib_record **entry);

void free_entry(struct grib_record **entry);

// File I/O functions
FILE* open_input_files(const char *filename);

FILE* open_output_wind_file(const char *filename);

void write_grid_file(FILE *fp, double *x_wind_grid, double *y_wind_grid, int nnodes, int date);

// Read the grib record from inventory file functions
void read_grib_record(struct grib_record *entry, char *rec);

void parse_entry(int *id_array, int *id_len, int id_max_len, const char *in_str);

void check_for_entry_errors(int *PDS, int PDS_length, int *GDS, int GDS_length);

// Read grib values function
void read_entry_values(FILE *fp, struct grib_record *entry);

// Determine the projection to plane from lat/lon
void compute_projection(double ref_lat, double ref_lon, struct grib_record *entry);

// Interpolates the grid
void interpolate_grid(struct grib_record *entry, double *out_grid, double *x, double *y, int nnodes);

// Helper function for projection
double mod_circle(double angle);

// === Functions below are for testing only ===
void write_projection_info(FILE *fp, struct grib_record *entry);

void test_data_read(int *data, int size, int DS_flag);

void print_grib_record(struct grib_record *entry);

void print_human_readable_date(int date);
// === Functions above are for testing only ===
#endif
