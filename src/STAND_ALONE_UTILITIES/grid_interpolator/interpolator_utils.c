#include "interpolator_utils.h"
#include "global_header.h"

void read_write_grib_records(char *adh_root, double *x, double* y, int nnodes, double ref_lat, double ref_lon) {

    FILE *inventory_file, *grid_file, *output_file, *projection_info_file;
    char inventory_filename[MAXLINE], grid_filename[MAXLINE];
    char output_filename[MAXLINE], projection_info_filename[MAXLINE];
    char rec[512], str[512];
    int i;
    struct grib_record *current_entry=NULL, *old_entry=NULL, *tmp_entry=NULL;
    double *x_wind_grid, *y_wind_grid;
    bool x_found, y_found;
   
    // Create filenames for input/output
    strcpy(inventory_filename,adh_root);
    strcat(inventory_filename,".inv");
    strcpy(grid_filename,adh_root);
    strcat(grid_filename,".grid");
    strcpy(output_filename,adh_root);
    strcat(output_filename,"_wind.grid");
    strcpy(projection_info_filename,adh_root);
    strcat(projection_info_filename,"_projection_info.txt");

    // Input files for the inventory and data
    inventory_file = io_fopen(inventory_filename, "r", TRUE);
    grid_file = io_fopen(grid_filename, "r", TRUE);
   
    // Output file
    output_file = io_fopen(output_filename, "w", TRUE);
    projection_info_file = io_fopen(projection_info_filename, "w", TRUE);

    // Write header for projection info file
    fprintf(projection_info_file, "Reference latitude and longitude: %f, %f\n", ref_lat, ref_lon);
    
    // Initialize entry data structs
    init_entry(&current_entry);
    init_entry(&old_entry);

    // Allocate memory for the output grid
    x_wind_grid = malloc(sizeof(double)*nnodes);
    y_wind_grid = malloc(sizeof(double)*nnodes);

    int rec_no = 0;
    x_found = false; y_found = false;

    // Loop over the entries of inventory file.
    // Each entry corresponds to a grid of values
    // for a single variable at a single time-step.
    // The spatial grid is assumed to be regular.
    while(fgets(rec, (512), inventory_file) != NULL){

        // Read in record
        rec_no = rec_no+1;

        // Read grib record inventory
        read_grib_record(current_entry, rec);

        // Check if record is wind velocity
        if(current_entry->var == 1 || current_entry->var == 2){
            // Entry is for x-wind velocity
            if(current_entry->var == 1){
                // Check if we have 2 consecutive x-winds with no y-wind.
                // If so, exit. We assume that x- and y-wind grids are
                // provided in alternating fashion.
                if(x_found == true){
                    printf("Exiting. Consecutive x-vel records were \
                            found without corresponding y-vel record. \
                            This code assumes that x and y wind velocities \
                            are present in pairs in the input files.\n"); 
                    exit(1);
                }
                // Interpolate and set x_found to be true
                else{
                    // Read in the grid values
                    read_entry_values(grid_file, current_entry);

                    // Interpolate grid values
                    // Compute the quantities describing the regular CPP projection
                    compute_projection(ref_lat, ref_lon, current_entry);

                    // Interpolate input grid to adh grid
                    interpolate_grid(current_entry, x_wind_grid, x, y, nnodes);
                    x_found = true;
                }
            }
            // Entry is for y-wind velocity
            if(current_entry->var == 2){
                // Check if we have 2 consecutive y-winds with no x-wind.
                // If so, exit. We assume that x- and y-wind grids are
                // provided in alternating fashion.
                if(y_found == true){
                    printf("Exiting. Consecutive y-vel records were \
                            found without corresponding x-vel record. \
                            This code assumes that x and y wind velocities \
                            are present in pairs in the input files.\n");
                    exit(1);
                }
                // Interpolate and set y_found to be true
                else{
                    // Read in the grid values
                    read_entry_values(grid_file, current_entry);

                    // Interpolate grid values
                    // Compute the quantities describing the regular CPP projection
                    compute_projection(ref_lat, ref_lon, current_entry);

                    // Interpolate input grid to adh grid
                    interpolate_grid(current_entry, y_wind_grid, x, y, nnodes);
                    y_found = true;
                }
            }

            // Check if we've already found corresponding wind velocity
            if(x_found && y_found){
                // Check entries for compatibility
                if(current_entry->date != old_entry->date){
                    printf("Error. Date mismatch for pair of x and \
                            y wind velocity wind grids.\n");
                    exit(1);
                }

                // Write to file
                write_grid_file(output_file, x_wind_grid, y_wind_grid, nnodes, current_entry->date);
                
                // Set flags to false
                x_found = false; y_found = false;
            }
            else{
                // We haven't found both pairs of wind velocities
                // Save the entry
                tmp_entry = old_entry;
                old_entry = current_entry;
                current_entry = tmp_entry;
            }
        }
        else{
            // Fast-forward values file
            for(i = 0; i < current_entry->num_pts+1; i++){
                if(fgets(str, 512, grid_file) == NULL){
                    break;
                }
            }
        }
    }

    // Close files
    fclose(inventory_file);
    fclose(grid_file);
    fclose(output_file);
    fclose(projection_info_file);

    // Free memory storing the entry struct
    free_entry(&current_entry);
    free_entry(&old_entry);
    tmp_entry = NULL;

    free(x_wind_grid);
    free(y_wind_grid);

    return;
}

// Entry struct functions ============================
void init_entry(struct grib_record **new_entry){

    int max_PDS_length, max_GDS_length;

    // Max length of the PDS/GDS arrays
    max_PDS_length = 100;
    max_GDS_length = 100;

    (*new_entry) = (struct grib_record*) malloc(sizeof(struct grib_record));

    (*new_entry)->vals = NULL;
    (*new_entry)->PDS  = NULL;
    (*new_entry)->GDS  = NULL;

    (*new_entry)->vals_length = 0;
    (*new_entry)->PDS_length  = 0;
    (*new_entry)->GDS_length  = 0;

    (*new_entry)->max_PDS_length = max_PDS_length;
    (*new_entry)->max_GDS_length = max_GDS_length;

    // Allocate arrays for storing the PDS and GDS numbers of each entry
    // These describe the grid (parameter type, lat/lon, number of points, etc.)
    (*new_entry)->PDS = malloc(max_PDS_length*sizeof(int));
    (*new_entry)->GDS = malloc(max_GDS_length*sizeof(int));
    
    return;
}

void free_entry(struct grib_record **entry){
    
    if((*entry)->vals != NULL) free((*entry)->vals);
    if((*entry)->PDS != NULL)  free((*entry)->PDS);
    if((*entry)->GDS != NULL)  free((*entry)->GDS);

    free(*entry);
    *entry = NULL;

    return;
}

// Probably don't need these file I/O functions anymore
// File I/O functions ================================
FILE* open_input_files(const char *filename){

    FILE *fp = NULL;
    if((fp = fopen(filename, "r"))){
        return fp;
    }
    printf("Error, file %s not found.\n",filename);
    exit(1);
}

FILE* open_output_wind_file(const char *filename){

    FILE *fp = NULL;
    if((fp = fopen(filename, "w"))){
        return fp;
    }else{
        printf("Warning! Problem opening output file!\n");
        exit(1);
    }
}

void write_grid_file(FILE *fp, double *x_wind_grid, double *y_wind_grid, int nnodes, int date){
    int i;

    fprintf(fp, "Date %i\n", date);

    for(i=0;i<nnodes;i++){
        fprintf(fp, "%E %E\n", x_wind_grid[i], y_wind_grid[i]);
    }

    return;
}

// Read the grib record from inv file functions ======
void read_grib_record(struct grib_record *entry, char *rec){

    char *PDS_str, *GDS_str;
    int PDS_size=0, GDS_size=0;
    int *PDS=NULL, *GDS=NULL;
    
    // Parse string to get PDS10 and GDS sections
    PDS_str = strstr(rec, "PDS");
    GDS_str = strstr(rec, "GDS");
    PDS_str = strtok(PDS_str, ":");
    
    // Check if we found reference to GDS
    if(GDS_str!=NULL){
        GDS_str = strtok(GDS_str, ":");
    }

    // Parse the PDS and GDS strings
    parse_entry(entry->PDS, &(entry->PDS_length), entry->max_PDS_length, PDS_str);
    parse_entry(entry->GDS, &(entry->GDS_length), entry->max_GDS_length, GDS_str);

    // Check for errors in the PDS/GDS numbers
    check_for_entry_errors(entry->PDS, entry->PDS_length, entry->GDS, entry->GDS_length);

    // Assign pointers for entry member arrays
    PDS = entry->PDS;
    GDS = entry->GDS;

    // Determine the data origin and process ID
    entry->table_version = entry->PDS[3];
    entry->center_id     = PDS[4];
    
    // Extract date from PDS
    entry->date = PDS[15]+100*PDS[14]+10000*PDS[13]+1000000*(PDS[12]+100*(PDS[24]-1));

    // Determine variable type (e.g. x- or y-component of wind velocity, rain, etc.)
    switch(PDS[8]){
        
        // u-component wind vel
        // Corresponds to x-direction
        case 165:
            entry->var = 1;
            break;

        // v-component wind vel
        // Corresponds to y-direction
        case 166:
            entry->var = 2;
            break;
        
        // Non-wind velocity variable
        default:
            entry->var = 0;
            //printf("Invalid variable number. Setting var to 0.\n");
    }

    // Read spatial grid data
    // Get latitude and longitude bounds
    entry->lat1     = ((GDS[11]<<8) + GDS[12]);
    entry->lon1     = ((GDS[14]<<8) + GDS[15]);
    entry->lat2     = ((GDS[18]<<8) + GDS[19]);
    entry->lon2     = ((GDS[21]<<8) + GDS[22]);
    // Account for the leftmost octet for min/max lat/long
    // min lat
    if(GDS[10] >= 128){
        entry->lat1 = -(entry->lat1+((GDS[10]-128)<<16))/1000.;
    }else{
        entry->lat1 = (entry->lat1+(GDS[10]<<16))/1000.;
    }
    // min long
    if(GDS[13] >= 128){
        entry->lon1 = -(entry->lon1+((GDS[13]-128)<<16))/1000.;
    }else{
        entry->lon1 = (entry->lon1+(GDS[13]<<16))/1000.;
    }
    // max lat
    if(GDS[17] >= 128){
        entry->lat2 = -(entry->lat2+((GDS[17]-128)<<16))/1000.;
    }else{
        entry->lat2 = (entry->lat2+(GDS[17]<<16))/1000.;
    }
    // max long
    if(GDS[20] >= 128){
        entry->lon2 = -(entry->lon2+((GDS[20]-128)<<16))/1000.;
    }else{
        entry->lon2 = (entry->lon2+(GDS[20]<<16))/1000.;
    }
    // Get number of latitude and longitude points
    // TODO: Double-check  GDS[6] and GDS[8].
    // I think these should be bit-shifted by 8.
    entry->ni       = (GDS[6]<<8) + GDS[7];
    entry->nj       = (GDS[8]<<8) + GDS[9];
    entry->num_pts  = (entry->ni)*(entry->nj);

    // Compute latitude and longitude resolution
    // TODO: This is already contained in the PDS. Should add code to
    // read this information from the PDS file and then do an error
    // check since first+last point + ni/nj + resolution contains one
    // piece of redundant information
    entry->lat_res  = fabs(entry->lat1 - entry->lat2)/(entry->nj-1.);
    entry->lon_res  = fabs(entry->lon1 - entry->lon2)/(entry->ni-1.);

    // Scan order
    entry->scan_order = GDS[27];

    // Grid type (lat/long, mercator projection, etc.)
    entry->data_rep_type = GDS[5];
    
    // Allocate memory for grid values array if needed
    if(entry->num_pts > entry->vals_length){
        if(entry->vals_length == 0){
            entry->vals = (double*) malloc(entry->num_pts*sizeof(double));
        }
        else{
            entry->vals = realloc((entry->vals), (entry->num_pts)*sizeof(double));
        }
        entry->vals_length = entry->num_pts;
    }

    return;
}

void parse_entry(int *id_array, int *id_len, int id_max_len, const char *in_str){
    
    char str[512], *token;
    int i;


    // Copy string
    // If NULL, then return with length 0
    if(in_str != NULL){
        strncpy(str, in_str, 512);
    }else{
        *id_len = 0;

        return;
    }

    // Split string at the equals sign
    token = strtok(str, "=");

    // Get first entry of PDS/GDS
    token = strtok(NULL, " ");

    // Loop over remaining space-separated entries of PDS/GDS
    i = 0;
    while(token!=NULL){
        if(i >= id_max_len){
            printf("Warning. PDS or GDS exceeds bounds of array.\n");
            break;
        }
        id_array[i] = atoi(token);
        token = strtok(NULL, " ");
        i = i + 1;
    }
    *id_len = i;

    return;
}

// TODO: This function checks the PDS and GDS for errors
// Currently it only checks that the PDS/GDS numbers are long enough
// and that the data originates from the ECMWF (PDS[4]=98).
// We may wish to add more checks in the future or allow for data
// from other sources.
// This function should evolve to reflect that.
void check_for_entry_errors(int *PDS, int PDS_length, int *GDS, int GDS_length){
    
    // Check that there are enough entries in PDS/GDS id's
    if(PDS_length < 25 || GDS_length < 28){
        printf("Error. PDS or GDS id not long enough. Exiting.\n");
        exit(1);
    } 
    
    // Check for correct parameter table version and originating center
    if( (PDS[3] != 128) || (PDS[4] != 98) ){
        printf("This utility is only set up to hand data from the ECMWF.\n");
        printf("It seems that the input file contains data from another source.\n");
        printf("Stopping PDS read and exiting.\n");
        exit(1);
    }

    return;
}

// Read grib values function ========================
void read_entry_values(FILE *fp, struct grib_record *entry){
    int i = 0, grid_ni, grid_nj;
    double *grid_vals = NULL;
    char str[512], *token;

    // Pull first line of the grid values entry
    fgets(str, 512, fp);

    // Check that our ni*nj values match that of the grid
    grid_ni = atoi(strtok(str, " "));
    grid_nj = atoi(strtok(NULL, " "));
    if(grid_ni*grid_nj != entry->num_pts){
        printf("Warning! Number of points in entry of grid file \
does not match number of points in inventory file.\n");
        exit(1);
    }

    // Pull the num_pts grid values
    while(i<entry->num_pts){
        fgets(str, 512, fp);
        entry->vals[i] = atof(str);
        i = i+1;
    }

    return;
}

// Determine projection to plane ====================
void compute_projection(double ref_lat, double ref_lon, struct grib_record *entry){
    
    const double Rearth = 6378206.4;
    double phi0, lambda0; // Reference latitude/longitude at center of projection
    double phi, lambda;   // Latitude/longitude of first grid point in radians
    double lambda_split;

    // Set values in radians
    phi0    = ref_lat*(M_PI/180.);                  // Latitude of reference point in radians
    lambda0 = mod_circle(ref_lon*(M_PI/180.));      // Longitude of reference point in radians
    
    // Location of split determines negative or positive x-values
    lambda_split = mod_circle(lambda0+M_PI);

    // Convert lat/lon starting point from degrees to radians
    phi     = entry->lat1*(M_PI/180.);              // Latitude of first grid point in radians
    lambda  = mod_circle(entry->lon1*(M_PI/180.));  // Longitude of first grid point in radians

    if(entry->scan_order == 0){
        entry->dx = (M_PI/180.) * Rearth * entry->lon_res * cos(phi0);
        entry->dy = (M_PI/180.) * Rearth * entry->lat_res;

        entry->x1 = Rearth * cos(phi0) * (mod_circle(lambda-lambda_split) - M_PI);
        entry->y1 = Rearth * (phi - phi0);
        
        entry->x2 = entry->x1 + (entry->ni-1)*entry->dx;
        entry->y2 = entry->y1 - (entry->nj-1)*entry->dy;
    }else{
        printf("Scan order %i not valid at this time\n",entry->scan_order);
        exit(1);
    }
    
    return;
}

// Interpolate to AdH grid
void interpolate_grid(struct grib_record *entry, double *out_grid, double *x, double *y, int nnodes){
    int i, nx, ny;
    double x1_loc, x2_loc, y1_loc, y2_loc;
    double val_11, val_12, val_21, val_22;
    double fy1, fy2;

    if(entry->scan_order == 0){
        for(i=0;i<nnodes;i++){
            // Find 4 point stencil in input grid
            if((x[i]-entry->x1<0) || (entry->y1-y[i]<0)){
                printf("Error on AdH point %f,%f. Point outside top left bounds of domain.\n",x[i],y[i]);
            }
            nx = (x[i]-entry->x1)/entry->dx;
            ny = (entry->y1-y[i])/entry->dy;
            if(nx>(entry->ni-2) || ny>(entry->nj-2)){
                printf("Error on AdH point %f,%f. One of nx=%i and ny=%i are \
invalid. Is your point inside the domain of the input grid?\n",x[i],y[i],nx,ny);
                //exit(1);
            }

            // Store locations
            x1_loc = entry->x1 +  nx    * entry->dx;
            x2_loc = entry->x1 + (nx+1) * entry->dx;
            y1_loc = entry->y1 - (ny+1) * entry->dy;
            y2_loc = entry->y1 -  ny    * entry->dy;

            // Get and store values
            val_11 = entry->vals[(ny+1)*entry->ni + nx];
            val_21 = entry->vals[(ny+1)*entry->ni + (nx+1)];
            val_12 = entry->vals[ny*entry->ni + nx];
            val_22 = entry->vals[ny*entry->ni + (nx+1)];

            // Compute intermediate values for bilinear interpolation
            fy1 = ((x2_loc-x[i])/(x2_loc-x1_loc))*val_11 + ((x[i]-x1_loc)/(x2_loc-x1_loc))*val_21;
            fy2 = ((x2_loc-x[i])/(x2_loc-x1_loc))*val_12 + ((x[i]-x1_loc)/(x2_loc-x1_loc))*val_22;

            // Store bilinear interpolated value
            out_grid[i] = ((y2_loc-y[i])/(y2_loc-y1_loc))*fy1 + ((y[i]-y1_loc)/(y2_loc-y1_loc))*fy2;

            //// Testing
            //printf("X: %f, %f, %f, %f, %i\n",x1_loc, x[i], x2_loc, entry->dx, nx);
            //printf("Y: %f, %f, %f, %f, %i\n",y1_loc, y[i], y2_loc, entry->dy, ny);
            ////printf("Locs are x1=%f, x2=%f, y1=%f, y2=%f\n",x1_loc,x2_loc,y1_loc,y2_loc);
            //printf("=================\n");
        }
    }
    else{
        printf("Scan order %i not valid at this time\n",entry->scan_order);
        exit(1);
    }
    return;
}

// Helper function for projection
double mod_circle(double angle){
    return angle-((2*M_PI)*floor(angle/(2*M_PI)));
}

// Functions for testing below ==================
void write_projection_info(FILE *fp, struct grib_record *entry){
    fprintf(fp, "Date: %i\n", entry->date);
    fprintf(fp, "Variable: %i\n", entry->var);
    fprintf(fp, "First point: %f, %f\n",entry->x1, entry->y1);
    fprintf(fp, "dx, dy: %f, %f\n",entry->dx, entry->dy);
    fprintf(fp, "nx, ny: %i, %i\n",entry->ni, entry->nj);
    fprintf(fp,"\n");
    
    return;
}

void test_data_read(int *data, int size, int DS_flag){
    int j = 0;
    
    if(DS_flag == 0){
        printf("PDS10= ");
    }else{
        printf("GDS10= ");
    }
    while(j<size){
        printf("%i ",data[j]);
        j = j+1;
    }
    printf("\n");
    return;
}

void print_grib_record(struct grib_record *entry){

    int i;
    
    printf("Latitude bounds are: %f, %f\n",entry->lat1,entry->lat2);
    printf("Number of point in North-South direction and latitude \
resolution is: %i, %f\n", entry->ni, entry->lat_res);
    printf("Longitude bounds are: %f, %f\n",entry->lon1,entry->lon2);
    printf("Number of point in East-West direction and longitude \
resolution is: %i, %f\n", entry->nj, entry->lat_res);
    printf("Scanning order is %i\n",entry->scan_order);

    printf("Raw date is: %i\n",entry->date);

    print_human_readable_date(entry->date);

    printf("Grid values:\n");
    for(i = 0; i<(entry->ni*entry->nj); i++){
        printf("%E\n",entry->vals[i]);
    }

    printf("XY grid data in km:\n");
    printf("x/y resolution: %f, %f\n", (entry->dx/1000.), (entry->dy/1000.));
    printf("min x/y: %f, %f\n",(entry->x1/1000.),(entry->y1/1000.));
    printf("max x/y: %f, %f\n",(entry->x2/1000.),(entry->y2/1000.));

    return;
}

void print_human_readable_date(int date){

    // Break up the date
    int year = (date/1000000);
    int month = ((date-year*1000000)/10000);
    int day = ((date-year*1000000-month*10000)/100);
    int hour = date-year*1000000-month*10000-day*100;
    
    printf("Date is: %i/%i/%i at hour %i\n", day, month, year, hour);

    return;
}

void test_interpolate_grid(double *out_grid, int nnodes, double *x, double *y, double x1, double y1, double x2, double y2, int nx, int ny, double m1, double m2){
    
    struct grib_record *test_entry;
    int i,j;
    
    FILE *fp;
    char filename[512];
    double eps = 0.001, err;
    
    test_entry = (struct grib_record*) malloc(sizeof(struct grib_record));

    strcpy(filename,"test_interpolation.txt");

    fp = fopen(filename, "w");

    // Create linearly increasing input grid
    test_entry->x1 = x1;
    test_entry->y1 = y2;
    test_entry->ni = nx;
    test_entry->nj = ny;
    test_entry->dx = (x2-x1)/(nx-1);
    test_entry->dy = (y2-y1)/(ny-1);
    test_entry->scan_order = 0;

    test_entry->num_pts = nx*ny;

    printf("Entry starting point: %f, %f\n",test_entry->x1, test_entry->y1);
    printf("Entry num points: %i, %i\n",test_entry->ni, test_entry->nj);
    printf("Entry dx, dy: %f, %f\n",test_entry->dx, test_entry->dy);

    // Create artificial grid
    test_entry->vals = (double*) malloc(test_entry->num_pts*sizeof(double));
    for(j=0;j<ny;j++){
        for(i=0;i<nx;i++){
            //test_entry->vals[j*nx+i] = i*test_entry->dx*m1 + (ny-j-1)*test_entry->dy*m2;
            test_entry->vals[j*nx+i] = (x1+i*test_entry->dx)*m1 + (y2-j*test_entry->dy)*m2;
        }
    }
    fprintf(fp,"Input grid values going from NE to SW:\n");
    for(i=0;i<test_entry->num_pts;i++){
        fprintf(fp,"%i, %f\n",i,test_entry->vals[i]);
    }

    interpolate_grid(test_entry, out_grid, x, y, nnodes);

    fprintf(fp,"Output grid values:\n");
    for(i=0;i<nnodes;i++){
        fprintf(fp,"%i: %f, %f: %f\n",i,x[i],y[i],out_grid[i]);
        err = fabs((x[i]+y[i])-out_grid[i]);
    }

    fclose(fp);

    return;
}

void test_compute_projection(double ref_lat, double ref_lon){
    
    struct grib_record *test_entry;
    int i,j;
    test_entry = (struct grib_record*) malloc(sizeof(struct grib_record));
    FILE *fp; char filename[512];

    strcpy(filename,"test_projection.txt");

    fp = fopen(filename, "w");

    // Set to some initial values
    test_entry->lat1 = 0;
    test_entry->lon1 = 0;
    test_entry->lat_res = 0.1;
    test_entry->lon_res = 0.1;
    test_entry->scan_order = 0;

    fprintf(fp,"Reference lat/lon: %f, %f\n\n",ref_lat,ref_lon);

    fprintf(fp,"Testing starting point projection:\n");
    fprintf(fp,"==============\n");
    for(i=0;i<37;i++){
        for(j=0;j<19;j++){
            test_entry->lat1 = -90.+10.*j;
            test_entry->lon1 = -180+10.*i;
            
            compute_projection(ref_lat,ref_lon,test_entry);

            fprintf(fp,"First point lat/lon: %f, %f\n",test_entry->lat1,test_entry->lon1);
            fprintf(fp,"First point X/Y: %f, %f\n",test_entry->x1,test_entry->y1);
            fprintf(fp,"==============\n");
        }
    }
    
    fprintf(fp,"==============\n");

    fprintf(fp,"Testing dx/dy projection:\n");
    fprintf(fp,"==============\n");
    for(i=0;i<10;i++){
        for(j=0;j<10;j++){
            test_entry->lat_res = 10./(2.*(j+1));
            test_entry->lon_res = 10./(2.*(i+1));
            
            compute_projection(ref_lat,ref_lon,test_entry);

            fprintf(fp,"Resolution lat/lon: %f, %f\n",test_entry->lat_res,test_entry->lon_res);
            fprintf(fp,"Resolution X/Y: %f, %f\n",test_entry->dx,test_entry->dy);
            fprintf(fp,"==============\n");
        }
    }
    fclose(fp);

    return;

}

