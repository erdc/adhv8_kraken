#include "global_header.h"
#include "grid_interpolator.h"

//FILE* utility_open_input_file(const char* filename, int is_bin);
//
//FILE* utility_open_output_file(const char* filename, int is_bin);
//
//void interpolate_grid(double *out_grid, SGRID *adh_grid, float *vals,
//        int nnodes, double x1, double x2, double y1, double y2, int ni, int nj);
//
//void utility_write_winds(FILE *fp, double *x_wind, double *y_wind, int nnodes, int time, int is_bin);

int main(int argc, char *argv[]) {

    // Check for correct number of command line arguments
    if (argc < 13 || argc > 13) {
        fprintf(stderr, "\n Incorrect number of arguments.\n");
        exit(0);
    }

    double x1, x2, y1, y2;
    int ni, nj, nnodes = 0;
    int u_first, num_times, i, j, binary_output;
    FILE *time_fp=NULL, *bin_fp=NULL, *out_fp=NULL;
    char *token, in_str[MAXLINE];
    char adh_root[MAXLINE], output_filename[MAXLINE];
    SGRID *adh_grid;

    // ==============================================
    // Read in the grid parameters and descriptors
    // ==============================================
    // Determine which wind direction is first
    u_first = atoi(argv[4]);

    // Get number of times (i.e. number of time steps at which
    // we have wind data)
    num_times = atoi(argv[5]);
    
    // Read in lat/lon values from argv
    x1    = atof(argv[6]);
    x2    = atof(argv[7]);
    y1    = atof(argv[8]);
    y2    = atof(argv[9]);
    ni      = atoi(argv[10]);
    nj      = atoi(argv[11]);

    // Determine output format (binary or text)
    binary_output = atoi(argv[12]);

    // Set number of ERA5 grid points and intialize
    // ERA5 array for storing the values
    int num_grid_pts = ni*nj;
    float era5_grid[num_grid_pts];

    // ==============================================
    // Open files
    // ==============================================
    // Read in AdH grid
    strcpy(adh_root, argv[1]);
    read_adh_grid(adh_root, &adh_grid);

    // 2D case
    if(adh_grid->ndim == 2){
        nnodes = adh_grid->nnodes;
    }
    else if((adh_grid->ndim == 3) && (adh_grid->type == COLUMNAR)){
        nnodes = adh_grid->nnodes_sur;
    }
    else{
        printf("Not ready for non-columnar 3D case.\n");
        exit(1);
    }

    // Open files
    time_fp = utility_open_input_file(argv[2], 0);
    bin_fp  = utility_open_input_file(argv[3], 1);
    strcpy(output_filename,adh_root);
    if(binary_output){
        strcat(output_filename,"_wind.bin");
    }
    else{
        strcat(output_filename,"_wind.txt");
    }
    out_fp  = utility_open_output_file(output_filename, binary_output);

    // Initialize grid for storing interpolated values
    double x_wind[nnodes], y_wind[nnodes];

    // Read in times from time file
    int times_arr[num_times];
    for(i = 0; i < num_times; i++){
        if(fgets(in_str,100,time_fp) == NULL){
            printf("Error. Not enough times in time file %s.\n",argv[2]);
        }
        times_arr[i] = atoi(in_str);
    }

    // Close times file
    fclose(time_fp);

    // ==============================================
    // Read ERA5 values, interpolate and write output to file
    // ==============================================
    int header;
    for(i = 0; i < num_times; i++){
        
        // First wind component (determined by u_first)
        // Read header which contains the size of the grid in bytes
        fread(&header, sizeof(int), 1, bin_fp);

        // Check that array allocated for correct size
        if(header/sizeof(float) != num_grid_pts){
            printf("Warning. ERA5 values array not allocated correctly\n");
            exit(1);
        }

        // Read in first ERA5 grid values
        fread(era5_grid, sizeof(float), num_grid_pts, bin_fp);

        // Read the trailing "header"
        fread(&header, sizeof(int), 1, bin_fp);
        
        // Interpolate to AdH grid
        if(u_first){
            interpolate_grid(x_wind, adh_grid, era5_grid,
                    nnodes, x1, x2, y1, y2, ni, nj);
        }
        else{
            interpolate_grid(y_wind, adh_grid, era5_grid,
                    nnodes, x1, x2, y1, y2, ni, nj);
        }
        
        // Second wind component (determined by u_first)
        // Read header which contains the size of the grid in bytes
        fread(&header, sizeof(int), 1, bin_fp);

        // Check that array allocated for correct size
        if(header/sizeof(float) != num_grid_pts){
            printf("Warning. ERA5 values array not allocated correctly\n");
            exit(1);
        }

        // Read in second ERA5 grid values
        fread(era5_grid, sizeof(float), num_grid_pts, bin_fp);

        // Read the trailing "header"
        fread(&header, sizeof(int), 1, bin_fp);
        
        // Interpolate to AdH grid (reversed from above)
        if(u_first){
            interpolate_grid(y_wind, adh_grid, era5_grid,
                    nnodes, x1, x2, y1, y2, ni, nj);
        }
        else{
            interpolate_grid(x_wind, adh_grid, era5_grid,
                    nnodes, x1, x2, y1, y2, ni, nj);
        }

        // Write output to file
        utility_write_winds(out_fp,x_wind,y_wind,nnodes,times_arr[i],binary_output);
    }

    // Free AdH grid
    sgrid_free(adh_grid);

    // Close remaining files
    fclose(bin_fp);
    fclose(out_fp);

    return 0;
}

FILE* utility_open_input_file(const char* filename, int is_bin){ 
   
    FILE *fp;
    
    if(is_bin){
        // Open binary file
        if((fp = fopen(filename, "rb")) != NULL){
            printf("Opened %s as binary\n",filename);
        }
        else{
            printf("Error, %s not found.\n",filename);
            exit(0);
        }
    }
    else{
        // Open text file
        if((fp = fopen(filename, "r")) != NULL){
            printf("Opened %s\n",filename);
        }
        else{
            printf("Error, %s not found.\n",filename);
            exit(0);
        }
    }

    return fp;
}

FILE* utility_open_output_file(const char* filename, int is_bin){
   
    FILE *fp;
    
    if(is_bin){
        // Open binary file
        if((fp = fopen(filename, "w")) != NULL){
            printf("Opened %s\n",filename);
        }
        else{
            printf("Error, %s not found.\n",filename);
            exit(0);
        }
    }
    else{
        // Open text file
        if((fp = fopen(filename, "w")) != NULL){
            printf("Opened %s\n",filename);
        }
        else{
            printf("Error, %s not found.\n",filename);
            exit(0);
        }
    }

    return fp;
}

// Interpolate to AdH grid
void interpolate_grid(double *out_grid, SGRID *adh_grid, float *vals,
        int nnodes, double x1, double x2, double y1, double y2, int ni, int nj){
    
    int i, nx, ny, id;
    double x1_loc, x2_loc, y1_loc, y2_loc;
    double val_11, val_12, val_21, val_22;
    double fy1, fy2;
    double dx, dy;
    double x, y;

    // Compute dx and dy; the spacing of the flattened ERA5 grid
    dx = (x2-x1)/(ni-1);
    dy = (y1-y2)/(nj-1);

    // Loop over AdH nodes
    for(i=0;i<nnodes;i++){

        // Find coordinates of (surface) node i
        if(adh_grid->ndim == 2){
            x = adh_grid->node[i].x;
            y = adh_grid->node[i].y;
        }
        else if((adh_grid->ndim == 3) && (adh_grid->type == COLUMNAR)){
            id = adh_grid->vertical_list[i]->id;
            x = adh_grid->node[id].x;
            y = adh_grid->node[id].y;
        }
        else{
            printf("Not ready for non-columnar case.\n");
            exit(1);
        }

        // Find 4 point stencil in input grid
        if((x-x1<0) || (y1-y<0)){
            printf("Error on AdH point %f,%f. Point outside domain.\n",x,y);
            exit(1);
        }
        nx = (x-x1)/dx;
        ny = (y1-y)/dy;
        if(nx>(ni-2) || ny>(nj-2)){
            printf("Error on AdH point %f,%f. Point outside domain.\n",x,y);
            exit(1);
        }

        // Store locations
        x1_loc = x1 +  nx    * dx;
        x2_loc = x1 + (nx+1) * dx;
        y1_loc = y1 - (ny+1) * dy;
        y2_loc = y1 -  ny    * dy;

        // Get and store values
        val_11 = vals[(ny+1)*ni + nx];
        val_21 = vals[(ny+1)*ni + (nx+1)];
        val_12 = vals[ny*ni + nx];
        val_22 = vals[ny*ni + (nx+1)];

        // Compute intermediate values for bilinear interpolation
        fy1 = ((x2_loc-x)/(x2_loc-x1_loc))*val_11 + ((x-x1_loc)/(x2_loc-x1_loc))*val_21;
        fy2 = ((x2_loc-x)/(x2_loc-x1_loc))*val_12 + ((x-x1_loc)/(x2_loc-x1_loc))*val_22;

        // Store bilinear interpolated value
        out_grid[i] = ((y2_loc-y)/(y2_loc-y1_loc))*fy1 + ((y-y1_loc)/(y2_loc-y1_loc))*fy2;

        //// Testing
        //printf("X: %f, %f, %f, %f, %i\n",x1_loc, x, x2_loc, dx, nx);
        //printf("Y: %f, %f, %f, %f, %i\n",y1_loc, y, y2_loc, dy, ny);
        //printf("Locs are x1=%f, x2=%f, y1=%f, y2=%f\n",x1_loc,x2_loc,y1_loc,y2_loc);
        //printf("val=%f\n",out_grid[i]);
        //printf("=================\n");
    }
    
    return;
}

void utility_write_winds(FILE *fp, double *x_wind, double *y_wind, int nnodes, int time, int is_bin){

    int j;
    
    // Write output to binary file
    if(is_bin){
        // Write time header
        fwrite(&time, sizeof(int), 1, fp);
        // Write x-component of wind velocity first
        fwrite(x_wind, sizeof(double), nnodes, fp);
        // Write y-component of wind velocity second
        fwrite(y_wind, sizeof(double), nnodes, fp);
    }
    // Write output to text file
    else{
        // write time header
        fprintf(fp,"TS 2 %i\n",time);

        // loop over the nodes of the grid and write x/y wind velocities
        for(j=0; j < nnodes; j++){
            fprintf(fp,"%f %f\n",x_wind[j],y_wind[j]);
        }
    }

    return;
}
