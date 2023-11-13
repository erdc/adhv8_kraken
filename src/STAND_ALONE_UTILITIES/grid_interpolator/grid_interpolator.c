#include "global_header.h"
#include "grid_interpolator.h"

int main(int argc, char *argv[]) {

    if (argc < 4 || argc > 4) {
        fprintf(stderr, "\n Wrong number of arguments given \n");
        exit(0);
    }

    char adh_root[MAXLINE];
    SGRID *adh_grid;
    int i;
    double *x, *y, ref_lat, ref_lon;
    int nnodes = 0;
    
    strcpy(adh_root, argv[1]);
    ref_lat = atof(argv[2]);
    ref_lon = atof(argv[3]);
   
    // Read in AdH grid to interpolate input values onto
    read_adh_grid(adh_root, &(adh_grid));

    // 2D case
    if(adh_grid->ndim == 2){

        printf("2D case\n");
        
        // Number of nodes in AdH grid
        nnodes = adh_grid->nnodes;

        // Get x,y coordinates of the nodes
        x = malloc(sizeof(double)*nnodes);
        y = malloc(sizeof(double)*nnodes);
        
        for(i=0;i<nnodes;i++){
            x[i] = adh_grid->node[i].x;
            y[i] = adh_grid->node[i].y;
        }
    }    
    // 3D columnar case
    else if((adh_grid->ndim == 3) && (adh_grid->type == COLUMNAR)){

        printf("3D columnar case\n");
        
        // Number of surface nodes in AdH grid
        nnodes = adh_grid->nnodes_sur;

        // Get x,y coordinates of the surface nodes
        x = malloc(sizeof(double)*nnodes);
        y = malloc(sizeof(double)*nnodes);

        ID_LIST_ITEM *ptr;

        for(i=0;i<nnodes;i++){
            // Get the first node in the vertical list (surface nodes)
            ptr = adh_grid->vertical_list[i];
            x[i] = adh_grid->node[ptr->id].x;
            y[i] = adh_grid->node[ptr->id].y;
        }
    }
    // General 3D non-columnar case
    else{
        printf("Not ready for general 3D non-columnar case.\n");
        return 1;
    }

    read_write_grib_records(adh_root, x, y, nnodes, ref_lat, ref_lon);
    
    free(x);
    free(y);

    return 0;
}

