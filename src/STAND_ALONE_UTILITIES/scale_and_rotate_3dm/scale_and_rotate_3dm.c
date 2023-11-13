#include "global_header.h"

int main(int argc, char *argv[]) {
 
    debug_initialize();

    int i,j;
    double x,y,z,xt,yt,zt;
    
    if (argc > 2 || argc < 1) {
        printf("Usage: ./scale_and_rotate_3dm [adh_file_base]\n");
        exit(0);
    }
    
    /* model base name */
    char adh_root[20];
    strcpy(adh_root, argv[1]);
    
    /* initialize i/o */
    SIO *io = (SIO *) tl_alloc(sizeof(SIO), 1);
    sio_init(io, adh_root);
    assert(io);
    
    /* opens the geo file */
    open_input_file(&(io->geo), "geometry file", 0);
    open_input_file(&(io->face), "boundary face file", 0);
    
    /* count the number of series and strings */
    //read_bc_prep(mod);
    
    /* allocate and initialize grids */
    SGRID *grid;
    sgrid_alloc_init(&grid, io, 0);
    
    // scaling coefficients
    double a_x = 1., b_x = 0.;
    double a_y = 1., b_y = 0.;
    double a_z = 1., b_z = 0.;
    
    //rotation angle
    double theta = 45. * (PI / 180.); // convert to radians
    
    // open new file
    FILE *fp_new_3dm = fopen("scaled_grid.3dm", "w");
    fprintf(fp_new_3dm,"MESH SCALED_AND_ROTATED\n");
    
    // write element connectivities
    for (i=0; i<grid->nelems3d; i++) {
        if (grid->elem3d[i].nnodes == 4) {
            fprintf(fp_new_3dm,"TET\t%d\t",i+1);
        } else if (grid->elem3d[i].nnodes == 6) {
            fprintf(fp_new_3dm,"PRISM\t%d\t",i+1);
        }
        for (j=0; j<grid->elem3d[i].nnodes; j++) {
            fprintf(fp_new_3dm,"%d\t",grid->elem3d[i].nodes[j]+1);
        }
        fprintf(fp_new_3dm,"%d\n",grid->elem3d[i].mat+1);
    }
    
    // write nodes
    for (i=0; i<grid->nnodes; i++) {
        
        // transform
        x = grid->node[i].x;
        y = grid->node[i].y;
        z = grid->node[i].z;
        x = a_x*x + b_x;
        y = a_y*y + b_y;
        z = a_z*z + b_z;
        xt = x*cos(theta) - y*sin(theta);
        yt = x*sin(theta) + y*cos(theta);
        
        // write
        fprintf(fp_new_3dm,"ND\t%d\t%20.10e\t%20.10e\t%20.10e\n",i+1,xt,yt,z);

    }

    return 1;
}
