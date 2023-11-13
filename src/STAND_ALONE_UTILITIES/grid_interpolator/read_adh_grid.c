#include "global_header.h"
#ifdef _ADH_HDF5
#include "adh_hdf5.h"
#endif

void read_adh_grid(char *adh_root, SGRID **grid) {

#ifdef _ADH_HDF5

#ifdef _MESSG
    debug_initialize(MPI_COMM_WORLD);
#else
    debug_initialize();
#endif

    xmlNodePtr node1 = NULL, node2 = NULL;
    FILE *fp;

    char buff1[MAXLINE] = "";
    char elem_type[20];
    char string[MAXLINE];

    int i;

    SMODEL *mod;
    HDF5 *hdf5 = &(mod->hdf5);
    mod = (SMODEL *) tl_alloc(sizeof(SMODEL), 1);
    /* standard model set up incase results need to get converted*/
    smodel_defaults(mod);
    mod->io = (SIO *) tl_alloc(sizeof(SIO), 1);
    sio_init(mod->io, adh_root);

    // check which geo file exists
    if (doesFileExist(mod->io->geo2d.filename)) {
        open_input_file(&(mod->io->geo2d), "2d geometry file", 0);
        mod->flag.SW2_FLOW = TRUE;
    } else if (doesFileExist(mod->io->geo3d.filename)) {
        open_input_file(&(mod->io->geo3d), "3d geometry file", 0); 
    } else {
        tl_error("2d or 3d geo file needed\n");
    }

    // open face file is 3D
    if (doesFileExist(mod->io->face.filename)) {
        ssw_3d_open_input(mod);
    }

    #ifdef _MESSG
        sgrid_alloc_init(&(mod->grid), mod->io, mod->proc_flag, mod->file_output, mod->model_comm);
    #else
        //sgrid_alloc_init(&(mod->grid), mod->io, mod->proc_flag, mod->file_output);
        sgrid_alloc_init(&(mod->grid), mod->io, 1, mod->file_output);
    #endif
    
    //if(mod->grid->ndim==3) mod->flag.SW3_FLOW=TRUE;
    //SGRID *grid = mod->grid;
    *grid = mod->grid;
    //SIO *io = mod->io;

    ////int nnodes = grid->nnodes;
    //nnodes = &(grid->nnodes);
    //int my_nnode = grid->my_nnodes;
    //int nelems3d = grid->nelems3d;
    //int nelems2d = grid->nelems2d;
    //int nelem;

    //if(grid->ndim == 3){
    //    printf("Error. Input grid is 3d and this code does not support that \
righ//t now. Exiting.\n");
    //    exit(1);
    //}
    //if(grid->ndim == 2){
    //    printf("Good, input grid is 2d.\n");
    //}

    //x = malloc(sizeof(double)*(*nnodes));
    //y = malloc(sizeof(double)*(*nnodes));

    //for(i = 0; i < *nnodes; i++){
    //    x[i] = grid->node[i].x;
    //    y[i] = grid->node[i].y;
    //}
    //for(i = 0; i < *nnodes; i++){
    //    printf("In fct, XY coordinates of node %i are: %f, %f\n",i,x[i],y[i]);
    //}

#endif
}



