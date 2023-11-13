#include "extrusion.h"

void read_zbins_nodal(SIO *io, SGRID *grid2d, double *head, char *filebase, zbins **bbin);
void read_nodal_bins(SIO *io, char *filebase, int *node_bin);
double create_mat_layers(SGRID *grid2d, int min_bed_layers, int *mat_layers, double *head, int *node_count);
void write_3d_geometry_adh(SIO *io, SGRID *geo3d, SGRID *geo2d);

//------------------------------------------------------------------
//------------------------------------------------------------------

static int mesh_type = TET_MESH;
static int *new_node_number2d;
static int *new_node_number3d;
static int *node_flags2d;
static int *node_flags3d;
static int which_flag;
static int flag_bins = OFF;
static int bin_flag = OFF;
static int flag_min_layer = OFF;
static int nbins_mat;

//------------------------------------------------------------------
//------------------------------------------------------------------

int extrudeAdH(int flag_bins,int flag_min_layer, int minimum_layers, char *adh_root, int mesh_type) {

    int i,j;

    printf("\n*****************************************************************\n");
    printf(  "*************** Beginning ADH 2D To 3D Conversion ***************\n\n");
    
    printf(">> 2D Project: %s\n",adh_root);
    
    if (mesh_type == MIXED_ELEMENT_MESH) {
        printf(">> Creating 3D mixed triangular prism/tetrahedral grid\n");
    } else {
        printf(">> Creating 3D tetrahedral grid\n");
    }
    
    //------------------------------------------------------------------
    //------------------------------------------------------------------
#ifdef _MESSG
    debug_initialize(MPI_COMM_WORLD);
#else    
    debug_initialize();
#endif
    
    int nsubModels = 1, ninterfaces = 0;
    SSUPER_MODEL *superModel = NULL;
    ssuperModel_alloc_init(&superModel, 1, &nsubModels, &ninterfaces);
    printf("\n*****************************************************");
    smodel_init(&(superModel[0].submodel[0]), adh_root);
    printf("*****************************************************\n\n");
    SMODEL *mod2d = &(superModel[0].submodel[0]);
    SGRID *grid2d = mod2d->grid;
    double *head2d = mod2d->sw->d2->head;
    
    //------------------------------------------------------------------
    //------------------------------------------------------------------
    
    zbins *bin = NULL;
    int *node_bin = NULL;
    
    int *mat_layers = (int *) tl_alloc(sizeof(int), grid2d->nelems2d);
    int *node_count = (int *) tl_alloc(sizeof(int), grid2d->nnodes);
    int *node_start = (int *) tl_alloc(sizeof(int), grid2d->nnodes);
    for (i=0; i<grid2d->nnodes; i++) {
        node_count[i] = -1;
        node_start[i] = -1;
    }
    for(i=0; i<grid2d->nmat; i++)  {
        mat_layers[i] = UNSET_INT;
    }
    
    if (flag_bins) {
        node_bin = (int *) tl_alloc(sizeof(int), grid2d->nnodes);
    }
    
    /* Find shallowest region of domain (for minimum # of layers) */
    double dz = 0.;
    if (flag_bins) {
        bin_flag = 1;
        read_zbins_nodal(mod2d->io, grid2d, head2d, adh_root, &bin);
        
        // set the default nodal bins to all 1
        for (i=0; i<grid2d->nnodes; i++) {
            node_bin[i] = 0;
        }
        
        // read if user variation given
        read_nodal_bins(mod2d->io, adh_root, node_bin);
    } else {
        dz = create_mat_layers(grid2d, minimum_layers, mat_layers, head2d, node_count);
    }
    
    /* count the number of nodes in each column */
    double zz, dz_fraction = 0.20; // Adjust so that no small layer exist
    for (i = 0; i < grid2d->nnodes; i++) {
        zz = head2d[i] + grid2d->node[i].z;
        node_count[i] = 1;
        for (j=0; j<10000; j++) {
            if (flag_bins) {
                dz = find_dz_bin(bin, zz, node_bin[i]);
            }
            if (zz - dz <= grid2d->node[i].z + dz_fraction * dz) {
                zz = grid2d->node[i].z;
                node_count[i]++;
                break;
            }
            zz -= dz;
            node_count[i]++;
        }
        //printf("node_count[%d]: %d \n",i,node_count[i]);
    }
    
    
    // create 3D grid
    SGRID *grid3d = (SGRID *) tl_alloc(sizeof(SGRID), 1);
#ifdef _MESSG
    sgrid_init(grid3d,0,NULL); // this should be ran serially 
#else
    sgrid_init(grid3d,0);
#endif
    create_3d_geometry_adh(grid2d, grid3d, node_count, node_start, node_bin, dz, head2d, bin, bin_flag, mesh_type);
    
    
    // write 3D grid
    write_3d_geometry_adh(mod2d->io, grid3d, grid2d);
    
    // create 3D bc
    create_3d_bc_adh(mod2d, grid3d, node_count, node_start, mesh_type, node_flags2d, node_flags3d, new_node_number2d, new_node_number3d);
    
    // create 3D hot
    create_sw3_hot(mod2d, grid3d, node_start, node_count);
    
    // write the output file in XDMF format
    char new_output_root[100];
    strcpy(new_output_root,adh_root);
    strcat(new_output_root, "_3d");
    geo_xdmf_write(new_output_root);
    
    return 0;
}

//------------------------------------------------------------------
//------------------------------------------------------------------
//------------------------------------------------------------------
//------------------------------------------------------------------

void read_zbins_nodal(SIO *io, SGRID *grid2d, double *head, char *filebase, zbins **bbin) {
    
    open_input_file(&(io->extrude_bin), "extrusion bin file", FALSE);
    assert(io && io->extrude_bin.fp);
    
    int i;
    CARD card, subCard;
    char line[MAXLINE];           /* the input line */
    char *data = NULL;
    
    // to make sure the first bin z values starts at top z
    /* Find shallowest region of domain (for minimum # of layers) */
    double max_z = -1e6;
    for (i=0; i<grid2d->nnodes; i++) {
        if (head[i] + grid2d->node[i].z > max_z) max_z = head[i] + grid2d->node[i].z;
    }
    
    // count bins for each nodal material
    nbins_mat = 0;
    int mat = UNSET_INT;
    while(fgets(line, MAXLINE, io->extrude_bin.fp) != NULL) {
        //printf("line: %s\n",line);
        io_save_line(io, io->extrude_bin.fp, io->extrude_bin.filename, line);
        if (strip_comments(line) <= 1) continue;
        card = parse_card(line, &data);
        //printf("card: %d\n",card);
        switch (card) {
            case CARD_BIN:
                //printf("line: %s\n",line);
                mat = read_int_field(*io,&data) - 1;
                if (mat != nbins_mat && mat != nbins_mat-1) {
                    printf("BIN ERROR: Bin nodal material IDS must be listed in ascending order :: last mat: %d new mat: %d\n",mat,nbins_mat);
                    exit(-1);
                }
                if (mat == nbins_mat) nbins_mat++;
                break;
            default:
                break;
        }
    }
    printf(">> total number of nodal bins: %d\n",nbins_mat);
    rewind(io->extrude_bin.fp);
    
    // allocate and initialize nodal bins
    (*bbin) = (zbins *)tl_alloc(sizeof(zbins), nbins_mat);
    zbins *bin = (*bbin); // alias
    for (i=0; i<nbins_mat; i++) {
        bin[i].n = 0;
    }
    
    while(fgets(line, MAXLINE, io->extrude_bin.fp) != NULL) {
        io_save_line(io, io->extrude_bin.fp, io->extrude_bin.filename, line);
        if (strip_comments(line) <= 1) continue;
        card = parse_card(line, &data);
        switch (card) {
            case CARD_BIN:
                mat = read_int_field(*io,&data) - 1;
                (bin[mat].n)++;
                break;
            default:
                break;
        }
    }
    rewind(io->extrude_bin.fp);
    
    for (i=0; i<nbins_mat; i++) {
        //printf("-- bin[%d] nbins: %d\n",i+1,bin[i].n);
        bin[i].dz = (double *)tl_alloc(sizeof(double), bin[i].n);
        bin[i].ztop = (double *)tl_alloc(sizeof(double), bin[i].n);
    }
    
    int bin_id = -1, last_bin_id = -1, old_mat = -1;
    while(fgets(line, MAXLINE, io->extrude_bin.fp) != NULL) {
        io_save_line(io, io->extrude_bin.fp, io->extrude_bin.filename, line);
        if (strip_comments(line) <= 1) continue;
        card = parse_card(line, &data);
        switch (card) {
            case CARD_BIN:
                
                mat = read_int_field(*io,&data) - 1;
                if (mat != old_mat) last_bin_id = -1;
                old_mat = mat;
                bin_id = read_int_field(*io,&data) - 1;
                if(bin_id <= last_bin_id) {
                    printf("BIN ERROR: Bins must be given in ascending order\n");
                    exit(-1);
                }
                
                bin[mat].ztop[bin_id] = read_dbl_field(*io,&data);
                if(bin_id != 0 && (bin[mat].ztop[bin_id] >= bin[mat].ztop[bin_id-1]) ) {
                    printf("BIN ERROR: The z-value for each bin but be descending with ascending ID\n");
                    exit(-1);
                }
                if(bin_id == 0) {
                    //printf("max_z: %10.5f \t bin->ztop[bin_id]: :%10.5f \n",max_z, bin->ztop[bin_id]); exit(-1);
                    bin[mat].ztop[bin_id] += 1e-6;
                    if (max_z > bin[mat].ztop[bin_id]) {
                        printf("BIN ERROR: The first bin top elevation should be zbed + depth = %20.10f :: your top bin (mat: %d): %20.10f\n",max_z,mat,bin[mat].ztop[bin_id]);
                        exit(-1);
                    }
                }
                
                bin[mat].dz[bin_id] = read_dbl_field(*io,&data);
                if (bin[mat].dz[bin_id] <= 0.0) {
                    printf("BIN ERROR: Bin dz values must be >= 0\n");
                    exit(-1);
                }
                
                printf(">> ---- bin[%d] bank: %d :: z: %10.5f dz: %10.5f\n",mat+1,bin_id+1,bin[mat].ztop[bin_id], bin[mat].dz[bin_id]);
                
                last_bin_id = bin_id;
                
                break;
            default:
                break;
        }
    }
    rewind(io->extrude_bin.fp);
    printf("*****************************************************\n\n");
}

//------------------------------------------------------------------
//------------------------------------------------------------------

void read_nodal_bins(SIO *io, char *filebase, int *node_bin) {
    
    open_input_file(&(io->extrude_node), "extrusion node file", FALSE);
    assert(io && io->extrude_node.fp);
    
    CARD card;
    char line[MAXLINE];           /* the input line */
    char *data = NULL;

    int number, bin_ID, max_bins = 0, max_nnodes = 0;

    while(fgets(line, MAXLINE, io->extrude_node.fp) != NULL) {
        io_save_line(io, io->extrude_node.fp, io->extrude_node.filename, line);
        if (strip_comments(line) <= 1) continue;
        //printf("line: %s\n",line);
        card = parse_card(line, &data);
        //printf("card: %d\n",card);
        switch (card) {
            case CARD_ND:
                number = read_int_field(*io,&data) - 1;
                bin_ID = read_int_field(*io,&data) - 1;
                node_bin[number] =  bin_ID;
                if (bin_ID > max_bins) bin_ID = max_bins;
                if (number > max_nnodes) max_nnodes = number;
                break;
            default:
                break;
        }
    }
    max_bins++;
    max_nnodes++;
    
    printf(">> total number of bins from node file: %d\n",max_bins);
    printf(">> total number of nodes from node file: %d\n",max_nnodes);
    
    int i, j;
    for (i=0; i<max_bins; i++) {
        int total_nodes_in_bin = 0;
        for (j=0; j<max_nnodes; j++) {
            if (node_bin[j] == i) total_nodes_in_bin++;
        }
        printf(">> # of nodes in bin %d: %d\n",i+1,total_nodes_in_bin);
    }
    printf("*****************************************************\n\n");
}

//------------------------------------------------------------------
//------------------------------------------------------------------

double find_dz_bin(zbins *bin, double z, int mat) {
    
    int ibin = -1;
    assert(z<=bin[mat].ztop[0]);
    double dz = bin[mat].dz[bin->n-1];
    for (ibin=0; ibin<(bin[mat].n)-1; ibin++) {
        if (z <= bin[mat].ztop[ibin] && z > bin[mat].ztop[ibin+1]) {
            dz = bin[mat].dz[ibin];
            break;
        }
    }
    return dz;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

double create_mat_layers(SGRID *grid2d, int min_bed_layers, int *mat_layers, double *head, int *node_count) {
    
    int i, j;
    
    /* Find shallowest region of domain (for minimum # of layers) */
    double min_depth = 1e6, depth = 0.;
    for (i=0; i<grid2d->nnodes; i++) {
        depth = head[i];
        if (depth < min_depth)  min_depth = depth;
    }
    double dz = min_depth/((double) min_bed_layers);
    
    /* count the number of nodes in each column */
    double zz = 0., surface_elevation = 0.;
    for (i = 0; i < grid2d->nnodes; i++) {
        surface_elevation = head[i] + grid2d->node[i].z;
        zz = surface_elevation;
        node_count[i] = 1;
        for (j=0; j<10000; j++) {
            if (zz - dz <= grid2d->node[i].z) {
                zz = grid2d->node[i].z;
                node_count[i]++;
                break;
            }
            zz -= dz;
            node_count[i]++;
        }
    }
    
    /* set the number of 2d elem bed layers to minimum over 3 nodes */
    int elem2d_layer_min = 1000;
    for (i = 0; i < grid2d->nelems2d; i++) {
        elem2d_layer_min = 1000;
        for (j=0; j<3; j++) {
            if (node_count[grid2d->elem2d[i].nodes[j]]-1 < elem2d_layer_min) elem2d_layer_min = node_count[grid2d->elem2d[i].nodes[j]]-1;
        }
        mat_layers[i] = elem2d_layer_min;
    }
    
    return dz;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

void write_3d_geometry_adh(SIO *io, SGRID *geo3d, SGRID *geo2d) {
    
    int ie,n;
    
    open_output_file(&(io->fout_grid), "3D grid file", FALSE); assert(io && io->fout_grid.fp);
    fprintf(io->fout_grid.fp, "MESH3D\n");
    fprintf(io->fout_grid.fp, "GRID NOCOL\n");
    for(ie = 0; ie < geo3d->nelems3d; ie++) {
        assert(geo3d->elem3d[ie].nnodes == 4 || geo3d->elem3d[ie].nnodes == 6);
        
        if (geo3d->elem3d[ie].nnodes == 4) {
            fprintf(io->fout_grid.fp,"TET  %8d  %8d  %8d  %8d  %8d  %3d  %6d\n",ie + 1,
                    geo3d->elem3d[ie].nodes[0] + 1,geo3d->elem3d[ie].nodes[1] + 1,
                    geo3d->elem3d[ie].nodes[2] + 1,geo3d->elem3d[ie].nodes[3] + 1,
                    geo3d->elem3d[ie].mat + 1, geo3d->elem3d[ie].elem2d_sur + 1);
        } else {
            fprintf(io->fout_grid.fp,"PRISM  %8d  %8d  %8d  %8d  %8d  %8d  %8d  %3d  %6d\n",ie + 1,
                    geo3d->elem3d[ie].nodes[0] + 1,geo3d->elem3d[ie].nodes[1] + 1,
                    geo3d->elem3d[ie].nodes[2] + 1,geo3d->elem3d[ie].nodes[3] + 1,
                    geo3d->elem3d[ie].nodes[4] + 1,geo3d->elem3d[ie].nodes[5] + 1,
                    geo3d->elem3d[ie].mat + 1, geo3d->elem3d[ie].elem2d_sur + 1);
        }
    }
    for(n = 0; n < geo3d->nnodes; n++) {
        fprintf(io->fout_grid.fp,"ND  %8d  %20.10lf  %20.10lf  %20.10lf  %6d\n", n + 1,
                geo3d->node[n].x, geo3d->node[n].y, geo3d->node[n].z, geo3d->node[n].global_surf_id + 1);
    }
    fclose(io->fout_grid.fp);
}


