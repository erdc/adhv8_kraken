/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates, initializes and reads in finite element geometry files.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] pgrid            (SGRID *)  a pointer to an AdH grid struture
 * @param[in]  io                   (SIO *) a pointer to an AdH structure for holding input/output files
 * @param[in]  flag               (int) a flag indicating =0 processors get divided in superfile read or =1 no division of processors necessary
 * @param[in]  file_output     (SFILE_OUTPUT) An AdH file output pointer
 * @param[in]  model_comm       (MPI_Comm)   an AdH model MPI communicator
 *
 * \note CJT \:: NOTE :: if flag = 0, then myid = -1 and npes = 0!
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// headers
#include "global_header.h"

// file scope variables
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;

// prototypes
int get_macro_nelems2d_bed(SGRID *grid);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void sgrid_alloc_init(SGRID **pgrid, SIO *io, int flag, SFILE_OUTPUT file_output
#ifdef _MESSG
                      ,MPI_Comm model_comm
#endif
) {
    
    int myid_model = UNSET_INT, myid_world = UNSET_INT, ierr = UNSET_INT;
    
#ifdef _DEBUG
    printf("\n-------------------------------------------------------\n");
#ifdef _MESSG
    ierr = MPI_Comm_rank(model_comm, &myid_model);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid_world);
    printf("CSTORM COMM PE: %d :: MODEL PE: %d :: FLAG: %d :: Allocating, initializing and reading model grid\n",myid_world,myid_model,flag);
#else
    printf("Allocating, initalizing and reading model grid\n");
#endif
    if (DEBUG_WITH_PICKETS == ON) DEBUG = ON;
#endif
    
    int inode, ie, i;
    double mesh_vol = 0.;
    SGRID *grid = NULL;
    char *filename = NULL;
    
    /* check inputs */
    assert(io);
    
    (*pgrid) = (SGRID *) tl_alloc(sizeof(SGRID), 1);
    grid = (*pgrid);  // alias
#ifdef _MESSG
    sgrid_init(grid,flag,model_comm);
#else
    sgrid_init(grid,flag);
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /* allocate and define edges maps ------------------------ */
    // cjt :: the ordering here has to agree with order of basis function order
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    /* the nodes on an edge of a 2d element */
    grid->nd_on_TriEdge = (int **) tl_alloc(sizeof(int *), 3);
    for (i=0; i<3; i++) {
        grid->nd_on_TriEdge[i] = (int *) tl_alloc(sizeof(int), 2);
    }
    grid->nd_on_TriEdge[0][0] = 0; grid->nd_on_TriEdge[0][1] = 1;
    grid->nd_on_TriEdge[1][0] = 1; grid->nd_on_TriEdge[1][1] = 2;
    grid->nd_on_TriEdge[2][0] = 2; grid->nd_on_TriEdge[2][1] = 0;
    
    grid->nd_on_QuadEdge = (int **) tl_alloc(sizeof(int *), 4);
    for (i=0; i<4; i++) {
        grid->nd_on_QuadEdge[i] = (int *) tl_alloc(sizeof(int), 2);
    }
    grid->nd_on_QuadEdge[0][0] = 0; grid->nd_on_QuadEdge[0][1] = 1;
    grid->nd_on_QuadEdge[1][0] = 1; grid->nd_on_QuadEdge[1][1] = 2;
    grid->nd_on_QuadEdge[2][0] = 2; grid->nd_on_QuadEdge[2][1] = 3;
    grid->nd_on_QuadEdge[3][0] = 3; grid->nd_on_QuadEdge[3][1] = 0;
    
    /* the nodes on an edge of a 3d element */
    grid->nd_on_PrismEdge = (int **) tl_alloc(sizeof(int *), 9);
    for (i=0; i<9; i++) {
        grid->nd_on_PrismEdge[i] = (int *) tl_alloc(sizeof(int), 2);
    }
    grid->nd_on_PrismEdge[0][0] = 0; grid->nd_on_PrismEdge[0][1] = 1;
    grid->nd_on_PrismEdge[1][0] = 1; grid->nd_on_PrismEdge[1][1] = 2;
    grid->nd_on_PrismEdge[2][0] = 2; grid->nd_on_PrismEdge[2][1] = 0;
    grid->nd_on_PrismEdge[6][0] = 3; grid->nd_on_PrismEdge[6][1] = 4;
    grid->nd_on_PrismEdge[7][0] = 4; grid->nd_on_PrismEdge[7][1] = 5;
    grid->nd_on_PrismEdge[8][0] = 5; grid->nd_on_PrismEdge[8][1] = 3;
    grid->nd_on_PrismEdge[3][0] = 0; grid->nd_on_PrismEdge[3][1] = 3;
    grid->nd_on_PrismEdge[4][0] = 1; grid->nd_on_PrismEdge[4][1] = 4;
    grid->nd_on_PrismEdge[5][0] = 2; grid->nd_on_PrismEdge[5][1] = 5;
    
    /* the nodes on an edge of a 3d element */
    grid->nd_on_TetEdge = (int **) tl_alloc(sizeof(int *), 6);
    for (i=0; i<6; i++) {
        grid->nd_on_TetEdge[i] = (int *) tl_alloc(sizeof(int), 2);
    }
    grid->nd_on_TetEdge[0][0] = 0; grid->nd_on_TetEdge[0][1] = 1;
    grid->nd_on_TetEdge[1][0] = 0; grid->nd_on_TetEdge[1][1] = 2;
    grid->nd_on_TetEdge[2][0] = 0; grid->nd_on_TetEdge[2][1] = 3;
    grid->nd_on_TetEdge[3][0] = 1; grid->nd_on_TetEdge[3][1] = 2;
    grid->nd_on_TetEdge[4][0] = 1; grid->nd_on_TetEdge[4][1] = 3;
    grid->nd_on_TetEdge[5][0] = 2; grid->nd_on_TetEdge[5][1] = 3;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
#ifdef _DEBUG
    if (io->geo2d.fp != NULL) {
        filename = io->geo2d.filename;
    } else if (io->geo3d.fp != NULL) {
        filename = io->geo3d.filename;
    }
#endif
    
    if (flag == 1) {
#ifdef _DEBUG
        if (DEBUG) printf("-- MPI_COMM_WORLD PE: %d :: MODEL PE: %d :: FLAG: %d :: FILENAME: %s CALLING READ_GEO\n",myid_world,myid_model,flag,filename);
#endif
        grid->ndim = read_geo(io, grid);
#ifdef _DEBUG
        if (DEBUG) printf("-- MPI_COMM_WORLD PE: %d :: MODEL PE: %d :: FLAG: %d :: FILENAME: %s DONE CALLING READ_GEO\n",myid_world,myid_model,flag,filename);
        if (DEBUG_WITH_PICKETS) tl_check_all_pickets(__FILE__, __LINE__);
#endif
    } else {
#ifdef _MESSG
#ifdef _DEBUG
        if (DEBUG)  printf("-- MPI_COMM_WORLD PE: %d :: MODEL PE: %d :: FLAG: %d :: FILENAME: %s CALLING READ INTERFACE\n",myid_world,myid_model,flag,filename);
#endif
        grid->ndim = read_interface_geo(io,grid);
        
#ifdef _DEBUG
        if (DEBUG) printf("-- MPI_COMM_WORLD PE: %d :: MODEL PE: %d :: FLAG: %d :: FILENAME: %s DONE CALLING READ INTERFACE\n",myid_world,myid_model,flag,filename);
        if (DEBUG_WITH_PICKETS) tl_check_all_pickets(__FILE__, __LINE__);
#endif
#endif
    }
    
#ifdef _MESSG
    grid->smpi->partition_info = (int *) tl_alloc(sizeof(int), grid->max_nnodes);
    for (i=0; i<grid->nnodes; i++){
        grid->smpi->partition_info[i] = grid->node[i].resident_pe;
    }
#endif
    
    /* check grid integrity */
    if (grid->nnodes <= 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> There are no nodes on this 2d grid.");
    }
    if (grid->ndim == 2 && grid->nelems2d <= 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> There are no 2d elements on this 2d grid.");
    }
    if (grid->ndim == 3 && grid->nelems3d <= 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> There are no 3d elements on this 3d grid.");
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /* initialize struct variables */
    grid->mesh_volume = 0.;
    grid->initial_nnodes = grid->nnodes;
    grid->orig_initial_nnodes = grid->initial_nnodes;
    if (grid->ndim == 2) {
        grid->initial_nelems = grid->nelems2d;
    }
    else {
        grid->initial_nelems = grid->nelems3d;
    }
    grid->orig_initial_nelems = grid->initial_nelems;
    grid->nnodes_matrix = 0;
    
    grid->max_nnodes = grid->nnodes;
    grid->max_nelems1d = grid->nelems1d;
    grid->max_nelems2d = grid->nelems2d;
    grid->max_nelems3d = grid->nelems3d;
    grid->nnodes_old = grid->nnodes;
    grid->nelems2d_old = grid->nelems2d;
    grid->nelems3d_old = grid->nelems3d;
    
    /* grid bounding rectangle */
    grid->nmat = 0;
    grid->x_min = 9.9e+99;
    grid->x_max = -9.9e+99;
    grid->y_min = 9.9e+99;
    grid->y_max = -9.9e+99;
    
    /* find bounding box */
    for (inode = 0; inode < grid->nnodes; inode++) {
        grid->x_min = MIN(grid->node[inode].x, grid->x_min);
        grid->x_max = MAX(grid->node[inode].x, grid->x_max);
        grid->y_min = MIN(grid->node[inode].y, grid->y_min);
        grid->y_max = MAX(grid->node[inode].y, grid->y_max);
    }
    
    /* flag a wet/dry element */
    if (grid->ndim == 2 && grid->nelems2d > 0) {
        grid->wd_flag = (int *) tl_alloc(sizeof(int), grid->nelems2d);
        for (i=0; i<grid->nelems2d; i++) grid->wd_flag[i] = 0;
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /* if the grid is 3d, build columns */
    if (grid->ndim == 3) {
        if (grid->type == COLUMNAR) {
            if (grid->smpi->myid == 0) printf("BUILDING COLUMNS IN COLUMNAR GRID\n");
            build_columns(grid, 1);
        } else {
            classify_2d_elements(pgrid);
        }
        
        // flag nodes as surface, bed or sidewall
        for (i=0; i<grid->nelems2d; i++) {
            if (grid->elem2d[i].bflag == 0) {
                grid->node[grid->elem2d[i].nodes[0]].bflag = 0;
                grid->node[grid->elem2d[i].nodes[1]].bflag = 0;
                grid->node[grid->elem2d[i].nodes[2]].bflag = 0;
            } else if (grid->elem2d[i].bflag == 1) {
                grid->node[grid->elem2d[i].nodes[0]].bflag = 1;
                grid->node[grid->elem2d[i].nodes[1]].bflag = 1;
                grid->node[grid->elem2d[i].nodes[2]].bflag = 1;
            } else {
                //grid->node[grid->elem2d[i].nodes[2]].bflag = 2;
                //grid->node[grid->elem2d[i].nodes[2]].bflag = 2;
                //grid->node[grid->elem2d[i].nodes[2]].bflag = 2;
            }
        }
         
        grid->nnodes_sur_orig = grid->nnodes_sur;
        grid->max_nnodes_sur = grid->nnodes_sur;
        grid->max_nnodes_bed = grid->nnodes_bed;
        grid->old_max_nnodes_sur = grid->nnodes_sur;
        grid->old_max_nnodes_bed = grid->nnodes_bed;
        
        // find the total number of 2d bed/surface elements on the global grid
        grid->macro_nelems2d_bed = get_macro_nelems2d_bed(grid);
        
#ifdef _MESSG
        if (grid->type == COLUMNAR) {
            grid->smpi->surface_partition_info = (int *)tl_alloc(sizeof(int), grid->nnodes_sur_orig);
            grid->old_global_surf = (int *)tl_alloc(sizeof(int), grid->nnodes_sur_orig);
            grid->old_global_bed = (int *)tl_alloc(sizeof(int), grid->nnodes_sur_orig);
            for (i=0; i <  grid->nnodes_sur; i++){
                grid->smpi->surface_partition_info[i] = grid->node[grid->nodeID_2d_to_3d_sur[i]].resident_pe;
                grid->node[grid->nodeID_2d_to_3d_bed[i]].global_bed_id = grid->node[grid->nodeID_2d_to_3d_sur[i]].global_surf_id;
                grid->old_global_surf[i]=grid->node[grid->nodeID_2d_to_3d_sur[i]].global_surf_id;
                grid->old_global_bed[i]=grid->node[grid->nodeID_2d_to_3d_bed[i]].global_bed_id;
            }
        }
        grid->orig_macro_nnodes = grid->macro_nnodes;
        grid->orig_macro_nnodes_sur = grid->macro_nnodes_sur;
        grid->orig_macro_nnodes_bed = grid->macro_nnodes_bed;
#else
        // These should be set in read_geo_mpi for multi processor runs
        grid->macro_nnodes_sur = grid->nnodes_sur;
        grid->macro_nnodes_bed = grid->nnodes_bed;
        grid->macro_nnodes = grid->nnodes;
        grid->orig_macro_nnodes = grid->macro_nnodes;
        grid->orig_macro_nnodes_sur = grid->macro_nnodes_sur;
        grid->orig_macro_nnodes_bed = grid->macro_nnodes_bed;
        if (grid->type == COLUMNAR) {
            grid->old_global_surf = (int *)tl_alloc(sizeof(int), grid->nnodes_sur_orig);
            grid->old_global_bed = (int *)tl_alloc(sizeof(int), grid->nnodes_sur_orig);
            /* Gajanan gkc - Shouldn't there be nnodes_bed_orig as well above? */
            /*mwf add separate loop for bed*/
            
            for(i=0; i <  grid->nnodes_sur; i++){
                grid->node[grid->nodeID_2d_to_3d_sur[i]].global_surf_id = i;
                grid->old_global_surf[i]=grid->node[grid->nodeID_2d_to_3d_sur[i]].global_surf_id;
            }
            for(i=0; i <  grid->nnodes_bed; i++){
                grid->node[grid->nodeID_2d_to_3d_bed[i]].global_bed_id = i;
                grid->old_global_bed[i]=grid->node[grid->nodeID_2d_to_3d_bed[i]].global_bed_id;
            }
        }
#endif
        /* sanity checking */
        if (grid->type == COLUMNAR) {
            if (grid->nelems2d_sur != grid->ncolumns) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> Number of 2d surace elements should be equal to number of columns.");
            }
            if (grid->nelems2d_bed != grid->ncolumns) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> Number of 2d bedom elements should be equal to number of columns.");
            }
        }
    }
    
#ifdef _DEBUG
        if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    if (grid->ndim == 2) {
        
        grid->elem_error = (double *) tl_alloc(sizeof(double), grid->nelems2d); // allocate errors
        
        /**********************************************************************/
        /**********************************************************************/
        /* assign elements to IDs on multicore runs to not overcalculate mass */
        /* cjt :: seems a bit convoluted to me */
        int **nds_on_edge;
        int flag1 = 1, iedge = UNSET_INT, nnodes2d = UNSET_INT;
        for (i = 0; i < grid->nelems2d; i++) {
            
            if (grid->elem2d[i].nnodes == NDONTRI) {
                nnodes2d = NDONTRI;
                nds_on_edge = grid->nd_on_TriEdge;
            } else {
                nnodes2d = NDONQUAD;
                nds_on_edge = grid->nd_on_QuadEdge;
            }
            
#ifdef _MESSG
	    /*mwf pe1,pe2 not defined?*/
            grid->elem2d[i].my_pe = 0; // cjt :: note, this is just a flag, not actualy PE id
            flag1 = 1;
            for (iedge=0; iedge<nnodes2d; iedge++) {
                int pe1 = grid->node[grid->elem2d[i].nodes[ nds_on_edge[iedge][0] ]].resident_pe;
                int pe2 = grid->node[grid->elem2d[i].nodes[ nds_on_edge[iedge][1] ]].resident_pe;
                if (pe1 == pe2) flag = 0;
            }
            if (flag1 == 1) {
                grid->elem2d[i].my_pe = 1;
                for (inode=0; inode<nnodes2d; inode++) {
                    if (grid->node[grid->elem2d[i].nodes[inode]].resident_pe > grid->smpi->myid) grid->elem2d[i].my_pe++;
                }
            } else {
                for (inode=0; inode<nnodes2d; inode++) {
                    if (grid->node[grid->elem2d[i].nodes[inode]].resident_pe == grid->smpi->myid) grid->elem2d[i].my_pe++;
                }
            }
#else
            grid->elem2d[i].my_pe = nnodes2d;
#endif
            /**********************************************************************/
            /**********************************************************************/
            
            grid->elem2d[i].id = i;
            grid->elem2d[i].id_orig = i;
            
            // get 2D element djacs, normals and basis function gradients (elemental constants)
            SVECT nds[nnodes2d];
            for (inode=0; inode<nnodes2d; inode++) {
                nds[inode].x = grid->node[ grid->elem2d[i].nodes[inode] ].x;
                nds[inode].y = grid->node[ grid->elem2d[i].nodes[inode] ].y;
                nds[inode].z = grid->node[ grid->elem2d[i].nodes[inode] ].z; // initial displacement should be added here
            }
            
            if (nnodes2d == NDONTRI) {
                get_triangle_linear_djac_nrml_gradPhi(&(grid->elem2d[i]), NULL, nds);
            } else {
                grid->elem2d[i].nrml = get_elem2d_normals(nds);
            }
            
            if (grid->elem2d[i].djac < SMALL6) {
                fprintf(stderr, "WARNING :: Improperly numbered triangle=%d on PE=%d || Nodes %d %d %d Jacobian is: %20.10f \n",
                        (i + 1), grid->smpi->myid, grid->elem2d[i].nodes[0], grid->elem2d[i].nodes[1], grid->elem2d[i].nodes[2],grid->elem2d[i].djac);
                
                snode_printScreen(grid->node[grid->elem2d[i].nodes[0]]);
                snode_printScreen(grid->node[grid->elem2d[i].nodes[1]]);
                snode_printScreen(grid->node[grid->elem2d[i].nodes[2]]);
                tl_error("ERROR: Improperly numbered triangle");
            }
            
            grid->mesh_volume += get_elem2d_area2d(&(grid->elem2d[i]),grid->node);
            grid->elem_error[i] = 0.;
        }
        
        /* computes the number material types */
        for (i = 0; i < grid->nelems2d; i++) {
            if (grid->elem2d[i].mat >= grid->nmat) {
                grid->nmat = grid->elem2d[i].mat + 1;
            }
        }
        
        /* these are for 3d, but they will keep us from if-checking constantly elsewhere */
        grid->my_nnodes_bed = grid->my_nnodes;
        grid->my_nnodes_sur = grid->my_nnodes;
        grid->nnodes_bed = grid->nnodes;
        grid->nnodes_sur = grid->nnodes;
        grid->initial_nnodes_bed = grid->nnodes;
        grid->initial_nelems_bed = grid->nelems2d;
#ifndef _MESSG
        /*These should be set in read_geo_mpi for multi processor runs */
        grid->macro_nnodes = grid->nnodes;
#endif
        grid->macro_nnodes_sur = grid->macro_nnodes;
        grid->macro_nnodes_bed = grid->macro_nnodes;
        grid->orig_macro_nnodes = grid->macro_nnodes;
        grid->orig_macro_nnodes_sur = grid->macro_nnodes_sur;
        grid->orig_macro_nnodes_bed = grid->macro_nnodes_bed;
        
    }
    
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    if (DEBUG)  printf("myid_model: %d || myid_world: %d || grid->ndim: %d grid->nelems2d: %d || grid->nelems3d: %d\n",grid->smpi->myid,myid_world,grid->ndim,grid->nelems2d,grid->nelems3d);
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    if (grid->ndim == 3) {
        
#ifdef _MESSG
#else
        for (i = 0; i < grid->nelems3d; i++) {
            grid->elem3d[i].my_pe = grid->elem3d[i].nnodes;
        }
        for (i = 0; i < grid->nelems2d; i++) {
            grid->elem2d[i].my_pe = grid->elem2d[i].nnodes;
        }
#endif
        
        // -------------------------------------------------------------//
        // build 2D element face info
        for (i = 0; i < grid->nelems2d; i++) {
            grid->elem2d[i].id = i;
            grid->elem2d[i].id_orig = i;
            
            // get 2D element djacs, normals and basis function gradients (elemental constants)
            SVECT nds[grid->elem2d[i].nnodes];
            for (inode=0; inode<grid->elem2d[i].nnodes; inode++) {
                nds[inode].x = grid->node[ grid->elem2d[i].nodes[inode] ].x;
                nds[inode].y = grid->node[ grid->elem2d[i].nodes[inode] ].y;
                nds[inode].z = grid->node[ grid->elem2d[i].nodes[inode] ].z;  // initial displacement should be added here
            }
            if (grid->elem2d[i].nnodes == NDONTRI) {
                get_triangle_linear_djac_nrml_gradPhi(&(grid->elem2d[i]), NULL, nds);
            } else {
                grid->elem2d[i].nrml = get_elem2d_normals(nds);
            }
            
            // save t=0 grid 2D jacobian in 3D space
            if (fabs(grid->elem2d[i].djac3d_fixed) < 1e-7) {
                grid->elem2d[i].djac3d_fixed = grid->elem2d[i].djac3d;
            }
            
            // cjt :: note :: here is what I learned for 3d bed faces due to backward numbering
            // cjt :: the basis gradients are ok (same as surface)
            // cjt :: the normals are flipped (as they should be)
            // cjt :: all jacobians are (-), which will cause issues ...
            if (grid->elem2d[i].bflag == 1) { // on bed
                grid->elem2d[i].djac = fabs(grid->elem2d[i].djac);
                grid->elem2d[i].djac3d = fabs(grid->elem2d[i].djac3d);
                grid->elem2d[i].djac3d_fixed = fabs(grid->elem2d[i].djac3d_fixed);
            }
            
            //selem2d_printScreen(grid->elem2d[i]); exit(-1);
            if (grid->elem2d[i].djac3d < SMALL6 || grid->elem2d[i].djac3d_fixed < SMALL6) {
                selem2d_printScreen(&grid->elem2d[i]);
                fprintf(stderr, "Improperly numbered triangle number %d :: djac:  %20.10e\n", (i + 1),grid->elem2d[i].djac);
                tl_error("ERROR: Improperly numbered triangle number\n");
            }
        }

        // -------------------------------------------------------------//
        // build 3D element info
        grid->elem_error = (double *) tl_alloc(sizeof(double), grid->nelems3d); // allocate errors
        grid->hyd_eddy = (double *) tl_alloc(sizeof(double), grid->nelems3d);
        grid->trn_diff = (double *) tl_alloc(sizeof(double), grid->nelems3d);
        for (i = 0; i < grid->nelems3d; i++) {
            grid->elem_error[i] = 0.;
            grid->hyd_eddy[i] = 0.;
            grid->trn_diff[i] = 0.;
            grid->elem3d[i].id = i;
            grid->elem3d[i].id_orig = i;
            // get 3D element djacs and basis function gradients
            if (grid->elem3d[i].nnodes == NDONTET) {
                SVECT nds[grid->elem3d[i].nnodes];
                for (inode=0; inode<grid->elem3d[i].nnodes; inode++) {
                    nds[inode].x = grid->node[ grid->elem3d[i].nodes[inode] ].x;
                    nds[inode].y = grid->node[ grid->elem3d[i].nodes[inode] ].y;
                    nds[inode].z = grid->node[ grid->elem3d[i].nodes[inode] ].z;  // initial displacement should be added here
                }
                get_tet_linear_djac_gradPhi(&(grid->elem3d[i]), NULL, nds);
                if (grid->elem3d[i].djac < 0.0) {
                    fprintf(stderr, "Improperly numbered tetrahedron number %d :: djac:  %20.10e\n", (i + 1),grid->elem3d[i].djac);
                    tl_error("ERROR: Improperly numbered tetrahedron number\n");
                }
            }
        }
        
        // assume initial dpl is zero!
        if (flag == 1) {
            double dpl[grid->nnodes]; sarray_init_dbl(dpl,grid->nnodes); // temporary
            grid->mesh_volume = tl_find_grid_mass_elem3d(1., grid, dpl);
        }
        
        /* computes the number material types */
        for (i = 0; i < grid->nelems3d; i++) {
            if (grid->elem3d[i].mat >= grid->nmat) {
                grid->nmat = grid->elem3d[i].mat + 1;
            }
        }
        
        // determine 3d elements col, bed and surface elements
        // note the original code, this is a loop over nelem3d, not columns, which leads to slightly different residuals
        if (grid->type == COLUMNAR) {
            int ie3d = UNSET_INT, icol=UNSET_INT;
            ID_LIST_ITEM *ptr;
            for(icol=0; icol<grid->ncolumns; icol++)  {
                ptr = grid->column_list[icol];
                while(ptr->next != NULL) {
                    ie3d = ptr->id;
                    grid->elem3d[ie3d].icol = icol;
                    grid->elem3d[ie3d].elem2d_sur = grid->elem2d_sur[icol];
                    grid->elem3d[ie3d].elem2d_bed = grid->elem2d_bed[icol];
                    ptr = ptr->next;
                }
            }
        }
        
        grid->initial_nnodes_bed = grid->nnodes_bed;
        grid->initial_nelems_bed = grid->nelems2d_bed;
    }
    grid->nnodes_prev = grid->nnodes;
    grid->nnodes_sur_prev = grid->nnodes_sur;
    grid->nnodes_bed_prev = grid->nnodes_bed;
    
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /* calculate quadrature points/weights */
    int quad_flag;
    if (grid->haveTets == TRUE) quad_flag = SQUAD_tetrahedron_alloc_init(&(grid->quad_tet));
    if (grid->haveQuads == TRUE) quad_flag = SQUAD_rectangle_alloc_init(&(grid->quad_rect));
    if (grid->havePrisms == TRUE) quad_flag = SQUAD_triprism_alloc_init(&(grid->quad_prism));
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /* write grid info to screen -- CJT :: do this after partition now! */
#ifdef _DEBUG
//    /*
    
    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
    if (grid->ndim == 2) {
        sgrid_printScreen(grid, io->geo2d.filename);
    } else if (grid->ndim == 3) {
        sgrid_printScreen(grid, io->geo3d.filename);
    }
    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
//     */
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /* create a child mesh */
    if (grid->ndim == 3 && file_output.grid2dm == ON && flag == 1) {
        //createChildMesh(grid, io);
    }

    /**mwf debug 
    for (i=0; i < grid->nelems3d; i++) {
      printf("sgrid_alloc_init elem=%d area=%g\n",i,grid->elem3d[i].djac);
    }
    */
#ifdef _DEBUG
#ifdef _MESSG
    printf("CSTORM COMM PE: %d :: MODEL PE: %d :: FLAG: %d :: FINISHED Allocating, initializing and reading model grid\n",myid_world,myid_model,flag);
#else
    printf("Finished allocating, initializing and reading model grid\n");
#endif
    if (myid_world == 0) printf("-------------------------------------------------------\n\n");
    
#endif
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void sgrid_printScreen(SGRID *grid, char *filename) {
    
    // if serial, print to screen, otherwise, print each PE to a file
    char fn[50+1];
    snprintf(fn, 50, "%s_mesh_info_pe%d.dat", filename,grid->smpi->myid);

    
    FILE *fp;
    if (grid->smpi->npes == 1) {
        fp = stdout;
    } else {
        fp = fopen(fn, "w");
    }
    
    fprintf(fp,"\n");
    fprintf(fp,"***********************************************************\n");
    fprintf(fp,"***********************************************************\n");
    fprintf(fp,"----- Grid: %s statistics\n", filename);
    fprintf(fp,"---------- initial_nnodes: %d \n", grid->initial_nnodes);
    fprintf(fp,"---------- max_nnodes: %d \n", grid->max_nnodes);
    fprintf(fp,"---------- macro_nnodes: %d \n", grid->macro_nnodes);
    fprintf(fp,"---------- macro_nelems1d: %d \n", grid->macro_nelems1d);
    fprintf(fp,"---------- macro_nelems2d: %d \n", grid->macro_nelems2d);
    fprintf(fp,"---------- macro_nelems3d: %d \n", grid->macro_nelems3d);
    fprintf(fp,"---------- total number of nodes: %d\n", grid->nnodes);
    fprintf(fp,"---------- total number of my_nodes: %d\n", grid->my_nnodes);
    fprintf(fp,"---------- total number of 1d elements: %d\n", grid->nelems1d);
    fprintf(fp,"---------- total number of 2d elements: %d\n", grid->nelems2d);
    if (grid->ndim == 3) fprintf(fp,"---------- total number of 3d elements: %d\n", grid->nelems3d);
    fprintf(fp,"---------- total number of materials: %d\n", grid->nmat);
    if (grid->ndim == 3 && grid->type == COLUMNAR) {
        fprintf(fp,"---------- surface statistics\n");
        fprintf(fp,"---------------- total number of surface nodes on global grid: %d\n", grid->macro_nnodes_sur);
        fprintf(fp,"---------------- total number of 2d surface elements on global grid: %d\n", grid->macro_nelems2d_bed);
        fprintf(fp,"---------------- total number of surface nodes: %d\n", grid->nnodes_sur);
        fprintf(fp,"---------------- total number of surface my_nodes: %d\n", grid->my_nnodes_sur);
        fprintf(fp,"---------------- total number of 2d & 3d surface elements (columns): %d\n", grid->nelems2d_sur);
        fprintf(fp,"---------- bed statistics\n");
        fprintf(fp,"---------------- total number of bed nodes on global grid: %d\n", grid->macro_nnodes_bed);
        fprintf(fp,"---------------- total number of 2d bed elements on global grid: %d\n", grid->macro_nelems2d_bed);
        fprintf(fp,"---------------- total number of bed nodes: %d\n", grid->nnodes_bed);
        fprintf(fp,"---------------- total number of 2d & 3d bed elements (columns): %d\n", grid->nelems2d_bed);
    }
    fprintf(fp,"---------- total mesh area: %14.6e\n", grid->mesh_volume);
    fprintf(fp,"***********************************************************\n");
    fprintf(fp,"***********************************************************\n");

    if (grid->smpi->npes != 1) {
        int i,icol;
        fprintf(fp,"NODES\n");
        for (i=0; i<grid->nnodes; i++) {
            //if (grid->node[i].gid == 6558) assert(fabs(grid->node[i].x - 36415.99923111)<1e-4 && fabs(grid->node[i].y - 33587.57210636)<1e-4);
            fprintf(fp,"%d %d %d %d %d ",grid->node[i].gid,grid->node[i].id,grid->node[i].myid,grid->node[i].resident_pe,grid->node[i].resident_id);
            if (grid->ndim == 3 && grid->type == COLUMNAR) {
                fprintf(fp,"%d %d",grid->nodeID_3d_to_2d_sur[i],grid->node[i].global_surf_id);
            }
            fprintf(fp,"\n");
        }
        
        
        fprintf(fp,"3d elements\n");
        for (i=0; i<grid->nelems3d; i++) {
            fprintf(fp,"%d %d %d %d %d %d %d \t nodes: %d %d %d %d \t nedges: %d nnodes: %d \t nnodes_quad: %d\n",
                    grid->elem3d[i].gid,grid->elem3d[i].id,grid->elem3d[i].id_orig,grid->elem3d[i].my_pe,
                    grid->elem3d[i].icol,grid->elem3d[i].elem2d_sur,grid->elem3d[i].elem2d_bed,
                    grid->elem3d[i].nodes[0],grid->elem3d[i].nodes[1],grid->elem3d[i].nodes[2],grid->elem3d[i].nodes[2],
                    grid->elem3d[i].nedges,grid->elem3d[i].nnodes,grid->elem3d[i].nnodes_quad);
        }
        fprintf(fp,"2d elements\n");
        for (i=0; i<grid->nelems2d; i++) {
            fprintf(fp,"%d %d %d %d %d %d \t nodes: %d %d %d \t nedges: %d nnodes: %d \t nnodes_quad: %d\n",
                    grid->elem2d[i].gid,grid->elem2d[i].id,grid->elem2d[i].id_orig,grid->elem2d[i].id_3d,
                    grid->elem2d[i].my_pe,grid->elem2d[i].resident_pe,
                    grid->elem2d[i].nodes[0],grid->elem2d[i].nodes[1],grid->elem2d[i].nodes[2],
                    grid->elem2d[i].nedges,grid->elem2d[i].nnodes,grid->elem2d[i].nnodes_quad);
        }
        
        if (grid->ndim == 3 && grid->type == COLUMNAR) {
            
            fprintf(fp,"elem2d_sur :: ncolumns: %d\n",grid->ncolumns);
            for (icol = 0; icol < grid->ncolumns; icol++) {
                fprintf(fp,"icol: %d \t gelem2d_sur[icol]: %d\n",icol,grid->elem2d_sur[icol]);
            }
            
            fprintf(fp,"column_list\n");
            ID_LIST_ITEM *ptr;
            for(icol=0; icol<grid->ncolumns; icol++)  {
                ptr = grid->column_list[icol];
                while(ptr->next != NULL) {
                    i = ptr->id;
                    ptr = ptr->next;
                    fprintf(fp,"icol: %d \t ie3d: %d \t elem3d.id: %d\n",icol,i,grid->elem3d[i].id);
                }
            }
            
            fprintf(fp,"midpoint list\n");
            MIDPT_LIST_ITEM *myptr;
            for (i=0; i<grid->num_midpts; i++)  {
                myptr = grid->midpt_list[i];
                fprintf(fp,"midpoint: %d \t top_node: %d \t bot_node: %d\n",i,myptr->node1,myptr->node2);
            }
        }
        

        fclose(fp);
    }
#ifdef _DEBUG
    if (DEBUG_WITH_PICKETS) tl_check_all_pickets(__FILE__, __LINE__);
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void sgrid_free(SGRID *grid) {
    
    
    int inode, ie, i, k;
    
    if (grid != NULL) {
        
        if (grid->nd_on_TriEdge != NULL) {
            for (i=0; i<3; i++) {
                grid->nd_on_TriEdge[i] = (int *) tl_free(sizeof(int), 2, grid->nd_on_TriEdge[i]);
            }
            grid->nd_on_TriEdge = (int **) tl_free(sizeof(int *), 3, grid->nd_on_TriEdge);
        }
        
        if (grid->nd_on_QuadEdge != NULL) {
            for (i=0; i<4; i++) {
                grid->nd_on_QuadEdge[i] = (int *) tl_free(sizeof(int), 2, grid->nd_on_QuadEdge[i]);
            }
            grid->nd_on_QuadEdge = (int **) tl_free(sizeof(int *), 4, grid->nd_on_QuadEdge);
        }
        
        if (grid->nd_on_PrismEdge != NULL) {
            for (i=0; i<9; i++) {
                grid->nd_on_PrismEdge[i] = (int *) tl_free(sizeof(int), 2, grid->nd_on_PrismEdge[i]);
            }
            grid->nd_on_PrismEdge = (int **) tl_free(sizeof(int *), 9, grid->nd_on_PrismEdge);
        }
        
        if (grid->nd_on_TetEdge != NULL) {
            for (i=0; i<6; i++) {
                grid->nd_on_TetEdge[i] = (int *) tl_free(sizeof(int), 2, grid->nd_on_TetEdge[i]);
            }
            grid->nd_on_TetEdge = (int **) tl_free(sizeof(int *), 6, grid->nd_on_TetEdge);
        }

        if (grid->elem1d != NULL && grid->nelems1d > 0) {
            selem1d_free_array(grid->elem1d, grid->max_nelems1d);
        }
        if (grid->elem2d != NULL && grid->nelems2d > 0) {
            //printf("pe: %d :: grid->nelems2d: %d grid->max_nelems2d: %d\n",grid->smpi->myid,grid->nelems2d,grid->max_nelems2d);
            selem2d_free_array(grid->elem2d, grid->max_nelems2d);
        }
        if (grid->elem3d != NULL && grid->nelems3d > 0) {
            //printf("pe: %d :: grid->nelems3d: %d grid->max_nelems3d: %d\n",grid->smpi->myid,grid->nelems3d,grid->max_nelems3d);
            selem3d_free_array(grid->elem3d, grid->max_nelems3d);
        }
        
        if (grid->nnodes > 0) {
            for (inode=0; inode<grid->max_nnodes; inode++) {
                snode_free(&(grid->node[inode]));
            }
            grid->node = (SNODE *) tl_free(sizeof(SNODE), grid->max_nnodes, grid->node);
        }
        
        if (grid->old_global_surf != NULL) grid->old_global_surf = (int *)tl_free(sizeof(int), grid->max_nnodes_sur,grid->old_global_surf);
        if (grid->old_global_bed  != NULL) grid->old_global_bed = (int *)tl_free(sizeof(int), grid->max_nnodes_bed,grid->old_global_bed);
        
        if (grid->ndim == 2 && grid->nelems2d > 0) {
            grid->wd_flag = (int *) tl_free(sizeof(int), grid->max_nelems2d, grid->wd_flag);
            grid->elem_error = (double *) tl_free(sizeof(double), grid->max_nelems2d, grid->elem_error);
        } else if (grid->nelems3d > 0) {
            grid->elem_error = (double *) tl_free(sizeof(double), grid->max_nelems3d, grid->elem_error);
            grid->hyd_eddy = (double *) tl_free(sizeof(double), grid->max_nelems3d, grid->hyd_eddy);
            grid->trn_diff = (double *) tl_free(sizeof(double), grid->max_nelems3d, grid->trn_diff);
        }
        
        // free column variables :: cjt :: if unstructured grid, these are NULL so ok to call
        // Gajanan gkc - Even unstructured GW model has surf and bed elements
        // which are only freed in this function. So we need this to always be
        // called. So I am removing this if condition: if (grid->ncolumns > 0)
        global_free_columns(grid);
        
        if (grid->start_node_ID !=NULL)
            grid->start_node_ID = (int *) tl_free(sizeof(int), grid->smpi->npes,grid->start_node_ID);
        if (grid->end_node_ID !=NULL)
            grid->end_node_ID = (int *) tl_free(sizeof(int), grid->smpi->npes, grid->end_node_ID);
        
        
        // free mpi stuff
        if(grid->smpi != NULL){
#ifdef _MESSG
            smpi_free(grid->smpi);
            
            if(grid->smpi->partition_info != NULL)
                grid->smpi->partition_info = (int *) tl_free(sizeof(int), grid->max_nnodes, grid->smpi->partition_info);
            
            if(grid->smpi->surface_partition_info != NULL)
                grid->smpi->surface_partition_info = (int *) tl_free(sizeof(int), grid->max_nnodes_sur, grid->smpi->surface_partition_info);
#endif
            if(grid->part_map != NULL){
                grid->part_map = (int *)tl_free(sizeof(int),grid->part_smpi->npes,grid->part_map);
            }
            if(grid->part_smpi != NULL){
                smpi_free(grid->part_smpi);
                grid->part_smpi = (SMPI *) tl_free(sizeof(SMPI), 1, grid->part_smpi);
            }
            grid->smpi = (SMPI *) tl_free(sizeof(SMPI), 1, grid->smpi);
        }
        
        /* free quadrature points */
        if (grid->haveQuads == TRUE) SQUAD_free(grid->quad_rect, NDONQUAD);
        if (grid->haveTets == TRUE) SQUAD_free(grid->quad_tet, NDONTET);
        if (grid->havePrisms == TRUE) SQUAD_free(grid->quad_prism, NDONPRISM);
        
        /* free grid */
        grid = (SGRID *) tl_free(sizeof(SGRID), 1, grid);
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void sgrid_init(SGRID *g, int flag
#ifdef _MESSG
                , MPI_Comm model_comm
#endif
) {
    
    g->type = UNSTRUCTURED;
    g->isGridMixedElement = FALSE;
    g->haveTets = FALSE;
    g->havePrisms = FALSE;
    g->haveTris = FALSE;
    g->haveQuads = FALSE;
    
    g->nd_on_TriEdge = NULL;
    g->nd_on_QuadEdge = NULL;
    g->nd_on_PrismEdge = NULL;
    g->nd_on_TetEdge = NULL;
    
    g->ndim = UNSET_INT;
    g->macro_nnodes = UNSET_INT;
    g->macro_nelems1d =  UNSET_INT;
    g->macro_nelems2d =  UNSET_INT;
    g->macro_nelems3d =  UNSET_INT;
    g->macro_nelems2d_bed = UNSET_INT;
    g->initial_nnodes = UNSET_INT;
    g->initial_nelems = UNSET_INT;
    g->initial_nnodes_bed = UNSET_INT;
    g->initial_nelems_bed = UNSET_INT;
    g->orig_initial_nnodes = UNSET_INT;
    g->orig_initial_nelems = UNSET_INT;
    g->nnodes = UNSET_INT;
    g->nelems3d = UNSET_INT;
    g->nelems2d = UNSET_INT;
    g->nelems1d = UNSET_INT;
    g->nelems1d2d = UNSET_INT;
    g->nedges = UNSET_INT;
    g->nmat = UNSET_INT;
    g->nnodes_bed_prev = UNSET_INT;
    g->nnodes_sur_prev = UNSET_INT;
    g->nnodes_matrix = UNSET_INT;
    g->max_nelems1d = UNSET_INT;
    g->my_nnodes = UNSET_INT;
    g->interface = 0;

    g->elem1d = NULL;
    g->elem2d = NULL;
    g->elem3d = NULL;
    g->node = NULL;
    
    g->x_min = 0.;
    g->x_max = 0.;
    g->y_min = 0.;
    g->y_max = 0.;
    g->xL = 0.;
    g->yL = 0.;
    
    g->nodeID_3d_to_2d_sur = NULL;
    g->nodeID_3d_to_2d_bed = NULL;
    g->nodeID_2d_to_3d_sur = NULL;
    g->nodeID_2d_to_3d_bed = NULL;
    g->wd_flag = NULL;
    g->elem_error = NULL;
    g->hyd_eddy = NULL;
    g->trn_diff = NULL;
    
    g->start_node_ID = NULL;
    g->end_node_ID = NULL;
    
    /* columns ***************************************************/
    g->nnodes_sur = UNSET_INT;
    g->nnodes_bed = UNSET_INT;
    g->macro_nnodes_sur = UNSET_INT;
    g->macro_nnodes_bed = UNSET_INT;
    g->num_midpts = UNSET_INT;
    g->num_vert_segments = UNSET_INT;
    g->nelems2d_sur = UNSET_INT;
    g->nelems2d_bed = UNSET_INT;
    g->nelems2d_sidewall = UNSET_INT;
    
    g->elem3d_sur = NULL;
    g->elem2d_sur = NULL;
    g->elem3d_bed = NULL;
    g->elem2d_bed = NULL;
    g->elem3d_sidewall = NULL;
    g->elem2d_sidewall = NULL;
    
    g->isize_sidewall_elems = UNSET_INT;
    g->isize_sur_elems = UNSET_INT;
    g->isize_bed_elems = UNSET_INT;
    g->isize3= 0; g->isize4 = 0; g->isize5 = 0;
    g->scale_factor = 0;
    g->shift_factor = 0;
    
    g->ncolumns = UNSET_INT;
    g->isize_ncolumns = 0;
    
    g->column_list = NULL;
    g->column_list2d = NULL;
    g->vertical_list = NULL;
    g->sidewall_list = NULL;
    g->midpt_list = NULL;
    
    g->hash_size = 0;
    g->is_allocated_column_hash = UNSET_INT;
    g->column_hash = NULL;
    g->is_allocated_column_hash2d = UNSET_INT;
    g->column_hash2d = NULL;
    g->is_allocated_midpt_hash = UNSET_INT;
    g->midpt_hash = NULL;
    g->is_allocated_vertical_hash = UNSET_INT;
    g->vertical_hash = NULL;
    
    g->nnodes_sur = 0;
    g->nnodes_bed = 0;
    g->max_nnodes_sur = 0;
    g->max_nnodes_bed = 0;
    g->num_midpts = 0;
    g->num_vert_segments = 0;
    g->nelems2d_sur = 0;
    g->nelems2d_bed = 0;
    g->nelems2d_sidewall = 0;
    g->isize_sidewall_elems = 0;
    g->isize_sur_elems = 0;
    g->isize_bed_elems = 0;
    g->isize3 = 0;
    g->isize4 = 0;
    g->isize5 = 0;
    g->scale_factor = 0;
    g->shift_factor = 0;
    g->ncolumns = 0;
    g->isize_ncolumns = 0;
    g->hash_size = 0;
    g->is_allocated_column_hash = 0;
    g->is_allocated_column_hash2d = 0;
    g->is_allocated_midpt_hash = 0;
    g->is_allocated_vertical_hash = 0;
    
    g->elem3d_sur = NULL;
    g->elem2d_sur = NULL;
    g->elem3d_bed = NULL;
    g->elem2d_bed = NULL;
    g->elem3d_sidewall = NULL;
    g->elem2d_sidewall = NULL;
    g->column_list = NULL;
    g->column_list2d = NULL;
    g->vertical_list = NULL;
    g->sidewall_list = NULL;
    g->midpt_list = NULL;
    g->column_hash = NULL;
    g->column_hash2d = NULL;
    g->midpt_hash = NULL;
    g->vertical_hash = NULL;
    
    g->old_global_surf = NULL;
    g->old_global_bed = NULL;
    
    g->smpi = (SMPI *) tl_alloc(sizeof(SMPI), 1);
    g->part_smpi = NULL;
    g->part_map = NULL;
    smpi_defaults(g->smpi);
    if(flag==1){
        smpi_init(g->smpi
#ifdef _MESSG
                  , model_comm
#endif
                  );
    }
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// cjt :: needs to be updated for more general element meshes

void sgrid_print_adapted_ts(SGRID *grid, SFILE sup, char *root_name, int it1, int it2, int **ndata, int my_nnode_max, int *my_nnode_ext) {
    int super = strlen(sup.filename);
    int is=UNSET_INT, ie=UNSET_INT, n=UNSET_INT, ierr=UNSET_INT;
    int ip;
    
    // define, initialize and open adapted grid file
    SFILE fout_geo_adapt;
    if(grid->smpi->myid<=0){
        init_adh_file(&(fout_geo_adapt));
        build_filename2(fout_geo_adapt.filename, MAXLINE, root_name, ".3dm-", it1,".", it2);
        open_output_file(&(fout_geo_adapt), "adapted 3dm file", super);
    }
    
    /* Find Maximum Array Sizes And Allocate Temporary Arrays */
    int nelem2d_max = grid->macro_nelems2d;
    int nelem3d_max = grid->macro_nelems3d;
    
    
    /* Determine Processor Owning Element (Processor Owning Minimum Node Number Of Element) */
    int itdim = nelem2d_max * 4;
    int *itemp = (int *) tl_alloc(sizeof(int), itdim);
    int mec = 0;
    for (ie = 0; ie < grid->nelems2d; ie++) {
#ifdef _MESSG
        if (grid->elem2d[ie].my_pe > 1) {
            for (n = 0; n < NDPRFC; n++) {
                itemp[mec++] = grid->node[grid->elem2d[ie].nodes[n]].gid + 1; /* element nodes */
            }
            itemp[mec++] = grid->elem2d[ie].mat + 1;    /* material type */
        }
#else
        for (n = 0; n < NDPRFC; n++) {
            itemp[mec++] = grid->node[grid->elem2d[ie].nodes[n]].gid + 1; /* element nodes */
        }
        itemp[mec++] = grid->elem2d[ie].mat + 1;    /* material type */
#endif
    }
    int nelem2d_loc = mec / 4;
#ifdef _MESSG
    ierr = MPI_Allreduce(&nelem2d_loc, &nelem2d_max, 1, MPI_INT, MPI_MAX, grid->smpi->ADH_COMM);
    if (ierr != MPI_SUCCESS) {
        messg_err(ierr);
    }
    ierr = MPI_Allreduce(&grid->nelems3d, &nelem3d_max, 1, MPI_INT, MPI_MAX, grid->smpi->ADH_COMM);
    if (ierr != MPI_SUCCESS) {
        messg_err(ierr);
    }
#endif
    /* Write Element Information */
#ifdef _MESSG
    MPI_Status *msg_status = grid->smpi->msg_status;
    if (grid->smpi->myid <= 0) {
#endif
        fprintf(fout_geo_adapt.fp, "MESH2D\n");
        ie = 1;
        for (n = 0, is = 0; n < nelem2d_loc; n++) {
            fprintf(fout_geo_adapt.fp, "E3T  %8d %8d %8d %8d %4d\n", ie++, itemp[is], itemp[is + 1], itemp[is + 2], itemp[is + 3]);
            is += 4;
        }
#ifdef _MESSG
        for (ip = 1; ip < grid->smpi->npes; ip++) {
            ierr = MPI_Recv(itemp, nelem2d_max * 4, MPI_INT, ip, 997, grid->smpi->ADH_COMM, msg_status);
            if (ierr != MPI_SUCCESS) {
                messg_err(ierr);
            }
            ierr = MPI_Get_count(msg_status, MPI_INT, &nelem2d_loc);
            if (ierr != MPI_SUCCESS) {
                messg_err(ierr);
            }
            nelem2d_loc = nelem2d_loc / 4;
            for (n = 0, is = 0; n < nelem2d_loc; n++) {
                fprintf(fout_geo_adapt.fp, "E3T  %8d %8d %8d %8d %4d\n", ie++, itemp[is], itemp[is + 1], itemp[is + 2], itemp[is + 3]);
                is += 4;
            }
        }
    }
    else {
        ierr = MPI_Send(itemp, mec, MPI_INT, 0, 997, grid->smpi->ADH_COMM);
        if (ierr != MPI_SUCCESS) {
            messg_err(ierr);
        }
    }
    messg_barrier(grid->smpi->ADH_COMM);
#endif
    itemp = (int *) tl_free(sizeof(int), itdim, itemp);
    
    /* Write Node Information */
    int i;
    int icnt=0, index=0;   /* external processor nodal information */
    double *node_data, *gdata;      /* local node and global data array */
    node_data = (double *) tl_alloc(sizeof(double), my_nnode_max*3);
    for (i=0; i<grid->my_nnodes; i++) {
        node_data[icnt++] = grid->node[i].x;
        node_data[icnt++] = grid->node[i].y;
        node_data[icnt++] = grid->node[i].z;
    }
    icnt=0;
    if(grid->smpi->myid==0){
        gdata = (double *) tl_alloc(sizeof(double), grid->macro_nnodes*3);
        
        /* load pe 0 nodes into gdata */
        for (i=0; i<grid->my_nnodes; i++) {
            icnt = grid->node[i].gid *3;
            gdata[icnt] = node_data[i*3];
            gdata[icnt+1] = node_data[(i*3)+1];
            gdata[icnt+2] = node_data[(i*3)+2];
        }
#ifdef _MESSG
        /*load external processor node_data into gdata */
        for (ip = 1; ip < grid->smpi->npes; ip++) {
            ierr = MPI_Recv(node_data, my_nnode_ext[ip]*3, MPI_DOUBLE, ip, 999, grid->smpi->ADH_COMM, msg_status);
            if (ierr != MPI_SUCCESS)
                messg_err(ierr);
            index=0;icnt=0;
            for (i=0;i<my_nnode_ext[ip];i++){
                index = ndata[ip][i]*3;
                gdata[index] = node_data[icnt++];
                gdata[index + 1] = node_data[icnt++];
                gdata[index + 2] = node_data[icnt++];
            }
        }
#endif
        index=0;
        for (i=0;i<grid->macro_nnodes;i++){
            index=i*3;
            fprintf(fout_geo_adapt.fp, "ND %d %16.8e %16.8e %16.8e\n",i+1,gdata[index],gdata[index+1],gdata[index+2]);
        }
        // added ENNDS and close file
        
        fclose(fout_geo_adapt.fp);
        fout_geo_adapt.fp = NULL;
        gdata = (double *) tl_free(sizeof(double), grid->macro_nnodes*3, gdata);
        
    }else{
#ifdef _MESSG
        ierr = MPI_Send(node_data, my_nnode_ext[grid->smpi->myid]*3, MPI_DOUBLE, 0, 999, grid->smpi->ADH_COMM);
        if (ierr != MPI_SUCCESS)
            messg_err(ierr);
#endif
    }
    node_data = (double *) tl_free(sizeof(double), my_nnode_max*3, node_data);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* This routine checks the integrity of the grid for hanging and missing nodes */
void sgrid_check(SGRID *grid, char *filename, int linenumber) {
    int i, j, k, inode, ierror=0;
    double sum, sum_check;
    int temp1, temp2, temp3, temp4;
#ifdef _MESSG
    printf("MYID %d checking grid at %s:%d ......", grid->smpi->myid, filename, linenumber);
#else
    printf("Checking grid at %s:%d ......", filename, linenumber);
#endif
    /* check the 3d element array */
    for (i = 0; i < grid->nelems3d; i++) {
        for(j=0; j<grid->elem3d[i].nnodes; j++){
            if(grid->elem3d[i].nodes[j] > grid->nnodes){ tl_error("3d Element has node outside of nnode range");}
#ifdef _MESSG
            if((grid->elem3d[i].nodes[j] < grid->my_nnodes) && (grid->node[grid->elem3d[i].nodes[j]].resident_pe != grid->smpi->myid)) {printf("Resident node %d has the wrong resident_pe", grid->elem3d[i].nodes[j] );ierror++;}
            if((grid->elem3d[i].nodes[j] >= grid->my_nnodes) && (grid->node[grid->elem3d[i].nodes[j]].resident_pe == grid->smpi->myid)) {printf("Ghost node %d is owned by me!",grid->elem3d[i].nodes[j] ); ierror++;}
#endif
        }
    }
    
    /* check the 2d element array */
    for (i = 0; i < grid->nelems2d; i++) {
        for(j=0; j<grid->elem2d[i].nnodes; j++){
            if(grid->elem2d[i].nodes[j] > grid->nnodes) tl_error("2D Element has node outside of nnode range");
#ifdef _MESSG
            if((grid->elem2d[i].nodes[j] < grid->my_nnodes) && (grid->node[grid->elem2d[i].nodes[j]].resident_pe != grid->smpi->myid)) tl_error("Resident node has the wrong resident_pe");
            if((grid->elem2d[i].nodes[j] >= grid->my_nnodes) && (grid->node[grid->elem2d[i].nodes[j]].resident_pe == grid->smpi->myid)) tl_error("Ghost node is owned by me!");
#endif
        }
    }
    
    /* check the 1d element array */
    for (i = 0; i < grid->nelems1d; i++) {
        for(j=0; j<grid->elem1d[i].nnodes; j++){
            if(grid->elem1d[i].nodes[j] > grid->nnodes) tl_error("1D Element has node outside of nnode range");
#ifdef _MESSG
            if((grid->elem1d[i].nodes[j] < grid->my_nnodes) && (grid->node[grid->elem1d[i].nodes[j]].resident_pe != grid->smpi->myid)) tl_error("Resident node has the wrong resident_pe");
            if((grid->elem1d[i].nodes[j] >= grid->my_nnodes) && (grid->node[grid->elem1d[i].nodes[j]].resident_pe == grid->smpi->myid)) tl_error("Ghost node is owned by me!");
#endif
        }
    }
    
    /* check for hanging nodes */
    for (k = 0; k < grid->nnodes; k++) {
        j=0;
        if (grid->ndim == 2) {
            for (i = 0; i < grid->nelems2d; i++) {
                for (inode=0; inode<grid->elem2d[i].nnodes; inode++) {
                    if (grid->elem2d[i].nodes[inode] == k) j = 1;
                }
            }
        } else {
            for (i = 0; i < grid->nelems3d; i++) {
                for (inode=0; inode<grid->elem3d[i].nnodes; inode++) {
                    if (grid->elem3d[i].nodes[inode] == k) j = 1;
                }
            }
        }
        if (j < 1) {
            printf("node: k %d\n", k);
            tl_error("Hanging node!");
        }
    }
    
    /* check that all nodes are used */
    sum = 0;
    for (k = 0; k < grid->my_nnodes; k++) {
        sum += grid->node[k].gid + 1;
    }
#ifdef _MESSG
    sum = messg_dsum(sum, grid->smpi->ADH_COMM);
#endif
    
    if(grid->smpi->npes>1){
        sum_check = grid->macro_nnodes * (grid->macro_nnodes + 1.0) * 0.5;
        
        if(fabs(sum - sum_check) > SMALL) {
            printf("Macro_nnodes %d sum_check %f sum %f \n", grid->my_nnodes, sum_check, sum);
            tl_error("Some nodes not owned!\n");
        }
    }
    /*check parents*/
    /*for (k = 0; k < grid->nnodes; k++) {
     if((grid->node[k].parent[0] < 0 || grid->node[k].parent[0] >= grid->nnodes || grid->node[k].parent[1] < 0 || grid->node[k].parent[1] >= grid->nnodes) && grid->node[k].original_id == UNSET_INT) {printf("MYID %d problem with parents for node %d of nnodes %d\n",grid->smpi->myid, k, grid->nnodes);snode_printScreen(grid->node[k]);ierror=1;}}*/
#ifdef _MESSG
    if(grid->smpi->npes>1) messg_barrier(grid->smpi->ADH_COMM);
#endif
    
    if(ierror>0) tl_error("Grid is NOT ok");
    printf("  Grid OK! \n");
    
    return;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// AdH surface and bed faces are always triangles
int get_macro_nelems2d_bed(SGRID *grid) {
    
    if (grid->type != COLUMNAR) return UNSET_INT; // unstructured grid
    
#ifdef _MESSG
    
    int myid = grid->smpi->myid;
    int npes = grid->smpi->npes;
    
    // now write element and conectivities
    int gid_0, gid_1, gid_2, total_nelems2d_bed = UNSET_INT;
    int ie, ie2d, nd0_3d, nd1_3d, nd2_3d, count=0;
    
    int kk = 0;
    for (ie = 0; ie < grid->nelems2d_sur; ie++) {
        ie2d = grid->elem2d_sur[ie]; // 2d element id in the complete list of 2d elements
        
        count = 0;
        nd0_3d = grid->elem2d[ie2d].nodes[0];
        nd1_3d = grid->elem2d[ie2d].nodes[1];
        nd2_3d = grid->elem2d[ie2d].nodes[2];
        
        if (grid->node[nd0_3d].resident_pe == myid) count++;
        if (grid->node[nd1_3d].resident_pe == myid) count++;
        if (grid->node[nd2_3d].resident_pe == myid) count++;
        
        // at least 2 nodes must be residential to own this element
        if (count > 1) {
            kk++;
        }
    }
    if(npes>1){
        MPI_Allreduce(&kk, &total_nelems2d_bed, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM);
    }else{
        total_nelems2d_bed = kk;
    }
    return total_nelems2d_bed;
#else
    return grid->nelems2d;
#endif
}



