/*! \file sw2_nb_test.c This file tests the sw2 engine */
#include "adh.h"
static double NEWTON_TEST_TOL = 1e-7;
static int NEWTON_TEST_NX = 101;//16;
static int NEWTON_TEST_NY = 6;//6;
static double write_testcase_error_wet_dry(SMODEL_SUPER *mod, double initial_grid_mass);
static void permute_array(double *arr,int *p, int n);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function tests a basic sw2 case where we have sloping beach and still
 *  conditions. Mass conservation is tested
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int sw2_nb_test(int argc, char **argv) {
	//create a grid
	SGRID *grid;
	grid = (SGRID *) tl_alloc(sizeof(SGRID), 1);
	//let's just do something simple
	//3x3 triangular element grid
	double xmin = 0.0;
	double xmax = 1000.0;
	double ymin = 0.0;
	double ymax = 50.0;
	int npx = NEWTON_TEST_NX;
	int npy = NEWTON_TEST_NY;
	double theta = 45.0;
	double dz = 1.0;
	double a0 = 0.0;
	double ax = 0.0015;
	double ax2 = 0.0;
	double ay = 0.0;
	double ay2 = 0.0;
	double axy = 0.0;
	double ax2y = 0.0;
	double axy2 = 0.0;
	double ax2y2 = 0.0;
	int flag3d =0;
	int err_code=-1;
    *grid = create_rectangular_grid(xmin, xmax, ymin, ymax, npx, npy,
 	theta, dz, a0, ax, ax2, ay, ay2, axy,
    ax2y, axy2, ax2y2, flag3d );

    //append 1d elements for BC enforcement
    //create a series and string for boundary data without file read
	//first add a string, this will come from smat physics file
	//! Downstream edges
	//EGS 1 2 2 
	//EGS 2 3 2 
	//EGS 3 4 2 
	//EGS 4 5 2 
	//EGS 5 6 2 
	//! Upstream edges
	//EGS 601 602 3 
	//EGS 602 603 3 
	//EGS 603 604 3 
	//EGS 604 605 3 
	//EGS 605 606 3
	grid->nelems1d = 10;
	grid->max_nelems1d = 10;
	grid->macro_nelems1d = grid->nelems1d;
	grid->orig_macro_nelems1d = grid->macro_nelems1d;
	grid->my_nelems1d = grid->nelems1d;
	printf("Attempting to add 1d elements: %d\n",grid->nelems1d);
	selem1d_init_alloc_array( &(grid->elem1d), grid->nelems1d);
	printf("Addied 1d elements\n");


	//hard coded, in reality this is given in geo now
	int nodes_segment[2];
	int nnodes_on_elem = 2;
	SVECT nds[nnodes_on_elem]; svect_init_array(nds, nnodes_on_elem);

	//Downstream Edges
	int node_ids[2] = {0,1};
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = grid->node[node_ids[i]].x;
    	nds[i].y = grid->node[node_ids[i]].y;
    	nds[i].z = grid->node[node_ids[i]].z;
    }
	selem1d_load(&(grid->elem1d[0]), 0, 0, 2, node_ids, 0, nds, 0);
	
	node_ids[0] = 1; node_ids[1] = 2;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = grid->node[node_ids[i]].x;
    	nds[i].y = grid->node[node_ids[i]].y;
    	nds[i].z = grid->node[node_ids[i]].z;
    }
	selem1d_load(&(grid->elem1d[1]), 0, 0, 2, node_ids, 0, nds, 0);
	
	node_ids[0] = 2; node_ids[1] = 3;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = grid->node[node_ids[i]].x;
    	nds[i].y = grid->node[node_ids[i]].y;
    	nds[i].z = grid->node[node_ids[i]].z;
    }
	selem1d_load(&(grid->elem1d[2]), 0, 0, 2, node_ids, 0, nds, 0);

	node_ids[0] = 3; node_ids[1] = 4;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = grid->node[node_ids[i]].x;
    	nds[i].y = grid->node[node_ids[i]].y;
    	nds[i].z = grid->node[node_ids[i]].z;
    }
	selem1d_load(&(grid->elem1d[3]), 0, 0, 2, node_ids, 0, nds, 0);

	node_ids[0] = 4; node_ids[1] = 5;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = grid->node[node_ids[i]].x;
    	nds[i].y = grid->node[node_ids[i]].y;
    	nds[i].z = grid->node[node_ids[i]].z;
    }
	selem1d_load(&(grid->elem1d[4]), 0, 0, 2, node_ids, 0, nds, 0);

	//Downstream Edges
	node_ids[0] = 600; node_ids[1] = 601;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = grid->node[node_ids[i]].x;
    	nds[i].y = grid->node[node_ids[i]].y;
    	nds[i].z = grid->node[node_ids[i]].z;
    }
	selem1d_load(&(grid->elem1d[5]), 0, 0, 2, node_ids, 0, nds, 0);
	
	node_ids[0] = 601; node_ids[1] = 602;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = grid->node[node_ids[i]].x;
    	nds[i].y = grid->node[node_ids[i]].y;
    	nds[i].z = grid->node[node_ids[i]].z;
    }
	selem1d_load(&(grid->elem1d[6]), 0, 0, 2, node_ids, 0, nds, 0);
	
	node_ids[0] = 602; node_ids[1] = 603;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = grid->node[node_ids[i]].x;
    	nds[i].y = grid->node[node_ids[i]].y;
    	nds[i].z = grid->node[node_ids[i]].z;
    }
	selem1d_load(&(grid->elem1d[7]), 0, 0, 2, node_ids, 0, nds, 0);

	node_ids[0] = 603; node_ids[1] = 604;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = grid->node[node_ids[i]].x;
    	nds[i].y = grid->node[node_ids[i]].y;
    	nds[i].z = grid->node[node_ids[i]].z;
    }
	selem1d_load(&(grid->elem1d[8]), 0, 0, 2, node_ids, 0, nds, 0);

	node_ids[0] = 604; node_ids[1] = 605;
	for (int i=0;i<nnodes_on_elem;i++){
		nds[i].x = grid->node[node_ids[i]].x;
    	nds[i].y = grid->node[node_ids[i]].y;
    	nds[i].z = grid->node[node_ids[i]].z;
    }
	selem1d_load(&(grid->elem1d[9]), 0, 0, 2, node_ids, 0, nds, 0);


	//reorder grid with Gibbs Poole
    int nnodes;
    nnodes = grid->nnodes;
    sgrid_reorder(grid,2);
    printf("Where are we\n");

	//create a supermodel with a grid in it
	//sm.grid = &grid;
    SMODEL_DESIGN dm;
	//specify elemental physics and other properties in super model
	double dt = 600.0;
	double t0 = 0.0;
	double tf = 864000.0;

	char elemVarCode[4]; 
	strcpy(&elemVarCode[0],"2");//SW2D
	strcpy(&elemVarCode[1],"0"); //GW
	strcpy(&elemVarCode[2],"0"); //Transport
	printf("GRID NELEMS2D = %d\n",grid->nelems2d);
	//smodel_super_no_read_simple(&sm, dt, t0, tf, 0 , 1, 0, elemVarCode);
	smodel_design_no_read_simple(&dm, dt, t0, tf,1, 1, 0, elemVarCode, grid);
	printf("NDOFS %d\n",dm.ndofs[0]);

	//use pointer for short hand
	SMODEL_SUPER *sm = &(dm.superModel[0]);

	//OVER WRITE TO SW2 and the BCs
	sm->elem2d_physics_mat[0].model[0].fe_resid = SW2;
	sm->elem2d_physics_mat[0].model[0].fe_init = SW2_INIT;
	sm->elem2d_physics_mat[0].model[0].nvar = 3;
	(dm.superModel[0].elem2d_physics_mat[0].model[0].physics_vars) = (int*) tl_realloc(sizeof(int), 3,1, dm.superModel[0].elem2d_physics_mat[0].model[0].physics_vars);
    sm->elem2d_physics_mat[0].model[0].physics_vars[0] = PERTURB_H;
    sm->elem2d_physics_mat[0].model[0].physics_vars[1] = PERTURB_U;
    sm->elem2d_physics_mat[0].model[0].physics_vars[2] = PERTURB_V;
    sm->LINEAR_PROBLEM = NO;
    sm->elem1d_physics_mat[0].model[0].fe_resid = SW2_BC_FLUX;
//    //need a skip routine (keep it null for now)
//	sm->elem2d_physics_mat[0].model[0].fe_init = SW2_INIT;
	sm->elem1d_physics_mat[0].model[0].nvar = 3;
	(dm.superModel[0].elem1d_physics_mat[0].model[0].physics_vars) = (int*) tl_realloc(sizeof(int), 3,1, dm.superModel[0].elem1d_physics_mat[0].model[0].physics_vars);
    sm->elem1d_physics_mat[0].model[0].physics_vars[0] = PERTURB_H;
    sm->elem1d_physics_mat[0].model[0].physics_vars[1] = PERTURB_U;
    sm->elem1d_physics_mat[0].model[0].physics_vars[2] = PERTURB_V;

    //hack together the sw2 structure
    //allocate
    sm->sw = (SSW*) tl_alloc(sizeof(SSW), 1);
    //use an alias
    SSW *sw = dm.superModel[0].sw; //alias for convenience
    //set it up
    printf("allocating ssw\n");
    

    //new way
    //ssw_alloc_init(sw, grid->nnodes, grid->nnodes, grid->nelems2d, 0, 0, 1);
    //Old way
    ssw_alloc_init(sw);
    //set up the dvar map
    printf("allocating sdvar\n");
    sdvar_alloc_init(&(sw->dvar), grid->nnodes, 0, 0, 1, grid->nnodes, grid->nelems2d);
    //hard code set index for wd flag
    sw->WD_FLAG = 0;


    //must also set up the dacont extra terms
    //could be clever hear and pick nnode on each element
    sw->elem_rhs_dacont_extra_terms = (double **) tl_alloc(sizeof(double *), grid->nelems2d);
    for (int i=0; i<grid->nelems2d; i++) {
        sw->elem_rhs_dacont_extra_terms[i] = (double *) tl_alloc(sizeof(double), grid->elem2d[i].nnodes);
        for (int j=0; j<grid->elem2d[i].nnodes; j++) {
        	sw->elem_rhs_dacont_extra_terms[i][j] = 0.;
        }
    }
    //printf("sw vals %d\n",sw->WD_FLAG);
    printf("dvar?? %d\n",sw->dvar.n_dvar_elem_int);
    //initial conditions and things
    // intialize dirichlet and old sol (initial guess)
	for (int local_index=0; local_index<dm.ndofs[0]; local_index++){
		sm->dirichlet_data[local_index] = 0.0;
		sm->sol_old[local_index] = 0.0;
		sm->sol[local_index] = 0.0;
		sm->lin_sys->dsol[local_index] = 0.0;
		sm->bc_mask[local_index] = NO;
	}

	//overwrite intial condition
	double x_coord, y_coord, z_coord;
	int id;
	for (int i=0; i<nnodes; i++){
		//mark the boundary only
		x_coord = grid->node[i].x;
		y_coord = grid->node[i].y;
		z_coord = grid->node[i].z;
		//id = grid->node[i].id;
		id=i;
		//need to set IC
		sm->sol[id*3] = 1.0 - z_coord;
		sm->sol_old[id*3] = 1.0 - z_coord;
		sm->sol_older[id*3] = 1.0 - z_coord;
		sm->dirichlet_data[id*3] = 1.0 - z_coord;
		//no dirichlet condition
		//if ( is_near(x_coord,xmin) || is_near(x_coord,xmax) || is_near(y_coord,ymin) || is_near(y_coord,ymax) ){
		//	continue;
		//}else{
		//	dm.superModel[0].bc_mask[id*3+1]=NO;
		//}
		//printf("Dirichlet data node[%d] = %f\n", i, sm.dirichlet_data[i*3+1]);
	}

	//create a series and string for boundary data without file read
	//first add a string, this will come from smat physics file
	//! Downstream edges
	//EGS 1 2 2 
	//EGS 2 3 2 
	//EGS 3 4 2 
	//EGS 4 5 2 
	//EGS 5 6 2 
	//! Upstream edges
	//EGS 601 602 3 
	//EGS 602 603 3 
	//EGS 603 604 3 
	//EGS 604 605 3 
	//EGS 605 606 3
	//add this to  



	//printf("Supermodel no read complete\n");


	//allocate linear system
	//doesn't currently work, need to go back and fix
	//fe_allocate_initialize_linear_system(&sm);
//	sm.cols_diag = NULL;
//	sm.vals_diag = NULL;
//	sm.cols_off_diag=NULL;
//	sm.vals_off_diag=NULL;
//	sm.bc_mask = NULL;
//	sm.nnz_diag_old=0;
//	sm.nnz_off_diag_old=0;
//	printf("Calling sparsity split CSR\n");
//	create_sparsity_split_CSR(&sm, sm.grid);
	//Screen_print_CSR(sm.indptr_diag, sm.cols_diag, sm.vals_diag, sm.ndofs);

	//do we want to stor nnz? it is stored in sm->indptr[nrows]
    //do we want to store local_size = local_range[1]-local_range[0]
//    printf("NNZ = %d = %d\n",sm.indptr_diag[sm.my_ndofs], sm.nnz_diag);
//    sarray_init_dbl(sm.vals_diag, sm.indptr_diag[sm.my_ndofs]);//

//	//manually set up some solver parameters
//	sm.macro_ndofs = sm.ndofs;
//	sm.dsol = (double*) tl_alloc(sizeof(double), sm.ndofs);
//	
//	sm.dirichlet_data = (double*) tl_alloc(sizeof(double), sm.ndofs);
//	sm.scale_vect = (double*) tl_alloc(sizeof(double), sm.ndofs);
////    for(int i=0;i<sm.ndofs;i++){
////    	sm.scale_vect[i] = 1.0;
////    }
//	sm.tol_nonlin = 1e-5;
//	sm.inc_nonlin = 1e-3;
//	sm.max_nonlin_linesearch_cuts = 5;
//	sm.it_count_nonlin_failed = 0;
//	sm.max_nonlin_it = 20;
//	sm.LINEAR_PROBLEM = NO;
//	sm.force_nonlin_it = NO;
//	sm.force_nonlin_it = NO;
//	sm.nonlinear_it_total = 0;
//	sm.nghost=0;
//	sm.ghosts = NULL;
//	sm.local_size = sm.ndofs;
//	sm.sol_old = (double*) tl_alloc(sizeof(double), sm.ndofs);
//	sm.sol_older = (double*) tl_alloc(sizeof(double), sm.ndofs);

	//set up bc mask

//	for (int i=0;i<sm.ndofs;i++){
//		printf("sm bc mask[%d] = %d\n",i,sm.bc_mask[i]);
//	}


	//see if it works
	//apply_Dirichlet_BC(&sm);
	//see if something is happening within newton loop or something else
//	initialize_system(&sm);
//	assemble_residual(&sm,sm.grid);
//	assemble_jacobian(&sm,sm.grid);
//	apply_Dirichlet_BC(&sm);
//	int status;
//	status = prep_umfpack(sm.indptr_diag,sm.cols_diag,sm.vals_diag, sm.dsol, sm.residual, sm.local_size);
//	status = solve_umfpack(sm.dsol, sm.indptr_diag, sm.cols_diag, sm.vals_diag, sm.residual, sm.local_size);
//	increment_function(&sm);
	//Screen_print_CSR(sm.indptr_diag, sm.cols_diag, sm.vals_diag, sm.ndofs);
	printf("Calling time loop\n");
	//set forward_step and call timeloop
	time_loop(&dm); 

	//get initial mass
	double *init_head;
    init_head = (double *) tl_alloc(sizeof(double), nnodes);
    SVECT2D *init_vel;
    init_vel = (SVECT2D *) tl_alloc(sizeof(SVECT2D), nnodes);
    //fill in temp arrays
    for (int inode=0; inode<nnodes; inode++){
    	init_head[inode] = dm.superModel[0].dirichlet_data[inode*3];
    	init_vel[inode].x = dm.superModel[0].dirichlet_data[inode*3+1];
    	init_vel[inode].y = dm.superModel[0].dirichlet_data[inode*3+2];
    }
    SFLAGS dummy;
	double initial_grid_mass = tl_find_grid_mass_elem2d(dm.superModel[0].density, NULL, NULL, init_head,dm.superModel[0].grid, dummy);
	//printf("Initial grid mass = %f\n",initial_grid_mass);
	//extract second variable here
	double total_error = write_testcase_error_wet_dry(&(dm.superModel[0]),initial_grid_mass); 
	printf("Final error: %6.4e\n", total_error);
	//plot grid in h5?
//    strcpy(sm.grid->filename, "residtest");
//    init_hdf5_file(sm.grid);
//    printf("hdf5 initialized\n");
//    sgrid_write_hdf5(sm.grid);
//    printf("hdf5 written\n");
//    sgrid_write_xdmf(sm.grid);
//    printf("xmf written\n");

	//return -1 if failed, 0 if good
	
	if(total_error < NEWTON_TEST_TOL){
		err_code=0;
	}
	//printf("Final error code %d\n",err_code);




	smodel_design_free(&dm);


	return err_code;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes analytic solution used to calculate error
 *  solution
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double write_testcase_error_wet_dry(SMODEL_SUPER *mod, double initial_grid_mass) {
    int inode;
    double error_vx = 0., error_vy = 0., error_h = 0., max_head = 0., max_u = 0., max_v = 0., total_time_mass_flux= 0.0;
    double error_vx_max, error_vy_max, error_h_max, max_head_grid, max_u_grid, max_v_grid;

    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        error_vx += fabs(mod->sol[inode*3+1]); // water initially at rest
        error_vy += fabs(mod->sol[inode*3+2]); // water initially at rest
        error_h  += fabs(mod->sol[inode*3] - mod->dirichlet_data[inode*3]);
        if (fabs(mod->dirichlet_data[inode*3]) > max_head) max_head = fabs(mod->dirichlet_data[inode*3]);
    }
    //NEED TO DO
    double *head;
    head = (double *) tl_alloc(sizeof(double), mod->grid->nnodes);
    SVECT2D *vel;
    vel = (SVECT2D *) tl_alloc(sizeof(SVECT2D), mod->grid->nnodes);
    //fill in temp arrays
    for (inode=0; inode<mod->grid->my_nnodes; inode++){
    	head[inode] = mod->sol[inode*3];
    	vel[inode].x = mod->sol[inode*3+1];
    	vel[inode].y = mod->sol[inode*3+2];
    }
    SFLAGS dummy;
    double grid_mass_error = tl_find_grid_mass_error_elem2d(mod->density, head, vel, mod->grid, dummy, initial_grid_mass, NULL, NULL, *(mod->dt), &total_time_mass_flux);
    
    printf("Grid mass error %6.4e\n", grid_mass_error);


    head = tl_free(sizeof(double),mod->grid->nnodes,head);
    vel = tl_free(sizeof(SVECT2D),mod->grid->nnodes,vel);

    return error_vx + error_vy + error_h;

//    if (max_head < 1e-6) max_head = 1.;
//    if (max_u < 1e-6) max_u = 1.;
//    if (max_v < 1e-6) max_v = 1.;
//#ifdef _MESSG
//    error_vx_max=messg_dsum(error_vx,mod->grid->smpi->ADH_COMM);
//    error_vy_max=messg_dsum(error_vy,mod->grid->smpi->ADH_COMM);
//    error_h_max=messg_dsum(error_h,mod->grid->smpi->ADH_COMM);
//    max_head_grid=messg_dmax(max_head,mod->grid->smpi->ADH_COMM);
//    max_u_grid=messg_dmax(max_u,mod->grid->smpi->ADH_COMM);
//    max_v_grid=messg_dmax(max_v,mod->grid->smpi->ADH_COMM);
//    error_vx=error_vx_max;
//    error_vy=error_vy_max;
//    error_h=error_h_max;
//    max_head=max_head_grid;
//    max_u=max_u_grid;
//    max_v=max_v_grid;
//#endif
//    if(mod->grid->smpi->myid==0){
//        FILE *fp;
//        fp = fopen("error.out", "w");
//        fprintf(fp,"x-velocity abs error: %30.20e :: relative error: %30.20e \n",error_vx/mod->grid->macro_nnodes,error_vx/mod->grid->macro_nnodes/max_u);
//        fprintf(fp,"y-velocity abs error: %30.20e :: relative error: %30.20e \n",error_vy/mod->grid->macro_nnodes,error_vy/mod->grid->macro_nnodes/max_v);
//        fprintf(fp,"head abs error: %30.20e :: relative error: %30.20e \n",error_h/mod->grid->macro_nnodes,error_h/mod->grid->macro_nnodes/max_head);
//        fprintf(fp,"grid_mass_error: %30.20e :: relative error: %30.20e \n",grid_mass_error, grid_mass_error/mod->initial_grid_mass);
//        fclose(fp);
//    }
}
void permute_array(double *arr,int *p, int n){
	double *temp;
	temp = (double *) tl_alloc(sizeof(double),n);
	for(int i =0;i<n;i++){
		temp[p[i]] = arr[i];
	}
	// Copy permuted elements back to the original array
    for (int i = 0; i < n; i++) {
        arr[i] = temp[i];
    }
	temp = (double *) tl_free(sizeof(double),n,temp);
}
