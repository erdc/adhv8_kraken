/*! \file sw2_wd_test.c This file tests the sw2 engine */
#include "adh.h"
static double NEWTON_TEST_TOL = 1e-7;
static int NEWTON_TEST_NX = 16;
static int NEWTON_TEST_NY = 6;
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
int sw2_wd_test(int argc, char **argv) {
	//create a grid
	SGRID *grid;
	grid = (SGRID *) tl_alloc(sizeof(SGRID), 1);
	//let's just do something simple
	//3x3 triangular element grid
	double xmin = 0.0;
	double xmax = 1000.0;
	double ymin = 0.0;
	double ymax = 500.0;
	int npx = NEWTON_TEST_NX;
	int npy = NEWTON_TEST_NY;
	double theta = 0.0;
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
    int nnodes;
    nnodes = grid->nnodes;
    sgrid_reorder(grid,2);
	//print coordinates
//  for(int local_index =0; local_index<grid.nnodes; local_index++){
//		printf("Node %d: (x,y) = {%f,%f}\n",grid.node[local_index].gid,grid.node[local_index].x,grid.node[local_index].y);
//	}
//	//print connectivity
//	for(int local_index =0; local_index<grid.nelems2d; local_index++){
//		printf("Element %d: (nd1,nd2,nd3) = {%d,%d,%d}\n",local_index ,grid.elem2d[local_index].nodes[0], grid.elem2d[local_index].nodes[1], grid.elem2d[local_index].nodes[2]);
//	}


	//create a supermodel with a grid in it
	//SMODEL_SUPER sm;
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
	smodel_design_no_read_simple(&dm, dt, t0, tf,0, 1, 0, elemVarCode, grid);
	printf("NDOFS %d\n",dm.ndofs[0]);

	//OVER WRITE TO SW2
	dm.superModel[0].elem2d_physics_mat[0].model[0].fe_resid = SW2;
	dm.superModel[0].elem2d_physics_mat[0].model[0].fe_init = SW2;
	dm.superModel[0].elem2d_physics_mat[0].model[0].nvar = 3;
    dm.superModel[0].elem2d_physics_mat[0].model[0].physics_vars[0] = PERTURB_H;
    dm.superModel[0].elem2d_physics_mat[0].model[0].physics_vars[1] = PERTURB_U;
    dm.superModel[0].elem2d_physics_mat[0].model[0].physics_vars[2] = PERTURB_V;
    dm.superModel[0].LINEAR_PROBLEM = NO;

    //hack together the sw2 structure
    //allocate
    dm.superModel[0].sw = (SSW*) tl_alloc(sizeof(SSW), 1);
    //use an alias
    SSW *sw = dm.superModel[0].sw; //alias for convenience
    //set it up
    printf("allocating ssw\n");
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
    printf("sw vals %d\n",sw->WD_FLAG);
    printf("dvar?? %d\n",sw->dvar.n_dvar_elem_int);
    //initial conditions and things
    // intialize dirichlet and old sol (initial guess)
	for (int local_index=0; local_index<dm.ndofs[0]; local_index++){
		dm.superModel[0].dirichlet_data[local_index] = 0.0;
		dm.superModel[0].sol_old[local_index] = 0.0;
		dm.superModel[0].sol[local_index] = 0.0;
		dm.superModel[0].lin_sys->dsol[local_index] = 0.0;
		dm.superModel[0].bc_mask[local_index] = NO;
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
		dm.superModel[0].sol[id*3] = 1.0 - z_coord;
		dm.superModel[0].sol_old[id*3] = 1.0 - z_coord;
		dm.superModel[0].sol_older[id*3] = 1.0 - z_coord;
		dm.superModel[0].dirichlet_data[id*3] = 1.0 - z_coord;
		//no dirichlet condition
		//if ( is_near(x_coord,xmin) || is_near(x_coord,xmax) || is_near(y_coord,ymin) || is_near(y_coord,ymax) ){
		//	continue;
		//}else{
		//	dm.superModel[0].bc_mask[id*3+1]=NO;
		//}
		//printf("Dirichlet data node[%d] = %f\n", i, sm.dirichlet_data[i*3+1]);
	}

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
