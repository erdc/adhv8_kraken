#include "global_header.h"
#include <complex.h>    /* Standard Library of Complex Numbers */

// SOME OF THESE SHOULD ONLY BE CALLED AT THE END OF THE SIMULATION
// CJT :: THESE REALLY NEED TO BE INSIDE THE MODEL STRUCT!!!

static int DEBUG = ON;
static double pi = 3.14159265359;

static int ALREADY_RECORDED = NO;
static int MY_NODE_PE = UNSET_INT; // cjt :: for MPI testing when a node is given

static double error_dpl_tmax = 0;
static double error_h_tmax = 0;
static double error_u_tmax = 0;
static double error_v_tmax = 0;
static double error_w_tmax = 0;
static double error_c_tmax = 0;

static double rmse_error_dpl_tmax = 0;
static double rmse_error_u_tmax = 0;
static double rmse_error_v_tmax = 0;
static double rmse_error_w_tmax = 0;

static double error_dpl_tmax_3d = 0;
static double error_u_tmax_3d = 0;
static double error_w_tmax_3d = 0;

static double rmse_error_dpl_tmax_3d = 0;
static double rmse_error_u_tmax_3d = 0;
static double rmse_error_w_tmax_3d = 0;

static double error_dpl_tmax_2d = 0;
static double error_u_tmax_2d = 0;

static double rmse_error_dpl_tmax_2d = 0;
static double rmse_error_u_tmax_2d = 0;

static double *initial_depth = NULL;
static int *nodeIDs_on_slice = NULL;
static int nnodes_Xsection = 0;
static int nnodes = 0;

static double total_time_mass_flux = 0.;

static int first_time_1 = YES;
static int first_time_2 = YES;
static int first_time_3 = YES;

/* prototypes ----------------------------------------------------*/
void write_testcase_error_ptest(SMODEL *);
void write_testcase_error_coriolis(SMODEL *);
void write_testcase_error_xcoriolis(SMODEL *);
void write_testcase_error_ycoriolis(SMODEL *);
void write_testcase_error_discharge(SMODEL *);
void write_testcase_error_slosh(SMODEL *);
void write_testcase_error_salt(SMODEL *);
void write_testcase_error_winds_Huang(SMODEL *);
void write_testcase_error_water_source(SMODEL *);
void write_testcase_error_nb(SMODEL *);
void write_testcase_error_lock(SMODEL *);
void write_testcase_error_winds_and_waves(SMODEL *);
void write_testcase_error_wet_dry(SMODEL *);
void write_testcase_error_dam2d(SMODEL *);
void write_testcase_error_steepfront(SMODEL *);
void write_testcase_error_circular_slosh(SMODEL *);
void write_testcase_error_outflow(SMODEL *);
void write_testcase_error_slosh2d3d(SSUPER_MODEL *); // gkc
void write_testcase_error_tide_2d(SMODEL *mod);
void write_testcase_error_tide_3d(SMODEL *);
void write_testcase_error_dwe_hunter(SMODEL *);

void testcase_prep_lock(SMODEL *);  // initializes salt concentration
void testcase_prep_slosh(SMODEL *); // initializes water displacement and velocities
void testcase_prep_xcoriolis(SMODEL *);
void testcase_prep_ycoriolis(SMODEL *);
void testcase_prep_winds_and_waves(SMODEL *);
void testcase_prep_wet_dry(SMODEL *);
void testcase_prep_dam2d(SMODEL *);
void testcase_prep_circular_slosh(SMODEL *);
void testcase_prep_slosh2d3d(SSUPER_MODEL *); // initializes water displacement and velocities // gkc
void testcase_prep_tide_2d(SMODEL *);
void testcase_prep_tide_3d(SMODEL *);
void testcase_prep_dwe_hunter(SMODEL *);

//****************************************************************//
//****************************************************************//

bool isPointOnLine(SVECT p1, SVECT p2, SVECT pointCheck) {
    double d1 = sqrt(pow(p1.x - pointCheck.x,2) + pow(p1.y - pointCheck.y,2));
    double d2 = sqrt(pow(p2.x - pointCheck.x,2) + pow(p2.y - pointCheck.y,2));
    double d3 = sqrt(pow(p1.x - p2.x,2) + pow(p1.y - p2.y,2));
    if (fabs((d1 + d2) - d3) < SMALL6) {return TRUE;}
    return FALSE;
}


//****************************************************************//
//****************************************************************//
// initialize (hot-start) test cases here so not changing input constantly
void testcase_prep(SSUPER_MODEL *sm) {
    int i,j,k;
    SMODEL *mod;
    
#ifdef _MESSG
    // CJT :: distribute testcase parameters to all the processors if 2d/3d coupling in HPC
    // This will need to be done if there are coupled subsmodels where one model calls a test case.
    MPI_Allreduce(&test_case_flag.slosh2d3d, &test_case_flag.slosh2d3d, 1, MPI_INT, MPI_MAX, cstorm_comm);
    MPI_Allreduce(&test_case_flag.a, &test_case_flag.a, 1, MPI_DOUBLE, MPI_MAX, cstorm_comm);
    MPI_Allreduce(&test_case_flag.L, &test_case_flag.L, 1, MPI_DOUBLE, MPI_MAX, cstorm_comm);
    MPI_Allreduce(&test_case_flag.B, &test_case_flag.B, 1, MPI_DOUBLE, MPI_MAX, cstorm_comm);
    MPI_Allreduce(&test_case_flag.H, &test_case_flag.H, 1, MPI_DOUBLE, MPI_MAX, cstorm_comm);
    
    // move test nodes to their local IDs
//    int nodeID_flag = 0, node1_flag = 0, node2_flag = 0, node2d_flag = 0, node3d_flag = 0;
//    for (j=0; j<sm->nsubmodels; j++) {
//        for (i=0; i<sm->submodel[j].grid->nnodes; i++) {
//            if (sm->submodel[j].grid->node[i].resident_pe == sm->submodel[j].grid->smpi->myid) {
//
//                if (sm->submodel[j].grid->node[i].gid == 6558) printf("test_case_flag.nodeID: %d x: %f %f i: %d myid: %d\n",
//                                                                      test_case_flag.node1,sm->submodel[j].grid->node[i].x,sm->submodel[j].grid->node[i].y
//                                                                      ,i,sm->submodel[j].grid->smpi->myid);
//
//                MY_NODE_PE = sm->submodel[j].grid->smpi->myid;
//                if (sm->submodel[j].grid->node[i].gid == test_case_flag.nodeID && nodeID_flag == 0)  {test_case_flag.nodeID = i; nodeID_flag = 1;};
//                if (sm->submodel[j].grid->node[i].gid == test_case_flag.node1 && node1_flag == 0)    {test_case_flag.node1  = i; node1_flag = 1;};
//                if (sm->submodel[j].grid->node[i].gid == test_case_flag.node2 && node2_flag == 0)    {test_case_flag.node2  = i; node2_flag = 1;};
//                if (sm->submodel[j].grid->ndim == 2 && sm->submodel[j].grid->node[i].gid == test_case_flag.nodeID_2d && node2d_flag == 0)  {
//                    test_case_flag.nodeID_2d = i;
//                    test_case_flag.node1_2d  = i;
//                    test_case_flag.node2_2d  = i;
//                    node2d_flag = 1;
//                }
//                if (sm->submodel[j].grid->ndim == 3 && sm->submodel[j].grid->node[i].gid == test_case_flag.nodeID_3d && node3d_flag == 0)  {
//                    test_case_flag.nodeID_3d = i;
//                    test_case_flag.node1_3d  = i;
//                    test_case_flag.node2_3d  = i;
//                    node3d_flag = 1;
//                }
//
//            }
//        }
//    }
//
//    for (j=0; j<sm->nsubmodels; j++) {
//        for (i=0; i<sm->submodel[j].grid->nnodes; i++) {
//            if (sm->submodel[j].grid->node[i].resident_pe == sm->submodel[j].grid->smpi->myid) {
//                if (sm->submodel[j].grid->node[i].gid == 6558) printf("test_case_flag.nodeID: %d x: %f %f i: %d myid: %d\n",
//                                                                      test_case_flag.node1,sm->submodel[j].grid->node[i].x,sm->submodel[j].grid->node[i].y
//                                                                      ,i,sm->submodel[j].grid->smpi->myid);
//            }
//        }
//    }
#endif
//    MPI_Barrier(MPI_COMM_WORLD);
//    tl_error("test");
    
//#ifdef _MESSG // cjt :: do this for coupled models!
//    MPI_Allreduce(&test_case_flag.tide2d,    &test_case_flag.tide2d,    1, MPI_INT, MPI_MAX, cstorm_comm);
//    MPI_Allreduce(&test_case_flag.nodeID_2d, &test_case_flag.nodeID_2d, 1, MPI_INT, MPI_MAX, cstorm_comm);
//    MPI_Allreduce(&test_case_flag.node1_2d,  &test_case_flag.node1_2d,  1, MPI_INT, MPI_MAX, cstorm_comm);
//    MPI_Allreduce(&test_case_flag.node2_2d,  &test_case_flag.node2_2d,  1, MPI_INT, MPI_MAX, cstorm_comm);
//
//    MPI_Allreduce(&test_case_flag.tide3d,    &test_case_flag.tide3d,    1, MPI_INT, MPI_MAX, cstorm_comm);
//    MPI_Allreduce(&test_case_flag.nodeID_3d, &test_case_flag.nodeID_3d, 1, MPI_INT, MPI_MAX, cstorm_comm);
//    MPI_Allreduce(&test_case_flag.node1_3d,  &test_case_flag.node1_3d,  1, MPI_INT, MPI_MAX, cstorm_comm);
//    MPI_Allreduce(&test_case_flag.node2_3d,  &test_case_flag.node2_3d,  1, MPI_INT, MPI_MAX, cstorm_comm);
//#endif
    
    if (test_case_flag.lock) {
        for (j=0; j<sm->nsubmodels; j++) testcase_prep_lock(&(sm->submodel[j]));
    } else if (test_case_flag.slosh) {
        for (j=0; j<sm->nsubmodels; j++) testcase_prep_slosh(&(sm->submodel[j]));
    } else if (test_case_flag.xcoriolis) {
        for (j=0; j<sm->nsubmodels; j++) testcase_prep_xcoriolis(&(sm->submodel[j]));
    } else if (test_case_flag.ycoriolis) {
        for (j=0; j<sm->nsubmodels; j++) testcase_prep_ycoriolis(&(sm->submodel[j]));
    } else if (test_case_flag.winds_and_waves) {
        for (j=0; j<sm->nsubmodels; j++) testcase_prep_winds_and_waves(&(sm->submodel[j]));
    } else if (test_case_flag.wet_dry) {
        for (j=0; j<sm->nsubmodels; j++) testcase_prep_wet_dry(&(sm->submodel[j]));
    } else if (test_case_flag.dam2d) {
        for (j=0; j<sm->nsubmodels; j++) testcase_prep_dam2d(&(sm->submodel[j]));
    } else if (test_case_flag.cslosh) {
        for (j=0; j<sm->nsubmodels; j++) testcase_prep_circular_slosh(&(sm->submodel[j]));
    } else if (test_case_flag.slosh2d3d) {
        testcase_prep_slosh2d3d(sm);
    } else if (test_case_flag.tide2d || test_case_flag.tide3d) {
        for (j=0; j<sm->nsubmodels; j++) {
            if (sm->submodel[j].grid->ndim == 3 && test_case_flag.tide3d) {
                testcase_prep_tide_3d(&(sm->submodel[j]));
            }
            if (sm->submodel[j].grid->ndim == 2 && test_case_flag.tide2d) {
                testcase_prep_tide_2d(&(sm->submodel[j]));
            }
        }
    } else if (test_case_flag.dwe_hunter) {
        for (j=0; j<sm->nsubmodels; j++) testcase_prep_dwe_hunter(&(sm->submodel[j]));
    }
    
    // calculate initial grid mass after prepping testcase, in case initial fields have changed
    for (j=0; j<sm->nsubmodels; j++) {
        mod = &(sm->submodel[j]);
        if(mod->proc_flag == 1){
            if (mod->flag.SW2_FLOW || mod->flag.DIFFUSIVE_WAVE) {
                mod->initial_grid_mass = tl_find_grid_mass_elem2d(mod->density, mod->str_values, mod->series_head, mod->sw->d2->head, mod->grid, mod->flag);
            } else if (mod->flag.SW3_FLOW) {
                mod->initial_grid_mass = tl_find_grid_mass_elem3d(mod->density, mod->grid, mod->sw->d3->displacement);
                tl_calculate_depthavgvel(mod->grid, mod->sw->d3);
            } else if (mod->flag.NS3_FLOW) {
                mod->initial_grid_mass = tl_find_grid_mass_elem3d(mod->density, mod->grid, mod->ns->d3->displacement);
            } else {
                tl_error("ERROR: Model must have a physics type associated with it!\n");
            }
        }
    }
}

//****************************************************************//
//****************************************************************//

void write_testcase_error(SSUPER_MODEL *sm) {
    int j;
    if (test_case_flag.xcoriolis) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_xcoriolis(&(sm->submodel[j]));
    } else if (test_case_flag.ycoriolis) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_ycoriolis(&(sm->submodel[j]));
    } else if (test_case_flag.coriolis) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_coriolis(&(sm->submodel[j]));
    } else if (test_case_flag.ptest) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_ptest(&(sm->submodel[j]));
    } else if (test_case_flag.discharge) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_discharge(&(sm->submodel[j]));
    } else if (test_case_flag.slosh) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_slosh(&(sm->submodel[j]));
    } else if (test_case_flag.salt) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_salt(&(sm->submodel[j]));
    } else if (test_case_flag.winds_Huang) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_winds_Huang(&(sm->submodel[j]));
    } else if (test_case_flag.water_source) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_water_source(&(sm->submodel[j]));
    } else if (test_case_flag.nb) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_nb(&(sm->submodel[j]));
    } else if (test_case_flag.lock && ALREADY_RECORDED == NO) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_lock(&(sm->submodel[j]));
    } else if (test_case_flag.winds_and_waves) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_winds_and_waves(&(sm->submodel[j]));
    } else if (test_case_flag.wet_dry) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_wet_dry(&(sm->submodel[j]));
    } else if (test_case_flag.dam2d) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_dam2d(&(sm->submodel[j]));
    } else if (test_case_flag.steep_front) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_steepfront(&(sm->submodel[j]));
    } else if (test_case_flag.cslosh) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_circular_slosh(&(sm->submodel[j]));
    } else if (test_case_flag.outflow) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_outflow(&(sm->submodel[j]));
    } else if (test_case_flag.slosh2d3d) {
        write_testcase_error_slosh2d3d(sm);
    } else if (test_case_flag.tide2d || test_case_flag.tide3d) {
        for (j=0; j<sm->nsubmodels; j++) {
            if (sm->submodel[j].grid->ndim == 3 && test_case_flag.tide3d) {
                write_testcase_error_tide_3d(&(sm->submodel[j]));
            }
            if (sm->submodel[j].grid->ndim == 2 && test_case_flag.tide2d) {
                write_testcase_error_tide_2d(&(sm->submodel[j]));
            }
        }
    } else if (test_case_flag.dwe_hunter) {
        for (j=0; j<sm->nsubmodels; j++) write_testcase_error_dwe_hunter(&(sm->submodel[j]));
    }
}

//****************************************************************//
//****************************************************************//

void testcase_clean() {
    if (nodeIDs_on_slice != NULL) {
        nodeIDs_on_slice = (int*) tl_free(sizeof(int), nnodes_Xsection, nodeIDs_on_slice);
    }
    if (initial_depth != NULL) {
        initial_depth = (double *) tl_free(sizeof(double), nnodes, initial_depth);
    }
}

//****************************************************************//
//****************************************************************//
void testcase_prep_wet_dry(SMODEL *mod) {
    initial_depth = (double *) tl_alloc(sizeof(double), mod->grid->nnodes);
    int inode;
    for (inode=0; inode<mod->grid->nnodes; inode++) {
        initial_depth[inode] = mod->sw->d2->head[inode];
    }
}

//****************************************************************//
void write_testcase_error_wet_dry(SMODEL *mod) {
    int inode;
    double error_vx = 0., error_vy = 0., error_h = 0., max_head = 0., max_u = 0., max_v = 0.;
    double error_vx_max, error_vy_max, error_h_max, max_head_grid, max_u_grid, max_v_grid;
    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        error_vx += fabs(mod->sw->d2->vel[inode].x); // water initially at rest
        error_vy += fabs(mod->sw->d2->vel[inode].y); // water initially at rest
        error_h  += fabs(mod->sw->d2->head[inode] - initial_depth[inode]);
        if (fabs(initial_depth[inode]) > max_head) max_head = fabs(initial_depth[inode]);
    }
    double grid_mass_error = tl_find_grid_mass_error_elem2d(mod->density, mod->sw->d2->head, mod->sw->d2->vel, mod->grid, mod->flag, mod->initial_grid_mass, mod->series_head, mod->str_values, mod->dt, &total_time_mass_flux);
    
    if (max_head < 1e-6) max_head = 1.;
    if (max_u < 1e-6) max_u = 1.;
    if (max_v < 1e-6) max_v = 1.;
#ifdef _MESSG
    error_vx_max=messg_dsum(error_vx,mod->grid->smpi->ADH_COMM);
    error_vy_max=messg_dsum(error_vy,mod->grid->smpi->ADH_COMM);
    error_h_max=messg_dsum(error_h,mod->grid->smpi->ADH_COMM);
    max_head_grid=messg_dmax(max_head,mod->grid->smpi->ADH_COMM);
    max_u_grid=messg_dmax(max_u,mod->grid->smpi->ADH_COMM);
    max_v_grid=messg_dmax(max_v,mod->grid->smpi->ADH_COMM);
    error_vx=error_vx_max;
    error_vy=error_vy_max;
    error_h=error_h_max;
    max_head=max_head_grid;
    max_u=max_u_grid;
    max_v=max_v_grid;
#endif
    if(mod->grid->smpi->myid==0){
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"x-velocity abs error: %30.20e :: relative error: %30.20e \n",error_vx/mod->grid->macro_nnodes,error_vx/mod->grid->macro_nnodes/max_u);
        fprintf(fp,"y-velocity abs error: %30.20e :: relative error: %30.20e \n",error_vy/mod->grid->macro_nnodes,error_vy/mod->grid->macro_nnodes/max_v);
        fprintf(fp,"head abs error: %30.20e :: relative error: %30.20e \n",error_h/mod->grid->macro_nnodes,error_h/mod->grid->macro_nnodes/max_head);
        fprintf(fp,"grid_mass_error: %30.20e :: relative error: %30.20e \n",grid_mass_error, grid_mass_error/mod->initial_grid_mass);
        fclose(fp);
    }
}


//****************************************************************//
//****************************************************************//
// note :: may get natural errors from water head not remaining flat at in flow edge
void write_testcase_error_nb(SMODEL *mod) {
    
    int inode;
    double convert_to_rads = 3.141592653589793 / 180.;
    
    // the inflow is 0.1 m/s at a 45 degree angle
    double analytic_vx = 0.1/sqrt(2.);
    double analytic_vy = analytic_vx;
    
    double error_vx = 0., error_vx_max;
    double error_vy = 0., error_vy_max;
    double error_h = 0., error_h_max;
    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        error_vx += fabs(mod->sw->d2->vel[inode].x - analytic_vx);
        error_vy += fabs(mod->sw->d2->vel[inode].y - analytic_vy);
        error_h  += fabs(mod->sw->d2->head[inode] - 1.);  // this is just to see how far from rigid lid
    }
    
    double grid_mass_error = tl_find_grid_mass_error_elem2d(mod->density, mod->sw->d2->head, mod->sw->d2->vel, mod->grid, mod->flag, mod->initial_grid_mass, mod->series_head, mod->str_values, mod->dt, &total_time_mass_flux);
    
#ifdef _MESSG
    error_vx_max=messg_dsum(error_vx,mod->grid->smpi->ADH_COMM);
    error_vy_max=messg_dsum(error_vy,mod->grid->smpi->ADH_COMM);
    error_h_max=messg_dsum(error_h,mod->grid->smpi->ADH_COMM);
    error_vx=error_vx_max;
    error_vy=error_vy_max;
    error_h=error_h_max;
#endif
    if(mod->grid->smpi->myid==0){
        
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"x-velocity abs error: %30.20e :: relative error: %30.20e \n",error_vx/mod->grid->macro_nnodes,error_vx/mod->grid->macro_nnodes/analytic_vx);
        fprintf(fp,"y-velocity abs error: %30.20e :: relative error: %30.20e \n",error_vy/mod->grid->macro_nnodes,error_vy/mod->grid->macro_nnodes/analytic_vy);
        fprintf(fp,"head abs error: %30.20e :: relative error: %30.20e \n",error_h/mod->grid->macro_nnodes,error_h/mod->grid->macro_nnodes/1.);
        if (error_vx/mod->grid->macro_nnodes > 1.e-6 || error_vx/mod->grid->macro_nnodes/analytic_vx > 1.e-6)
            fprintf(fp,"VelX FAIL\n");
        else
            fprintf(fp, "PASS\n");
        
        if (error_vy/mod->grid->macro_nnodes > 1.e-6 || error_vy/mod->grid->macro_nnodes/analytic_vy > 1.e-6)
            fprintf(fp,"VelY FAIL\n");
        else
            fprintf(fp, "PASS\n");
        
        fclose(fp);
    }
}

//********************************************************************************/
//********************************************************************************/
// CORIOLIS VERIFICATION
//****************************************************************//
//****************************************************************//

//****************************************************************//
//****************************************************************//
void testcase_prep_xcoriolis(SMODEL *mod) {
    
    int i, id;
    
    // find nodes on user-defined cross-section
    if (mod->grid->ndim == 3) {
        //nnodes_Xsection = tl_find_surface_nodes_on_sliceX(mod->grid, test_case_flag.xslice, &nodeIDs_on_slice);
        nnodes_Xsection = tl_find_surface_nodes_on_sliceX(mod->grid, test_case_flag.node1, test_case_flag.node2, &nodeIDs_on_slice);
    } else {
        //nnodes_Xsection = tl_find_nodes_on_sliceX(mod->grid, test_case_flag.xslice, &nodeIDs_on_slice, 0.);
        nnodes_Xsection = tl_find_nodes_on_sliceX(mod->grid, test_case_flag.node1, test_case_flag.node2, &nodeIDs_on_slice);
    }
    
    printf("\n** prepping angle coriolis test -- angling x direction flume in AdH **");
    printf("** analytic results using the following X-section nodes: **\n");
    if(nnodes_Xsection>0){
        for(i=0; i<nnodes_Xsection; i++) {
            id = nodeIDs_on_slice[i];
            printf("nodeID: %d {x,y,z} = {%20.10e,%20.10e,%20.10e} \n",id+1,mod->grid->node[id].x,mod->grid->node[id].y,mod->grid->node[id].z);
        }
    }
}

//****************************************************************//
//****************************************************************//
void testcase_prep_ycoriolis(SMODEL *mod) {
    
    int i, id;
    
    // find nodes on user-defined cross-section
    if (mod->grid->ndim == 3) {
        nnodes_Xsection = tl_find_surface_nodes_on_sliceY(mod->grid, test_case_flag.node1, test_case_flag.node2, &nodeIDs_on_slice);
    } else {
        nnodes_Xsection = tl_find_nodes_on_sliceY(mod->grid, test_case_flag.node1, test_case_flag.node2, &nodeIDs_on_slice);
    }
    
    printf("\n** analytic results using the following X-section nodes: **\n");
    if(nnodes_Xsection>0){
        for(i=0; i<nnodes_Xsection; i++) {
            id = nodeIDs_on_slice[i];
            printf("nodeID: %d {x,y,z} = {%20.10e,%20.10e,%20.10e} \n",id+1,mod->grid->node[id].x,mod->grid->node[id].y,mod->grid->node[id].z);
        }
    }
}

//****************************************************************//
//****************************************************************//
void write_testcase_error_xcoriolis(SMODEL *mod) {
    
    // this test case finds the water surace slope at an x-transect at y = 50,000
    // if the grid changes, so must this routines!!!
    double convert_to_rads = 3.141592653589793 / 180.;
    double rads = convert_to_rads * mod->mat[0].sw->coriolis;  // assume all materials have same coriolis
    double angular_speed = 2. * EARTH_ROTAT * sin(rads);
    
    int i;
    int nd1 = 0;
    int nd2 = 0;
    int nnodes_Xsection_total;
    double error_max = 0.;
    double error = 0., dhdy = 0., slope = 0.;
    printf("\n");
    if(nnodes_Xsection > 1){
        for (i=0; i<nnodes_Xsection-2; i++) {
            nd1 = nodeIDs_on_slice[i+1];
            nd2 = nodeIDs_on_slice[i];
            if (mod->grid->ndim == 3) {
                slope = - (angular_speed * (mod->sw->d3->vel[nd1].x + mod->sw->d3->vel[nd2].x)/2.)/mod->gravity;
                dhdy = (mod->sw->d3->displacement[nd1] - mod->sw->d3->displacement[nd2])/(mod->grid->node[nd1].y - mod->grid->node[nd2].y);
            } else {
                slope = - (angular_speed * (mod->sw->d2->vel[nd1].x + mod->sw->d2->vel[nd2].x)/2.)/mod->gravity;
                dhdy = (mod->sw->d2->head[nd1] - mod->sw->d2->head[nd2])/(mod->grid->node[nd1].y - mod->grid->node[nd2].y);
            }
            error += fabs(dhdy - slope);
            printf("X-Coriolis time-step results :: dh/dy: %30.20e rhs: %30.20e :: Error: %30.20e :: dy: %20.10f :: coriolis angle: %20.4f \n",
                   dhdy, slope, dhdy-slope, mod->grid->node[nd1].y - mod->grid->node[nd2].y,mod->mat[0].sw->coriolis);
        }
    }
#ifdef _MESSG
    error_max = messg_dsum(error,mod->grid->smpi->ADH_COMM);
    nnodes_Xsection_total = messg_isum(nnodes_Xsection, mod->grid->smpi->ADH_COMM);
#else
    error_max=error;
    nnodes_Xsection_total = nnodes_Xsection;
#endif
    if(mod->grid->smpi->myid==0){
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"average transect error:  %30.20e\n",error_max/(nnodes_Xsection_total-1));
        if (error/10. > 1.e-6)
            fprintf(fp,"FAIL\n");
        else
            fprintf(fp, "PASS\n");
        
        fclose(fp);
    }
}

//****************************************************************//
//****************************************************************//

void write_testcase_error_ycoriolis(SMODEL *mod) {
    
    // this test case finds the water surace slope at an x-transect at y = 50,000
    // if the grid changes, so must this routines!!!
    double convert_to_rads = 3.141592653589793 / 180.;
    double rads = convert_to_rads * mod->mat[0].sw->coriolis;  // assume all materials have same coriolis
    double angular_speed = 2. * EARTH_ROTAT * sin(rads);
    
    int i;
    int nd1 = 0;
    int nd2 = 0;
    int nnodes_Xsection_total;
    double error_max = 0.;
    double error = 0., dhdx = 0., slope = 0.;
    printf("\n");
    if(nnodes_Xsection>1){
        for (i=0; i<nnodes_Xsection-1; i++) {
            nd1 = nodeIDs_on_slice[i+1];
            nd2 = nodeIDs_on_slice[i];
            if (mod->grid->ndim == 3) {
                slope = - (angular_speed * (mod->sw->d3->vel[nd1].y + mod->sw->d3->vel[nd2].y)/2.)/mod->gravity;
                dhdx = -(mod->sw->d3->displacement[nd1] - mod->sw->d3->displacement[nd2])/(mod->grid->node[nd1].x - mod->grid->node[nd2].x);
            } else {
                slope = - (angular_speed * (mod->sw->d2->vel[nd1].y + mod->sw->d2->vel[nd2].y)/2.)/mod->gravity;
                dhdx = -(mod->sw->d2->head[nd1] - mod->sw->d2->head[nd2])/(mod->grid->node[nd1].x - mod->grid->node[nd2].x);
            }
            error += fabs(dhdx - slope);
            printf("Y-Coriolis time-step results :: dh/dx: %30.20e rhs: %30.20e :: Error: %30.20e :: dx: %10.10f :: coriolis angle: %10.4f\n",
                   dhdx, slope, dhdx-slope, mod->grid->node[nd1].x - mod->grid->node[nd2].x,mod->mat[0].sw->coriolis);
        }
    }
#ifdef _MESSG
    error_max = messg_dsum(error,mod->grid->smpi->ADH_COMM);
    nnodes_Xsection_total = messg_isum(nnodes_Xsection, mod->grid->smpi->ADH_COMM);
#else
    error_max=error;
    nnodes_Xsection_total = nnodes_Xsection;
#endif
    if(mod->grid->smpi->myid==0){
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"average transect error:  %30.20e\n",error_max/(nnodes_Xsection_total-1));
        if (error/10. > 1.e-6)
            fprintf(fp,"FAIL\n");
        else
            fprintf(fp, "PASS\n");
        
        fclose(fp);
    }
    
}

//****************************************************************//
//****************************************************************//
// cjt :: coriolis X angled at 45 degrees
// assumes initial elevation is 0
// because this grid is angled, we want the a force resolved diagonal coriolis and diagonal derivative balance
// problem with resolving in orthogonal coordinates, x, y, for angled grid, convection cannot be neglected

void write_testcase_error_coriolis(SMODEL *mod) {
    
    printf("\n");
    
    double g = mod->gravity;
    double convert_to_rads = pi / 180.;
    double rads = convert_to_rads * mod->mat[0].sw->coriolis;  // assume all materials have same coriolis
    double angular_speed = 2. * EARTH_ROTAT * sin(rads);
    
    // if the grid changes, so must this routines!!! Test node must be centered.
    int nz = 6;
    int ny = 11;
    
    // find local node ID for testing node
     int i, nd = UNSET_INT, nd1_x = UNSET_INT, nd2_x = UNSET_INT, nd1_y = UNSET_INT, nd2_y = UNSET_INT;
#ifdef _MESSG
    for (i=0; i<mod->grid->my_nnodes; i++) {
        //if (mod->grid->node[i].resident_pe == mod->grid->smpi->myid && mod->grid->node[i].gid == test_case_flag.node1) {
        if (mod->grid->node[i].gid == test_case_flag.node1) {
            nd = i;
            break;
        }
    }
    if (nd == UNSET_INT) return; // test node not on this PE
    assert(mod->grid->nodeID_3d_to_2d_sur[nd] != UNSET_INT);
    for (i=0; i<mod->grid->nnodes; i++) {
        //if (mod->grid->node[i].resident_pe == mod->grid->smpi->myid && mod->grid->node[i].gid == test_case_flag.node1) {
        if (mod->grid->node[i].gid == test_case_flag.node1 + nz*ny) {
            nd1_x = i;
            break;
        }
    }
    for (i=0; i<mod->grid->nnodes; i++) {
        //if (mod->grid->node[i].resident_pe == mod->grid->smpi->myid && mod->grid->node[i].gid == test_case_flag.node1) {
        if (mod->grid->node[i].gid == test_case_flag.node1 - nz*ny) {
            nd2_x = i;
            break;
        }
    }
    for (i=0; i<mod->grid->nnodes; i++) {
        //if (mod->grid->node[i].resident_pe == mod->grid->smpi->myid && mod->grid->node[i].gid == test_case_flag.node1) {
        if (mod->grid->node[i].gid == test_case_flag.node1 - nz) {
            nd1_y = i;
            break;
        }
    }
    for (i=0; i<mod->grid->nnodes; i++) {
        //if (mod->grid->node[i].resident_pe == mod->grid->smpi->myid && mod->grid->node[i].gid == test_case_flag.node1) {
        if (mod->grid->node[i].gid == test_case_flag.node1 + nz) {
            nd2_y = i;
            break;
        }
    }
#else
    nd = test_case_flag.node1;
    nd1_x = nd + nz*ny; assert(mod->grid->nodeID_3d_to_2d_sur[nd1_x] != UNSET_INT); // + dx
    nd2_x = nd - nz*ny; assert(mod->grid->nodeID_3d_to_2d_sur[nd2_x] != UNSET_INT); // - dx
    nd1_y = nd - nz; assert(mod->grid->nodeID_3d_to_2d_sur[nd1_y] != UNSET_INT);    // + dy
    nd2_y = nd + nz; assert(mod->grid->nodeID_3d_to_2d_sur[nd2_y] != UNSET_INT);    // - dy
#endif
    
    
    // test must be surface and centered
    //if (mod->grid->node[test_case_flag.node1].resident_pe != mod->grid->smpi->myid) return;
    //if (MY_NODE_PE != mod->grid->smpi->myid) return;
    //printf("myid: %d nd: %d nnodes: %d mod->grid->nodeID_3d_to_2d_sur[nd]: %d \n",mod->grid->smpi->myid,nd,mod->grid->nnodes,mod->grid->nodeID_3d_to_2d_sur[nd]);
    //int nd = test_case_flag.node1;
//    snode_printScreen(mod->grid->node[nd]);
//    snode_printScreen(mod->grid->node[nd1_x]);
//    snode_printScreen(mod->grid->node[nd2_x]);
    
    //  sanity checks on grid
    assert(fabs(mod->grid->node[nd1_x].x - mod->grid->node[nd2_y].x)<1e-6 );
    assert(fabs(mod->grid->node[nd2_x].x - mod->grid->node[nd1_y].x)<1e-6 );
    
    
    double dist_x = sqrt(pow(mod->grid->node[nd1_x].x - mod->grid->node[nd2_x].x,2) + pow(mod->grid->node[nd1_x].y - mod->grid->node[nd2_x].y,2));
    double dist_y = sqrt(pow(mod->grid->node[nd1_y].x - mod->grid->node[nd2_y].x,2) + pow(mod->grid->node[nd1_y].y - mod->grid->node[nd2_y].y,2));
    
    //    printf("%10.5f \t %10.5f \t %10.5f\n",mod->grid->node[nd].x,mod->grid->node[nd].y,mod->grid->node[nd].z);
    //    printf("%10.5f \t %10.5f \t %10.5f\n",mod->grid->node[nd1_x].x,mod->grid->node[nd1_x].y,mod->grid->node[nd1_x].z);
    //    printf("%10.5f \t %10.5f \t %10.5f\n",mod->grid->node[nd2_x].x,mod->grid->node[nd2_x].y,mod->grid->node[nd2_x].z);
    //    printf("%10.5f \t %10.5f \t %10.5f\n",mod->grid->node[nd1_y].x,mod->grid->node[nd1_y].y,mod->grid->node[nd1_y].z);
    //    printf("%10.5f \t %10.5f \t %10.5f\n",mod->grid->node[nd2_y].x,mod->grid->node[nd2_y].y,mod->grid->node[nd2_y].z);
    //    printf("dist_x: %10.5e \t dist_y: %10.5e\n",dist_x,dist_y);
    
    double term_x, term_y, dhdx, dhdy, dhdx2, dhdy2, error, error_x, error_y;
    if (mod->grid->ndim == 3) {
        term_x = -angular_speed * mod->sw->d3->vel[nd].y;
        term_y =  angular_speed * mod->sw->d3->vel[nd].x;
        //dhdx = (mod->sw->d3->displacement[nd1_x] - mod->sw->d3->displacement[nd2_x])/dist_x;
        dhdy = (mod->sw->d3->displacement[nd2_y] - mod->sw->d3->displacement[nd1_y])/dist_y; // derivative along nd1,2_y string
        
        term_x = sqrt( term_x*term_x + term_y*term_y); // get force vector along nd1,2_y direction
        
        error_x = fabs(dhdy - term_x/g);
        error_y = 0.0;  //fabs(dhdy - term_y/g);
        error = error_x + error_y;
    }
    printf("dh: %10.5e    g*term_x: %10.5e    error_x: %10.5e\n",dhdy,term_x/g,error_x);
    //printf("dh/dy: %10.5e    g*term_y: %10.5e    error_y: %10.5e\n",dhdy,term_y/g,error_y);
    printf("total absolute error: %10.5e = %10.5f%%\n",error,100*error/(dhdy));
    
    //if(mod->grid->smpi->myid==0){
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"absolute error:  %10.5e\n",error);
        fprintf(fp,"absolute relative error:  %10.5e%%\n",100*error/(dhdx + dhdy));
        fclose(fp);
    //}
}

//********************************************************************************/
//********************************************************************************/
//****************************************************************//
//****************************************************************//

void write_testcase_error_ptest(SMODEL *mod) {
    
    // CJT :: just check surface displacements and velocities
    
    int i=0, nd=0, inode=0;
    ID_LIST_ITEM *ptr;
    double elevation = 0.;
    
    double error_dpl = 0., error_dpl_max = 0.;
    double error_velx = 0., error_velx_max = 0.;
    double error_vely = 0., error_vely_max = 0.;
    double error_velz = 0., error_velz_max = 0.;
    
    int count = 0, count_max = 0.;
    if (mod->flag.NS3_FLOW) {
        for (inode=0; inode<mod->grid->nnodes; inode++) {
            if (mod->grid->node[inode].resident_pe == mod->grid->smpi->myid && mod->grid->node[inode].bflag == 0) {
                count++;
                elevation = mod->grid->node[inode].z + mod->ns->d3->displacement[inode];
                error_dpl += fabs(0.05 - elevation);
                error_velx += fabs(mod->ns->d3->vel[inode].x);
                error_vely += fabs(mod->ns->d3->vel[inode].y);
                error_velz += fabs(mod->ns->d3->vel[inode].z);
            }
            
        }
    } else if (mod->flag.SW3_FLOW) {
        for (inode=0; inode<mod->grid->nnodes_sur; inode++) {
            ptr = mod->grid->vertical_list[inode];
            nd = ptr->id;
            if (mod->grid->node[nd].resident_pe == mod->grid->smpi->myid){
                count ++;
                elevation = mod->grid->node[nd].z + mod->sw->d3->displacement[nd];
                error_dpl += fabs(0.05 - elevation);
                error_velx += fabs(mod->sw->d3->vel[nd].x);
                error_vely += fabs(mod->sw->d3->vel[nd].y);
                error_velz += fabs(mod->sw->d3->vel[nd].z);
            }
        }
    }
    
#ifdef _MESSG
    count_max=messg_isum(count,mod->grid->smpi->ADH_COMM);
    error_velx_max=messg_dsum(error_velx,mod->grid->smpi->ADH_COMM);
    error_vely_max=messg_dsum(error_vely,mod->grid->smpi->ADH_COMM);
    error_vely_max=messg_dsum(error_velz,mod->grid->smpi->ADH_COMM);
    error_dpl_max=messg_dsum(error_dpl,mod->grid->smpi->ADH_COMM);
    count = count_max;
    error_dpl=error_dpl_max;
    error_velx=error_velx_max;
    error_vely=error_vely_max;
#endif
    
    if(mod->grid->smpi->myid==0){
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"displacement error:  %30.20e\n",error_dpl/count);
        fprintf(fp,"x-velocity error:  %30.20e\n",error_velx/count);
        fprintf(fp,"y-velocity error:  %30.20e\n",error_vely/count);
        fprintf(fp,"z-velocity error:  %30.20e\n",error_velz/count);
        if (error_dpl/count > 1.e-6){
            fprintf(fp,"Displacment FAIL\n");
        }
        if (error_velx/count > 1.e-6)
            fprintf(fp,"VelX FAIL\n");
        else
            fprintf(fp, "PASS\n");
        if (error_vely/count > 1.e-6)
            fprintf(fp,"VelY FAIL\n");
        else
            fprintf(fp, "PASS\n");
        if (error_velz/count > 1.e-6)
            fprintf(fp,"VelZ FAIL\n");
        else
            fprintf(fp, "PASS\n");
        fclose(fp);
    }
}

//****************************************************************//
//****************************************************************//

// the whole grid should be at vx=vy=0.003536
void write_testcase_error_discharge(SMODEL *mod) {
    int i = 0, inode = 0;
    double error_vx = 0., error_vy = 0., error_vx_max = 0., error_vy_max = 0.;
    double discharge = 0.353553;
    
    // include include downstream half of grid
    int node_start = 7855;
#ifdef _MESSG
    node_start = mod->grid->my_nnodes+1;
    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        if (mod->grid->node[inode].gid >= 7855){
            node_start = inode;
            break;
        }
    }
#endif
    for (inode=node_start; inode<mod->grid->my_nnodes; inode++) {
        error_vx += fabs(mod->sw->d3->vel[inode].x - discharge);
        error_vy += fabs(mod->sw->d3->vel[inode].y - discharge);
    }
#ifdef _MESSG
    error_vx_max=messg_dsum(error_vx,mod->grid->smpi->ADH_COMM);
    error_vy_max=messg_dsum(error_vy,mod->grid->smpi->ADH_COMM);
    error_vx=error_vx_max;
    error_vy=error_vy_max;
#endif
    
    
    double new_grid_mass = 0.;
    double mass_error = tl_find_3d_grid_mass_error(mod->str_values, mod->series_head, mod->initial_grid_mass, mod->density, mod->grid, mod->sw->d3->vel, mod->sw->d3->displacement,  mod->sw->d3->old_displacement, mod->sw->d3->older_displacement, mod->dt, &new_grid_mass, &total_time_mass_flux);
    
    if(mod->grid->smpi->myid == 0) {
        printf("\ninitial grid mass: %-8.4e(kg) \t grid mass: %-8.4e(kg) \t total_time_mass_flux: %-8.4e(kg) \t mass_error: %-8.4e(kg) \t relative mass eror: %-8.4e \n",mod->initial_grid_mass,new_grid_mass,total_time_mass_flux,mass_error,100*(mass_error/(mod->initial_grid_mass)));
        
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"x-velocity abs error: %30.20e :: relative error: %30.20e \n",error_vx/mod->grid->macro_nnodes,error_vx/mod->grid->macro_nnodes/discharge);
        fprintf(fp,"y-velocity abs error: %30.20e :: relative error: %30.20e \n",error_vy/mod->grid->macro_nnodes,error_vy/mod->grid->macro_nnodes/discharge);
        if (error_vx/mod->grid->nnodes > 1.e-6 || error_vx/mod->grid->macro_nnodes/discharge > 1.e-6)
            fprintf(fp,"VelX FAIL\n");
        else
            fprintf(fp, "PASS\n");
        
        if (error_vy/mod->grid->nnodes > 1.e-6 || error_vy/mod->grid->macro_nnodes/discharge > 1.e-6)
            fprintf(fp,"VelY FAIL\n");
        else
            fprintf(fp, "PASS\n");
        
        fclose(fp);
    }
    
}

//****************************************************************//
//****************************************************************//

// NOTE:  (cjt) w is linear along the z-axis, so one layer should recover this

void write_testcase_error_slosh(SMODEL *mod) {
    
    double a = test_case_flag.a;
    double L = test_case_flag.L;
    double B = test_case_flag.B;
    double H = test_case_flag.H;
    
#ifdef _MESSG
    int i;
    for (i=0;i<mod->grid->my_nnodes;i++){
        if(mod->grid->node[i].gid==test_case_flag.nodeID){
            ALREADY_RECORDED = 1;
            test_case_flag.nodeID = i;
            break;
        }
    }
#else
    int ALREADY_RECORDED=1;
#endif
    
    int nodeID = test_case_flag.nodeID;
    double x = mod->grid->node[nodeID].x;
    double y = mod->grid->node[nodeID].y;
    double z = mod->grid->node[nodeID].z;
    
    double time = mod->t_prev + mod->dt;
    
    double term1 = pi/L;
    double term2 = term1 * sqrt(mod->gravity * H);
    
    double analytic_dpl = a * cos(term1 * x) * cos(term2 * time);
    double analytic_u = a * sqrt(mod->gravity * H) / H * sin(term1 * x) * sin(term2 * time);
    double analytic_w = -a * term2 * ((z + H)/ H) * cos(term1 * x) * sin(term2 * time);
    double error_dpl,error_u,error_w;
    
    if(ALREADY_RECORDED > 0){
        error_dpl = fabs(mod->sw->d3->displacement[nodeID]  - analytic_dpl);
        error_u = fabs(mod->sw->d3->vel[nodeID].x - analytic_u);
        error_w = fabs(mod->sw->d3->vel[nodeID].z - analytic_w);
        
        
        printf("\nnode: %d :: (x,y,z): (%20.10f, %20.10f, %20.10f)\n",nodeID,x,y,z);
        printf("dpl :: analyic: %20.10f \t numeric: %20.10f \t error: %20.10f \n",analytic_dpl,  mod->sw->d3->displacement[nodeID], error_dpl);
        printf("u :: analyic: %20.10f \t numeric: %20.10f \t error: %20.10f \n",analytic_u, mod->sw->d3->vel[nodeID].x,error_u);
        printf("w :: analyic: %20.10f \t numeric: %20.10f \t error: %20.10f \n",analytic_w, mod->sw->d3->vel[nodeID].z,error_w);
        
        if (fabs(error_dpl) > error_dpl_tmax) error_dpl_tmax = error_dpl;
        if (fabs(error_u) > error_u_tmax) error_u_tmax = error_u;
        if (fabs(error_w) > error_w_tmax) error_w_tmax = error_w;
    }
    
    // calculate errors over entire grid (may include wacky boundary corners)
    int inode;
    double analytic_grid_dpl = 0., analytic_grid_u = 0, analytic_grid_w = 0.;
    double sum_error_dpl = 0., sum_error_u = 0., sum_error_w = 0.;
    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        x = mod->grid->node[inode].x;
        y = mod->grid->node[inode].y;
        z = mod->grid->node[inode].z;
        analytic_grid_dpl = a * cos(term1 * x) * cos(term2 * time);
        analytic_grid_u = a * sqrt(mod->gravity * H) / H * sin(term1 * x) * sin(term2 * time);
        analytic_grid_w = -a * term2 * ((z + H)/ H) * cos(term1 * x) * sin(term2 * time);
        sum_error_dpl += pow(mod->sw->d3->displacement[inode]  - analytic_grid_dpl,2);
        sum_error_u += pow(mod->sw->d3->vel[inode].x - analytic_grid_u,2);
        sum_error_w += pow(mod->sw->d3->vel[inode].z - analytic_grid_w,2);
    }
    
#ifdef _MESSG
    
    sum_error_dpl =  messg_dsum(sum_error_dpl, mod->grid->smpi->ADH_COMM);
    sum_error_u =  messg_dsum(sum_error_u, mod->grid->smpi->ADH_COMM);
    sum_error_w =  messg_dsum(sum_error_w, mod->grid->smpi->ADH_COMM);
#endif
    
    double rmse_dpl = sqrt(sum_error_dpl/(float)mod->grid->macro_nnodes);
    double rmse_u = sqrt(sum_error_u/(float)mod->grid->macro_nnodes);
    double rmse_w = sqrt(sum_error_w/(float)mod->grid->macro_nnodes);
    
    if (fabs(rmse_dpl) > rmse_error_dpl_tmax) rmse_error_dpl_tmax = rmse_dpl;
    if (fabs(rmse_u) > rmse_error_u_tmax) rmse_error_u_tmax = rmse_u;
    if (fabs(rmse_w) > rmse_error_w_tmax) rmse_error_w_tmax = rmse_w;
    
    FILE *fp;
    
    // plot the error
    //fp = fopen("plot.out","a+");
    //fprintf(fp,"%20.10f %20.10f\n",time,error_dpl);
    //fclose(fp);
    if (ALREADY_RECORDED >0){
        fp = fopen("error.out", "w");
        
        //fprintf(fp,"node %d max error last time-step displacment: %30.20e :: relative error: %30.20e \n",nodeID+1,error_dpl,error_dpl/analytic_dpl);
        //fprintf(fp,"node %d max error last time-step x-velocity: %30.20e :: relative error: %30.20e \n",nodeID+1,error_u,error_u/analytic_u);
        //fprintf(fp,"node %d max error last time-step z-velocity: %30.20e :: relative error: %30.20e \n",nodeID+1,error_w,error_w/analytic_w);
        //fprintf(fp,"\n");
        fprintf(fp,"node %d error max over time displacment: %30.20e \n",nodeID+1,error_dpl_tmax);
        fprintf(fp,"node %d error max over time x-velocity: %30.20e \n",nodeID+1,error_u_tmax);
        fprintf(fp,"node %d error max over time z-velocity: %30.20e \n",nodeID+1,error_w_tmax);
        fprintf(fp,"\n");
        fprintf(fp,"grid RMSE error max over time displacment: %30.20e \n",rmse_error_dpl_tmax);
        fprintf(fp,"grid RMSE error max over time x-velocity: %30.20e \n",rmse_error_u_tmax);
        fprintf(fp,"grid RMSE error max over time z-velocity: %30.20e \n",rmse_error_w_tmax);
        
        if (error_dpl > 1.e-6 || error_dpl/analytic_dpl > 1.e-6 || error_dpl_tmax > 1.e-6)
            fprintf(fp,"Displacement FAIL\n");
        else
            fprintf(fp, "PASS\n");
        
        if (error_u > 1.e-6 || error_u/analytic_u > 1.e-6 || error_u_tmax > 1.e-6)
            fprintf(fp,"XVel FAIL\n");
        else
            fprintf(fp, "PASS\n");
        
        if (error_w > 1.e-6 || error_w/analytic_w > 1.e-6 || error_w_tmax > 1.e-6)
            fprintf(fp,"ZVel FAIL \n");
        else
            fprintf(fp, "PASS\n");
        
        fclose(fp);
    }
    
    
}

//----------------------------------------------------------------//
//****************************************************************//

void testcase_prep_slosh(SMODEL *mod) {
    double a = test_case_flag.a;
    double L = test_case_flag.L;
    double B = test_case_flag.B;
    double H = test_case_flag.H;
    
    double time = mod->t_prev;
    double time_old = mod->t_prev - mod->dt;
    double time_older = mod->t_prev - 2*mod->dt;
    
    
    //printf("a: %20.10f L: %20.10f B: %20.10f H:%20.10f time: %20.10f time_old: %20.10f time_older: %20.10f\n",a,L,B,H,time,time_old,time_older);
    //exit(-1);
    
    if (time < 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> T0 < 0.\n");
    }
    if (time_old < 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> For slosh model error calculation, please set initial time to T0 + 2 * dt (for hotstart).\n");
    }
    if (time_older < 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> For slosh model error calculation, please set initial time to T0 + 2 * dt (for hotstart).\n");
    }
    
    double term1 = pi/L;
    double term2 = term1 * sqrt(mod->gravity * H);
    
    int inode;
    double x = 0., y=0., z=0.;
    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        x = mod->grid->node[inode].x;
        y = mod->grid->node[inode].y;
        z = mod->grid->node[inode].z;
        
        mod->sw->d3->displacement[inode] = a * cos(term1 * x) * cos(term2 * time);
        mod->sw->d3->old_displacement[inode] = a * cos(term1 * x) * cos(term2 * time_old);
        mod->sw->d3->older_displacement[inode] = a * cos(term1 * x) * cos(term2 * time_older);
        
        mod->sw->d3->vel[inode].x =  a * (sqrt(mod->gravity * H) / H) * sin(term1 * x) * sin(term2 * time);
        mod->sw->d3->old_vel[inode].x =  a * (sqrt(mod->gravity * H) / H) * sin(term1 * x) * sin(term2 * time_old);
        mod->sw->d3->older_vel[inode].x =  a * (sqrt(mod->gravity * H) / H) * sin(term1 * x) * sin(term2 * time_older);
        
        mod->sw->d3->vel[inode].y = 0.;
        mod->sw->d3->old_vel[inode].y = 0.;
        mod->sw->d3->older_vel[inode].y = 0.;
        
        mod->sw->d3->vel[inode].z = -a * term2 * ((z + H)/ H) * cos(term1 * x) * sin(term2 * time);
        mod->sw->d3->old_vel[inode].z = -a * term2 * ((z + H)/ H) * cos(term1 * x) * sin(term2 * time_old);
        mod->sw->d3->older_vel[inode].z = -a * term2 * ((z + H)/ H) * cos(term1 * x) * sin(term2 * time_older);
    }
    
    // now distribute displacement down the column
    tl_vertical_adapt(mod->grid, mod->sw->d3->displacement);
    tl_vertical_adapt(mod->grid, mod->sw->d3->old_displacement);
    tl_vertical_adapt(mod->grid, mod->sw->d3->older_displacement);
}

//****************************************************************//
//****************************************************************//

// NOTE:  (cjt) w is linear along the z-axis, so one layer should recover this

void write_testcase_error_slosh2d3d(SSUPER_MODEL *sm) {
    
    int i,j;
    SMODEL *mod;
    SGRID *grid;
    
    //    int myid,npes,ierr_code;
    //    for (i=0;i<1000000;i++) {}
    //    ierr_code = MPI_Comm_rank(cstorm_comm, &myid);
    //    ierr_code = MPI_Comm_size(cstorm_comm, &npes);
    //    printf("world pe %d of %d in write_testcase_error_slosh2d3d\n",myid,npes);
    //    fflush(stdout); messg_barrier(cstorm_comm); //tl_error("fornow\n");
    
    double a = test_case_flag.a;
    double L = test_case_flag.L;
    double B = test_case_flag.B;
    double H = test_case_flag.H;
    
    double x, y, z, time, term1, term2;
    double analytic_dpl=0., analytic_u=0., analytic_w=0.;
    double error_dpl=0., error_u=0., error_w=0.;
    
    int nodeID, inode;
    
    //messg_barrier(cstorm_comm); tag();
#ifdef _MESSG
    for (i=0; i<sm->nsubmodels; i++){
        mod = &(sm->submodel[i]);
        if (mod->proc_flag == 1) {
            grid = mod->grid;
            for (j=0;j<mod->grid->my_nnodes;j++){
                if(mod->grid->node[j].gid==test_case_flag.nodeID){
                    ALREADY_RECORDED = 1;
                    test_case_flag.nodeID = j;
                    break;
                }
            }
        }
    }
#else
    ALREADY_RECORDED=1;
#endif
    //messg_barrier(cstorm_comm); tag();
    
    // printf("world pe %d of %d :: test_case_flag.model_id_2d3d: %d\n",myid,npes,test_case_flag.model_id_2d3d);
    mod = &(sm->submodel[test_case_flag.model_id_2d3d]);
    if (mod->proc_flag == 1) {
        
        grid = mod->grid;
        //printf("world pe %d of %d  :: model grid ndim: %d proc_flag: %d\n",myid,npes,grid->ndim,mod->proc_flag);
        nodeID = test_case_flag.nodeID;
        x = grid->node[nodeID].x;
        y = grid->node[nodeID].y;
        z = grid->node[nodeID].z;
        
        time = mod->t_prev + mod->dt;
        
        term1 = pi/L;
        term2 = term1 * sqrt(mod->gravity * H);
        
        analytic_dpl = a * cos(term1 * x) * cos(term2 * time);
        analytic_u = a * sqrt(mod->gravity * H) / H * sin(term1 * x) * sin(term2 * time);
        if (mod->flag.SW3_FLOW){
            analytic_w = -a * term2 * ((z + H)/ H) * cos(term1 * x) * sin(term2 * time);
        }
        
        if(ALREADY_RECORDED > 0){
            if (mod->flag.SW2_FLOW){
                error_dpl = fabs((mod->sw->d2->head[nodeID] - H) - analytic_dpl);
                error_u = fabs(mod->sw->d2->vel[nodeID].x - analytic_u);
                printf("\nnode: %d :: (x,y): (%20.10f, %20.10f)\n",nodeID,x,y);
                printf("dpl ::\tanalytic: %20.10e \n\t numeric: %20.10e \n\t   error: %20.10e \n",analytic_dpl,  mod->sw->d2->head[nodeID]-H, error_dpl);
                printf("u   ::\tanalytic: %20.10e \n\t numeric: %20.10e \n\t   error: %20.10e \n",analytic_u, mod->sw->d2->vel[nodeID].x,error_u);
                
                if (fabs(error_dpl) > error_dpl_tmax) error_dpl_tmax = error_dpl;
                if (fabs(error_u) > error_u_tmax) error_u_tmax = error_u;
            }
            else if (mod->flag.SW3_FLOW){
                error_dpl = fabs(mod->sw->d3->displacement[nodeID]  - analytic_dpl);
                error_u = fabs(mod->sw->d3->vel[nodeID].x - analytic_u);
                error_w = fabs(mod->sw->d3->vel[nodeID].z - analytic_w);
                
                printf("\nnode: %d :: (x,y,z): (%20.10f, %20.10f, %20.10f)\n",nodeID,x,y,z);
                printf("dpl ::\tanalytic: %20.10e \n\t numeric: %20.10e \n\t   error: %20.10e \n",analytic_dpl,  mod->sw->d3->displacement[nodeID], error_dpl);
                printf("u   ::\tanalytic: %20.10e \n\t numeric: %20.10e \n\t   error: %20.10e \n",analytic_u, mod->sw->d3->vel[nodeID].x,error_u);
                printf("w   ::\tanalytic: %20.10e \n\t numeric: %20.10e \n\t   error: %20.10e \n",analytic_w, mod->sw->d3->vel[nodeID].z,error_w);
                
                if (fabs(error_dpl) > error_dpl_tmax) error_dpl_tmax = error_dpl;
                if (fabs(error_u) > error_u_tmax) error_u_tmax = error_u;
                if (fabs(error_w) > error_w_tmax) error_w_tmax = error_w;
            }
        }
    }
    
    char filename[30];
#ifdef _MESSG
#else
    FILE *depanalyticfp;
    FILE *velanalyticfp;
    FILE *deperrfp;
    FILE *velerrfp;
#endif
    
    // calculate errors over entire grid (may include wacky boundary corners)
    double analytic_grid_dpl = 0., analytic_grid_u = 0, analytic_grid_w = 0.;
    double sum_error_dpl_2d  = 0., sum_error_u_2d  = 0.;
    double sum_error_dpl_3d  = 0., sum_error_u_3d  = 0., sum_error_w_3d = 0.;
    int macro_nnodes_2d = 0, macro_nnodes_3d_depth = 0, macro_nnodes_3d_vel = 0;
    //tag();
    for (i=0; i<sm->nsubmodels; i++){
        mod = &(sm->submodel[i]);
        if(mod->proc_flag==1){
            grid = mod->grid;
            
#ifdef _MESSG
#else
            sprintf(filename,"%s_analytic_dep.dat",mod->io->proj_name);
            depanalyticfp = fopen(filename,"a");
            sprintf(filename,"%s_analytic_vel.dat",mod->io->proj_name);
            velanalyticfp = fopen(filename,"a");
            
            sprintf(filename,"%s_error_dep.dat",mod->io->proj_name);
            deperrfp = fopen(filename,"a");
            sprintf(filename,"%s_error_vel.dat",mod->io->proj_name);
            velerrfp = fopen(filename,"a");
            
            fprintf(depanalyticfp, "TS 0 %15.8e\n", time);
            fprintf(velanalyticfp, "TS 0 %15.8e\n", time);
            fprintf(deperrfp, "TS 0 %15.8e\n", time);
            fprintf(velerrfp, "TS 0 %15.8e\n", time);
#endif
            
            //tag();
            if (mod->flag.SW2_FLOW){
                //tag();
                for (inode=0; inode<grid->my_nnodes; inode++) {
                    x = grid->node[inode].x;
                    y = grid->node[inode].y;
                    z = grid->node[inode].z;
                    analytic_grid_dpl = a * cos(term1 * x) * cos(term2 * time);
                    analytic_grid_u = a * sqrt(mod->gravity * H) / H * sin(term1 * x) * sin(term2 * time);
                    sum_error_dpl_2d += pow((mod->sw->d2->head[inode]-H) - analytic_grid_dpl,2);
                    sum_error_u_2d += pow(mod->sw->d2->vel[inode].x - analytic_grid_u,2);
#ifdef _MESSG
#else
                    fprintf(depanalyticfp,"%16.8e\n",H+analytic_grid_dpl);
                    fprintf(velanalyticfp,"%16.8e %16.8e %16.8e\n",analytic_grid_u, 0.0, 0.0);
                    fprintf(deperrfp,"%16.8e\n",(mod->sw->d2->head[inode]-H) - analytic_grid_dpl);
                    fprintf(velerrfp,"%16.8e %16.8e %16.8e\n"
                            ,mod->sw->d2->vel[inode].x - analytic_grid_u, mod->sw->d2->vel[inode].y, 0.0);
#endif
                }
                macro_nnodes_2d+=grid->macro_nnodes;
                //tag();
            }
            else if (mod->flag.SW3_FLOW){
                //tag();
                /* Loop over the surface nodes */
                for (j=0; j<grid->my_nnodes_sur; j++) {
                    inode = grid->vertical_list[j]->id;
                    x = grid->node[inode].x;
                    y = grid->node[inode].y;
                    z = grid->node[inode].z;
                    analytic_grid_dpl = a * cos(term1 * x) * cos(term2 * time);
                    sum_error_dpl_3d += pow(mod->sw->d3->displacement[inode]  - analytic_grid_dpl,2);
#ifdef _MESSG
#else
                    fprintf(depanalyticfp,"%16.8e\n", H + analytic_grid_dpl);
                    fprintf(deperrfp,"%16.8e\n",mod->sw->d3->displacement[inode] - analytic_grid_dpl);
#endif
                }
                macro_nnodes_3d_depth +=grid->macro_nnodes_sur;
                for (inode=0;inode<grid->my_nnodes;inode++){
                    x = grid->node[inode].x;
                    y = grid->node[inode].y;
                    z = grid->node[inode].z;
                    analytic_grid_u = a * sqrt(mod->gravity * H) / H * sin(term1 * x) * sin(term2 * time);
                    analytic_grid_w = -a * term2 * ((z + H)/ H) * cos(term1 * x) * sin(term2 * time);
                    sum_error_u_3d += pow(mod->sw->d3->vel[inode].x - analytic_grid_u,2);
                    sum_error_w_3d += pow(mod->sw->d3->vel[inode].z - analytic_grid_w,2);
#ifdef _MESSG
#else
                    fprintf(velanalyticfp,"%15.8e %15.8e %15.8e\n"
                            ,analytic_grid_u
                            ,0.0
                            ,analytic_grid_w);
                    fprintf(velerrfp,"%15.8e %15.8e %15.8e\n"
                            ,mod->sw->d3->vel[inode].x - analytic_grid_u
                            ,mod->sw->d3->vel[inode].y
                            ,mod->sw->d3->vel[inode].z - analytic_grid_w);
#endif
                }
                macro_nnodes_3d_vel+=grid->macro_nnodes;
                //tag();
            }
            
#ifdef _MESSG
#else
            fclose(depanalyticfp);
            fclose(velanalyticfp);
            fclose(deperrfp);
            fclose(velerrfp);
#endif
            //tag();
            
#ifdef _MESSG
            if (mod->flag.SW2_FLOW){
                //tag();
                //if (mod->proc_flag == 1) {
                sum_error_dpl_2d =  messg_dsum(sum_error_dpl_2d, mod->grid->smpi->ADH_COMM);
                sum_error_u_2d =  messg_dsum(sum_error_u_2d, mod->grid->smpi->ADH_COMM);
                //}
                //tag();
            }
            else if (mod->flag.SW3_FLOW){
                //tag();
                //if (mod->proc_flag == 1) {
                sum_error_dpl_3d =  messg_dsum(sum_error_dpl_3d, mod->grid->smpi->ADH_COMM);
                sum_error_u_3d =  messg_dsum(sum_error_u_3d, mod->grid->smpi->ADH_COMM);
                sum_error_w_3d =  messg_dsum(sum_error_w_3d, mod->grid->smpi->ADH_COMM);
                //}
                //tag();
            }
#endif
        }
    }
    //tag();
    
    double rmse_dpl = sqrt(sum_error_dpl_2d/(double)macro_nnodes_2d);
    double rmse_u = sqrt(sum_error_u_2d/(double)macro_nnodes_2d);
    
    if (fabs(rmse_dpl) > rmse_error_dpl_tmax_2d) rmse_error_dpl_tmax_2d = rmse_dpl;
    if (fabs(rmse_u) > rmse_error_u_tmax_2d) rmse_error_u_tmax_2d = rmse_u;
    
    rmse_dpl = sqrt(sum_error_dpl_3d/(double)macro_nnodes_3d_depth);
    rmse_u = sqrt(sum_error_u_3d/(double)macro_nnodes_3d_vel);
    double rmse_w = sqrt(sum_error_w_3d/(double)macro_nnodes_3d_vel);
    
    if (fabs(rmse_dpl) > rmse_error_dpl_tmax_3d) rmse_error_dpl_tmax_3d = rmse_dpl;
    if (fabs(rmse_u) > rmse_error_u_tmax_3d) rmse_error_u_tmax_3d = rmse_u;
    if (fabs(rmse_w) > rmse_error_w_tmax_3d) rmse_error_w_tmax_3d = rmse_w;
    
    mod = &(sm->submodel[test_case_flag.model_id_2d3d]);
    FILE *fp;
    
    /*    //plot the error
     if (mod->flag.SW2_FLOW){
     sprintf(filename,"plot_error_dpl_node_%i_2d.out",nodeID+1);
     }
     else if (mod->flag.SW3_FLOW){
     sprintf(filename,"plot_error_dpl_node_%i_3d.out",nodeID+1);
     }
     fp = fopen(filename,"a");
     fprintf(fp,"%20.10e    %20.10e\n",time,error_dpl);
     fclose(fp);
     
     //plot the error
     if (mod->flag.SW2_FLOW){
     sprintf(filename,"plot_error_xvel_node_%i_2d.out",nodeID+1);
     }
     else if (mod->flag.SW3_FLOW){
     sprintf(filename,"plot_error_xvel_node_%i_3d.out",nodeID+1);
     }
     fp = fopen(filename,"a");
     fprintf(fp,"%20.10e    %20.10e\n",time,error_u);
     fclose(fp);
     */
    if (ALREADY_RECORDED >0){
        if (mod->flag.SW2_FLOW){
            sprintf(filename,"error_node_%i_2d.out",nodeID+1);
        }
        else if (mod->flag.SW3_FLOW){
            sprintf(filename,"error_node_%i_3d.out",nodeID+1);
        }
        fp = fopen(filename, "w");
        
        fprintf(fp,"%s node %d max error last time-step displacment: %20.10e :: relative error: %20.10e \n",mod->io->proj_name,nodeID+1,error_dpl,error_dpl/analytic_dpl);
        fprintf(fp,"%s node %d max error last time-step x-velocity:  %20.10e :: relative error: %20.10e \n",mod->io->proj_name,nodeID+1,error_u,error_u/analytic_u);
        if (mod->flag.SW3_FLOW){
            fprintf(fp,"%s node %d max error last time-step z-velocity:  %20.10e :: relative error: %20.10e \n",mod->io->proj_name,nodeID+1,error_w,error_w/analytic_w);
        }
        fprintf(fp,"\n");
        fprintf(fp,"%s node %d error max over time displacment: %20.10e \n",mod->io->proj_name,nodeID+1,error_dpl_tmax);
        fprintf(fp,"%s node %d error max over time x-velocity:  %20.10e \n",mod->io->proj_name,nodeID+1,error_u_tmax);
        if (mod->flag.SW3_FLOW){
            fprintf(fp,"%s node %d error max over time z-velocity:  %20.10e \n",mod->io->proj_name,nodeID+1,error_w_tmax);
        }
        fprintf(fp,"\n");
        fprintf(fp,"2D submodels: grid RMSE error max over time displacment: %20.10e \n",rmse_error_dpl_tmax_2d);
        fprintf(fp,"2D submodels: grid RMSE error max over time x-velocity:  %20.10e \n",rmse_error_u_tmax_2d);
        fprintf(fp,"\n");
        fprintf(fp,"3D submodels: grid RMSE error max over time displacment: %20.10e \n",rmse_error_dpl_tmax_3d);
        fprintf(fp,"3D submodels: grid RMSE error max over time x-velocity:  %20.10e \n",rmse_error_u_tmax_3d);
        fprintf(fp,"3D submodels: grid RMSE error max over time z-velocity:  %20.10e \n",rmse_error_w_tmax_3d);
        fprintf(fp,"\n");
        
        if (error_dpl > 1.e-6 || error_dpl/analytic_dpl > 1.e-6 || error_dpl_tmax > 1.e-6)
            fprintf(fp,"Displacement FAIL\n");
        else
            fprintf(fp, "Displacement PASS\n");
        
        if (error_u > 1.e-6 || error_u/analytic_u > 1.e-6 || error_u_tmax > 1.e-6)
            fprintf(fp,"XVel FAIL\n");
        else
            fprintf(fp, "XVel PASS\n");
        
        if (mod->flag.SW3_FLOW){
            if (error_w > 1.e-6 || error_w/analytic_w > 1.e-6 || error_w_tmax > 1.e-6)
                fprintf(fp,"ZVel FAIL \n");
            else
                fprintf(fp, "ZVel PASS\n");
        }
        
        fclose(fp);
    }
}

//----------------------------------------------------------------//
//****************************************************************//

void testcase_prep_slosh2d3d(SSUPER_MODEL *sm) {
    
    int i;
    // int mod_id = test_case_flag.model_id_2d3d;
    SMODEL *mod;
    SGRID *grid;
    
    double a = test_case_flag.a;
    double L = test_case_flag.L;
    double B = test_case_flag.B;
    double H = test_case_flag.H;
    
#ifdef _MESSG
#else
    FILE *depanalyticfp;
    FILE *velanalyticfp;
    FILE *deperrfp;
    FILE *velerrfp;
#endif
    
    for (i=0; i<sm->nsubmodels; i++){
        mod = &(sm->submodel[i]);
        if (mod->proc_flag == 1) {
            grid = mod->grid;
            double time = mod->t_prev;
            double time_old = mod->t_prev - mod->dt;
            double time_older = mod->t_prev - 2*mod->dt;
            
            if (time < 0) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> T0 < 0.\n");
            }
            if (time_old < 0) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> For slosh model error calculation, please set initial time to T0 + 2 * dt (for hotstart).\n");
            }
            if (time_older < 0) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> For slosh model error calculation, please set initial time to T0 + 2 * dt (for hotstart).\n");
            }
            
#ifdef _MESSG
#else
            char meshdim[15],filename[50];
            
            if (mod->flag.SW2_FLOW) sprintf(meshdim,"mesh2d");
            else if (mod->flag.SW3_FLOW) sprintf(meshdim,"mesh3d");
            sprintf(filename,"%s_analytic_dep.dat",mod->io->proj_name);
            depanalyticfp = fopen(filename,"w");
            sprintf(filename,"%s_analytic_vel.dat",mod->io->proj_name);
            velanalyticfp = fopen(filename,"w");
            
            sprintf(filename,"%s_error_dep.dat",mod->io->proj_name);
            deperrfp = fopen(filename,"w");
            sprintf(filename,"%s_error_vel.dat",mod->io->proj_name);
            velerrfp = fopen(filename,"w");
            
            
            /* Headers */
            fprintf(depanalyticfp, "DATASET\n");
            fprintf(velanalyticfp, "DATASET\n");
            fprintf(deperrfp,      "DATASET\n");
            fprintf(velerrfp,      "DATASET\n");
            pdata(mod->io->proj_name, "", depanalyticfp, 0, "Analytic Depth",    "", "BEGSCL", "mesh2d", grid->macro_nnodes_bed, grid->initial_nelems_bed);
            pdata(mod->io->proj_name, "", velanalyticfp, 0, "Analytic Velocity", "", "BEGVEC", meshdim, grid->initial_nnodes, grid->initial_nelems);
            pdata(mod->io->proj_name, "", deperrfp,      0, "Depth error",       "", "BEGSCL", "mesh2d", grid->macro_nnodes_bed, grid->initial_nelems_bed);
            pdata(mod->io->proj_name, "", velerrfp,      0, "Velocity error",    "", "BEGVEC", meshdim, grid->initial_nnodes, grid->initial_nelems);
            
            tc_timeunits(depanalyticfp,1.0);
            tc_timeunits(velanalyticfp,1.0);
            tc_timeunits(deperrfp,     1.0);
            tc_timeunits(velerrfp,     1.0);
            fprintf(depanalyticfp, "TS 0 %15.8e\n", time);
            fprintf(velanalyticfp, "TS 0 %15.8e\n", time);
            fprintf(deperrfp, "TS 0 %15.8e\n", time);
            fprintf(velerrfp, "TS 0 %15.8e\n", time);
#endif
            
            double term1 = pi/L;
            double term2 = term1 * sqrt(mod->gravity * H);
            
            int inode;
            double x = 0., y=0., z=0.;
            if (mod->flag.SW2_FLOW){
#ifdef _MESSG
#else
                sprintf(meshdim,"mesh2d");
#endif
                for (inode=0; inode<grid->nnodes; inode++) {
                    x = grid->node[inode].x;
                    y = grid->node[inode].y;
                    z = grid->node[inode].z;
                    
                    mod->sw->d2->head[inode] = H + a * cos(term1 * x) * cos(term2 * time);
                    mod->sw->d2->old_head[inode] = H + a * cos(term1 * x) * cos(term2 * time_old);
                    mod->sw->d2->older_head[inode] = H + a * cos(term1 * x) * cos(term2 * time_older);
                    
                    mod->sw->d2->vel[inode].x =  a * (sqrt(mod->gravity * H) / H) * sin(term1 * x) * sin(term2 * time);
                    mod->sw->d2->old_vel[inode].x =  a * (sqrt(mod->gravity * H) / H) * sin(term1 * x) * sin(term2 * time_old);
                    mod->sw->d2->older_vel[inode].x =  a * (sqrt(mod->gravity * H) / H) * sin(term1 * x) * sin(term2 * time_older);
                    
                    mod->sw->d2->vel[inode].y = 0.;
                    mod->sw->d2->old_vel[inode].y = 0.;
                    mod->sw->d2->older_vel[inode].y = 0.;
                    
#ifdef _MESSG
#else
                    fprintf(depanalyticfp,"%16.8e\n",mod->sw->d2->head[inode]);
                    fprintf(velanalyticfp,"%16.8e %16.8e %16.8e\n",mod->sw->d2->vel[inode].x, 0.0, 0.0);
                    fprintf(deperrfp,"%16.8e\n",0.0);
                    fprintf(velerrfp,"%16.8e %16.8e %16.8e\n",0.0,0.0,0.0);
#endif
                }
            }
            else if (mod->flag.SW3_FLOW){
#ifdef _MESSG
#else
                sprintf(meshdim,"mesh3d");
#endif
                for (inode=0; inode<grid->nnodes; inode++) {
                    x = grid->node[inode].x;
                    y = grid->node[inode].y;
                    z = grid->node[inode].z;
                    
                    mod->sw->d3->displacement[inode] = a * cos(term1 * x) * cos(term2 * time);
                    mod->sw->d3->old_displacement[inode] = a * cos(term1 * x) * cos(term2 * time_old);
                    mod->sw->d3->older_displacement[inode] = a * cos(term1 * x) * cos(term2 * time_older);
                    
                    mod->sw->d3->vel[inode].x =  a * (sqrt(mod->gravity * H) / H) * sin(term1 * x) * sin(term2 * time);
                    mod->sw->d3->old_vel[inode].x =  a * (sqrt(mod->gravity * H) / H) * sin(term1 * x) * sin(term2 * time_old);
                    mod->sw->d3->older_vel[inode].x =  a * (sqrt(mod->gravity * H) / H) * sin(term1 * x) * sin(term2 * time_older);
                    
                    mod->sw->d3->vel[inode].y = 0.;
                    mod->sw->d3->old_vel[inode].y = 0.;
                    mod->sw->d3->older_vel[inode].y = 0.;
                    
                    mod->sw->d3->vel[inode].z = -a * term2 * ((z + H)/ H) * cos(term1 * x) * sin(term2 * time);
                    mod->sw->d3->old_vel[inode].z = -a * term2 * ((z + H)/ H) * cos(term1 * x) * sin(term2 * time_old);
                    mod->sw->d3->older_vel[inode].z = -a * term2 * ((z + H)/ H) * cos(term1 * x) * sin(term2 * time_older);
#ifdef _MESSG
#else
                    fprintf(velanalyticfp,"%15.8e %15.8e %15.8e\n"
                            ,mod->sw->d3->vel[inode].x
                            ,0.0
                            ,mod->sw->d3->vel[inode].z);
                    fprintf(velerrfp,"%15.8e %15.8e %15.8e\n",0.0,0.0,0.0);
#endif
                }
                
#ifdef _MESSG
#else
                int j;
                for (j=0; j<grid->my_nnodes_sur; j++) {
                    inode = grid->vertical_list[j]->id;
                    fprintf(depanalyticfp,"%16.8e\n", H + mod->sw->d3->displacement[inode]);
                    fprintf(deperrfp,"%16.8e\n",0.0);
                }
#endif
                
                // now distribute displacement down the column
                tl_vertical_adapt(mod->grid, mod->sw->d3->displacement);
                tl_vertical_adapt(mod->grid, mod->sw->d3->old_displacement);
                tl_vertical_adapt(mod->grid, mod->sw->d3->older_displacement);
            }
#ifdef _MESSG
#else
            fclose(depanalyticfp);
            fclose(velanalyticfp);
            fclose(deperrfp);
            fclose(velerrfp);
#endif
        }
    }
}

//****************************************************************//
//****************************************************************//

void write_testcase_error_salt(SMODEL *mod) {
    int nnodes = 0;
    double error = 0., error_conc = 0., error_max = 0;
    
    nnodes = mod->grid->my_nnodes;
    
    int inode;
    for (inode=0; inode<nnodes; inode++) {
        error = fabs(17.5 - mod->con[0].concentration[inode]);
        if (error > error_max) error_max = error;
        if (error > error_c_tmax) error_c_tmax = error;
        error_conc += error;
    }
#ifdef _MESSG
    error_max = messg_dmax(error_max, mod->grid->smpi->ADH_COMM);
    if(mod->grid->smpi->myid==0){
#endif
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"angle_salt abs concentration error: %30.20e :: relative error: %30.20e \n",error_max,error_max/17.5);
        if (error_max > 1.e-6 || error_max/17.5 > 1.e-6) {
            fprintf(fp,"FAIL\n");
        } else {
            fprintf(fp, "PASS\n");
        }
        
        fclose(fp);
#ifdef _MESSG
    }
#endif
}

//****************************************************************//
//****************************************************************//
// lock exchange test case
void testcase_prep_lock(SMODEL *mod) {
    int inode;
    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        if (mod->grid->node[inode].x <= 1.0) {
            mod->con[0].concentration[inode] = 35.;
            mod->con[0].old_concentration[inode] = 35.;
            mod->con[0].older_concentration[inode] = 35.;
        } else {
            mod->con[0].concentration[inode] = 0.;
            mod->con[0].old_concentration[inode] = 0.;
            mod->con[0].older_concentration[inode] = 0.;
        }
    }
}

void write_testcase_error_lock(SMODEL *mod) {
    int nnodes = 0;
    double froude_number = 0.5;
    double analytic_wedge_speed = froude_number * sqrt(mod->gravity * 0.2 * (1-0.997));
    double time = mod->t_prev + mod->dt;
    double numeric_wedge_speed = 1./time; // half the grid is lx = 1
    nnodes = mod->grid->my_nnodes;
    
    int inode;
    for (inode=0; inode<nnodes; inode++) {
        // the concentration has risen at all on the edge of the grid ...
        if (fabs(mod->grid->node[inode].x - 2.0) < 1e-6 && mod->con[0].concentration[inode] > 1e-1) {
            FILE *fp;
            fp = fopen("error.out", "w");
            fprintf(fp,"time: %20.10e concentration: %20.10e analytic_wedge speed: %20.10e numeric_wedge_speed: %20.10e  error: %30.20e \n",time, mod->con[0].concentration[inode], analytic_wedge_speed, numeric_wedge_speed, fabs(numeric_wedge_speed - analytic_wedge_speed));
            fclose(fp);
            ALREADY_RECORDED = YES;
            exit(-1);  // exit run
        }
    }
    
}


/******************************************************************/
// Error computation for the steep front test case

void write_testcase_error_steepfront(SMODEL *mod) {
    
    int nodeID = test_case_flag.nodeID;
    double theta = -45. * (PI / 180.);
    
    int inode;
#ifdef _MESSG
    int i, imit = 0;
    for (i=0;i<mod->grid->my_nnodes;i++){
        if(mod->grid->node[i].gid==nodeID){
            imit = 1;
            nodeID = i;
            break;
        }
    }
#else
    int imit=1;
#endif
    double x,y,z,xt,yt;
    double sum_error_conc = 0.;
    double rmse_conc = 0.;
    double conc_init = test_case_flag.init_conc;
    double u_vel = sqrt(debug.u_vel*debug.u_vel + debug.v_vel*debug.v_vel); // the flow along the flume direction
    
    // angled coordinates
    xt = mod->grid->node[nodeID].x;
    yt = mod->grid->node[nodeID].y;
    z = mod->grid->node[nodeID].z;
    
    // orthogonal coordinates
    x = xt*cos(theta) - yt*sin(theta);
    y = xt*sin(theta) + yt*cos(theta);
    
    double time = mod->t_prev + mod->dt;
    double analytic_conc;
    double model_conc = mod->con[0].concentration[nodeID];
    double error;
    double erfc1, erfc2, erfc3;
    double diffusion = mod->mat[0].trn[0].d_m;
    FILE *fp;
    erfc1 = 0.;
    erfc2 = 0.;
    erfc3 = 0.;
    if(imit == 1) {
        erfc1 = erfc((x-u_vel*time)/pow(4.*diffusion*time, 0.5));
        erfc2 = erfc((x+u_vel*time)/pow(4.*diffusion*time, 0.5));
        if (fabs(erfc2) > 1e-6)
            erfc3 = erfc2 * exp(x*u_vel/diffusion);
        else
            erfc3 = 0.;
        
        analytic_conc = (conc_init/2.)*(erfc1 +  erfc3);
        error = analytic_conc - model_conc;
    }
    // calculate errors over entire grid (may include wacky boundary corners)
    double Ganalytic_conc = 0.;
    double Gmodel_conc = 0.;
    
    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        
        // angled coordinates
        xt = mod->grid->node[inode].x;
        yt = mod->grid->node[inode].y;
        z = mod->grid->node[inode].z;
        
        // orthogonal coordinates
        x = xt*cos(theta) - yt*sin(theta);
        y = xt*sin(theta) + yt*cos(theta);
        
        Gmodel_conc = mod->con[0].concentration[inode];
        erfc1 = erfc((x-u_vel*time)/pow(4.*diffusion*time, 0.5));
        erfc2 = erfc((x+u_vel*time)/pow(4.*diffusion*time, 0.5));
        if (fabs(erfc2) > 1e-6)
            erfc3 = erfc2 * exp (x*u_vel/diffusion);
        else
            erfc3 = 0.;
        Ganalytic_conc = (conc_init/2.)*(erfc1 +  erfc3);
        sum_error_conc += pow(Gmodel_conc  - Ganalytic_conc,2);
        
    }
    rmse_conc = sqrt(sum_error_conc/(double)mod->grid->my_nnodes);
    
#ifdef _MESSG
    double rmse_conc_tot = 0.;
    rmse_conc_tot = messg_dsum(rmse_conc, mod->grid->smpi->ADH_COMM);
    rmse_conc = rmse_conc_tot;
#endif
    if(mod->grid->smpi->myid==0){
        if (mod->t_prev == mod->t_init){
            fp = fopen("error.out", "w");
            fprintf(fp,"Node \t Time \t Analytic Concentration \t Model Computed Concentration \t Node Error \t GRID RMS ERROR\n");
        } else {
            fp = fopen("error.out", "a");
        }
        printf("\nnode: %d\t time: %20.10f\t analytic_c: %20.10f\t model_c: %20.10f\t error: %20.10f\t grid_rms_error: %20.10f\n",nodeID+1, time, analytic_conc, model_conc, error,rmse_conc);
        fprintf(fp," %d\t  %20.10f \t  %20.10f \t %20.10f \t  %20.10f \t %20.10f\n",nodeID+1, time, analytic_conc, model_conc, error,rmse_conc);
        fflush(fp);
        fclose(fp);
    }
    
}

//****************************************************************//
//****************************************************************//
// partial slip wind test case
// error is averaged over a verticle line dropped in the center of the grid
void write_testcase_error_winds_Huang(SMODEL *mod) {
    
    double k = 0; // linearalized bottom friction coeff -- set to 0 since AdH does not have linear bottom friction
    double tau = 0.1;  // surface stress mag
    double h = 40.;  // water depth
    
    if (mod->flag.WIND) {
        tau = sqrt( pow(mod->sw->d3->winds[0].stress.x,2) + pow(mod->sw->d3->winds[0].stress.y,2));
        //tau /= 1000.;
    } else if (mod->flag.WAVE) {
        tau = sqrt( pow(1000*mod->sw->d3->waves[0].stress.x,2) + pow(1000*mod->sw->d3->waves[0].stress.y,2));
    }
    //printf("tau: %20.10f \n",tau); //exit(-1);
    
    int inode, nodeID;
    double z, delta, eta_slope, analytic_u, vel_mag, error_u;
    double viscosity = mod->mat[0].sw->ev.xz + mod->viscosity;
    
    if (fabs(mod->mat[0].sw->ev.xz - mod->mat[0].sw->ev.yz) > 1e-6) {
        printf("ev.xz: %20.10f ev.yz: %20.10f\n",mod->mat[0].sw->ev.xz, mod->mat[0].sw->ev.yz);
        tl_error("ev.xz != ev.yz :: please fix wind/wave test case bc");
    }
    
    int count = 0;
    error_u = 0.;
    printf("\n");
    
#ifdef _MESSG
    int imit = 0;
    int i;
    for (i=0;i<mod->grid->nnodes_sur;i++){
        if((mod->grid->node[mod->grid->nodeID_2d_to_3d_sur[i]].gid==test_case_flag.nodeID) &&
           (mod->grid->node[mod->grid->nodeID_2d_to_3d_sur[i]].resident_pe == mod->grid->smpi->myid)){
            imit = 1;
            inode = i;
            break;
        }
    }
#else
    inode=test_case_flag.nodeID;
    int imit=1;
#endif
    
    if(imit > 0){
        
        ID_LIST_ITEM *ptr = mod->grid->vertical_list[inode];
        double surf_dpl = mod->sw->d3->displacement[ptr->id];
        while(ptr->next != NULL) {
            count ++;
            nodeID = ptr->id;
            z = mod->grid->node[nodeID].z;
            delta = (z+h)/h;
            // the actual solutions are on grid + dpl (matters on ALE grid most)
            delta = (z + mod->sw->d3->displacement[nodeID] + h)/(h + surf_dpl);
            eta_slope = (3./2.) * (tau/(mod->density * mod->gravity * h)) * ((2 * viscosity) + k * h)/(((3 * viscosity) + k * h));
            analytic_u = (1./(6 * viscosity)) * mod->gravity * eta_slope * pow(h,2) * (3 * pow(delta - 1,2) - 1) + tau * h * (2*delta - 1) / (2 * mod->density * viscosity);
            vel_mag = (mod->sw->d3->vel[nodeID].x/fabs(mod->sw->d3->vel[nodeID].x)) * sqrt( pow(mod->sw->d3->vel[nodeID].x,2) + pow(mod->sw->d3->vel[nodeID].y,2) );
            error_u += fabs(vel_mag - analytic_u);
            printf("node: %d [ %10.2f, %10.2f, %10.2f ] :: analytic: %20.10f  numeric: %20.10f  error: %20.10f\n",nodeID,mod->grid->node[nodeID].x,mod->grid->node[nodeID].y, z+ mod->sw->d3->displacement[nodeID], analytic_u, vel_mag,fabs(vel_mag - analytic_u));
            ptr = ptr->next;
        }
        error_u /= ((double) count);
        printf("average column u error: %20.10f\n",error_u);
        
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"x-velocity abs error last time: %30.20e \n",error_u);
        if (error_u > 1.e-6) {
            fprintf(fp,"XVel FAIL \n");
        } else {
            fprintf(fp, "PASS\n");
        }
        fclose(fp);
    }
    
    // check mass
    double new_grid_mass = 0.;
    double mass_error = tl_find_3d_grid_mass_error(mod->str_values, mod->series_head, mod->initial_grid_mass, mod->density, mod->grid, mod->sw->d3->vel, mod->sw->d3->displacement,  mod->sw->d3->old_displacement,  mod->sw->d3->older_displacement, mod->dt, &new_grid_mass, &total_time_mass_flux);
    if (mod->grid->smpi->myid == 0){
        printf("\ninitial grid mass: %-12.10e(kg) \t grid mass: %-12.10e(kg) total_time_mass_flux: %-12.10e(kg) \t mass_error: %-12.10e(kg) \t relative mass eror: %-12.10e%%\n",
               mod->initial_grid_mass,
               new_grid_mass,
               total_time_mass_flux,
               mass_error,
               100*(mass_error/(mod->initial_grid_mass)));
    }
    
}



//****************************************************************//
//****************************************************************//
// Tests water sourcing (rain)
// Depth = Depth(t=0) + ((1/4*Omega1*RI_omega1 + 1/4*Omega2*RI_omega2) / Total_Grid_Area) * t
void write_testcase_error_water_source(SMODEL *mod) {
    int i = 0, inode = 0;
    
    double omega1 = 2500;   // area of material 2 - 1/4 total grid area
    double omega2 = 2500;   // area of material 3 - 1/4 total grid area
    double RI1 = 0.001;     // rain intensity of material 2
    double RI2 = 0.003;     // rain intensity of material 3
    double total_area = 10000;  // total grid area
    double head0 = 100.;  // initial water depth (t=0)
    double time = mod->t_prev + mod->dt;
    
    double grid_mass = 0.;
    if (mod->flag.SW2_FLOW) {
        grid_mass = tl_find_grid_mass_elem2d(mod->density, mod->str_values, mod->series_head, mod->sw->d2->head, mod->grid, mod->flag);
    } else if (mod->flag.SW3_FLOW) {
        grid_mass = tl_find_grid_mass_elem3d(mod->density, mod->grid, mod->sw->d3->displacement);
    }
    
//#ifdef _MESSG
//    grid_mass = messg_dsum(grid_mass, mod->grid->smpi->ADH_COMM);
//#endif
    
    double analytic_depth = head0 + ((omega1 * RI1 + omega2 * RI2) / total_area) * time;
    double analytic_mass = analytic_depth * total_area * mod->density;
    double error_mass = fabs(grid_mass - analytic_mass);
    
    if(mod->grid->smpi->myid==0){
        printf("\ngrid_mass: %20.10e \t analytic_grid_mass: %20.10e \t absolute grid mass error: %20.10e\n",grid_mass, analytic_mass, error_mass);
        
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"mass abs error: %30.20e :: mass relative error: %30.20e \n",error_mass, error_mass/grid_mass);
        fclose(fp);
    }
}

//********************************************************************************/
//********************************************************************************/
// 2D WIND & WIND WAVE VERIFICTION
//****************************************************************//
//****************************************************************//
void testcase_prep_winds_and_waves(SMODEL *mod) {
    
    int i, id;
    
    // find nodes on user-defined cross-section
    if (mod->grid->ndim == 3) {
        nnodes_Xsection = tl_find_surface_nodes_on_sliceX(mod->grid, test_case_flag.node1, test_case_flag.node2, &nodeIDs_on_slice);
    } else {
        nnodes_Xsection = tl_find_nodes_on_sliceX(mod->grid, test_case_flag.node1, test_case_flag.node2, &nodeIDs_on_slice);
    }
    
    printf("\n** analytic results using the following X-section nodes: **\n");
    for(i=0; i<nnodes_Xsection; i++) {
        id = nodeIDs_on_slice[i];
        printf("nodeID: %d {x,y,z} = {%20.10e,%20.10e,%20.10e} \n",id+1,mod->grid->node[id].x,mod->grid->node[id].y,mod->grid->node[id].z);
    }
}
//****************************************************************//
//****************************************************************//
void write_testcase_error_winds_and_waves(SMODEL *mod) {
    
    int i;
    int nd1 = 0;
    int nd2 = 0;
    double error = 0., dhdx = 0., slope = 0., stress_mag1=0., stress_mag2=0., dist = 0.;
    double error_max;
    int nnodes_Xsection_max;
    printf("\n");
    if(nnodes_Xsection>1){
        for (i=0; i<nnodes_Xsection-1; i++) {
            nd1 = nodeIDs_on_slice[i+1];
            nd2 = nodeIDs_on_slice[i];
            if (mod->flag.WIND) {
                slope = sqrt( pow(mod->sw->d2->winds[nd1].stress.x,2) + pow(mod->sw->d2->winds[nd1].stress.y,2));
            } else if (mod->flag.WAVE) {
                slope = sqrt( pow(mod->sw->d2->waves[nd1].stress.x,2) + pow(mod->sw->d2->waves[nd1].stress.y,2));
                slope *= mod->density;
                //stress_mag1 = sqrt( pow(mod->sw->d2->waves[nd1].stress.x,2) + pow(mod->sw->d2->waves[nd1].stress.y,2));
                //stress_mag2 = sqrt( pow(mod->sw->d2->waves[nd2].stress.x,2) + pow(mod->sw->d2->waves[nd2].stress.y,2));
                //dist = sqrt(pow(mod->grid->node[nd1].x - mod->grid->node[nd2].x,2) + pow(mod->grid->node[nd1].y - mod->grid->node[nd2].y,2));
                //slope = (stress_mag1 - stress_mag2)/dist;
            }
            dist = sqrt(pow(mod->grid->node[nd1].x - mod->grid->node[nd2].x,2) + pow(mod->grid->node[nd1].y - mod->grid->node[nd2].y,2));
            dhdx = (mod->sw->d2->head[nd1] - mod->sw->d2->head[nd2])/dist;
            error += fabs(mod->gravity * 0.5*(mod->sw->d2->head[nd1]+mod->sw->d2->head[nd2]) * dhdx - slope/mod->density);
            
            printf("Surface Stress Time-Step Results :: dh/dx: %30.20e stress term: %30.20e :: Error: %30.20e\n",
                   mod->gravity * 0.5*(mod->sw->d2->head[nd1]+mod->sw->d2->head[nd2]) * dhdx, slope/mod->density,
                   fabs(mod->gravity * 0.5*(mod->sw->d2->head[nd1]+mod->sw->d2->head[nd2]) * dhdx - slope/mod->density));
        }
    }
#ifdef _MESSG
    error_max = messg_dsum(error,mod->grid->smpi->ADH_COMM);
    nnodes_Xsection_max = messg_isum(nnodes_Xsection, mod->grid->smpi->ADH_COMM);
#else
    error_max = error;
    nnodes_Xsection_max = nnodes_Xsection;
#endif
    if(mod->grid->smpi->myid==0){
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"average transect error:  %30.20e nnodes_Xsection_max %d \n",error_max/(nnodes_Xsection_max-1), nnodes_Xsection_max);
        if (error_max/10. > 1.e-6)
            fprintf(fp,"FAIL\n");
        else
            fprintf(fp, "PASS\n");
        
        fclose(fp);
    }
    
}

//********************************************************************************/
//********************************************************************************/
// 2D DAM BREAK ANALYTIC TEST CASE
//****************************************************************//
//****************************************************************//
void testcase_prep_dam2d(SMODEL *mod) {
    double x = 0.;
    int inode;
    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        x = mod->grid->node[inode].x;
        if (x <= test_case_flag.dam_location ) {
            mod->sw->d2->head[inode] = test_case_flag.hl;
        } else {
            mod->sw->d2->head[inode] = test_case_flag.hr;
        }
    }
}
//****************************************************************//
//****************************************************************//
void write_testcase_error_dam2d(SMODEL *mod) {
    double c0  = sqrt(mod->gravity * test_case_flag.hl);
    double u_analytic = 0., h_analytic = 0.;
    double x = 0.;
    double time = mod->t_prev + mod->dt;
    double xA = 0., xB = 0., xold = 1e10 ;
    
    int inode;
    if(mod->grid->smpi->myid==0){printf("\n");}
    double error_u = 0., error_v = 0., error_h = 0., max_head = 0., max_u = 0., max_v = 0.;
    double error_u_max, error_v_max, error_h_max, max_head_max, max_u_max, max_v_max;
    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        x = mod->grid->node[inode].x;
        xA = test_case_flag.dam_location - time*c0;
        xB = test_case_flag.dam_location + 2*time*c0;
        if (x <= xA) {
            h_analytic = test_case_flag.hl;
            u_analytic = 0.;
        } else if (x > xA && x <= xB) {
            h_analytic = (4./(9. * mod->gravity)) * pow(c0 - ((x-test_case_flag.dam_location)/(2*time)),2);
            u_analytic = (2./3.)*(((x-test_case_flag.dam_location)/mod->t_prev) + c0);
        } else {
            h_analytic = test_case_flag.hr;
            u_analytic = 0.;
        }
        error_u += fabs(mod->sw->d2->vel[inode].x - u_analytic);
        error_h  += fabs(mod->sw->d2->head[inode] - h_analytic);
        if (fabs(h_analytic) > max_head) max_head = fabs(h_analytic);
        if (fabs(u_analytic) > max_head) max_u= fabs(u_analytic);
#ifndef _MESSG
        if (fabs(x - xold) > 1e-6) {
            printf("timestep abs error: x: %13.8e \t analytic head: %13.8e \t numeric head: %13.8e \t error: %13.8e\n",x,h_analytic,mod->sw->d2->head[inode],h_analytic-mod->sw->d2->head[inode]);
            xold = x;
        }
#endif
    }
    
    
    if(mod->grid->smpi->myid==0){printf("time-step average absolute error: h: %20.10e \t u: %20.10e \n",error_h/mod->grid->nnodes,error_u/mod->grid->nnodes);}
    
    
    double grid_mass_error = tl_find_grid_mass_error_elem2d(mod->density, mod->sw->d2->head, mod->sw->d2->vel, mod->grid, mod->flag, mod->initial_grid_mass, mod->series_head, mod->str_values, mod->dt, &total_time_mass_flux);
    
#ifdef _MESSG
    error_v_max=messg_dsum(error_v,mod->grid->smpi->ADH_COMM);
    error_u_max=messg_dsum(error_u,mod->grid->smpi->ADH_COMM);
    error_h_max=messg_dsum(error_h,mod->grid->smpi->ADH_COMM);
    max_v_max=messg_dsum(max_v,mod->grid->smpi->ADH_COMM);
    max_u_max=messg_dsum(max_u,mod->grid->smpi->ADH_COMM);
    max_head_max=messg_dsum(max_head,mod->grid->smpi->ADH_COMM);
    error_v=error_v_max;
    error_u=error_u_max;
    error_h=error_h_max;
    max_head=max_head_max;
    max_v=max_v_max;
    max_u=max_u_max;
#endif
    
    if (max_head < 1e-6) max_head = 1.;
    if (max_u < 1e-6) max_u = 1.;
    if (max_v < 1e-6) max_v = 1.;
    
    error_h /= mod->grid->macro_nnodes;
    error_u /= mod->grid->macro_nnodes;
    
    if (error_h > error_h_tmax) error_h_tmax = error_h;
    if (error_u > error_u_tmax) error_u_tmax = error_u;
    if(mod->grid->smpi->myid==0){
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"last time-step average absoluted errors:\n");
        fprintf(fp,"h: %30.20e :: relative error: %30.20e \n",error_h,error_h/max_head);
        fprintf(fp,"u: %30.20e :: relative error: %30.20e \n",error_u,error_u/max_u);
        //fprintf(fp,"y-velocity abs error: %30.20e :: relative error: %30.20e \n",error_vy/mod->grid->nnodes,error_vy/mod->grid->nnodes/max_v);
        
        fprintf(fp,"\n");
        fprintf(fp,"max average absolute errors over all time: \n");
        fprintf(fp,"h: %20.10e\n",error_h_tmax);
        fprintf(fp,"u: %20.10e\n",error_u_tmax);
        fprintf(fp,"\n");
        fprintf(fp,"grid_mass_error: %30.20e :: relative error: %30.20e \n",grid_mass_error, grid_mass_error/mod->initial_grid_mass);
        fclose(fp);
    }
}


//****************************************************************//
//****************************************************************//
// circular slosh test case
// Bessel function (first kind, zero order) used with one unique nodal circle (s=0)//
// Start at time = 0 //

void testcase_prep_circular_slosh(SMODEL *mod) {
    
    // turn advection off, since analytic solutions assume this is 0
    debug.no_advection = ON;
    
    // since there is no advection, turn SUPG off
    mod->tau_pg = 0.;
    
    double a, R, H, k;
    double x, y, z, r;
    double J10, J11, DDJ0, sigma, g, theta, theta0, depth, depth_old, depth_older;
    int inode;
    
    
    a = test_case_flag.a;
    R = test_case_flag.R;
    H = test_case_flag.H;
    k = 3.83171 / R;
    
    g = mod->gravity;
    sigma = sqrt(g * H) * k;
    
    double time = mod->t_prev;
    double time_old = mod->t_prev - mod->dt;
    double time_older = mod->t_prev - 2*mod->dt;
    
    //printf("time:  %20.10f time_old: %20.10f time_older: %20.10f\n",time, time_old, time_older);
    
    if (time < 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> T0 < 0.\n");
    }
    if (time_old < 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> For slosh model error calculation, please set initial time to T0 + 2 * dt (for hotstart).\n");
    }
    if (time_older < 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> For slosh model error calculation, please set initial time to T0 + 2 * dt (for hotstart).\n");
    }
    
    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        x = mod->grid->node[inode].x;
        y = mod->grid->node[inode].y;
        z = mod->grid->node[inode].z;
        
        r = sqrt(pow(x,2)+pow(y,2));
        theta0 = atan2(y, x);
        if (x >= 0 && y >= 0) theta = pi / 2 - theta0; /* quad 1 */
        if (x <= 0 && y >= 0) theta = 2 * pi - theta0 + pi / 2; /* quad 2 */
        if (x < 0 && y < 0) theta = - theta0 + pi / 2; /* quad 3 */
        if (x > 0 && y < 0) theta = - theta0 + pi / 2; /* quad 4 */
        //printf("theta0: %6.2f theta: %6.2f\n",theta0 * 180/3.1418,theta * 180/3.1418);
        
        mod->sw->d3->displacement[inode] = a * jn(0,k*r) * cos(sigma * time);
        mod->sw->d3->old_displacement[inode] = a * jn(0,k*r) * cos(sigma * time_old);
        mod->sw->d3->older_displacement[inode] = a * jn(0,k*r) * cos(sigma * time_older);
        
        depth = H + mod->sw->d3->displacement[inode];
        depth_old = H + mod->sw->d3->old_displacement[inode];
        depth_older = H + mod->sw->d3->older_displacement[inode];
        
        if (fabs(r) > 1e-6) {
            mod->sw->d3->vel[inode].x = a*g*k*x*jn(1, r*k)*sin(sigma*time)/(r*sigma);
            mod->sw->d3->old_vel[inode].x = a*g*k*x*jn(1, r*k)*sin(sigma*time_old)/(r*sigma);
            mod->sw->d3->older_vel[inode].x = a*g*k*x*jn(1, r*k)*sin(sigma*time_older)/(r*sigma);
            
            mod->sw->d3->vel[inode].y = a*g*k*y*jn(1, r*k)*sin(sigma*time)/(r*sigma);
            mod->sw->d3->old_vel[inode].y = a*g*k*y*jn(1, r*k)*sin(sigma*time_old)/(r*sigma);
            mod->sw->d3->older_vel[inode].y = a*g*k*y*jn(1, r*k)*sin(sigma*time_older)/(r*sigma);
            
            mod->sw->d3->vel[inode].z = -0.5*(-r*k*depth*(jn(2,r*k) - jn(0,r*k)) + 2*depth*jn(1,r*k))*a*g*k*sin(sigma*time)/(r*sigma);
            mod->sw->d3->old_vel[inode].z = -0.5*(-r*k*depth_old*(jn(2,r*k) - jn(0,r*k)) + 2*depth_old*jn(1,r*k))*a*g*k*sin(sigma*time_old)/(r*sigma);
            mod->sw->d3->older_vel[inode].z = -0.5*(-r*k*depth_older*(jn(2,r*k) - jn(0,r*k)) + 2*depth_older*jn(1,r*k))*a*g*k*sin(sigma*time_older)/(r*sigma);
        } else {
            mod->sw->d3->vel[inode].x = 0;
            mod->sw->d3->old_vel[inode].x = 0;
            mod->sw->d3->older_vel[inode].x = 0;
            
            mod->sw->d3->vel[inode].y = 0;
            mod->sw->d3->old_vel[inode].y = 0;
            mod->sw->d3->older_vel[inode].y = 0;
            
            mod->sw->d3->vel[inode].z = 0;
            mod->sw->d3->old_vel[inode].z = 0;
            mod->sw->d3->older_vel[inode].z = 0;
        }
    }
    
    // now distribute displacement down the column
    tl_vertical_adapt(mod->grid, mod->sw->d3->displacement);
    tl_vertical_adapt(mod->grid, mod->sw->d3->old_displacement);
    tl_vertical_adapt(mod->grid, mod->sw->d3->older_displacement);
    
}
//****************************************************************//
//****************************************************************//
// circular slosh test case
// Bessel function (first kind, zero order) used with one unique nodal circle (s=0)//

void write_testcase_error_circular_slosh(SMODEL *mod) {
    
    double a, R, H, k, g;
    double J10, J11, analytic_dpl, analytic_u, analytic_v, analytic_w, sigma, time;
    int nodeID, inode;
    double x, y, z, r, theta0, theta;
    double error_dpl, error_u, error_v, error_w;
    double analytic_grid_dpl, analytic_grid_u, analytic_grid_v, analytic_grid_w;
    double sum_error_dpl, sum_error_u, sum_error_v, sum_error_w, rmse_dpl, rmse_u, rmse_v, rmse_w;
    
    a = test_case_flag.a;
    R = test_case_flag.R;
    H = test_case_flag.H;
    k = 3.83171/R;
    g = mod->gravity;
    
    nodeID = test_case_flag.nodeID;
    x = mod->grid->node[nodeID].x;
    y = mod->grid->node[nodeID].y;
    z = mod->grid->node[nodeID].z;
    r = sqrt(pow(x,2)+pow(y,2));
    
    theta0 = atan2(y,x);
    if (x >= 0 && y >= 0) theta = pi / 2 - theta0; /* quad 1 */
    if (x <= 0 && y >= 0) theta = 2 * pi - theta0 + pi / 2; /* quad 2 */
    if (x < 0 && y < 0) theta = - theta0 + pi / 2; /* quad 3 */
    if (x > 0 && y < 0) theta = - theta0 + pi / 2; /* quad 4 */
    
    sigma = sqrt(g * H) * k;
    time = mod->t_prev + mod->dt; // this routine called before time is updated
    
    // new stuff
    analytic_dpl = a * jn(0,k*r) * cos(sigma * time);
    double depth = H + analytic_dpl;
    analytic_u = a*g*k*x*jn(1,r*k)*sin(sigma*time)/(r*sigma);
    analytic_v = a*g*k*y*jn(1,r*k)*sin(sigma*time)/(r*sigma);
    analytic_w = -0.5*(-r*k*depth*(jn(2,r*k) - jn(0,r*k)) + 2*depth*jn(1,r*k))*a*g*k*sin(sigma*time)/(r*sigma);
    
    double analytic_u_max = fabs(a*g*k*x*jn(1,r*k)/(r*sigma));
    double analytic_v_max = fabs(a*g*k*y*jn(1,r*k)/(r*sigma));
    double analytic_w_max = fabs(-0.5*(-r*k*depth*(jn(2,r*k) - jn(0,r*k)) + 2*depth*jn(1,r*k))*a*g*k/(r*sigma));
    
    // l_infinity errors at user node
    error_dpl = fabs(mod->sw->d3->displacement[nodeID]  - analytic_dpl);
    error_u = fabs(mod->sw->d3->vel[nodeID].x - analytic_u);
    error_v = fabs(mod->sw->d3->vel[nodeID].y - analytic_v);
    error_w = fabs(mod->sw->d3->vel[nodeID].z - analytic_w);
    
    printf("\nnode: %d :: (x,y,z): (%20.10f, %20.10f, %20.10f)\n",nodeID+1,x,y,z);
    printf("max values :: dpl: %10.5e \t u: %10.5e \t v: %10.5e \t w: %10.5e \n",a,analytic_u_max,analytic_v_max,analytic_w_max);
    printf("dpl :: analytic: %10.5e \t numeric: %10.5e \t |error|: %10.5e = %10.5f%%\n",analytic_dpl,  mod->sw->d3->displacement[nodeID], error_dpl, 100*(error_dpl/a));
    printf("u   :: analytic: %10.5e \t numeric: %10.5e \t |error|: %10.5e = %10.5f%%\n",analytic_u, mod->sw->d3->vel[nodeID].x,error_u, 100*(error_u/analytic_u_max));
    printf("v   :: analytic: %10.5e \t numeric: %10.5e \t |error|: %10.5e = %10.5f%%\n",analytic_v, mod->sw->d3->vel[nodeID].y,error_v, 100*(error_v/analytic_v_max));
    printf("w   :: analytic: %10.5e \t numeric: %10.5e \t |error|: %10.5e = %10.5f%%\n",analytic_w, mod->sw->d3->vel[nodeID].z,error_w, 100*(error_w/analytic_w_max));
    
    FILE *fp;  /*temp*/
    fp = fopen("station.out", "a+");  /*temp*/
    fprintf(fp," %30.20e %30.20e %30.20e %30.20e %30.20e %30.20e %30.20e %30.20e %30.20e\n",time,
            mod->sw->d3->displacement[nodeID], analytic_dpl,
            mod->sw->d3->vel[nodeID].x, analytic_u,
            mod->sw->d3->vel[nodeID].y, analytic_v,
            mod->sw->d3->vel[nodeID].z, analytic_w);
    fclose(fp);
    
    // max l-infinity errors over time
    if (fabs(error_dpl) > error_dpl_tmax) error_dpl_tmax = error_dpl;
    if (fabs(error_u) > error_u_tmax) error_u_tmax = error_u;
    if (fabs(error_v) > error_v_tmax) error_v_tmax = error_v;
    if (fabs(error_w) > error_w_tmax) error_w_tmax = error_w;
    
    // calculate errors over entire grid (may include wacky boundary corners)
    analytic_grid_dpl = 0.;
    analytic_grid_u = 0.;
    analytic_grid_v = 0.;
    sum_error_dpl = 0.;
    sum_error_u = 0.;
    sum_error_v = 0.;
    for (inode=0; inode<mod->grid->nnodes; inode++) {
        x = mod->grid->node[inode].x;
        y = mod->grid->node[inode].y;
        z = mod->grid->node[inode].z;
        
        r = sqrt(pow(x,2)+pow(y,2));
        theta0 = atan2(y, x);
        if (x >= 0 && y >= 0) theta = pi / 2 - theta0; /* quad 1 */
        if (x <= 0 && y >= 0) theta = 2 * pi - theta0 + pi / 2; /* quad 2 */
        if (x < 0 && y < 0) theta = - theta0 + pi / 2; /* quad 3 */
        if (x > 0 && y < 0) theta = - theta0 + pi / 2; /* quad 4 */
        
        analytic_dpl= a * jn(0,k*r) * cos(sigma * time);
        depth = H + analytic_dpl;
        
        if (fabs(r) > 1e-6) {
            analytic_grid_u = a*g*k*x*jn(1, r*k)*sin(sigma*time)/(r*sigma);
            analytic_grid_v = a*g*k*y*jn(1, r*k)*sin(sigma*time)/(r*sigma);
            analytic_grid_w = -0.5*(-r*k*depth*(jn(2,r*k) - jn(0,r*k)) + 2*depth*jn(1,r*k))*a*g*k*sin(sigma*time)/(r*sigma);
        } else {
            analytic_grid_u = 0;
            analytic_grid_v = 0;
            analytic_grid_w = 0;
        }
        
        sum_error_dpl += pow(mod->sw->d3->displacement[inode]  - analytic_grid_dpl,2);
        sum_error_u += pow(mod->sw->d3->vel[inode].x - analytic_grid_u,2);
        sum_error_v += pow(mod->sw->d3->vel[inode].y - analytic_grid_v,2);
        sum_error_w += pow(mod->sw->d3->vel[inode].z - analytic_grid_w,2);
    }
    rmse_dpl = sqrt(sum_error_dpl/(float)mod->grid->nnodes);
    rmse_u = sqrt(sum_error_u/(float)mod->grid->nnodes);
    rmse_v = sqrt(sum_error_v/(float)mod->grid->nnodes);
    rmse_w = sqrt(sum_error_w/(float)mod->grid->nnodes);
    
    if (fabs(rmse_dpl) > rmse_error_dpl_tmax) rmse_error_dpl_tmax = rmse_dpl;
    if (fabs(rmse_u) > rmse_error_u_tmax) rmse_error_u_tmax = rmse_u;
    if (fabs(rmse_v) > rmse_error_v_tmax) rmse_error_v_tmax = rmse_v;
    if (fabs(rmse_w) > rmse_error_w_tmax) rmse_error_w_tmax = rmse_w;
    
    
    fp = fopen("error.out", "w");
    fprintf(fp,"node %d max error @ last time-step :: displacment: %10.7e :: relative error: %10.7e%% \n",nodeID+1,error_dpl,100*error_dpl/a);
    fprintf(fp,"node %d max error @ last time-step :: x-velocity: %10.7e :: relative error: %10.7e%% \n",nodeID+1,error_u,100*error_u/analytic_u_max);
    fprintf(fp,"node %d max error @ last time-step :: y-velocity: %10.7e :: relative error: %10.7e%% \n",nodeID+1,error_v,100*error_v/analytic_v_max);
    fprintf(fp,"node %d max error @ last time-step :: w-velocity: %10.7e :: relative error: %10.7e%% \n",nodeID+1,error_w,100*error_w/analytic_w_max);
    
    fprintf(fp,"\n");
    fprintf(fp,"node %d max error over time:: displacment: %10.7e :: relative error: %10.7e%% \n",nodeID+1,error_dpl_tmax,100*error_dpl_tmax/a);
    fprintf(fp,"node %d max error over time :: x-velocity: %10.7e :: relative error: %10.7e%% \n",nodeID+1,error_u_tmax,100*error_u_tmax/analytic_u_max);
    fprintf(fp,"node %d max error over time :: y-velocity: %10.7e :: relative error: %10.7e%% \n",nodeID+1,error_v_tmax,100*error_v_tmax/analytic_v_max);
    fprintf(fp,"node %d max error over time :: w-velocity: %10.7e :: relative error: %10.7e%% \n",nodeID+1,error_w_tmax,100*error_w_tmax/analytic_w_max);
    
    
    fprintf(fp,"\n");
    fprintf(fp,"grid RMSE max error over time:: displacment: %10.7e :: relative error: %10.7e%% \n",rmse_error_dpl_tmax,100*rmse_error_dpl_tmax/a);
    fprintf(fp,"grid RMSE max error over time :: x-velocity: %10.7e :: relative error: %10.7e%% \n",rmse_error_u_tmax,100*rmse_error_u_tmax/analytic_u_max);
    fprintf(fp,"grid RMSE max error over time :: y-velocity: %10.7e :: relative error: %10.7e%% \n",rmse_error_v_tmax,100*rmse_error_v_tmax/analytic_v_max);
    fprintf(fp,"grid RMSE max error over time :: w-velocity: %10.7e :: relative error: %10.7e%% \n",rmse_error_w_tmax,100*rmse_error_w_tmax/analytic_w_max);
    fclose(fp);
    
    
}

//********************************************************************************
//********************************************************************************
// OUTFLOW TEST CASE
//****************************************************************
//****************************************************************
// Note: at first the errors will grow, as it takes some finite amount of time for
// water to travel the length of the grid.  The final time errors are all the matters here.
void write_testcase_error_outflow(SMODEL *mod) {
    
    int i = 0, inode = 0;
    double new_grid_mass = 0.;
    double mass_error = tl_find_3d_grid_mass_error(mod->str_values, mod->series_head, mod->initial_grid_mass, mod->density, mod->grid, mod->sw->d3->vel, mod->sw->d3->displacement,  mod->sw->d3->old_displacement,  mod->sw->d3->older_displacement, mod->dt, &new_grid_mass, &total_time_mass_flux);
    if (mod->grid->smpi->myid == 0){
        printf("\ninitial grid mass: %-12.10e(kg) \t grid mass: %-12.10e(kg) total_time_mass_flux: %-12.10e(kg) \t mass_error: %-12.10e(kg) \t relative mass eror: %-12.10e%%\n",
               mod->initial_grid_mass,
               new_grid_mass,
               total_time_mass_flux,
               mass_error,
               100*(mass_error/(mod->initial_grid_mass)));
        
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"final time mass absolute error: %-12.10e(kg) :: relative error: %-12.10e%% \n",mass_error, 100*(mass_error/mod->initial_grid_mass));
        fclose(fp);
    }
}

//********************************************************************************
//********************************************************************************
// TIDAL TEST CASE (ASCE_35)
// cjt :: for a rotated tidal grid with ASCE tidal inputs

void get_tide_analytic(double x, double t, double depth, int flag, double *dpl, double *u, double *v, double *w) {
    
    //    // solutions in orthogonal coordinates
    //    *dpl = ((-132.1868227653527*cos(0.00014075235903181*t) + (6.66133814775094e-15)*sin(0.00014075235903181*t))*cos(1.330997269696963*log(x)) + 3.297864281790462*cos(0.00014075235903181*t)*sin(1.330997269696963*log(x)))/sqrt(x);
    //
    //    if (flag==1) {
    //        depth += *dpl;
    //    }
    //
    //    *v= 0;
    //
    //    *u = ((-(2.31900602934853e-10)*cos(0.00014075235903181*t) - 4907427.699112926*sin(0.00014075235903181*t) + 2.31900602934853e-10)*cos(1.330997269696963*log(x)) + (-(6.17318138694737e-10)*cos(0.00014075235903181*t) - 12135181.3826516*sin(0.00014075235903181*t) + 6.17318138694737e-10)*sin(1.330997269696963*log(x)))/pow(x,1.5);
    //
    //    *w = (((4.63801205869706e-19)*cos(0.00014075235903181*t) + 0.009814855398225852*sin(0.00014075235903181*t) - 4.63801205869706e-19)/sqrt(x) + ((4.73797852734827e-10)*depth*cos(0.00014075235903181*t) + 8790751.73891731*depth*sin(0.00014075235903181*t) - (4.73797852734827e-10)*depth)/pow(x,2.5))*cos(1.330997269696963*log(x)) + (((1.23463627738947e-18)*cos(0.00014075235903181*t) + 0.0242703627653032*sin(0.00014075235903181*t) - 1.23463627738947e-18)/sqrt(x) + (-(1.23463627738948e-09)*depth*cos(0.00014075235903181*t) - 24734544.942732*depth*sin(0.00014075235903181*t) + (1.23463627738948e-09)*depth)/pow(x,2.5))*sin(1.330997269696963*log(x));
    
    
    
    // solutions in orthogonal coordinates
    double tau = 0; // linear friction
    double pi = 3.14159265359;
    double g = 9.8000000000000000000;
    double Ho = 1.0000000000000e-9;
    double x1 = 100000.00000000000000;
    double x2 = 200000.0000000000000;
    double a = 0.250000000000;
    double wT = (2.000000000*pi)/(12.400000000000 * 3600.00000000000000);
    double complex B2 = (wT*wT - I*wT*tau)/(g*Ho);
    double complex s1 = -(1./2.) + csqrt((1./4.) - B2);
    double complex s2 = -(1./2.) - csqrt((1./4.) - B2);
    double complex A =  (a * s2 * cpow(x1,s2))/(s2*cpow(x2,s1)*cpow(x1,s2) - s1*cpow(x1,s1)*cpow(x2,s2));
    double complex B = -(a * s1 * cpow(x1,s1))/(s2*cpow(x2,s1)*cpow(x1,s2) - s1*cpow(x1,s1)*cpow(x2,s2));
    
    //    *dpl =  (cos(t*wT)*creal(A) - cimag(A)*sin(t*wT))*pow(x,creal(s1))*cos(cimag(s1)*log(x)) + (cos(t*wT)*creal(B) - cimag(B)*sin(t*wT))*pow(x,creal(s2))*cos(cimag(s2)*log(x)) - (cos(t*wT)*cimag(A) + creal(A)*sin(t*wT))*pow(x,creal(s1))*sin(cimag(s1)*log(x)) - (cos(t*wT)*cimag(B) + creal(B)*sin(t*wT))*pow(x,creal(s2))*sin(cimag(s2)*log(x));
    //
    //    //depth += *dpl;
    //
    //
    //    *u =  ((cimag(A)*cimag(s1) - creal(A)*creal(s1))*pow(x,creal(s1))*cos(cimag(s1)*log(x)) + (cimag(B)*cimag(s2) - creal(B)*creal(s2))*pow(x,creal(s2))*cos(cimag(s2)*log(x)) + (cimag(s1)*creal(A) + cimag(A)*creal(s1))*pow(x,creal(s1))*sin(cimag(s1)*log(x)) + (cimag(s2)*creal(B) + cimag(B)*creal(s2))*pow(x,creal(s2))*sin(cimag(s2)*log(x)))*g*sin(t*wT)/(wT*x);
    //    *v =  0;
    //    *w =  ((cimag(A)*pow(cimag(s1),2) - cimag(A)*pow(creal(s1),2) - (cimag(A)*pow(cimag(s1),2) - cimag(A)*pow(creal(s1),2) + cimag(s1)*creal(A) - (2*cimag(s1)*creal(A) - cimag(A))*creal(s1))*cos(t*wT) + cimag(s1)*creal(A) - (2*cimag(s1)*creal(A) - cimag(A))*creal(s1) - (pow(cimag(s1),2)*creal(A) - creal(A)*pow(creal(s1),2) - cimag(A)*cimag(s1) + (2*cimag(A)*cimag(s1) + creal(A))*creal(s1))*sin(t*wT))*pow(x,creal(s1))*cos(cimag(s1)*log(x)) + (cimag(B)*pow(cimag(s2),2) - cimag(B)*pow(creal(s2),2) - (cimag(B)*pow(cimag(s2),2) - cimag(B)*pow(creal(s2),2) + cimag(s2)*creal(B) - (2*cimag(s2)*creal(B) - cimag(B))*creal(s2))*cos(t*wT) + cimag(s2)*creal(B) - (2*cimag(s2)*creal(B) - cimag(B))*creal(s2) - (pow(cimag(s2),2)*creal(B) - creal(B)*pow(creal(s2),2) - cimag(B)*cimag(s2) + (2*cimag(B)*cimag(s2) + creal(B))*creal(s2))*sin(t*wT))*pow(x,creal(s2))*cos(cimag(s2)*log(x)) + (pow(cimag(s1),2)*creal(A) - creal(A)*pow(creal(s1),2) - (pow(cimag(s1),2)*creal(A) - creal(A)*pow(creal(s1),2) - cimag(A)*cimag(s1) + (2*cimag(A)*cimag(s1) + creal(A))*creal(s1))*cos(t*wT) - cimag(A)*cimag(s1) + (2*cimag(A)*cimag(s1) + creal(A))*creal(s1) + (cimag(A)*pow(cimag(s1),2) - cimag(A)*pow(creal(s1),2) + cimag(s1)*creal(A) - (2*cimag(s1)*creal(A) - cimag(A))*creal(s1))*sin(t*wT))*pow(x,creal(s1))*sin(cimag(s1)*log(x)) + (pow(cimag(s2),2)*creal(B) - creal(B)*pow(creal(s2),2) - (pow(cimag(s2),2)*creal(B) - creal(B)*pow(creal(s2),2) - cimag(B)*cimag(s2) + (2*cimag(B)*cimag(s2) + creal(B))*creal(s2))*cos(t*wT) - cimag(B)*cimag(s2) + (2*cimag(B)*cimag(s2) + creal(B))*creal(s2) + (cimag(B)*pow(cimag(s2),2) - cimag(B)*pow(creal(s2),2) + cimag(s2)*creal(B) - (2*cimag(s2)*creal(B) - cimag(B))*creal(s2))*sin(t*wT))*pow(x,creal(s2))*sin(cimag(s2)*log(x)))*depth*g/(wT*x*x);
    //
    //    double wb = 2*(((cimag(s1)*creal(A) + cimag(A)*creal(s1))*cos(t*wT) - cimag(s1)*creal(A) - cimag(A)*creal(s1) - (cimag(A)*cimag(s1) - creal(A)*creal(s1))*sin(t*wT))*pow(x,creal(s1))*cos(cimag(s1)*log(x)) + ((cimag(s2)*creal(B) + cimag(B)*creal(s2))*cos(t*wT) - cimag(s2)*creal(B) - cimag(B)*creal(s2) - (cimag(B)*cimag(s2) - creal(B)*creal(s2))*sin(t*wT))*pow(x,creal(s2))*cos(cimag(s2)*log(x)) - ((cimag(A)*cimag(s1) - creal(A)*creal(s1))*cos(t*wT) - cimag(A)*cimag(s1) + creal(A)*creal(s1) + (cimag(s1)*creal(A) + cimag(A)*creal(s1))*sin(t*wT))*pow(x,creal(s1))*sin(cimag(s1)*log(x)) - ((cimag(B)*cimag(s2) - creal(B)*creal(s2))*cos(t*wT) - cimag(B)*cimag(s2) + creal(B)*creal(s2) + (cimag(s2)*creal(B) + cimag(B)*creal(s2))*sin(t*wT))*pow(x,creal(s2))*sin(cimag(s2)*log(x)))*Ho*g/wT;
    //
    //    *w += wb;
    
    *dpl =  (cos(t*wT)*creal(A) - cimag(A)*sin(t*wT))*pow(x,creal(s1))*cos(cimag(s1)*log(x)) + (cos(t*wT)*creal(B) - cimag(B)*sin(t*wT))*pow(x,creal(s2))*cos(cimag(s2)*log(x)) - (cos(t*wT)*cimag(A) + creal(A)*sin(t*wT))*pow(x,creal(s1))*sin(cimag(s1)*log(x)) - (cos(t*wT)*cimag(B) + creal(B)*sin(t*wT))*pow(x,creal(s2))*sin(cimag(s2)*log(x));
    
    //depth += *dpl;
    
    
    *u =  ((cimag(A)*cimag(s1) - creal(A)*creal(s1))*pow(x,creal(s1))*cos(cimag(s1)*log(x)) + (cimag(B)*cimag(s2) - creal(B)*creal(s2))*pow(x,creal(s2))*cos(cimag(s2)*log(x)) + (cimag(s1)*creal(A) + cimag(A)*creal(s1))*pow(x,creal(s1))*sin(cimag(s1)*log(x)) + (cimag(s2)*creal(B) + cimag(B)*creal(s2))*pow(x,creal(s2))*sin(cimag(s2)*log(x)))*g*sin(t*wT)/(wT*x);
    *v =  0;
    *w =  ((cimag(A)*pow(cimag(s1),2) - cimag(A)*pow(creal(s1),2) - (cimag(A)*pow(cimag(s1),2) - cimag(A)*pow(creal(s1),2) + cimag(s1)*creal(A) - (2*cimag(s1)*creal(A) - cimag(A))*creal(s1))*cos(t*wT) + cimag(s1)*creal(A) - (2*cimag(s1)*creal(A) - cimag(A))*creal(s1) - (pow(cimag(s1),2)*creal(A) - creal(A)*pow(creal(s1),2) - cimag(A)*cimag(s1) + (2*cimag(A)*cimag(s1) + creal(A))*creal(s1))*sin(t*wT))*pow(x,creal(s1))*cos(cimag(s1)*log(x)) + (cimag(B)*pow(cimag(s2),2) - cimag(B)*pow(creal(s2),2) - (cimag(B)*pow(cimag(s2),2) - cimag(B)*pow(creal(s2),2) + cimag(s2)*creal(B) - (2*cimag(s2)*creal(B) - cimag(B))*creal(s2))*cos(t*wT) + cimag(s2)*creal(B) - (2*cimag(s2)*creal(B) - cimag(B))*creal(s2) - (pow(cimag(s2),2)*creal(B) - creal(B)*pow(creal(s2),2) - cimag(B)*cimag(s2) + (2*cimag(B)*cimag(s2) + creal(B))*creal(s2))*sin(t*wT))*pow(x,creal(s2))*cos(cimag(s2)*log(x)) + (pow(cimag(s1),2)*creal(A) - creal(A)*pow(creal(s1),2) - (pow(cimag(s1),2)*creal(A) - creal(A)*pow(creal(s1),2) - cimag(A)*cimag(s1) + (2*cimag(A)*cimag(s1) + creal(A))*creal(s1))*cos(t*wT) - cimag(A)*cimag(s1) + (2*cimag(A)*cimag(s1) + creal(A))*creal(s1) + (cimag(A)*pow(cimag(s1),2) - cimag(A)*pow(creal(s1),2) + cimag(s1)*creal(A) - (2*cimag(s1)*creal(A) - cimag(A))*creal(s1))*sin(t*wT))*pow(x,creal(s1))*sin(cimag(s1)*log(x)) + (pow(cimag(s2),2)*creal(B) - creal(B)*pow(creal(s2),2) - (pow(cimag(s2),2)*creal(B) - creal(B)*pow(creal(s2),2) - cimag(B)*cimag(s2) + (2*cimag(B)*cimag(s2) + creal(B))*creal(s2))*cos(t*wT) - cimag(B)*cimag(s2) + (2*cimag(B)*cimag(s2) + creal(B))*creal(s2) + (cimag(B)*pow(cimag(s2),2) - cimag(B)*pow(creal(s2),2) + cimag(s2)*creal(B) - (2*cimag(s2)*creal(B) - cimag(B))*creal(s2))*sin(t*wT))*pow(x,creal(s2))*sin(cimag(s2)*log(x)))*depth*g/(wT*x*x);
    
    double wb = 2*(((cimag(s1)*creal(A) + cimag(A)*creal(s1))*cos(t*wT) - cimag(s1)*creal(A) - cimag(A)*creal(s1) - (cimag(A)*cimag(s1) - creal(A)*creal(s1))*sin(t*wT))*pow(x,creal(s1))*cos(cimag(s1)*log(x)) + ((cimag(s2)*creal(B) + cimag(B)*creal(s2))*cos(t*wT) - cimag(s2)*creal(B) - cimag(B)*creal(s2) - (cimag(B)*cimag(s2) - creal(B)*creal(s2))*sin(t*wT))*pow(x,creal(s2))*cos(cimag(s2)*log(x)) - ((cimag(A)*cimag(s1) - creal(A)*creal(s1))*cos(t*wT) - cimag(A)*cimag(s1) + creal(A)*creal(s1) + (cimag(s1)*creal(A) + cimag(A)*creal(s1))*sin(t*wT))*pow(x,creal(s1))*sin(cimag(s1)*log(x)) - ((cimag(B)*cimag(s2) - creal(B)*creal(s2))*cos(t*wT) - cimag(B)*cimag(s2) + creal(B)*creal(s2) + (cimag(s2)*creal(B) + cimag(B)*creal(s2))*sin(t*wT))*pow(x,creal(s2))*sin(cimag(s2)*log(x)))*Ho*g/wT;
    
    *w += wb;
    
    
    //    printf("A: %20.10f %20.10f\n",creal(A),cimag(A));
    //    printf("b: %20.10f %20.10f\n",creal(B),cimag(B));
    //    printf("s1: %20.10f %20.10f\n",creal(s1),cimag(s1));
    //    printf("s1: %20.10f %20.10f\n",creal(s1),cimag(s2));
    //    printf("t: %20.10f x: %20.10f dpl: %20.10f\n",t,x,*dpl);
    //    printf("u: %20.10f ugrid: %20.10e\n",*u,sqrt(2.)/2.*(*u));
    //    exit(-1);
}

//****************************************************************//
//****************************************************************//
void testcase_prep_tide_2d(SMODEL *mod) {
    
    //if (mod->proc_flag != 1) return;
    
    //debug.no_advection = ON; // turn advection off, since analytic solutions assume this is 0
    //mod->tau_pg = 0.; // since there is no advection, turn SUPG off
    
    int i, nodeID = UNSET_INT, nd1 = UNSET_INT, nd2 = UNSET_INT;
    
    int iseg, surface_node, inode;
    double g = mod->gravity;
    double x,y,z, depth, analytic_elevation, analytic_u, analytic_v, analytic_w;
    double time = mod->t_prev;
    double time_old = mod->t_prev - mod->dt;
    double time_older = mod->t_prev - 2*mod->dt;
    
    if (time < 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> T0 < 0.\n");
    }
    if (time_old < 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> For slosh model error calculation, please set initial time to T0 + 2 * dt (for hotstart).\n");
    }
    if (time_older < 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> For slosh model error calculation, please set initial time to T0 + 2 * dt (for hotstart).\n");
    }
    
    // make sure only the owner of the test node is calculating errors
    int have_node=0,have_node1=0, have_node2=0;
#ifdef _MESSG
    for(i=0;i<mod->grid->my_nnodes;i++){
        if(test_case_flag.node1_2d==mod->grid->node[i].gid){
            have_node1=1;
            nd1=i;
        }
        if(test_case_flag.node2_2d==mod->grid->node[i].gid){
            have_node2=1;
            nd2=i;
        }
        if(test_case_flag.nodeID_2d==mod->grid->node[i].gid){
            have_node=1;
            nodeID=i;
        }
    }
#else
    nodeID = test_case_flag.nodeID_2d;
    nd1    = test_case_flag.node1_2d;
    nd2    = test_case_flag.node2_2d;
    have_node1=1;
    have_node2=1;
    have_node=1;
#endif
    
    // for plotting points along a cross-section
    SVECT p1, p2, pointCheck;
    if((have_node1==1)&&(have_node2==1)){
        p1.x = mod->grid->node[nd1].x; p1.y = mod->grid->node[nd1].y;
        p2.x = mod->grid->node[nd2].x; p2.y = mod->grid->node[nd2].y;
    }
    
    FILE *fp=NULL, *fp_surface=NULL, *fp_bed=NULL;
    if(have_node==1) fp = fopen("station2d.out", "w");
    if(have_node1==1 && have_node2==1) fp_surface = fopen("surfaceXsection2d.out", "w");
    if(have_node1==1 && have_node2==1) fp_bed = fopen("bedXsection2d.out", "w");

    double ut,vt,xt;
    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        
        // rotated coordinates by 45 degrees
        x = mod->grid->node[inode].x;  if (fabs(x)<1) tl_error("x cannot be zero\n");
        y = mod->grid->node[inode].y;
        z = mod->grid->node[inode].z; // bed node
        
        depth = -z; // 2D doesn't use this, since it's only needed for w
        
        // convert back to orthogonal coordinates
        xt = (x+y)/sqrt(2.);
        
        if (fabs(xt)<1) tl_error("x cannot be zero\n");
        get_tide_analytic(xt, time_older, depth, 1, &analytic_elevation, &analytic_u, &analytic_v, &analytic_w);
        ut = analytic_u;
        vt = analytic_v;
        analytic_u = (sqrt(2.)/2.) * (ut - vt); // in rotated coordinates
        analytic_v = (sqrt(2.)/2.) * (ut + vt); // in rotated coordinates
        mod->sw->d2->older_head[inode]  = analytic_elevation - z;
        mod->sw->d2->older_vel[inode].x = analytic_u; // in rotated coordinates
        mod->sw->d2->older_vel[inode].y = analytic_v; // in rotated coordinates
        if ((mod->grid->node[inode].gid == nodeID)&&(have_node==1))
            fprintf(fp,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
                    time_older,analytic_elevation,analytic_elevation,analytic_u,analytic_u,analytic_v,analytic_v,analytic_w,analytic_w);
        
        get_tide_analytic(xt, time_old, depth, 1, &analytic_elevation, &analytic_u, &analytic_v, &analytic_w);
        ut = analytic_u;
        vt = analytic_v;
        analytic_u = (sqrt(2.)/2.) * (ut - vt); // in rotated coordinates
        analytic_v = (sqrt(2.)/2.) * (ut + vt); // in rotated coordinates
        mod->sw->d2->old_head[inode]  = analytic_elevation - z;
        mod->sw->d2->old_vel[inode].x = analytic_u; // in rotated coordinates
        mod->sw->d2->old_vel[inode].y = analytic_v; // in rotated coordinates
        if ((mod->grid->node[inode].gid == nodeID)&&(have_node==1))
            fprintf(fp,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
                    time_old,analytic_elevation,analytic_elevation,analytic_u,analytic_u,analytic_v,analytic_v,analytic_w,analytic_w);
        
        get_tide_analytic(xt, time, depth, 1, &analytic_elevation, &analytic_u, &analytic_v, &analytic_w);
        ut = analytic_u;
        vt = analytic_v;
        analytic_u = (sqrt(2.)/2.) * (ut - vt); // in rotated coordinates
        analytic_v = (sqrt(2.)/2.) * (ut + vt); // in rotated coordinates
        mod->sw->d2->head[inode]  = analytic_elevation - z;
        mod->sw->d2->vel[inode].x = analytic_u; // in rotated coordinates
        mod->sw->d2->vel[inode].y = analytic_v; // in rotated coordinates
        if ((mod->grid->node[inode].gid == nodeID)&&(have_node==1))
            fprintf(fp,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
                    time,analytic_elevation,analytic_elevation,analytic_u,analytic_u,analytic_v,analytic_v,analytic_w,analytic_w);
    }
    if(have_node==1) fclose(fp);
    if((have_node1==1)&&(have_node2==1)) {
        fclose(fp_surface);
        fclose(fp_bed);
    }
    
}

//****************************************************************//
//****************************************************************//
void testcase_prep_tide_3d(SMODEL *mod) {
    
    //if (mod->proc_flag != 1) return;
    
    //debug.no_advection = ON; // turn advection off, since analytic solutions assume this is 0
    //mod->tau_pg = 0.; // since there is no advection, turn SUPG off
    
    int i, nodeID = UNSET_INT, nd1 = UNSET_INT, nd2 = UNSET_INT;
    
    int iseg, surface_node, bed_node, inode;
    double g = mod->gravity;
    double x,y,z, depth, analytic_elevation, analytic_u, analytic_v, analytic_w;
    double time = mod->t_prev;
    double time_old = mod->t_prev - mod->dt;
    double time_older = mod->t_prev - 2*mod->dt;
    
    if (time < 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> T0 < 0.\n");
    }
    if (time_old < 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> For slosh model error calculation, please set initial time to T0 + 2 * dt (for hotstart).\n");
    }
    if (time_older < 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> For slosh model error calculation, please set initial time to T0 + 2 * dt (for hotstart).\n");
    }
    
    // update total depths just in case
    tl_calculate_depthavgvel(mod->grid, mod->sw->d3);
    
    // make sure only the owner of the test node is calculating errors
    int have_node=0,have_node1=0, have_node2=0;
#ifdef _MESSG
    for(i=0;i<mod->grid->my_nnodes;i++){
        if(test_case_flag.node1_3d==mod->grid->node[i].gid){
            have_node1=1;
            nd1=i;
        }
        if(test_case_flag.node2_3d==mod->grid->node[i].gid){
            have_node2=1;
            nd2=i;
        }
        if(test_case_flag.nodeID_3d==mod->grid->node[i].gid){
            have_node=1;
            nodeID=i;
        }
    }
#else
    nodeID = test_case_flag.nodeID_3d;
    nd1    = test_case_flag.node1_3d;
    nd2    = test_case_flag.node2_3d;
    have_node1=1;
    have_node2=1;
    have_node=1;
#endif
    
    // for cross-sectional printing
    SVECT p1, p2, pointCheck;
    if((have_node1==1)&&(have_node2==1)){
        p1.x = mod->grid->node[nd1].x; p1.y = mod->grid->node[nd1].y;
        p2.x = mod->grid->node[nd2].x; p2.y = mod->grid->node[nd2].y;
    }
    
    // open all files
    FILE *fp=NULL, *fp_surface=NULL, *fp_bed=NULL;
    if(have_node==1) fp = fopen("station3d.out", "w");
    if(have_node1==1 && have_node2==1) fp_surface = fopen("surfaceXsection3d.out", "w");
    if(have_node1==1 && have_node2==1) fp_bed = fopen("bedXsection3d.out", "w");
    
    double ut, vt, xt, xx;
    for (inode=0; inode<mod->grid->my_nnodes; inode++) {
        
        // rotated coordinates by 45 degrees
        x = mod->grid->node[inode].x;  if (fabs(x)<1) tl_error("x cannot be zero\n");
        y = mod->grid->node[inode].y;
        z = mod->grid->node[inode].z;
        
        // convert back to orthogonal coordinates
        xt = (x+y)/sqrt(2.);
        
        iseg = find_vertical_segment(mod->grid, inode, mod->grid->vertical_hash);
        surface_node = mod->grid->nodeID_3d_to_2d_sur[mod->grid->vertical_list[iseg]->id];
        depth = mod->sw->d3->depth[surface_node];
        
        if (fabs(xt)<1) tl_error("x cannot be zero\n");
        get_tide_analytic(xt, time_older, depth, 1, &analytic_elevation, &analytic_u, &analytic_v, &analytic_w);
        ut = analytic_u;
        vt = analytic_v;
        analytic_u = (sqrt(2.)/2.) * (ut - vt); // in rotated coordinates
        analytic_v = (sqrt(2.)/2.) * (ut + vt); // in rotated coordinates
        mod->sw->d3->older_displacement[inode]  = analytic_elevation;
        mod->sw->d3->older_vel[inode].x = analytic_u; // in rotated coordinates
        mod->sw->d3->older_vel[inode].y = analytic_v; // in rotated coordinates
        mod->sw->d3->older_vel[inode].z = analytic_w;
        if (inode == nodeID && have_node==1)
            fprintf(fp,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
                    time_older,analytic_elevation,analytic_elevation,analytic_u,analytic_u,analytic_v,analytic_v,analytic_w,analytic_w);
        
        get_tide_analytic(xt, time_old, depth, 1, &analytic_elevation, &analytic_u, &analytic_v, &analytic_w);
        ut = analytic_u;
        vt = analytic_v;
        analytic_u = (sqrt(2.)/2.) * (ut - vt); // in rotated coordinates
        analytic_v = (sqrt(2.)/2.) * (ut + vt); // in rotated coordinates
        mod->sw->d3->old_displacement[inode]  = analytic_elevation;
        mod->sw->d3->old_vel[inode].x = analytic_u; // in rotated coordinates
        mod->sw->d3->old_vel[inode].y = analytic_v; // in rotated coordinates
        mod->sw->d3->old_vel[inode].z = analytic_w;
        if (inode == nodeID && have_node==1)
            fprintf(fp,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
                    time_old,analytic_elevation,analytic_elevation,analytic_u,analytic_u,analytic_v,analytic_v,analytic_w,analytic_w);
        
        get_tide_analytic(xt, time, depth, 1, &analytic_elevation, &analytic_u, &analytic_v, &analytic_w);
        ut = analytic_u;
        vt = analytic_v;
        analytic_u = (sqrt(2.)/2.) * (ut - vt); // in rotated coordinates
        analytic_v = (sqrt(2.)/2.) * (ut + vt); // in rotated coordinates
        mod->sw->d3->displacement[inode]  = analytic_elevation;
        mod->sw->d3->vel[inode].x = analytic_u; // in rotated coordinates
        mod->sw->d3->vel[inode].y = analytic_v; // in rotated coordinates
        mod->sw->d3->vel[inode].z = analytic_w;
        if (inode == nodeID  && have_node==1)
            fprintf(fp,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
                    time,analytic_elevation,analytic_elevation,analytic_u,analytic_u,analytic_v,analytic_v,analytic_w,analytic_w);
    }
    if(have_node==1) fclose(fp);
    
    // now distribute displacement down the column
    tl_vertical_adapt(mod->grid, mod->sw->d3->displacement);
    tl_vertical_adapt(mod->grid, mod->sw->d3->old_displacement);
    tl_vertical_adapt(mod->grid, mod->sw->d3->older_displacement);
    
    // print initial cross-section reults
    if((have_node1==1)&&(have_node2==1)){
        double adh_elevation, adh_u, adh_v, adh_w;
        for (inode=0; inode<mod->grid->my_nnodes; inode++) {
            
            pointCheck.x = mod->grid->node[inode].x; pointCheck.y = mod->grid->node[inode].y;
            if (isPointOnLine(p1,p2,pointCheck) == TRUE)  {
                
                iseg = find_vertical_segment(mod->grid, inode, mod->grid->vertical_hash);
                surface_node = mod->grid->vertical_list[iseg]->id; // this time, the 3D node ID of surface node
                bed_node = (mod->grid->vertical_list[iseg]->prev)->id;
                
                x = mod->grid->node[inode].x;
                y = mod->grid->node[inode].y;
                z = mod->grid->node[inode].z;
                
                depth = mod->sw->d3->depth[mod->grid->nodeID_3d_to_2d_sur[surface_node]];
                
                adh_elevation = mod->grid->node[surface_node].z + mod->sw->d3->older_displacement[surface_node];
                adh_u = mod->sw->d3->older_vel[inode].x;
                adh_v = mod->sw->d3->older_vel[inode].y;
                adh_w = mod->sw->d3->older_vel[inode].z;
                
                xx = (x+y)/sqrt(2.); // convert to orthognal coordinates
                get_tide_analytic(xx, time_older, depth, 0, &analytic_elevation, &analytic_u, &analytic_v, &analytic_w);
                ut = analytic_u;
                vt = analytic_v;
                analytic_u = (sqrt(2.)/2.) * (ut - vt); // in rotated coordinates
                analytic_v = (sqrt(2.)/2.) * (ut + vt); // in rotated coordinates
                
                
                if (inode == surface_node) {
                    //printf("node: %d x: %20.10f y: %20.10f z: %20.10f\n",inode+1,mod->grid->node[inode].x,mod->grid->node[inode].y,mod->grid->node[inode].z);
                    fprintf(fp_surface,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
                            time_older,x,y,analytic_elevation,adh_elevation,analytic_u,adh_u,analytic_v,adh_v,analytic_w,adh_w);
                }
                if (inode == bed_node) {
                    //printf("node: %d x: %20.10f y: %20.10f z: %20.10f\n",inode+1,mod->grid->node[inode].x,mod->grid->node[inode].y,mod->grid->node[inode].z);
                    fprintf(fp_bed,   "%20.10e %20.10e%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
                            time_older,x,y,analytic_elevation,adh_elevation,analytic_u,adh_u,analytic_v,adh_v,analytic_w,adh_w);
                }
            }
        }
        fprintf(fp_surface,"\n\n");
        fprintf(fp_bed,"\n\n");
        fclose(fp_surface);
        fclose(fp_bed);
    }
}

//****************************************************************
//****************************************************************
void write_testcase_error_tide_2d(SMODEL *mod) {
    
    if (mod->proc_flag != 1) return;
    
    assert(test_case_flag.nodeID_2d >= 0 && test_case_flag.nodeID_2d < mod->grid->macro_nnodes);
    assert(test_case_flag.node1_2d  >= 0 && test_case_flag.node1_2d < mod->grid->macro_nnodes);
    assert(test_case_flag.node2_2d  >= 0 && test_case_flag.node2_2d < mod->grid->macro_nnodes);
    
    double g = mod->gravity;
    double time = mod->t_prev + mod->dt; // this routine called before time is updated
    SVECT p1, p2, pointCheck;
    int inode, i, iseg, myid = UNSET_INT;
    double x,y,z,depth,xx,analytic_elevation,analytic_u,analytic_v,analytic_w,ut,vt;
    double adh_elevation,adh_u,adh_v,abs_error_elevation,abs_error_u,abs_error_v;
    FILE *fp;
    
    // node to calculate error must be at surface
    int nodeID = UNSET_INT, nd1 = UNSET_INT, nd2 = UNSET_INT;
    
    // make sure only the owner of the test node is calculating errors
    int have_node=0,have_node1=0, have_node2=0;
#ifdef _MESSG
    for(i=0;i<mod->grid->my_nnodes;i++){
        if(test_case_flag.node1_2d==mod->grid->node[i].gid){
            have_node1=1;
            nd1=i;
        }
        if(test_case_flag.node2_2d==mod->grid->node[i].gid){
            have_node2=1;
            nd2=i;
        }
        if(test_case_flag.nodeID_2d==mod->grid->node[i].gid){
            have_node=1;
            nodeID=i;
        }
    }
#else
    nodeID = test_case_flag.nodeID_2d;
    nd1    = test_case_flag.node1_2d;
    nd2    = test_case_flag.node2_2d;
    have_node1=1;
    have_node2=1;
    have_node=1;
#endif
    
    FILE *fp_station=NULL, *fp_error=NULL, *fp_surface=NULL, *fp_bed=NULL;
    if(have_node==1) {
        fp_error = fopen("error2d.out", "w");
        fp_station = fopen("station2d.out", "a");
    }
    if(have_node1==1 && have_node2==1) {
        fp_surface = fopen("surfaceXsection2d.out", "a");
        fp_bed = fopen("bedXsection2d.out", "a");
    }
    
    // for cross section plot files
    if((have_node1==1)&&(have_node2==1)){
        p1.x = mod->grid->node[nd1].x; p1.y = mod->grid->node[nd1].y;
        p2.x = mod->grid->node[nd2].x; p2.y = mod->grid->node[nd2].y;
    }
    
    if(have_node==1){
        x = mod->grid->node[nodeID].x;
        y = mod->grid->node[nodeID].y;
        z = mod->grid->node[nodeID].z;
        
        // calculate total depth at test node
        depth = -z;
        
        xx = (x+y)/sqrt(2.); // convert to orthognal coordinates
        get_tide_analytic(xx, time, depth, 0, &analytic_elevation, &analytic_u, &analytic_v, &analytic_w);
        ut = analytic_u;
        vt = analytic_v;
        analytic_u = (sqrt(2.)/2.) * (ut - vt); // in rotated coordinates; // convert to rotated
        analytic_v = (sqrt(2.)/2.) * (ut + vt); // in rotated coordinates; // convert to rotated
        
        // calculate AdH elevation
        adh_elevation = z + mod->sw->d2->head[nodeID];
        adh_u = mod->sw->d2->vel[nodeID].x;
        adh_v = mod->sw->d2->vel[nodeID].y;
        
        // absolute errors
        abs_error_elevation = fabs(analytic_elevation - adh_elevation);
        abs_error_u = fabs(analytic_u - adh_u);
        abs_error_v = fabs(analytic_v - adh_v);
        
        // max l-infinity errors over time
        if (fabs(abs_error_elevation) > error_dpl_tmax) error_dpl_tmax = abs_error_elevation;
        if (fabs(abs_error_u) > error_u_tmax) error_u_tmax = abs_error_u;
        if (fabs(abs_error_v) > error_v_tmax) error_v_tmax = abs_error_v;

        // print errors to screen
        printf("-----------------------------------------------------------\n");
        printf("SW 2D Tidal Errors ----------------------------------------\n");
        printf("Errors at node location: {%10.5e, %10.5e, %10.5e}\n",x,y,z);
        printf("elevation || analytic: %10.5e \t numeric: %10.5e \t error: %10.5e = %10.5f%%\n",
               analytic_elevation,adh_elevation,abs_error_elevation,100*abs_error_elevation/analytic_elevation);
        printf("        u || analytic: %10.5e \t numeric: %10.5e \t error: %10.5e = %10.5f%%\n",
               analytic_u, adh_u,abs_error_u,100*abs_error_u/analytic_u);
        printf("        v || analytic: %10.5e \t numeric: %10.5e \t error: %10.5e = %10.5f%%\n",
               analytic_v, adh_v,abs_error_v,100*abs_error_v/analytic_v);

        
        // print node/station errors to GNUplot file
        fprintf(fp_station,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
                time,adh_elevation,analytic_elevation,adh_u,analytic_u,adh_v,analytic_v,analytic_w,analytic_w);
        fclose(fp_station);
        
        // print node errors to error file
        fprintf(fp_error,"node %d max error @ last time-step :: displacment: %10.7e\n",nodeID+1,abs_error_elevation);
        fprintf(fp_error,"node %d max error @ last time-step :: x-velocity: %10.7e\n",nodeID+1,abs_error_u);
        fprintf(fp_error,"node %d max error @ last time-step :: y-velocity: %10.7e\n",nodeID+1,abs_error_v);
        
        fprintf(fp_error,"\n");
        fprintf(fp_error,"node %d max error over time:: displacment: %10.7e\n",nodeID+1,error_dpl_tmax);
        fprintf(fp_error,"node %d max error over time :: x-velocity: %10.7e\n",nodeID+1,error_u_tmax);
        fprintf(fp_error,"node %d max error over time :: y-velocity: %10.7e\n",nodeID+1,error_v_tmax);
    }
    
    // calculate errors over entire grid (may include wacky boundary corners)
    int surface_node;
    double sum_error_dpl = 0., sum_error_u = 0., sum_error_v = 0., sum_error_w = 0.;
    double rmse_dpl = 0., rmse_u = 0., rmse_v = 0., rmse_w = 0.;
    for (inode=0; inode<mod->grid->nnodes; inode++) {
        x = mod->grid->node[inode].x;
        y = mod->grid->node[inode].y;
        z = mod->grid->node[inode].z;
        depth = mod->sw->d2->head[nodeID];
        adh_elevation = z + mod->sw->d2->head[inode];
        adh_u = mod->sw->d2->vel[inode].x;
        adh_v = mod->sw->d2->vel[inode].y;
        
        xx = (x+y)/sqrt(2.); // convert to orthognal coordinates
        get_tide_analytic(xx, time, depth, 0, &analytic_elevation, &analytic_u, &analytic_v, &analytic_w);
        ut = analytic_u;
        vt = analytic_v;
        analytic_u = (sqrt(2.)/2.) * (ut - vt); // in rotated coordinates
        analytic_v = (sqrt(2.)/2.) * (ut + vt); // in rotated coordinates
        
        sum_error_dpl += pow(adh_elevation - analytic_elevation,2);
        sum_error_u += pow(adh_u - analytic_u,2);
        sum_error_v += pow(adh_v - analytic_v,2);
    }
#ifdef _MESSG
    sum_error_dpl = messg_dsum(sum_error_dpl, mod->grid->smpi->ADH_COMM);
    sum_error_u   = messg_dsum(sum_error_u,   mod->grid->smpi->ADH_COMM);
    sum_error_v   = messg_dsum(sum_error_v,   mod->grid->smpi->ADH_COMM);
#endif
    rmse_dpl = sqrt(sum_error_dpl/(float)mod->grid->macro_nnodes);
    rmse_u = sqrt(sum_error_u/(float)mod->grid->macro_nnodes);
    rmse_v = sqrt(sum_error_v/(float)mod->grid->macro_nnodes);
    
    if (fabs(rmse_dpl) > rmse_error_dpl_tmax) rmse_error_dpl_tmax = rmse_dpl;
    if (fabs(rmse_u) > rmse_error_u_tmax) rmse_error_u_tmax = rmse_u;
    if (fabs(rmse_v) > rmse_error_v_tmax) rmse_error_v_tmax = rmse_v;
    
    // print errors
    if(have_node==1) {
        // print grid errors to screen
        printf("2d grid RMSE max error over time:: absolute d: %10.7e\n",rmse_error_dpl_tmax);
        printf("2d grid RMSE max error over time:: absolute u: %10.7e\n",rmse_error_u_tmax);
        printf("2d grid RMSE max error over time:: absolute v: %10.7e\n",rmse_error_v_tmax);
        printf("-----------------------------------------------------------\n");
        fflush(stdout);
        
        // print grid errors to file
        fprintf(fp_error,"\n");
        fprintf(fp_error,"2d grid RMSE max error over time:: absolute d: %10.7e\n",rmse_error_dpl_tmax);
        fprintf(fp_error,"2d grid RMSE max error over time:: absolute u: %10.7e\n",rmse_error_u_tmax);
        fprintf(fp_error,"2d grid RMSE max error over time:: absolute v: %10.7e\n",rmse_error_v_tmax);
        fclose(fp_error);
    }
    
    if(have_node1==1 && have_node2==1) {
        fclose(fp_surface);
        fclose(fp_bed);
    }
    
}

//****************************************************************
//****************************************************************
void write_testcase_error_tide_3d(SMODEL *mod) {
    
    if (mod->proc_flag != 1) return;
    
    assert(test_case_flag.nodeID_3d >= 0 && test_case_flag.nodeID_3d < mod->grid->macro_nnodes);
    assert(test_case_flag.node1_3d  >= 0 &&  test_case_flag.node1_3d < mod->grid->macro_nnodes);
    assert(test_case_flag.node2_3d  >= 0 &&  test_case_flag.node2_3d < mod->grid->macro_nnodes);
    
    SVECT p1, p2, pointCheck;
    int i, inode, iseg;
    double x,y,z,xx,analytic_elevation,analytic_u,analytic_v,analytic_w,depth=0.;
    double adh_elevation, adh_u, adh_v, adh_w, ut, vt;
    double rmse_dpl = 0., rmse_u = 0., rmse_v = 0., rmse_w = 0.;
    double abs_error_elevation = 0., abs_error_u = 0., abs_error_v = 0., abs_error_w = 0.;
    double g = mod->gravity;
    double time = mod->t_prev + mod->dt; // this routine called before time is updated
    FILE *fp = NULL;
    
    // node to calculate error must be at surface
    int nodeID = UNSET_INT, nd1 = UNSET_INT, nd2 = UNSET_INT;
    
    // make sure only the owner of the test node is calculating errors
    int have_node=0,have_node1=0, have_node2=0;
#ifdef _MESSG
    for(i=0;i<mod->grid->my_nnodes;i++){
        if(test_case_flag.node1_3d==mod->grid->node[i].gid){
            have_node1=1;
            nd1=i;
        }
        if(test_case_flag.node2_3d==mod->grid->node[i].gid){
            have_node2=1;
            nd2=i;
        }
        if(test_case_flag.nodeID_3d==mod->grid->node[i].gid){
            have_node=1;
            nodeID=i;
        }
    }
#else
    nodeID = test_case_flag.nodeID_3d;
    nd1    = test_case_flag.node1_3d;
    nd2    = test_case_flag.node2_3d;
    have_node1=1;
    have_node2=1;
    have_node=1;
#endif
    
    FILE *fp_station=NULL, *fp_error=NULL, *fp_surface=NULL, *fp_bed=NULL;
    if(have_node==1) {
        fp_error = fopen("error3d.out", "w");
        fp_station = fopen("station3d.out", "a"); 
    }
    if(have_node1==1 && have_node2==1) {
        fp_surface = fopen("surfaceXsection3d.out", "a");
        fp_bed = fopen("bedXsection3d.out", "a");
    }
    
    // for cross section plot files
    if((have_node1==1)&&(have_node2==1)){
        p1.x = mod->grid->node[nd1].x; p1.y = mod->grid->node[nd1].y;
        p2.x = mod->grid->node[nd2].x; p2.y = mod->grid->node[nd2].y;
    }
    
    // calculate total depth at test node
    if (have_node==1) {
        x = mod->grid->node[nodeID].x;
        y = mod->grid->node[nodeID].y;
        z = mod->grid->node[nodeID].z;
        
        iseg = find_vertical_segment(mod->grid, nodeID, mod->grid->vertical_hash);
        int test_node_surface = mod->grid->vertical_list[iseg]->id;
        assert(test_node_surface >= 0 && test_node_surface < mod->grid->nnodes);
        
        tl_calculate_depthavgvel(mod->grid, mod->sw->d3);
        depth = mod->sw->d3->depth[mod->grid->nodeID_3d_to_2d_sur[test_node_surface]];
        
        xx = (x+y)/sqrt(2.); // convert to orthognal coordinates
        get_tide_analytic(xx, time, depth, 0, &analytic_elevation, &analytic_u, &analytic_v, &analytic_w);
        ut = analytic_u;
        vt = analytic_v;
        analytic_u = (sqrt(2.)/2.) * (ut - vt); // in rotated coordinates; // convert to rotated
        analytic_v = (sqrt(2.)/2.) * (ut + vt); // in rotated coordinates; // convert to rotated
        
        // calculate AdH elevation
        adh_elevation = mod->grid->node[test_node_surface].z + mod->sw->d3->displacement[test_node_surface];
        adh_u = mod->sw->d3->vel[nodeID].x;
        adh_v = mod->sw->d3->vel[nodeID].y;
        adh_w = mod->sw->d3->vel[nodeID].z;
        
        // absolute errors
        abs_error_elevation = fabs(analytic_elevation - adh_elevation);
        abs_error_u = fabs(analytic_u - adh_u);
        abs_error_v = fabs(analytic_v - adh_v);
        abs_error_w = fabs(analytic_w - adh_w);
        
        // max l-infinity errors over time
        if (fabs(abs_error_elevation) > error_dpl_tmax) error_dpl_tmax = abs_error_elevation;
        if (fabs(abs_error_u) > error_u_tmax) error_u_tmax = abs_error_u;
        if (fabs(abs_error_v) > error_v_tmax) error_v_tmax = abs_error_v;
        if (fabs(abs_error_w) > error_w_tmax) error_w_tmax = abs_error_w;

        // print errors to screen
        printf("-----------------------------------------------------------\n");
        printf("SW 3D Tidal Errors ----------------------------------------\n");
        printf("Errors at node location: {%10.5e, %10.5e, %10.5e}\n",x,y,z);
        printf("elevation || analytic: %10.5e \t numeric: %10.5e \t error: %10.5e = %10.5f%%\n",
               analytic_elevation,adh_elevation,abs_error_elevation,100*abs_error_elevation/analytic_elevation);
        printf("        u || analytic: %10.5e \t numeric: %10.5e \t error: %10.5e = %10.5f%%\n",
               analytic_u, adh_u,abs_error_u,100*abs_error_u/analytic_u);
        printf("        v || analytic: %10.5e \t numeric: %10.5e \t error: %10.5e = %10.5f%%\n",
               analytic_v, adh_v,abs_error_v,100*abs_error_v/analytic_v);
        printf("        w || analytic: %10.5e \t numeric: %10.5e \t error: %10.5e = %10.5f%%\n",
               analytic_w, adh_w,abs_error_w,100*abs_error_w/analytic_w);
        
        
        // print node/station GNUplot file
        fprintf(fp_station,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
                time,adh_elevation,analytic_elevation,adh_u,analytic_u,adh_v,analytic_v,analytic_w,analytic_w);
        fclose(fp_station);
        
        // print node errors to error file :: replace each call
        fprintf(fp_error,"node %d max error @ last time-step :: displacment: %10.7e\n",nodeID+1,abs_error_elevation);
        fprintf(fp_error,"node %d max error @ last time-step :: x-velocity: %10.7e\n",nodeID+1,abs_error_u);
        fprintf(fp_error,"node %d max error @ last time-step :: y-velocity: %10.7e\n",nodeID+1,abs_error_v);
        fprintf(fp_error,"node %d max error @ last time-step :: w-velocity: %10.7e\n",nodeID+1,abs_error_w);
        
        fprintf(fp_error,"\n");
        fprintf(fp_error,"node %d max error over time:: displacment: %10.7e\n",nodeID+1,error_dpl_tmax);
        fprintf(fp_error,"node %d max error over time :: x-velocity: %10.7e\n",nodeID+1,error_u_tmax);
        fprintf(fp_error,"node %d max error over time :: y-velocity: %10.7e\n",nodeID+1,error_v_tmax);
        fprintf(fp_error,"node %d max error over time :: w-velocity: %10.7e\n",nodeID+1,error_w_tmax);
    }
    
    // calculate errors over entire grid (may include wacky boundary corners)
    int surface_node, bed_node;
    double sum_error_dpl = 0., sum_error_u = 0., sum_error_v = 0., sum_error_w = 0.;
    
    tl_calculate_depthavgvel(mod->grid, mod->sw->d3);
    
    for (inode=0; inode<mod->grid->nnodes; inode++) {
        x = mod->grid->node[inode].x;
        y = mod->grid->node[inode].y;
        z = mod->grid->node[inode].z;
        
        iseg = find_vertical_segment(mod->grid, inode, mod->grid->vertical_hash);
        surface_node = mod->grid->vertical_list[iseg]->id; // this time, the 3D node ID of surface node
        bed_node = (mod->grid->vertical_list[iseg]->prev)->id;
        
        assert(surface_node >= 0 && surface_node < mod->grid->nnodes);
        assert(bed_node >= 0 && bed_node < mod->grid->nnodes);
        
        depth = mod->sw->d3->depth[mod->grid->nodeID_3d_to_2d_sur[mod->grid->vertical_list[iseg]->id]];
        adh_elevation = mod->grid->node[surface_node].z + mod->sw->d3->displacement[surface_node];
        adh_u = mod->sw->d3->vel[inode].x;
        adh_v = mod->sw->d3->vel[inode].y;
        adh_w = mod->sw->d3->vel[inode].z;
        
        xx = (x+y)/sqrt(2.); // convert to orthognal coordinates
        get_tide_analytic(xx, time, depth, 0, &analytic_elevation, &analytic_u, &analytic_v, &analytic_w);
        ut = analytic_u;
        vt = analytic_v;
        analytic_u = (sqrt(2.)/2.) * (ut - vt); // in rotated coordinates
        analytic_v = (sqrt(2.)/2.) * (ut + vt); // in rotated coordinates
        
        sum_error_dpl += pow(adh_elevation - analytic_elevation,2);
        sum_error_u += pow(mod->sw->d3->vel[inode].x - analytic_u,2);
        sum_error_v += pow(mod->sw->d3->vel[inode].y - analytic_v,2);
        sum_error_w += pow(mod->sw->d3->vel[inode].z - analytic_w,2);
        
        
        if((have_node1==1)&&(have_node2==1)){
            pointCheck.x = x; pointCheck.y = y;
            if (isPointOnLine(p1,p2,pointCheck) == TRUE)  {
                if (inode == surface_node) {
                    //printf("node: %d x: %20.10f y: %20.10f z: %20.10f\n",inode+1,mod->grid->node[inode].x,mod->grid->node[inode].y,mod->grid->node[inode].z);
                    fprintf(fp_surface,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
                            time,x,y,analytic_elevation,adh_elevation,analytic_u,adh_u,analytic_v,adh_v,analytic_w,adh_w);
                }
                if (inode == bed_node) {
                    //printf("node: %d x: %20.10f y: %20.10f z: %20.10f\n",inode+1,mod->grid->node[inode].x,mod->grid->node[inode].y,mod->grid->node[inode].z);
                    fprintf(fp_bed,   "%20.10e %20.10e%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
                            time,x,y,analytic_elevation,adh_elevation,analytic_u,adh_u,analytic_v,adh_v,analytic_w,adh_w);
                }
            }
        }
        
    }
    
#ifdef _MESSG
    sum_error_dpl = messg_dsum(sum_error_dpl,mod->grid->smpi->ADH_COMM);
    sum_error_u   = messg_dsum(sum_error_u,  mod->grid->smpi->ADH_COMM);
    sum_error_v   = messg_dsum(sum_error_v,  mod->grid->smpi->ADH_COMM);
    sum_error_w   = messg_dsum(sum_error_w,  mod->grid->smpi->ADH_COMM);
#endif
    
    rmse_dpl = sqrt(sum_error_dpl/(float)mod->grid->macro_nnodes);
    rmse_u = sqrt(sum_error_u/(float)mod->grid->macro_nnodes);
    rmse_v = sqrt(sum_error_v/(float)mod->grid->macro_nnodes);
    rmse_w = sqrt(sum_error_w/(float)mod->grid->macro_nnodes);
    
    if(have_node1==1 && have_node2==1) {
        fprintf(fp_surface,"\n\n");
        fprintf(fp_bed,"\n\n");
        fclose(fp_surface); fclose(fp_bed);
    }
    
    if (fabs(rmse_dpl) > rmse_error_dpl_tmax) rmse_error_dpl_tmax = rmse_dpl;
    if (fabs(rmse_u) > rmse_error_u_tmax) rmse_error_u_tmax = rmse_u;
    if (fabs(rmse_v) > rmse_error_v_tmax) rmse_error_v_tmax = rmse_v;
    if (fabs(rmse_w) > rmse_error_w_tmax) rmse_error_w_tmax = rmse_w;
    
    // print grid errors
    if(have_node==1) {
        // to screen
        printf("3d grid RMSE max error over time:: absolute d: %10.7e\n",rmse_error_dpl_tmax);
        printf("3d grid RMSE max error over time:: absolute u: %10.7e\n",rmse_error_u_tmax);
        printf("3d grid RMSE max error over time:: absolute v: %10.7e\n",rmse_error_v_tmax);
        printf("3d grid RMSE max error over time:: absolute w: %10.7e\n",rmse_error_w_tmax);
        printf("-----------------------------------------------------------\n");
        fflush(stdout);
        
        // to file
        fprintf(fp_error,"\n");
        fprintf(fp_error,"3d grid RMSE max error over time:: absolute d: %10.7e\n",rmse_error_dpl_tmax);
        fprintf(fp_error,"3d grid RMSE max error over time:: absolute u: %10.7e\n",rmse_error_u_tmax);
        fprintf(fp_error,"3d grid RMSE max error over time:: absolute v: %10.7e\n",rmse_error_v_tmax);
        fprintf(fp_error,"3d grid RMSE max error over time:: absolute w: %10.7e\n",rmse_error_w_tmax);
        fclose(fp_error);
    }
}

//********************************************************************************/
//********************************************************************************/
// DIFFUSIVE WAVE VERIFICATION
//****************************************************************//
//****************************************************************//

//****************************************************************//
//****************************************************************//
double get_hunter_depth(double c, double u, double n, double x, double t) {
    // return the depth of flow
    // c - constant of integration
    // n - mannings coefficient
    // u - velocity magnitude
    // t - time
    //
    return ( pow((7.0/3.0) * (c - n*n*u*u*u*(x - u*t)),(3./7.)) );
}

void testcase_prep_dwe_hunter(SMODEL *mod) {
}
void write_testcase_error_dwe_hunter(SMODEL *mod) {
    int inode;
    double t = mod->t_prev + mod->dt;
    double c = 0.0;
    double n = 0.01;
    double u = 1;
    
    double sum_error_h = 0., h = 0., h_adh = 0., x = 0, max_h = 0., max_h_adh = 0.;
    for (inode=0; inode<mod->grid->nnodes; inode++) {
        x = mod->grid->node[inode].x;
        h = 0.;
        if (x < u * t) h = get_hunter_depth(c, u, n, x, t);
        h_adh = mod->sw->d2->head[inode];
        
        if (h > max_h) max_h = h;
        if (h_adh > max_h_adh) max_h_adh = h_adh;
        
        sum_error_h += pow(h_adh - h,2);
    }
    sum_error_h = sqrt(sum_error_h/(float)mod->grid->nnodes);
    if (sum_error_h > rmse_error_dpl_tmax) rmse_error_dpl_tmax = sum_error_h;
    printf("analytic max h: %-20.10e \t adh max h: %-20.10e \t mse @ %-20.10e\n",max_h,max_h_adh,sum_error_h);
    
    // grid mass error
    //double grid_mass_error = tl_find_grid_mass_error_elem2d(mod->density, mod->sw->d2->head, mod->sw->d2->vel, mod->grid, mod->flag, mod->initial_grid_mass, mod->series_head, mod->str_values, mod->dt, &total_time_mass_flux);
    
    
    FILE *fp;
    fp = fopen("error.out", "w");
    fprintf(fp,"mse: %-12.10e\n",rmse_error_dpl_tmax);
    //fprintf(fp,"grid mass error: %-12.10e\n",grid_mass_error);
    fclose(fp);
    
}


