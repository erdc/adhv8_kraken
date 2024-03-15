/* This is the main routine for stand-alone adh */
#include "global_header.h"

static time_t time1, time2, time_start;        /* for run-time calculation */

int petsc_test(int argc,char **argv);

/*******************************************************/
/*******************************************************/
/*******************************************************/

int main(int argc, char *argv[]) {
    
    int i,j;
    
    
    /* Check the arguments */
    if (argc < 2 || argc > 5) {
        fprintf(stderr, "\n Missing argument.\n Usage:\n" "   To run a simulation:\n     adh file_base -s (optional)\n" "   Or,\n   For version information:\n      adh -v\n");
        exit(0);
    }
    else if (strcmp(argv[1], "-v") == AGREE) {
        print_build_info();
        printf("\n");
        exit(0);
    }
    
    int input_test = OFF;
    double run_time = -3;
    int nsupermodels = 0;
    int nsuperinterfaces = 0;
    SSUPER_MODEL *superModel = NULL;
    SSUPER_INTERFACE *superinterface = NULL;
    
    /* start calculation time */
    time(&time_start);
    
    /*******************************************************/
    /*******************************************************/
    int ierr = 0;
    
    // AdH initialization
    time(&time1);
#ifdef _MESSG
    MPI_Comm comm_world = MPI_COMM_WORLD;
    ierr = adh_init_func_(&superModel, &superinterface, &nsupermodels, &nsuperinterfaces, &argc, &argv, &comm_world);
#else
    ierr = adh_init_func_(&superModel, &superinterface, &nsupermodels, &nsuperinterfaces, &argc, &argv);
#endif
    time(&time2);
    double time_initial = difftime(time2, time1);
    
    assert(&(superModel[0]) != NULL);
    
    // AdH execution
    time(&time1);
    ierr = adh_run_func_(superModel, superinterface, nsupermodels, nsuperinterfaces, &run_time);
    time(&time2);
    double time_run = difftime(time2, time1);
    
    int ierr_code, myid=0;
#ifdef _MESSG
    ierr_code = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
    if (myid <= 0) {
        printf("\n");
#ifdef _DEBUG
        printf("**********************************************************************\n");
        printf("Timings: \n");
        printf("-- Total initialization time: %lf seconds\n", time_initial);
        for (i=0; i<nsupermodels; i++) {
            for (j=0; j<superModel[i].nsubmodels; j++) {
                if (superModel[i].submodel[j].flag.SW2_FLOW) {
                    printf("SW 2D TIMINGS\n");
                    printf("---- HYD RESID runtime: %lf seconds\n", TIME_IN_2D_SW_RESID);
                    printf("---- HYD LOAD runtime: %lf seconds\n", TIME_IN_2D_SW_LOAD);
                    printf("------ HYD 2D ELEM RESID runtime: %lf seconds\n", TIME_IN_2D_SW_BODY_RESID);
                    printf("------ HYD 1D ELEM RESID runtime: %lf seconds\n", TIME_IN_2D_SW_BOUNDARY_RESID);
                } else if (superModel[i].submodel[j].flag.SW3_FLOW) {
                    printf("SW 3D TIMINGS\n");
                    printf("---- HVEL RESID runtime: %lf seconds\n", TIME_IN_HVEL_RESID);
                    printf("---- HVEL LOAD runtime: %lf seconds\n", TIME_IN_HVEL_LOAD);
                    printf("------ HVEL 3D ELEM RESID runtime: %lf seconds\n", TIME_IN_HVEL_BODY_RESID);
                    printf("------ HVEL 2D ELEM RESID runtime: %lf seconds\n", TIME_IN_HVEL_BOUNDARY_RESID);
                    printf("---- WVEL RESID runtime: %lf seconds\n", TIME_IN_WVEL_RESID);
                    printf("---- WVEL LOAD runtime: %lf seconds\n", TIME_IN_WVEL_LOAD);
                    printf("------ WVEL 3D ELEM RESID runtime: %lf seconds\n", TIME_IN_WVEL_BODY_RESID);
                    printf("------ WVEL 2D ELEM RESID runtime: %lf seconds\n", TIME_IN_WVEL_BOUNDARY_RESID);
                } else if (superModel[i].submodel[j].flag.GW_FLOW) {
                    printf("---- GW RESID runtime: %lf seconds\n", TIME_IN_GW_RESID);
                    printf("---- GW LOAD runtime: %lf seconds\n", TIME_IN_GW_LOAD);
                    printf("------ GW 3D ELEM RESID runtime: %lf seconds\n", TIME_IN_GW_BODY_RESID);
                    printf("------ GW 2D ELEM RESID runtime: %lf seconds\n", TIME_IN_GW_BOUNDARY_RESID);
                    
                }
            }
        }
#endif
    }
    
    // AdH finalization
    time(&time1);
    ierr = adh_finalize_func_(&superModel, &superinterface, nsupermodels, nsuperinterfaces);
    time(&time2);
    double time_final = difftime(time2, time1);
    
    /*******************************************************/
    /*******************************************************/
    /* getting calculation time */
    time(&time2);
    double t2 = difftime(time2, time_start);
    
    if (myid <= 0) {
#ifdef _DEBUG
        printf("-- Total execution time: %lf seconds\n", time_run);
        printf("-- Total finalization time: %lf seconds\n", time_final);
#endif
        printf("Total simulation runtime is %lf seconds\n", t2);
    }

#ifdef _PETSC
    ierr = PetscFinalize();
#endif
#ifdef _MESSG
MPI_Finalize();
#endif
return ierr;
}

/*******************************************************/
/*******************************************************/
/*******************************************************/
