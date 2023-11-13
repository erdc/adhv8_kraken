#include "global_header.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main(int argc, char *argv[]) {
    
    int ierror = 0;
    
#ifdef _MPI
    ierror = MPI_Init(&argc, &argv);
#endif
    
    /* Check the arguments */
    if (argc < 2 || argc > 3) {
        fprintf(stderr, "\n Missing argument.\n Usage:\n" "   To run a simulation:\n     pt file_base -s (optional)\n" "   or,\n   to test:\n      pt -t\n");
        exit(0);
    } else if (strcmp(argv[1], "-t") == 0) {
        pt_test_searchEngine(); // test :: geometry and search engine
        pt_verification();      // test :: verification
        exit(-1);
    }
    
    RUN_VERBOSE = true;
    
    // allocate the model
    SMODEL *mod = (SMODEL *) malloc(sizeof(SMODEL));
    
    // initialize project base filename
    strcpy(mod->file_base,argv[1]);
    
    // intialize test cases
    mod->get_analytic_p = NULL;
    mod->get_analytic_v = NULL;
    mod->get_analytic_p = get_analytic_p_2d_uniform_oscillating;
    mod->get_analytic_v = get_vel_2d_uniform_oscillating;
    
    // initialize the model
    ierror = smodel_init(mod);
    
    // run the model
    ierror = smodel_run_lockstep(mod);
    
    // finalize the model
    ierror = smodel_finalize(mod);

#ifdef _MPI
    ierror = MPI_Finalize();
#endif
    
    // finalize the model
    return ierror;
}

