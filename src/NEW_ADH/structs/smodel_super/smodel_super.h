// An AdH SuperModel
#ifndef H_SMODEL_SUPER_
#define H_SMODEL_SUPER_



/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {

    SIO *io;  /* file names associated with this model application */
    SFLAGS flag;
    int o_flag;
    SGRID *grid; // just a pointer to the grid in the design model

    //Mark added type flag, isSimple -> superModel has only one physics type on whole grid
    int isSimple;

    //just pointers to dm
    double *dt, *old_dt, *dt_err;
    double* dt_prev;
    double inc_nonlin;
    double tol_nonlin;
    double* t_init;
    double* t_prev;
    double* t_final;
    int* t_adpt_flag;
    double tau_temporal; //can change from super model to super model
    double gravity; //idk where this should be? probably design, this can just be pointer
    double density;
    int vorticity_id;

    //an additional integer if we want to subtimestep off the main dt
    //nsubstep of 1 would be dt = dt of super model
    int nsubsteps;

    int nseries;              // the number of series in this model 
    int itrns;

    //Mark adding other things that used to be part of solver
    int max_nonlin_linesearch_cuts;
    int max_nonlin_it;
    int it_count_nonlin;
    int force_nonlin_it;
    int nonlinear_it_total;
    int LINEAR_PROBLEM;
    int it_count_nonlin_failed;

    //Linear solver
    SLIN_SYS *lin_sys; //pointer to the design model's linear system
    double *sol; //solution variable, stored in each super model [ndofs]
    double *sol_old;
    double *sol_older;


    int *dof_map_local;    // a local map from the local node ID to the
                           // local equation number for building the
                           // FE residual and matrix [nnodes] long

    int meshcode; //need a way if supermodel is defined on entire mesh (0), surface(1), or floor(2)


    //Mark adding sizes for convenience
    int nphysics_mat_1d;
    int nphysics_mat_2d;
    int nphysics_mat_3d;

    SMAT_PHYSICS *elem1d_physics_mat; //[nphysics_mat_1d]
    SMAT_PHYSICS *elem2d_physics_mat; //[nphysics_mat_2d]
    SMAT_PHYSICS *elem3d_physics_mat; //[nphysics_mat_3d]

    //Mark, elements need integer to store the physics mat id
    int *elem1d_physics_mat_id; //[nelem1d]
    int *elem2d_physics_mat_id; //[nelem2d]
    int *elem3d_physics_mat_id; //[nelem3d]

    //Only necessary for CG:
    //Mark proposes swapping nodal vars above to node-based material, this will cut down on memory
    //but may be challenging to form. This won't be set by user but implicitly built
    //at run time
    //int nphysics_mat_node;
    //int *node_physics_mat_id; //[nnode] ? local vs what idk
    //SMAT_PHYSICS *node_physics_mat; // [nphysics_mat_node]
    //Mark, NOTE: CHANGE TO THIS
    SMAT_PHYSICS **node_physics_mat; //[nnode] array of pointers to SMAT_PHYSICS on elements

    //Mark a pointer to array of function pointers, using model_codes.h
    //void *forward_step;
    int forward_step;
    
    //void (*printFunc)(struct MyStruct *); // Function pointer
    //Mark, has moved to SMAT_PHYSICS
    //this should be based on number of physics materials, no longer numbr of elements
    //int *nSubMods1d;                // [nphysics_mat_1d] the total number of physics modules on each 1D element
    //int *nSubMods2d;                // [nphysics_mat_2d] the total number of physics modules on each 2D element
    //int *nSubMods3d;                // [nphysics_mat_3d] the total number of physics modules on each 3D element
    
    //Now part of SMAT_PHYSICS
    //SMODEL **elem1d_physics;  // [nphysics_mat_1d][nsubmods_1d] the fe routines for each type of physics on each 1D element
    //SMODEL **elem2d_physics;  // [nphysics_mat_2d][nsubmods_2d] the fe routines for each type of physics on each 2D element
    //SMODEL **elem3d_physics;  // [nphysics_mat_3d][nsubmods_3d] the fe routines for each type of physics on each 3D element

    //Mark added local to process, this should be same as local_range[1]-local_range[0]
    //are these redundant in any way?
    int *my_ndofs; //pointers to design model, not arrays
    int *my_ndofs_old;
    int *ndofs; // local number of degrees of freedom
    int *ndofs_old; //local numer of solution variables the processor is in charge of
    int *macro_ndofs;
    int *macro_ndofs_old;






    /* boundary conditions mask */
    //maybe we can get rid of this through weak enforcement
    int *bc_mask;
    //instead of mask, allocate array of ints that is dof of each dirichlet dof
    //int *dirchlet_dofs;
    double *dirichlet_data;
    //Mark, havent looked into this yet. Need to iron out
    STR_VALUE *str_values;    /* strings */
    SSERIES *series_head;
    SSERIES *series_curr;
    SSERIES *series_dt;       /* time-step series */
    SSERIES *series_out;      /* output series */
    SSERIES *series_wind_head, *series_wind_curr;
    SSERIES *series_wave_head, *series_wave_curr;
//    //stripping out models and putting info in here
//    // SOLUTION VARIABLES //
//    // DERIVATIVE QUANTITIES DO NOT GO HERE //
//    //any physics specific quantities will be mapped from sol vector
//    //any time we want a new physical variable it will need to be added here
//    //this is listed in order of how dofs will be sorted, water depth
//    //do we do this or separate SW2/SW3/...
//    int nhead;
//    double *head;    // present in SW, GW, NS, ...
//    double *old_head;    /* pressure from the previous time step at time t_{n-1} */
//    double *older_head;  /* pressure from the time step before last t_{n-2} */
//    //3d sw
//    int ndisplacement;
//    double *displacement;
//    double *old_displacement;
//    double *older_displacement;
//    //navier stokes
//    double *prs;
//    double *old_prs;
//    double *older_prs;
//    //2d velocity
//    //present in SW2D, NS2D
//    int nvel2d;
//    SVECT2D *vel2d;
//    SVECT2D *old_vel2d;
//    SVECT2D *older_vel2d;
//    //3d velocity
//    int nvel3d;
//    SVECT *vel3d;
//    SVECT *old_vel3d;
//    SVECT *older_vel3d;
//    //concentrations for transport (do we need more than one for sediment too?)
//    int *nconcentration;
//    double *concentration;
//    double *old_concentration;
//    double *older_concentration;
//    //sediment solution variable, maybe find better name than c
//    int nc;
//    double *c;          // concentration
//    double *old_c;      // old concentration
//    double *older_c;    // older concentration
//    // END OF SOLUTION VARIABLES
//    //structures to contain data from other models OTHER than solution variables
//    //maybe just add scon and ssed
//#ifdef _SEDIMENT
//    SSED *sed;
//#endif
    //potentially surround these with different ifdefs
//    SSW_2D *ssw2d;
//    SSW_3D *ssw3d;
//    SGW *sgw;
//    SCON *con;
//    SNS_2D *sns2d;
//    SNS_3D *sns3d;
    SSW *sw; //pointer or actual structure?
} SMODEL_SUPER;



/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void smodel_super_alloc_init(SMODEL_SUPER *sm);
void smodel_super_alloc_init_array(SMODEL_SUPER **smod, int nSuperModels);


void smodel_super_free_array(SMODEL_SUPER *sm, int nSuper);
void smodel_super_free(SMODEL_SUPER *sm);

void smodel_super_read(SMODEL_SUPER *smod, FILE *fp);
void smodel_super_printScreen(SMODEL_SUPER *smod);

//Mark added for testing
void smodel_super_no_read_simple(SMODEL_SUPER *sm, double* dt_in, double* t_init, double* t_prev,
    double* t_final, int nphysics_mat_1d, int nphysics_mat_2d, int nphysics_mat_3d, 
    char elemVarCode[4], int isSimple, SGRID *grid, SLIN_SYS *sys);
//from fe_newton.c
//int fe_newton(struct SMODEL_SUPER* sm);

//TODOD, write this up!
int smodel_super_forward_step(SMODEL_SUPER* sm, int (*ts_fnctn)(SMODEL_SUPER*));

int smodel_super_resid(SMODEL_SUPER* sm, double *rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG, int (*fe_resid)(SMODEL_SUPER *, double *, int, double, int, int, int, int));

void smodel_super_update_dirichlet_data(SMODEL_SUPER *sm);
void smodel_super_prep_sol(SMODEL_SUPER *sm);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif

