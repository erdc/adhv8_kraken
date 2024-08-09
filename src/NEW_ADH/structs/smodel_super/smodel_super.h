// An AdH SuperModel
#ifndef H_SMODEL_SUPER_
#define H_SMODEL_SUPER_

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
typedef struct {

    SIO *io;  /* file names associated with this model application */
    SFLAGS flag;
    int o_flag;
    SGRID *grid; // just a copy of the grid 

    double dt, old_dt, dt_err;
    double dt_prev;
    double inc_nonlin;
    double tol_nonlin;
    double t_init;
    double t_prev;
    double t_final;
    int t_adpt_flag;

    int nseries;              // the number of series in this model 
    int itrns;


    //FE_MATRIX *matrix;  // stores matrix
#ifdef _PETSC
        Mat         A;
        Mat         P;                      // Preallocator matrix for A
        int         PREALLOC_FLAG;          // Determines whether or not we preallocate space for the matrix
        int         Istart;                 // Start of PETSc matrix rows on pe
        int         Iend;                   // End of PETSc matrix rows on pe
        const int   *ownership_range;       // Gives the range of rows owned by each pe
        Vec         residual;
        Vec         sol;
        KSP         ksp;
        int         *bc_mask;
        int         old_bc_mask_size;
        double      *scale_vect; // Can I remove this too? - SAM
        // TODO: Remove two members below - SAM
        double *diagonal;
        SPARSE_VECT *matrix;
#else
        int    *bc_mask;
        double *residual;
        double *sol;
        double *scale_vect;
        double *diagonal;
        //SPARSE_VECT *matrix;
#endif


    int *physics_mat_code; // an code for each physics material that
                           // determinds which equations will be
                           // solved on that material
    int *dof_map_local;    // a local map from the local node ID to the
                           // local equation number for building the
                           // FE residual and matrix 
    int *dof_map_global;   // for ghost nodes - a map from the local
                           // node ID to the global equation number
                           // for building the FE matrix 

    int meshcode; //need a way if supermodel is defined on entire mesh (0), surface(1), or floor(2)



    int *nSubMods1d;                // [nelems1d] the total number of physics modules on each 1D element
    int *nSubMods2d;                // [nelems2d] the total number of physics modules on each 2D element
    int *nSubMods3d;                // [nelems3d] the total number of physics modules on each 3D element
    SELEM_PHYSICS **elem1d_physics;  // [nelems1d][nsubmods_1d] the fe routines for each type of physics on each 1D element
    SELEM_PHYSICS **elem2d_physics;  // [nelems2d][nsubmods_2d] the fe routines for each type of physics on each 2D element
    SELEM_PHYSICS **elem3d_physics;  // [nelems3d][nsubmods_3d] the fe routines for each type of physics on each 3D element

    
    int ndofs; // local number of degrees of freedom
    int dofs_old; //local numer of solution variables the processor is in charge of

    int **elem_nvars;
    int ***elem_vars;
    int *node_nvars;
    int **node_vars;

    /* boundary conditions mask */
    //bcmask *bc_mask;

    STR_VALUE *str_values;    /* strings */
    SSERIES *series_head;
    SSERIES *series_curr;
    SSERIES *series_dt;       /* time-step series */
    SSERIES *series_out;      /* output series */
    SSERIES *series_wind_head, *series_wind_curr;
    SSERIES *series_wave_head, *series_wave_curr;

    //stripping out models and putting info in here

    // SOLUTION VARIABLES //

    // DERIVATIVE QUANTITIES DO NOT GO HERE //

    //any physics specific quantities will be mapped from sol vector
    //any time we want a new physical variable it will need to be added here

    //this is listed in order of how dofs will be sorted, water depth
    int nhead;
    double *head;    // present in SW, GW, NS, ...
    double *old_head;    /* pressure from the previous time step at time t_{n-1} */
    double *older_head;  /* pressure from the time step before last t_{n-2} */

    //3d sw
    int ndisplacement;
    double *displacement;
    double *old_displacement;
    double *older_displacement;

    //navier stokes
    double *prs;
    double *old_prs;
    double *older_prs;

    //2d velocity
    //present in SW2D, NS2D
    int nvel2d;
    SVECT2D *vel2d;
    SVECT2D *old_vel2d;
    SVECT2D *older_vel2d;

    //3d velocity
    int nvel3d;
    SVECT *vel3d;
    SVECT *old_vel3d;
    SVECT *older_vel3d;

    //concentrations for transport (do we need more than one for sediment too?)
    int *nconcentration;
    double *concentration;
    double *old_concentration;
    double *older_concentration;


    //sediment solution variable, maybe find better name than c
    int nc;
    double *c;          // concentration
    double *old_c;      // old concentration
    double *older_c;    // older concentration



    // END OF SOLUTION VARIABLES



    //structures to contain data from other models OTHER than solution variables
    //maybe just add scon and ssed
#ifdef _SEDIMENT
    SSED *sed;
#endif

    //potentially surround these with different ifdefs
//    SSW_2D *ssw2d;
//    SSW_3D *ssw3d;
//    SGW *sgw;
//    SCON *con;
//    SNS_2D *sns2d;
//    SNS_3D *sns3d;

    



} SMODEL_SUPER;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void smodel_super_alloc_init(SMODEL_SUPER **smod, int nSuperModels);
void smodel_super_free(SMODEL_SUPER *smod, int nSuperModels);
void smodel_super_read(SMODEL_SUPER *smod, FILE *fp);
void smodel_super_printScreen(SMODEL_SUPER *smod);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif

