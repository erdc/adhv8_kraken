#ifndef H_SSUPER_MODEL_
#define H_SSUPER_MODEL_


/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {

    int nsubmodels;
    SMODEL *submodel;

    // gkc below
    int NumInterfaces;           // Total number of coupled interfaces in the supermodel
    SINTERFACE *interface;       // Length to be set to NumInterfaces
    //gkc above

    int **fmap;                  // forward mapping from submodel nodes to supermodel matrix
    int *fmap_nodes;

    /* WVEL coupling */
    int **fmap_wvel;
    int *fmap_wvel_nodes;
    int save_wvel_firstcall;

    int nsys;
    int nsys_sq;
    int max_nsys;
    int max_nsys_sq;


    int macro_nnodes_original;   // Initial size of the global matrix composed of all submodels
    int macro_nnodes;            // Size of the global matrix composed of all submodels, possibly after grid adaption
    int my_nnodes;               // Sum of # of residential nodes of all submodels on myid pe. ?? maybe (residential - interface considerations)?
    int nnodes;                  // Matrix size on myid pe (residential + ghost - interface considerations)
    int nnodes_matrix;           // Previous matrix size on myid pe
    int total_macro_nnodes;            // Blind sum of all submodel->grid->nnodes       // Used for (re)allocating matrix memory
    int total_nnodes;            // Blind sum of all submodel->grid->nnodes       // Used for (re)allocating matrix memory
    int total_my_nnodes;         // Blind sum of all submodel->grid->my_nnodes    // Wondering what this could be used for.
    int total_nnodes_matrix;     // Blind sum of all submodel->grid_nnodes_matrix // Used for (re)allocating matrix memory
    int wvel_nnodes;             // WVEL Matrix size on myid pe (residential + ghost - interface considerations)
    int wvel_my_nnodes;          // Resident WVEL Matrix size on myid pe (residential - interface considerations)
    int wvel_macro_nnodes;       // Resident WVEL Matrix size, Allgathered (residential - interface considerations)

    int nonlinear_it_total;
    int nonlinear_it_total_hvel;
    int nonlinear_it_total_wvel;
    int max_nonlin_it;
    int nblock;
    int nalloc_inc;
    double inc_nonlin;
    double tol_nonlin;
    double dt;
    double dt_prev;

    /* Transport */
    int ntransport;
    // int itrns;
    int *con_type;

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
    SPARSE_VECT *matrix;
#endif
    Solver_Info solver_info; //SAM - Maybe put this in non-PETSc part too

    SMPI *supersmpi;      /* supermodel SMPI used only for updating HVEL solver variables. */
    SMPI *supersmpi_wvel; /* supermodel SMPI used only for updating WVEL solver variables. */

    int **proc_map; /* maps supermodel procesor id to submodel ids */ 
    int **proc_map_wvel;
    int proc_flag; /* determines if a procesor should work on a supermodel */

} SSUPER_MODEL;

/*********************************************************/
/*********************************************************/
/* struct methods -------------------------------------- */

void ssuperModel_alloc_init(SSUPER_MODEL **superModel_ptr, int nsupermodels, int *nsubmodels, int *nInterfaces);
void ssuperModel_free(SSUPER_MODEL **superModel_ptr, SSUPER_INTERFACE **superInterface_ptr, int nsupermodels, int nsuperinterfaces);
void ssuperModel_printScreen(SSUPER_MODEL *sm, SSUPER_INTERFACE *si, int nsupermodels, int nsuperinterfaces);
void ssuperModel_readSuperfile(SSUPER_MODEL **superModel_ptr, SSUPER_INTERFACE **superInterface_ptr, int *nsupmods, int *nsupifcs);
void ssuperModel_printFmap(SSUPER_MODEL *sm);
void ssuperModel_create_interface_data(SSUPER_MODEL *);
void ssuperModel_forward_map(SSUPER_MODEL *);
void ssuperModel_forward_map_wvel(SSUPER_MODEL *);

void generate_2d_2d_interface(SSUPER_MODEL *, SINTERFACE *, int, int, int);
void generate_3d_3d_interface(SSUPER_MODEL *, SINTERFACE *, int, int, int);
void generate_2d_3d_interface(SSUPER_MODEL *, SINTERFACE *, int, int, int);
void fix_interface_ownership(SSUPER_MODEL *, int);

#endif
