#ifndef H_SSW_3D_
#define H_SSW_3D_

// dependencies :: SVECT, SVECT2D

/***********************************************************/
/***********************************************************/
/***********************************************************/

// NOTES:
//  - since surface and bed displacements change the jacobians, it's easier to keep them in 3d for elem transfer

typedef struct {
    
    /* 3d grid solution arrays */
    double *displacement;
    double *old_displacement;
    double *older_displacement;
    double *grid_speed;
    double *old_grid_speed;
    double *density;
    double *dpl_perturbation;
    double *prs;
    double *prs_plus;
    double *prs_minus;
    double *error;
    double *vertical_node_flux;
    double *hyd_viscosity;  /* GSAVANT */
    double *trn_diffusivity; /* GSAVANT */
    SVECT2D *grad_bed;        // bed gradient
    double *darray;
    int *iarray;
    
    SVECT *vel;
    SVECT *old_vel;
    SVECT *older_vel;
    SVECT2D *tanvec;
#ifdef _SEDIMENT
    double *bed_displacement;      // displacement from sediment transport
    double *old_bed_displacement;  // displacement from sediment transport
    double *older_bed_displacement;// displacement from sediment transport
#endif
    
    // cjt :: used to store SUPG rhs for later use (reset every dt - no adaption needed, just reallocation)
    double **elem_rhs_supg_dacont;
    double **elem_rhs_supg_cont;
    int elem_rhs_realloc;
    
    
    /* 2d grid solution arrays */
    double *depth;
    double *old_depth;
    double *bed_elevation;          // node.z + bed displacement from sediment
    SVECT2D *depth_avg_vel;
    SVECT2D *old_depth_avg_vel;
    SVECT2D *surface_vel;
    SVECT2D *bottom_vel;
    SWIND *winds;
    SWAVE *waves;
    
    /* general use array to keep from reallocating after adaption */
    void *vwork;
    int vwork_size;
    
} SSW_3D;

/*********************************************************/
/* struct methods -------------------------------------- */

//void ssw_3d_alloc_init(SSW_3D **, SGRID *, SIO *, SFLAGS);
void ssw_3d_realloc_init(SSW_3D *, int, int);
void ssw_3d_realloc_init_surface(SSW_3D *, int, int);
void ssw_3d_realloc_init_bed(SSW_3D *, int, int);
void ssw_3d_renumber(SSW_3D *, int, int *, int *, double *, SVECT2D *, SVECT *);
void ssw_3d_renumber_surface(SSW_3D *, int, int *, int *, double *, SVECT2D *);
//void ssw_3d_init(SSW_3D *, SGRID *);
//void ssw_3d_free(SSW_3D *, SGRID *, SFLAGS);
void ssw_3d_checkall(SSW_3D *, int, int);
//void ssw_3d_print_ts(SSW_3D *, SIO *, SGRID *, double, int, SFILE, char *, int, int, int **, int, int *, int **, int, int *, int);
//void ssw_3d_open_output(SMODEL *);
//void ssw_3d_open_input(SMODEL *);
//void ssw_3d_calculate_elem_error(SSW_3D *, SGRID *, SMAT *, double);
void ssw_3d_realloc_init(SSW_3D *, int , int );
void ssw_3d_realloc_init_surface(SSW_3D *, int , int );
void ssw_3d_renumber(SSW_3D *, int , int *, int *, double *, SVECT2D *, SVECT *);
void ssw_3d_renumber_surface(SSW_3D *, int , int *, int *, double *, SVECT2D *);
//void ssw_3d_node_avg(SSW_3D *, int, int, int, SGRID *);
//void ssw_3d_node_avg_sur(SSW_3D *, int, int, int, SGRID *);
//void ssw_3d_node_avg_bed(SSW_3D *, int, int, int, SGRID *);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
