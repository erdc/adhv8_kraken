#ifndef H_SNS_3D_
#define H_SNS_3D_

// dependencies :: SVECT, SVECT2D

/***********************************************************/
/***********************************************************/
/***********************************************************/

// NOTES:
//  - since surface and bed displacements change the jacobians, it's easier to keep them in 3d for elem transfer

typedef struct {
    
    /* 3d grid solution arrays */
    //migrated up to super model
    //SVECT *vel;
    //SVECT *old_vel;
    //SVECT *older_vel;
    //double *prs;
    //double *old_prs;
    //double *older_prs;



    double *displacement;
    double *old_displacement;
    double *older_displacement;
    double *grid_speed;
    double *old_grid_speed;
    double *density;
    double *dpl_perturbation;

    double *error;
    double *vertical_node_flux;
    double *hyd_viscosity;  /* GSAVANT */
    double *trn_diffusivity; /* GSAVANT */
    
    double *darray;
    int *iarray;
    

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
    
} SNS_3D;

/*********************************************************/
/* struct methods -------------------------------------- */

void sns_3d_realloc_init(SNS_3D *, int, int);
void sns_3d_realloc_init_surface(SNS_3D *, int, int);
void sns_3d_realloc_init_bed(SNS_3D *, int, int);
void sns_3d_renumber(SNS_3D *, int, int *, int *, double *, SVECT2D *, SVECT *);
void sns_3d_renumber_surface(SNS_3D *, int, int *, int *, double *, SVECT2D *);
void sns_3d_checkall(SNS_3D *, int, int);
void sns_3d_realloc_init(SNS_3D *, int , int );
void sns_3d_realloc_init_surface(SNS_3D *, int , int );
void sns_3d_renumber(SNS_3D *, int , int *, int *, double *, SVECT2D *, SVECT *);
void sns_3d_renumber_surface(SNS_3D *, int , int *, int *, double *, SVECT2D *);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
