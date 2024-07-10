#ifndef H_SSW_2D_
#define H_SSW_2D_

// dependencies :: SVECT2D

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {

    /* solution arrays */
    double *head;
    double *old_head;
    double *older_head;  
    SVECT2D *vel;
    SVECT2D *old_vel;
    SVECT2D *older_vel;     
    double *density;
    double *bed_displacement;   /* bed displacement from sediment */ 
    double *bed_elevation;      /* bed elevation (z + dpl) */
    double *dacontResid;        // store the nodal da-continuity residual for 3d model coupling

    /* meteoroligic surface stresses */
    SWIND *winds;
    SWAVE *waves;
    
    /* nodal errors */
    double *error;              /* nodal errors for shallow water continuity only */

    /* general use array to keep from reallocating after adaption */
    double *darray;
    int *iarray;
    void *vwork;
    int vwork_size;

    /* store these to add to other equations */
    double **elem_rhs_dacont_extra_terms;
    int elem_rhs_realloc;

} SSW_2D;

/*********************************************************/
/* struct methods -------------------------------------- */

void ssw_2d_init(SSW_2D *, int, int);
void ssw_2d_alloc_init(SSW_2D **, SGRID *, SIO *, SFLAGS);
void ssw_2d_realloc_init(SSW_2D *, int, int);
void ssw_2d_node_avg(SSW_2D *, int, int, int);
void ssw_2d_renumber(SSW_2D *, int, int *, int *, double *, SVECT2D *);
void ssw_2d_checkall(SSW_2D, int);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
