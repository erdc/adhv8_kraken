#ifndef H_SNS_2D_
#define H_SNS_2D_

// dependencies :: SVECT2D

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {

    /* solution arrays */
    //migrated up to super model
    //pressure? same or different than head?
    //double *head;
    //double *old_head;
    //double *older_head;  
    //SVECT2D *vel;
    //SVECT2D *old_vel;
    //SVECT2D *older_vel;     
    
    
    double *density;
    double *bed_displacement;   /* bed displacement from sediment */ 
    double *bed_elevation;      /* bed elevation (z + dpl) */

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

} SNS_2D;

/*********************************************************/
/* struct methods -------------------------------------- */

void sns_2d_init(SNS_2D *, int, int);
void sns_2d_alloc_init(SNS_2D **, SGRID *, SIO *, SFLAGS);
void sns_2d_realloc_init(SNS_2D *, int, int);
void sns_2d_node_avg(SNS_2D *, int, int, int);
void sns_2d_renumber(SNS_2D *, int, int *, int *, double *, SVECT2D *);
void sns_2d_checkall(SNS_2D, int);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
