#ifndef H_SCON_
#define H_SCON_

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {

    int type;                           /* type of transported constituent */
    double *property;                   /* properties like particle diameter, etc. */
    double *concentration;
    double *old_concentration;
    double *older_concentration;
    double *sink;
    double *source;
    double *nodal_decay_coefficient;    /* multiplies the concentration */
    double *error;
    double *mfcf;
    SVECT2D *vcf;

} SCON;

/*********************************************************/
/* struct methods -------------------------------------- */

void scon_alloc_init(int, int, int, SCON **);
void scon_realloc_init(SCON *, int, int, int);
void scon_init(SCON *, int, int, int);
void scon_realloc_init_nodes(SCON *, int, int, int);
void scon_node_avg(SCON *, int, int, int, int);
void scon_renumber(SCON *, int, int, int *, int *, double *, SVECT2D *);
void scon_free(int, int, int, SCON *);
//void scon_print_ts(int, SCON *, SIO *, SGRID *, double, int, SFLAGS, SFILE, char *, int, int, int **, int, int *, int);
//void scon_print_adapted_ts(int, double, SCON *, SGRID *, double, SFLAGS, SFILE, char *, int, int);
//void scon_open_output(SMODEL *);
//void scon_calculate_elem2d_error(int, SCON *, SSW_2D *, SGRID *, SMAT *, double);
//void scon_calculate_elem3d_error(int, SCON *, SSW_3D *, SGRID *, SMAT *, double);

/***********************************************************/
/***********************************************************/
/***********************************************************/


#endif
