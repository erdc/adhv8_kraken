#ifndef H_SSEDIMENT_
#define H_SSEDIMENT_

// dependencies :: SVECT, SVECT2D

/* active and bed are all layers .. */


/* allocate an array of SEDIMENT types */
/* sed[inode].depth_avg_flag */
/* sed[inode].layer[ilayer].sedflume_flag */

/* SEDIMENT ALGORITHM CONSTANTS FOR READING PARAMETERS FROM BC FILE */
#include "sedlib_header.h"

#define COHSET 4
#define WINDWAVE 10
#define NONCOENT 1 /* (cjt) -- mdw SW3_FLOW svn merge */
#define DIVCOEF 3


/***********************************************************************/
/***********************************************************************/
/* struct for suspended load */
typedef struct {
    
    double *c;          // concentration
    double *old_c;      // old concentration
    double *older_c;    // older concentration
    double *source;     // sediment source
    double *sink;       // sediment sink
    double *rouse_coef; // rouse profile coefficient
    double *mfcf;       // mass flux correction factor
    double *error;      // suspended load nodal error
    SVECT2D *vcf;       // velocity correction factor
    SVECT2D *vor_vel;   // the radial near-bed velocity due to vorticity
    
} SSUSLOAD;

/***********************************************************************/
/***********************************************************************/
/* struct for bed load */
typedef struct {
    
    double *c;           // concentration
    double *old_c;       // old concentration
    double *older_c;     // older concentration
    double *thick;       // thickness
    double *old_thick;   // old thickness
    double *older_thick; // older thickness
    double *source;      // sedlib source
    double *sink;        // sediment sink
    double *error;       // bed load node error
    double *shear;       // bed load shear stress
    SVECT2D *v;          // bedload velocity vector
    SVECT2D *flux;       // bedload flux vector
    
} SBEDLOAD;

/***********************************************************************/
/***********************************************************************/
typedef struct {

    double *shear_stress;
    double *erosion_rate;

} SSEDFLUME;

/***********************************************************************/
/***********************************************************************/
typedef struct {

	int type;
    int cbed_flag;
	int sedflume_flag;
	SSEDFLUME *sedflume;
    double *bulk_density;
    double *porosity;
	double *critical_erosion_shear;
	double *erosion_rate_constant;
	double *erosion_rate_exponent;
	double *thickness;
	double **distribution;		/* for each bed layer, of size [ngrains][nnodes] */

} SBED_LAYER;

/***********************************************************************/
/***********************************************************************/
typedef struct {

    int sed_to_clay;
    int clay_to_sed;
    double settling_velocity;
    double tau_ce;
    double tau_cd;
    double erode_const;

} SCLAY;

/***********************************************************************/
/***********************************************************************/
typedef struct {

    int sed_to_snd;
    int snd_to_sed;
    double settling_velocity;
    double tau_ce;
    double tau_cd;
    double erode_const;

} SSAND;

/***********************************************************************/
/***********************************************************************/
typedef struct {

	int type;
    double reference_c;
	double specific_gravity;
	double diameter;
	double porosity;
	SCLAY clay;
    SSAND sand;

} SGRAIN;

/***********************************************************************/
/***********************************************************************/
typedef struct {

    /* flags */
    int wind_wave_flag;         /* flag for wind wave algorithm  */
    int cohesive_settling_flag; /* flag for cohesive settling algorithm  */
    int noncoh_ent_flag;        /* flag for noncohesive sediment entrainment algorithm */
    int noncoh_ble_flag;        /* flag for noncohesive bedload entrainment */
    int hid_fact_flag;          /* flag for selection of the hiding factor */
    int bed_type_flag;          /* ?? */

    /* parameters (int) */
    int nnodes_bed;                 /* the number of nodes for bedload transport and sedlib arrays */
    int nnodes_sus;                 /* the number of nodes for suspended transport */
    int nse_params;                 /* The number of suspended entrainment process parameters */
    int nbe_params;                 /* The number of bedload flux process parameters */
    int nsand;                      /* number of sand constituents being transported */
    int nclay;                      /* number of clay constituents being transported */
    int nsilt;                      /* number of silt constituents being transported */
    int nsed;                       /* the total of nsand, nsilt, nclay */
    int nlayers;                    /* the number of sediment bed layers allowed */
    int nsfssi;                     /* the number of Sedflume shear stress intervals allowed */
    int ncti;                       /* the number of consolidation time intervals allowed */
    
    /* sediment materials */
    SMAT_TRN *mat;

    /* parameters (double) */
    double *cohesive_settling_const;/* constants for the cohesive settling algorithm */
    double *wind_wave_const;        /* constants for the wind wave algorithm */
    double *noncoh_ent_const;       /* constants for the sediment entrainment algorithm */
    double *noncoh_ble_const;       /* constants for the sediment bedload algorithm (cjt) */
    double *div_coef;               /* constant to define the sediment diversion specifications */
    
    int *call_flag;
    int *no_displacement_flag;
    int *friction_flag;
    int *cbed_flag;               /* flag that tells if cohesive bed properties are assigned */
    double *friction_coef;
    double *as_ceiling, *old_as_ceiling;      /* active stratum ceiling */
    double *bed_displacement, *old_bed_displacement;
    double *bed_shear_stress;
    SVECT2D *bedload_vector, *susload_vector;
    SVECT2D *bed_gradient;
    
    SGRAIN *grain; /* grain properties (array of size nsed = nsand + nclay + nsilt */
    SBEDLOAD *bedload;
    SSUSLOAD *susload;
    SBED_LAYER *active_layer;
    SBED_LAYER *old_active_layer;
    SBED_LAYER *bed_layer;      /* there are some number of bed layers at each node */
    SBED_LAYER *old_bed_layer;
    double ***consolidation_bed_properties;
    
    // The sedlib struct
    SSEDLIB *sedlib;


} SSED;

/*********************************************************/
/* struct methods -------------------------------------- */

#ifdef _SEDIMENT
void ssusload_alloc(SSUSLOAD **, int, int, int);
void ssusload_realloc(SSUSLOAD *, int, int, int, int, int);
void ssusload_init_range(SSUSLOAD *, int, int, int, int, int);
void ssusload_node_avg(SSUSLOAD *, int, int, int, int);
void ssusload_renumber(SSUSLOAD *, int, int, int *, int *, double *, SVECT2D *);
void ssusload_free(SSUSLOAD *, int, int, int);
void sbedload_alloc(SBEDLOAD **, int, int);
void sbedload_realloc(SBEDLOAD *, int, int, int);
void sbedload_init_range(SBEDLOAD *, int, int, int);
void sbedload_node_avg(SBEDLOAD *, int, int, int, int);
void sbedload_renumber(SBEDLOAD *, int, int, int *, int *, double *, SVECT2D *);
void sbedload_free(SBEDLOAD *, int, int);
void sbed_layer_alloc(SBED_LAYER **, int, int, int);
void sbed_layer_realloc(SBED_LAYER *, int, int, int, int);
void sbed_layer_init_range(SBED_LAYER *, int, int, int, int);
void sbed_layer_node_avg(SBED_LAYER *, int, int, int, int, int, double, double, int);
void sbed_layer_renumber(SBED_LAYER *, int, int, int, int *, int *, double *);
void sbed_layer_free(SBED_LAYER *, int, int, int);
void sbed_layer_copy(SBED_LAYER *, SBED_LAYER *, int, int, int);
void sgrain_init(SGRAIN *);
void sgrain_copy(SGRAIN *, SGRAIN *);
void sclay_init(SCLAY *);
void ssand_init(SSAND *);
void ssedflume_alloc_init(SSEDFLUME **, int);
void ssedflume_realloc_init(SSEDFLUME *, int, int);
void ssedflume_init_range(SSEDFLUME *, int, int);
void ssedflume_node_avg(SSEDFLUME *, int, int, int);
void ssedflume_renumber(SSEDFLUME *, int, int *, int *, double *);
void ssedflume_free(SSEDFLUME *, int);
void ssediment_alloc_init(int, int, int, int, int, int, int, int, int, int, SSED **);
void ssediment_realloc_init(SSED *, int, int, int, int);
void ssediment_realloc_init_sus(SSED *, int, int);
void ssediment_realloc_init_bed(SSED *, int, int);
void ssediment_init_range_bed(SSED *, int, int);
void ssediment_init_range_sus(SSED *, int, int);
void ssediment_renumber(SSED *, int, int *, int *, int *, double *, SVECT2D *);
void ssediment_renumber_3d_sus(SSED *, int, int *, int *, int *, double *, SVECT2D *);
void ssediment_renumber_3d_bed(SSED *, int, int *, int *, int *, double *, SVECT2D *);
void ssediment_free(SSED *);
//void ssediment_node_avg(SSED *, int, int, int, SGRID *);
//void ssediment_prep(SMODEL *, SSED *, int);
//void ssediment_open_output(SMODEL *);
//void ssediment_print_ts(SGRID *, SIO *, SSED *, double, double, int, SFILE, char *, int, int, int **, int, int *, int **, int, int *, int);
//void ssediment_print_adapted_ts(SGRID *, SSED *, double, double, int, SFILE, char *, int, int);
//void ssediment_read_hot(SMODEL *);
//void ssediment_calculate_2d_error(SSED *, SSW_2D *, SGRID *, SMAT *, double);
//void ssediment_calculate_3d_error(SSED *, SSW_3D *, SGRID *, SMAT *, double);
#endif

#endif

