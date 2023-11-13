#ifndef H_SSEDLIB_
#define H_SSEDLIB_

/************************************************************/
/************************************************************/
typedef struct {
    int ibreak;
    double height;
    double period;
    double angle;
    double number;
    double speed;
    double edb;
    double edf;
    double sef;
    double seb;
    double wcf;
    double beta;
    double gama;
} SSEDLIB_WAVE;
/************************************************************/
/************************************************************/
typedef struct {
    int type;
    int sedflume_flag;
    int fully_eroded_layer;

    double por;
    double ces;
    double erc;
    double ere;
    double thick;

} SSEDLIB_BED_LAYER;

/************************************************************/
/************************************************************/
/* allocate n_sed + 1 */
typedef struct {

    int type;
    double bed_mass_change;
    double net_mass_change;
    double diameter;
    double settling_vel;
    double coh_free_settling_vel;
    double por, ces, cds, erc, ere;
    double critical_shear_stress;
    double critical_shear_stress_plane_bed;
    double sg;
    double al_dist;
    double as_dist;
    double bl_con;
    double bl_thick;
    double sl_con;
    double sl_bcon;
    double sl_sg;
    double bl_erosion_flux;
    double bl_deposition_flux;
    double sl_erosion_flux;
    double sl_deposition_flux;
    double sl_erosion_flux_nbl;
    double bl_erosion_flux_nbl;
    double bl_deposition_flux_nbl;
    double bl_net_flux;
    double sl_net_flux;
    double sl_nbf;
    double old_al_dist;


} SSEDLIB_GRAIN;
/************************************************************/
/************************************************************/
typedef struct {
    int WAVE_SED;       /* flag set to TRUE when waves are implemented within ADH */
    int n_blay;
    int n_sed;
    int n_sp1;
    int n_sedfl;
    int n_conti;
    int roughness_type;
    int oald_resolve_flag;
    int layer_being_eroded;
    int consolidate_flag;
    int coh_settling_flag;
    int wave_current_flag;
    int noncoh_sedent_flag;
    int noncoh_bedload_flag;
    int hiding_factor_flag;
    int depth_averaged_flag;

    double al_por, old_al_por, as_por;
    double al_ces, old_al_ces, as_ces;
    double al_erc, old_al_erc, as_erc;
    double al_ere, old_al_ere, as_ere;
    double al_thick;
    double as_ceiling, old_as_ceiling;
    double disp, old_disp;
    double depth_old, depth_new;
    double kin_visc;
    double grav;
    double rho;
    double ice_roughness;
    double roughness_coefficient;
    double total_erodable_thickness;
    double minimum_erodable_thickness;
    double strata_dep_coef;
    double total_shear_stress_mag;
    double grain_shear_stress_mag;
    double near_bed_velocity_mag;
    double time_step;
    double ero_limit_ts;
    double sl_mfcf;
    double esrh, tsrh;

    SVECT2D grain_shear_stress;
    SVECT2D v_da;
    SVECT2D v_nb;
    SVECT2D v_wind;
    SVECT2D v_vort;
    SVECT2D bl_vel;
    SVECT2D bed_slope;
    SVECT2D sl_vcf;

    double *noncoh_sedent_constants; // ??
    double *coh_settling_constants; //[4];
    double *al_coh_prop; // [5];
    double **sedfl_shear;
    double **sedfl_erate;
    double **consolidate_table;
    double **bl_dist;

    SSEDLIB_GRAIN *grain;
    SSEDLIB_BED_LAYER *layer;
    SSEDLIB_WAVE wave;
} SSEDLIB;


#endif
