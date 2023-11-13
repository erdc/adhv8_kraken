/* This is the include file for */
/* SedLib */
/* V1.0 */
/* macros */
#define MAX(a,b) (a > b ? a : b)
#define MIN(a,b) (a < b ? a : b)

/**************************************************************************/
/* PROTOTYPES *************************************************************/

void ssedlib_alloc(SSEDLIB **, int, int, int, int);
void ssedlib_init(SSEDLIB *);
void ssedlib_free(SSEDLIB *);
void ssedlib_wave_init(SSEDLIB_WAVE *);
void ssedlib_bed_layer_init(SSEDLIB_BED_LAYER *, int);
void ssedlib_grain_init(SSEDLIB_GRAIN *, int);
void ssedlib_printScreen(SSEDLIB);
void ssedlib_grain_printScreen(SSEDLIB_GRAIN);
void ssedlib_bed_layer_printScreen(SSEDLIB_BED_LAYER);
void ssedlib_wave_printScreen(SSEDLIB_WAVE);

void check_sedlib_revision(void);
void sedlib_core_active_layer_properties(SSEDLIB* , double);
void sedlib_core_active_strata_source (SSEDLIB*, double, double, double, double);
void sedlib_core_bed(SSEDLIB *);
void sedlib_core_bed_layer_update(SSEDLIB *);
void sedlib_core_bed_shear_stress(SSEDLIB *);
void sedlib_core_bed_sorting_and_displacement(SSEDLIB*, double, int);
void sedlib_core_bedload_transport(SSEDLIB *, int);
void sedlib_core_coh_deposition_flux(SSEDLIB *);
void sedlib_core_critical_shear_stress(SSEDLIB *);
void sedlib_core_limit_erosion(SSEDLIB *);
void sedlib_core_mixed_sed_erosion_flux(SSEDLIB *);
void sedlib_core_noncoh_erosion_deposition_flux(SSEDLIB *);
void sedlib_core_suspended_transport(SSEDLIB *, int);
double sedlib_proc_active_layer_thickness(SSEDLIB *);
double sedlib_proc_find_percent_finer(SSEDLIB *, double, int);
double sedlib_proc_noncoh_fall_velocity_cheng(double, double, double, double);
double sedlib_proc_coh_fall_velocity_hwang_and_mehta(double, double, double*);
double sedlib_proc_sl_near_bed_correction_brown(int, double, double, double, double, double, double);
double sedlib_proc_bedshear_current_only_christensen(double, double);
double sedlib_proc_bedshear_under_ice_savant_brown(double, double, double, double);
double sedlib_proc_critical_shear_stress_shields_van_rijn(double, double, double, double, double);
double sedlib_proc_critical_shear_correction(double, double, double, double);
void sedlib_proc_bed_load_direction_correction(double, double, double, double, double, double, double, double, double*, double*);
double sedlib_proc_noncoh_entrainment_rate_garcia_parker(double, double, double, double, double, double, double, double);
double sedlib_proc_noncoh_entrainment_rate_wright_parker(double, double, double, double, double, double, double, double, double, double, double, double);

double sedlib_proc_noncoh_entrainment_rate_van_rijn(
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double
                             );

double sedlib_proc_noncoh_entrainment_rate_wave(
							double,
							double,
							double,
							double,
							double,
							double,
							double,
							double,
							double
);

double sedlib_proc_bed_load_flux_meyer_peter_mueller(
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double
                             );

double sedlib_proc_bed_load_flux_wong_parker(
                             double,
                             double,
                             double,
                             double,
                             double
                             );

double sedlib_proc_bed_load_flux_wave(
							double,
							double,
							double,
							double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double *,
                            double *
);

double sedlib_proc_hiding_factor_karim_holly_yang(
                             double,
                             double
                             );

double sedlib_proc_hiding_factor_egiazaroff(
                             double,
                             double
                             );

double sedlib_proc_hiding_factor_wu_wang_jia(
                             int,
                             SSEDLIB*
                             );

double sedlib_proc_transport_mode_allocation_van_rijn(
                             double,
                             double,
                             double,
                             double
                             );

double sedlib_proc_bed_load_flux_van_rijn(
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double
                             );

double sedlib_proc_bed_load_adaption_length_jain(
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double
                             );

double sedlib_proc_bedload_velocity_phillips_and_sutherland(
                             double,
                             double,
                             double
                             );

double sedlib_proc_coh_deposition_krone(
                             double,
                             double,
                             double
                             );

double sedlib_proc_coh_erosion_alishahi_and_krone(
                             double,
                             double,
                             double,
                             double
                             );

double sedlib_proc_sedflume_erosion(
                             SSEDLIB*
                             );

void sedlib_proc_bed_layer_consolidate(
                             int,
                             double,
                             double,
                             double,
                             double**,
                             double*,
                             double*,
                             double*,
                             double*,
                             double*
                             );

void sedlib_proc_remove_negative_concentration(
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double*,
                             double*,
                             double*
                             );


void sedlib_proc_transport_correction_factors(
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double*,
                             double*,
                             double*
                             );

double sedlib_proc_wave_current_shear_teeter(
                             double,
                             double,
                             double,
                             double,
                             double
                             );

double sedlib_proc_wave_current_shear_grant_madsen(
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
                             double,
			     double
                             );

