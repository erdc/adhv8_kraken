#ifndef H_SMAT_SW_
#define H_SMAT_SW_

typedef struct {
  int EEVF;                     /* flag that EEV is being used */
  int EVSF;                     /* flfat that EVS is being used */
  double coriolis;
  STENSOR ev;                    /* the eddy viscosities for a 2d shallow water model */
  int bed_disp_flag;            /* flag inidicating Y/N for sediment bed displacement */
  int vor_flag;                 /* flag indicating Y/N for vorticity calculations */
  double eev_coef;              /* coefficent used in viscosity calculation if eev_flag = ON */
  double fraction;              /* Used with FRC flag GSAVANT */
  int eev_mode;                 /* EEV Computation Flag */
  int d_flag;                   /* INdicate whether density coupling is off */
  int wind_flag;                /* flag for wind variable calculations (cjt) */
  double windatt;               /* wind attenuation */
  double smag_coeff;            /* turbulent smagorinski coefficient */
  int turbulence_model_xy;      /* lateral turbulence model */
  int turbulence_model_z;       /* vertical turbulence model */
  int supression_func;          /* Suppression function */
  int wall_func;                /* turbulence wall function */
  double min_tke;
  double min_tds;
  double len_max;               /* Maximum mixing length */
  double hyd_conductivity;      /* Green and Ampt Infiltration */
  double psi;
  double rooting_depth;
} SMAT_SW;

// Methods
void smat_sw_alloc_init(SMAT_SW **mat_sw);
void smat_sw_free(SMAT_SW *mat);

#endif
