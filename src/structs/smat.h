#ifndef H_SMAT_H_
#define H_SMAT_H_

/*******************************************************************/
/*******************************************************************/

typedef struct {
  int EEVF;                     /* flag that EEV is being used */
  int EVSF;                     /* flfat that EVS is being used */
  double coriolis;
  int max_lev;                  /* the maximum number of refinement levels in the material */
  STENSOR ev;                    /* the eddy viscosities for a 2d shallow water model */
  double refine_tolerance;      /* the refinement tolerances */
  double unrefine_tolerance;    /* the unrefinement tolerances */
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
  int supression_func;      /* Suppression function */
  int wall_func;      /* turbulence wall function */
  double min_tke;
  double min_tds;
  double len_max;   /* Maximum mixing length */
  double hyd_conductivity;            /* Green and Ampt Infiltration */
  double psi;
  double rooting_depth;
} SMAT_SW;

/*******************************************************************/
/*******************************************************************/

typedef struct {
  int EEVF;                     /* flag that EEV is being used */
  int EVSF;                     /* flfat that EVS is being used */
  double coriolis;
  int max_lev;                  /* the maximum number of refinement levels in the material */
  STENSOR3D ev;                    /* the eddy viscosities for a 2d shallow water model */
  double refine_tolerance;      /* the refinement tolerances */
  double unrefine_tolerance;    /* the unrefinement tolerances */
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
  int supression_func;      /* Suppression function */
  int wall_func;      /* turbulence wall function */
  double min_tke;
  double min_tds;
  double len_max;   /* Maximum mixing length */
  double hyd_conductivity;            /* Green and Ampt Infiltration */
  double psi;
  double rooting_depth;
} SMAT_NS;


/*******************************************************************/
/*******************************************************************/

typedef struct {
  int DIFF_FLAG;                /* flag the consituent diffusion is used */
  int max_lev;                  /* the maximum number of refinement levels in the material */
  double d_m;                   /* the turbulent diffusion coefficient */
  double source;                /* the strength of the volumetric source */
  double d_l, d_t;              /* the longitudinal and transverse dispersion coefficients */
  double *react;                /* linear reaction coefficient */
  double rd;                    /* retardation coefficient */
  double tortuosity;            /* tortuosity */
  double refine_tolerance;      /* the refinement tolerance */

} SMAT_TRN;

/*******************************************************************/
typedef struct
{
  double alpha0;
  double alpha1;
  double alpha2;
  double alpha3;
  double alpha4;
  double alpha5;
  double alpha6;
  double beta1;
  double beta2;
  double beta3;
  double beta4;
  double k1;
  double k2;
  double k3;
  double k4;
  double mu;
  double rho;
  double sigma1;
  double sigma2;
  double sigma3;
  double sigma4;
  double sigma5;
  double kl;
  double kn;
  double kp;
  double pn;
  double lambda0;
  double lambda1;
  double lambda2;
  double min_do_conc;
} SMAT_NSM;

/*******************************************************************/
#ifdef _ADH_GROUNDWATER
typedef struct {
  int max_lev;                  /* the maximum number of refinement levels in the material */
  double refine_tolerance;      /* the refinement tolerances */
  double unrefine_tolerance;    /* the unrefinement tolerances */
  STENSOR k;                     /* the hydraulic conductivities */
  double s_s;                   /* the storage coefficient */
  double porosity;              /* porosity */
  double water_vol;             /* the volume of water of a volumetric source */
  double residual_sat;          /* residual saturation */
  double vangen_alpha;          /* alpha from vangenuchten equations */
  double vangen_max_cp;         /* vangenuchten maximum capillary pressure head */
  double vangen_n;              /* vangenuchten exponent (assume mualem) */
  int vangen_num_xy;            /* vangenuchten number of entries in the xy series */
  double brooks_lambda;         /* brooks-corey pore size index exponent */
  double brooks_pd;             /* brooks-corey air entry pressure head */
  double brooks_max_cp;         /* brooks-corey maximum capillary pressure head */
  int brooks_num_xy;            /* brooks-corey number of entries in the xy series */
  int ikr;                      /* index to the pressure - relative conductivity series */
  int isat;                     /* index to the pressure - saturation curve */
  double tortuosity;            /* tortuosity (only used for diffusion transport) */
  double d_l;                   /* longitudinal dispersivity */
  double d_t;                   /* transverse dispersivity */
  double ss_area;               /* specific surface area */
  double bulk_density;          /* bulk density */
  int itype;                    /* ? */

} SMAT_GW;
#endif
/*******************************************************************/

/*******************************************************************/
/*******************************************************************/

typedef struct {
        SMAT_SW *sw;        // only ever allocate one of these
        SMAT_NS *ns;
        SMAT_TRN *trn;      // allocate ntransport of these (no grains here)
        SMAT_TRN *sed;      // allocate nsed of these
        SMAT_NSM *wnsm;     // NSM parameters
#ifdef _ADH_GROUNDWATER
        SMAT_GW *gw;
#endif   
} SMAT;

/*******************************************************************/
/*******************************************************************/
#endif 
