#ifndef H_SMAT_GW
#define H_SMAT_GW

typedef struct {
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

// methods
void smat_gw_alloc_init(SMAT_GW **mat_gw);
void smat_gw_free(SMAT_GW *mat_gw);

#endif 
