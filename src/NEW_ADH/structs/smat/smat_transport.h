#ifndef H_SMAT_TRANSPORT_
#define H_SMAT_TRANSPORT_

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

// methods
void smat_trn_alloc_init(SMAT_TRN **mat_trn);
void smat_trn_free(SMAT_TRN *mat);
#endif 
