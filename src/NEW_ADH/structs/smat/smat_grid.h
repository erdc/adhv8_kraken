#ifndef H_SMAT_GRID_
#define H_SMAT_GRID_

typedef struct {
  int max_lev;                  /* the maximum number of refinement levels in the material */
  double refine_tolerance;      /* the refinement tolerances */
  double unrefine_tolerance;    /* the unrefinement tolerances */
} SMAT_GRID;

// Methods
void smat_grid_alloc_init(SMAT_GRID **mat_grid);
void smat_grid_free(SMAT_GRID *mat_grid);

#endif
