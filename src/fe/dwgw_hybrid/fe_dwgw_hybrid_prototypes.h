#ifndef H_FE_DWGW_HYBRID_PROTOTYPES_
#define H_FE_DWGW_HYBRID_PROTOTYPES_

/***********************************************************/
/***********************************************************/
/***********************************************************/

void fe_dwgw_hybrid_inc(SSUPER_MODEL *sm, int imod);
void fe_dwgw_hybrid_init(SSUPER_MODEL *sm, int imod);
void fe_dwgw_hybrid_load(SSUPER_MODEL *sm, int imod);
void fe_dwgw_hybrid_resid(SSUPER_MODEL *sm, int imod);
void fe_dwgw_hybrid_update(SSUPER_MODEL *sm, int imod);
int fe_dwgw_hybrid_solve(SSUPER_MODEL *sm, int imod);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
