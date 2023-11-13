#ifndef H_FE_SW_HYBRID_PROTOTYPES_
#define H_FE_SW_HYBRID_PROTOTYPES_

/***********************************************************/
/***********************************************************/
/***********************************************************/

void fe_sw_hybrid_inc(SSUPER_MODEL *sm, int imod);
void fe_sw_hybrid_init(SSUPER_MODEL *sm, int imod);
void fe_sw_hybrid_load(SSUPER_MODEL *sm, int imod);
void fe_sw_hybrid_resid(SSUPER_MODEL *sm, int imod);
void fe_sw_hybrid_update(SSUPER_MODEL *sm, int imod);
int fe_sw_hybrid_solve(SSUPER_MODEL *sm, int imod);

void fe_sw_hybrid_wvel_inc(SSUPER_MODEL *sm, int imod);
void fe_sw_hybrid_wvel_init(SSUPER_MODEL *sm, int imod);
void fe_sw_hybrid_wvel_load(SSUPER_MODEL *sm, int imod);
void fe_sw_hybrid_wvel_resid(SSUPER_MODEL *sm, int imod);
void fe_sw_hybrid_wvel_update(SSUPER_MODEL *sm, int imod);


/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
