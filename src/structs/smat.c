#include "global_header.h"

/**************************************************************/
/**************************************************************/

void smat_init(SMODEL *mod, SMAT **mat, int nmat) {
    
    (*mat) = (SMAT *) tl_alloc(sizeof(SMAT), nmat);
    SMAT *mt = *mat; // alias
    
    int imat, itrns;
    for (imat=0; imat<nmat; imat++) {
        mt[imat].sw = NULL;
        mt[imat].trn = NULL;
        mt[imat].sed = NULL;
        
#ifdef _ADH_GROUNDWATER
        mt[imat].gw = NULL;
#endif
        if (mod->flag.SW_FLOW) {
            mt[imat].sw = (SMAT_SW *) tl_alloc(sizeof(SMAT_SW), 1);
            smat_sw_init(mod, mt[imat].sw);
        }
        
        if (mod->flag.NS_FLOW) {
            mt[imat].ns = (SMAT_NS *) tl_alloc(sizeof(SMAT_NS), 1);
            smat_ns_init(mod, mt[imat].ns);
        }
        
#ifdef _ADH_GROUNDWATER
        if (mod->flag.GW_FLOW) {
            mt[imat].gw = (SMAT_GW *) tl_alloc(sizeof(SMAT_GW), 1);
            smat_gw_init(mod, mt[imat].gw);
        }
#endif

        if (mod->flag.TRANSPORT) {
            mt[imat].trn = (SMAT_TRN *) tl_alloc(sizeof(SMAT_TRN), mod->ntransport);
            for (itrns=0; itrns<mod->ntransport; itrns++) {
                smat_trn_init(mod, &mt[imat].trn[itrns]);
            }
        }
        
        if (mod->flag.SEDIMENT) {
            mt[imat].sed = (SMAT_TRN *) tl_alloc(sizeof(SMAT_TRN), mod->nsed);
            for (itrns=0; itrns<mod->nsed; itrns++) {
                smat_trn_init(mod, &mt[imat].sed[itrns]);
            }
        }
        
        if (mod->flag.NSM) {
            mt[imat].wnsm = (SMAT_NSM *) tl_alloc(sizeof(SMAT_NSM), 1);
            smat_wnsm_init(mod, mt[imat].wnsm);
        }
    }
    
}

/**************************************************************/
/**************************************************************/

void smat_free(SMODEL *mod, SMAT *mat, int nmat) {
    int imat;
    
    if (mat != NULL) {
        for (imat=0; imat<nmat; imat++) {
            if (mod->flag.SW_FLOW) {
                mat[imat].sw = (SMAT_SW *) tl_free(sizeof(SMAT_SW), 1, mat[imat].sw);
            }
            
            if (mod->flag.NS_FLOW) {
                mat[imat].ns = (SMAT_NS *) tl_free(sizeof(SMAT_NS), 1, mat[imat].ns);
            }
            
#ifdef _ADH_GROUNDWATER
            if (mod->flag.GW_FLOW) {
                mat[imat].gw = (SMAT_GW *) tl_free(sizeof(SMAT_GW), 1, mat[imat].gw);
            }
#endif

            if (mod->flag.TRANSPORT) {
                mat[imat].trn = (SMAT_TRN *) tl_free(sizeof(SMAT_TRN), mod->ntransport, mat[imat].trn);
            }
            
            if (mod->flag.SEDIMENT) {
                mat[imat].sed = (SMAT_TRN *) tl_free(sizeof(SMAT_TRN), mod->nsed, mat[imat].sed);
            }
        }
        mat = (SMAT *) tl_free(sizeof(SMAT), nmat, mat);
    }
}

/**************************************************************/
/**************************************************************/

void smat_sw_check(SMODEL *mod, SMAT_SW mat, int imat) {
    
    if ( (mat.EEVF == ON) && (mat.EVSF == ON) ) {
        printf("ERROR: MATERIAL #: %d\n",imat);
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> eev and evs cards cannot be used together!");
    }
    
    if ( (mat.EEVF == OFF) && (mat.EVSF == OFF) ) {
        if (mod->flag.DIFFUSIVE_WAVE == OFF) {
            printf("ERROR: MATERIAL #: %d\n",imat);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> no eev or evs card defined.");
        }
    }
    
    if ( mod->flag.CORIOLIS && mat.coriolis == UNSET_FLT ) {
        printf("ERROR: MATERIAL #: %d\n",imat);
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> COR card not set for all materials");
    }
    
    if ( mat.EEVF == ON) {
        if ( (mat.EEVF < 1) || (mat.EEVF > 3) ) {
            printf("ERROR: MATERIAL #: %d\n",imat);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> EEVF mode should be 1,2 or 3.");
        }
    }
}
/**************************************************************/
/**************************************************************/

void smat_ns_check(SMODEL *mod, SMAT_NS mat, int imat) {
    
    if ( (mat.EEVF == ON) && (mat.EVSF == ON) ) {
        printf("ERROR: MATERIAL #: %d\n",imat);
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> eev and evs cards cannot be used together!");
    }
    
    if ( (mat.EEVF == OFF) && (mat.EVSF == OFF) ) {
        printf("ERROR: MATERIAL #: %d\n",imat);
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> no eev or evs card defined.");
    }
    
    if ( mod->flag.CORIOLIS && mat.coriolis == UNSET_FLT ) {
        printf("ERROR: MATERIAL #: %d\n",imat);
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> COR card not set for all materials");
    }
    
    if ( mat.EEVF == ON) {
        if ( (mat.EEVF < 1) || (mat.EEVF > 3) ) {
            printf("ERROR: MATERIAL #: %d\n",imat);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> EEVF mode should be 1,2 or 3.");
        }
    }
}
/**************************************************************/
/**************************************************************/

void smat_trn_check(SMODEL *mod, int EEVF, int EVSF, int DIFF_FLAG, int imat, int id) {
    
    // check that either DIFF_FLAG or EEV card is used
    if ( (EEVF == OFF) && (DIFF_FLAG == OFF) ) {
        printf("ERROR: MATERIAL #: %d :: CONSTITUENT: %d \n",imat,id);
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> must use either EEV or DIFF_FLAG card for constituent transport.");
    }
    
    // check if both cards are used - then warn of over-ride
    if ( (EVSF == ON) && (DIFF_FLAG == ON) ) {
        printf("\nWARNING **********************\n");
        printf("The specified value for DIFF_FLAG will over-ride EVS value for transport constituent id %d\n",id);
    }
    
}
/**************************************************************/
/**************************************************************/
#ifdef _ADH_GROUNDWATER
void smat_gw_check(SMODEL *mod, SMAT_GW mat, int imat) {

  /** mwf just some debugging output for now **/
  assert(mod);
  assert(mod->flag.GW_FLOW);
  /** 
  printf("*** begin material %d *** \n",imat);
  printf("k: %g %g %g %g %g %g \n",mat.k.xx,mat.k.yy,mat.k.zz,mat.k.xy,mat.k.xz,mat.k.yz);
  printf("s_s porosity water_vol: %g %g %g\n",mat.s_s,mat.porosity,mat.water_vol);
  printf("S_r vga vgn: %g %g %g\n",mat.residual_sat,mat.vangen_alpha,mat.vangen_n);
  printf("vg_cp vg_nxy bc_lam: %g %d %g\n",mat.vangen_max_cp,mat.vangen_num_xy,mat.brooks_lambda);
  printf("bc_pd bc_cp bc_nxy: %g %g %d\n",mat.brooks_pd,mat.brooks_max_cp,mat.brooks_num_xy);
  printf("ikr isat tor: %d %d %g\n",mat.ikr,mat.isat,mat.tortuosity);
  printf("d_l d_t ss_area: %g %g %g\n",mat.d_l,mat.d_t,mat.ss_area);
  printf("rho_b itype: %g %d \n",mat.bulk_density,mat.itype);
  printf("*** end material %d *** \n",imat);
  */
}
#endif
/**************************************************************/
/**************************************************************/

void smat_checkall(SMODEL *mod) {
    
    int imat=0;
    if (mod->flag.SW_FLOW) {
        for (imat = 0; imat<mod->nmat; imat++) {
            smat_sw_check(mod, *(mod->mat[imat].sw), imat);
        }
    }
    
    if (mod->flag.NS_FLOW) {
        for (imat = 0; imat<mod->nmat; imat++) {
            smat_ns_check(mod, *(mod->mat[imat].ns), imat);
        }
    }
#ifdef _ADH_GROUNDWATER
    if (mod->flag.GW_FLOW) {
      for (imat = 0; imat<mod->nmat; imat++) {
	smat_gw_check(mod, *(mod->mat[imat].gw), imat);
      }
    }
#endif
    int itrns=0;
    if (mod->flag.TRANSPORT) {
        for (imat = 0; imat<mod->nmat; imat++) {
            for (itrns=0; itrns<mod->ntransport; itrns++) {
                if (mod->flag.SW_FLOW) {
                    smat_trn_check(mod, mod->mat[imat].sw->EEVF, mod->mat[imat].sw->EVSF, mod->mat[imat].trn->DIFF_FLAG, imat, itrns);
                } else if (mod->flag.NS_FLOW) {
                    smat_trn_check(mod, mod->mat[imat].ns->EEVF, mod->mat[imat].ns->EVSF, mod->mat[imat].trn->DIFF_FLAG, imat, itrns);
                }
            }
        }
    }
    
#ifdef _SEDIMENT
    int ised=0;
    if (mod->flag.SEDIMENT) {
        for (imat = 0; imat<mod->nmat; imat++) {
            for (ised=0; ised<mod->nsed; ised++) {
                smat_trn_check(mod, *(mod->mat[imat].sw), mod->mat[imat].sed[ised], imat, ised);
            }
        }
    }
#endif
    
}

/**************************************************************/
/**************************************************************/

void smat_sw_init(SMODEL *mod, SMAT_SW *mat) {
    
    mat->EEVF = OFF;
    mat->EVSF = OFF;
    mat->max_lev = 0;
    mat->ev.xx = 0.;
    mat->ev.yy = 0.;
    mat->ev.xy = 0.;
    mat->ev.xz = 0.;
    mat->ev.yz = 0.;
    mat->ev.zz = 0.;
    mat->refine_tolerance = 1; //UNSET_FLT;
    mat->unrefine_tolerance = UNSET_FLT;
    mat->windatt = 1.0;
    mat->eev_mode = UNSET_INT;
    mat->fraction = 0.0;
    mat->turbulence_model_xy = UNSET_INT;
    mat->turbulence_model_z = UNSET_INT;
    if (mod->flag.CORIOLIS == OFF) {
        mat->coriolis = 0.0;
    } else {
        mat->coriolis = UNSET_FLT;
    }
}

/**************************************************************/
/**************************************************************/

void smat_ns_init(SMODEL *mod, SMAT_NS *mat) {
    
    mat->EEVF = OFF;
    mat->EVSF = OFF;
    mat->max_lev = 0;
    mat->ev.xx = 0.;
    mat->ev.yy = 0.;
    mat->ev.xy = 0.;
    mat->ev.xz = 0.;
    mat->ev.yz = 0.;
    mat->ev.zz = 0.;
    mat->refine_tolerance = 1; //UNSET_FLT;
    mat->unrefine_tolerance = UNSET_FLT;
    mat->windatt = 1.0;
    mat->eev_mode = UNSET_INT;
    mat->fraction = 0.0;
    mat->turbulence_model_xy = UNSET_INT;
    mat->turbulence_model_z = UNSET_INT;
    if (mod->flag.CORIOLIS == OFF) {
        mat->coriolis = 0.0;
    } else {
        mat->coriolis = UNSET_FLT;
    }
}

/**************************************************************/
/**************************************************************/

void smat_trn_init(SMODEL *mod, SMAT_TRN *mat) {
    
    mat->DIFF_FLAG = OFF;
    mat->react = NULL;
    mat->max_lev = UNSET_INT;
    mat->d_m = 0.;
    mat->d_l = 0.;
    mat->d_t = 0.;
    mat->source = 0.;
    mat->rd = 1.;
    mat->tortuosity = 1.;
    mat->refine_tolerance = 1; //UNSET_FLT;
}

/**************************************************************/
/**************************************************************/

void smat_wnsm_init(SMODEL *mod, SMAT_NSM *mat) {
#ifdef _NSM
    mat->alpha0 = 0.;
    mat->alpha1 = 0.;
    mat->alpha2 = 0.;
    mat->alpha3 = 0.;
    mat->alpha4 = 0.;
    mat->alpha5 = 0.;
    mat->alpha6 = 0.;
    mat->beta1 = 0.;
    mat->beta2 = 0.;
    mat->beta3 = 0.;
    mat->beta4 = 0.;
    mat->k1 = 0.;
    mat->k2 = 0.;
    mat->k3 = 0.;
    mat->k4 = 0.;
    mat->mu = 0.;
    mat->rho = 0.;
    mat->sigma1 = 0.;
    mat->sigma2 = 0.;
    mat->sigma3 = 0.;
    mat->sigma4 = 0.;
    mat->sigma5 = 0.;
    mat->kl = 0.;
    mat->kn = 0.;
    mat->kp = 0.;
    mat->pn = 0.;
    mat->lambda0 = 0.;
    mat->lambda1 = 0.;
    mat->min_do_conc = 0.5;
#endif
}

#ifdef _ADH_GROUNDWATER
void smat_gw_init(SMODEL *mod, SMAT_GW *mat) {
            
  mat->max_lev = 0;
  mat->refine_tolerance = 1; //UNSET_FLT;
  mat->unrefine_tolerance = UNSET_FLT;
  mat->k.xx = UNSET_FLT;
  mat->k.yy = UNSET_FLT;
  mat->k.xy = UNSET_FLT;
  mat->k.xz =UNSET_FLT;
  mat->k.yz =UNSET_FLT;
  mat->k.zz =UNSET_FLT;
  mat->s_s = UNSET_FLT;
  mat->water_vol = 0.0;
  mat->porosity = UNSET_FLT;
  mat->residual_sat = 0.0;
  mat->vangen_alpha = UNSET_FLT;
  mat->vangen_max_cp = 100.0;
  mat->vangen_n = UNSET_FLT;
  mat->vangen_num_xy = 400;
  mat->brooks_lambda = UNSET_FLT;
  mat->brooks_max_cp = 100.0;
  mat->brooks_pd = UNSET_FLT;
  mat->brooks_num_xy = 400;
  mat->ikr = UNSET_INT;
  mat->isat = UNSET_INT;
  mat->tortuosity = UNSET_FLT;
  mat->d_l = UNSET_FLT;
  mat->d_t = UNSET_FLT;
  mat->ss_area=UNSET_FLT;
  mat->bulk_density=UNSET_FLT;
  mat->itype = UNSET_INT;

}

#endif 
