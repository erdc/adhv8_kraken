#include "global_header.h"
void error_main(SMODEL *mod) {
    if (mod->flag.GRID_ADAPTION) {
        
        int ie;
        
        // CJT :: FORCE ADAPTIVE REFINEMENT!
        //mod->grid->elem_error[2609] = 0;
        
        // CJT :: FORCE ADAPTIVE UNREFINEMENT OF A REFINED ELEMENT
        //for (ie = 0; ie < mod->grid->nelems3d; ie++) {
        //    if (elem3d_level(mod->grid->elem3d[ie]) > 0) mod->grid->elem_error[ie] = 0;
        //}
        //printf("\nmod->grid->elem3d[ie].mat: %d UNREF_TOL: %20.10f \t elem3d_level(mod->grid->elem3d[ie]): %d\n",mod->grid->elem3d[2609].mat,UNREF_TOL,elem3d_level(mod->grid->elem3d[2609]));

        mod->flag.GRID_REFINED = NO;
        mod->flag.GRID_UNREFINED = NO;
        
        if (mod->flag.SW3_FLOW) {
            for (ie = 0; ie < mod->grid->nelems3d; ie++) {
                if (mod->grid->elem3d[ie].mat != UNSET_INT && mod->grid->elem_error[ie] > REF_TOL && elem3d_level(mod->grid->elem3d[ie]) < mod->mat[mod->grid->elem3d[ie].mat].sw->max_lev && mod->grid->elem3d[ie].interface==0) {
                    mod->flag.GRID_REFINED = YES;
                }
                if (mod->grid->elem3d[ie].mat != UNSET_INT && mod->grid->elem_error[ie] < UNREF_TOL && elem3d_level(mod->grid->elem3d[ie]) > 0){
                    mod->flag.GRID_UNREFINED = YES;
                }
            }
        } else if (mod->flag.NS3_FLOW) {
            for (ie = 0; ie < mod->grid->nelems3d; ie++) {
                if (mod->grid->elem3d[ie].mat != UNSET_INT && mod->grid->elem_error[ie] > REF_TOL && elem3d_level(mod->grid->elem3d[ie]) < mod->mat[mod->grid->elem3d[ie].mat].ns->max_lev && mod->grid->elem3d[ie].interface==0) {
                    mod->flag.GRID_REFINED = YES;
                }
                if (mod->grid->elem3d[ie].mat != UNSET_INT && mod->grid->elem_error[ie] < UNREF_TOL && elem3d_level(mod->grid->elem3d[ie]) > 0){
                    //printf("3D ie: %d mat: %d \t error: %20.10e level: %d \t UNREF_TOL: %20.10e \t node levls: %d %d %d %d\n",
                    //        ie,mod->grid->elem3d[ie].mat,mod->grid->elem_error[ie],elem3d_level(mod->grid->elem3d[ie]),UNREF_TOL,
                    //        mod->grid->elem3d[ie].levels[0],mod->grid->elem3d[ie].levels[1],mod->grid->elem3d[ie].levels[2],mod->grid->elem3d[ie].levels[3]);
                    mod->flag.GRID_UNREFINED = YES;
                }
            }
#ifdef _ADH_GROUNDWATER
        } else if (mod->flag.GW_FLOW) {
            for (ie = 0; ie < mod->grid->nelems3d; ie++) {
                if (mod->grid->elem3d[ie].mat != UNSET_INT && mod->grid->elem_error[ie] > REF_TOL && elem3d_level(mod->grid->elem3d[ie]) < mod->mat[mod->grid->elem3d[ie].mat].gw->max_lev && mod->grid->elem3d[ie].interface==0) {
                    mod->flag.GRID_REFINED = YES;
                }
                if (mod->grid->elem3d[ie].mat != UNSET_INT && mod->grid->elem_error[ie] < UNREF_TOL && elem3d_level(mod->grid->elem3d[ie]) > 0){
                    //printf("3D ie: %d mat: %d \t error: %20.10e level: %d \t UNREF_TOL: %20.10e \t node levls: %d %d %d %d\n",
                    //        ie,mod->grid->elem3d[ie].mat,mod->grid->elem_error[ie],elem3d_level(mod->grid->elem3d[ie]),UNREF_TOL,
                    //        mod->grid->elem3d[ie].levels[0],mod->grid->elem3d[ie].levels[1],mod->grid->elem3d[ie].levels[2],mod->grid->elem3d[ie].levels[3]);
                    mod->flag.GRID_UNREFINED = YES;
                }
            }
#endif
        } else {
            for (ie = 0; ie < mod->grid->nelems2d; ie++) {
                if (mod->grid->elem2d[ie].mat != UNSET_INT && mod->grid->elem_error[ie] > REF_TOL && elem2d_level(mod->grid->elem2d[ie]) < mod->mat[mod->grid->elem2d[ie].mat].sw->max_lev && mod->grid->elem2d[ie].interface==0) {
                    mod->flag.GRID_REFINED = YES;
                }
                if (mod->grid->elem2d[ie].mat != UNSET_INT && mod->grid->elem_error[ie] < UNREF_TOL && elem2d_level(mod->grid->elem2d[ie]) > 0){
                    mod->flag.GRID_UNREFINED = YES;
                    //printf("2D ie: %d mat: %d \t error: %20.10e level: %d \t UNREF_TOL: %20.10e\n",ie,mod->grid->elem3d[ie].mat,mod->grid->elem_error[ie],elem3d_level(mod->grid->elem3d[ie]),UNREF_TOL);
                }
            }
        }
    }
}
