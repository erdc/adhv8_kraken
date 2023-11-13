/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Initializes 2D and 3D transport solutions and boundary conditions for a Newton iterate.
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Gary Brown, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout]  mod (SMODEL *) pointer to the model struct
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_transport_init(SSUPER_MODEL *sm, int imod) {
    
    SMODEL *mod = &(sm->submodel[imod]);

    int DEBUG = OFF;
    
    double *c, *old_c;
    double c_inv;                 /* 1./con[itrn].property[0] */
    int i, ie;                    /* loop counter */
    int istr;                     /* the string number for evaluating boundary conditions */
    int isers;                    /* the series number for evaluating boundary conditions */
    int index = UNSET_INT;
    
    // is this being ran with sediment or as a general constituent?
    if (mod->is_sediment_running) {
#ifdef _SEDIMENT
        index = mod->ised;
        c_inv = 1.E-6 / mod->sed->grain[index].reference_c;
        c = mod->sed->susload[index].c;
        old_c = mod->sed->susload[index].old_c;
#endif
    } else {
        index = mod->itrns;
        c_inv = 1. / mod->con[index].property[0];
        c = mod->con[index].concentration;
        old_c = mod->con[index].old_concentration;
    }
    
    /* copies the old solution to the current solution */
    for (i = 0; i < mod->grid->nnodes; i++) {
        if (solv_isnan(c[i]) || solv_isinf(c[i])) {
            c[i] = 0.0;
        }
        
        mod->bc_mask[i] = NO;
        if (mod->grid->node[i].string > NORMAL) {
            istr = mod->grid->node[i].string;
            
            if (mod->is_sediment_running) {
#ifdef _SEDIMENT
                if (mod->str_values[istr].sed[index].bc_flag == BCT_CEQ) {
                    mod->bc_mask[i] = YES;
                }
#endif
            } else {
                if (mod->str_values[istr].trans[index].bc_flag == BCT_CEQ) {
                    mod->bc_mask[i] = YES;
                }
            }
        }
        
        if (mod->bc_mask[i] == NO) {
            c[i] = old_c[i];
        }
        
        if (mod->flag.SW2_FLOW && mod->bc_mask[i] == NO && mod->sw->d2->head[i] <= NOT_QUITE_SMALL) {
            c[i] = 0.0;
        }
        
    }
    
    
    /* enforce the Dirichlet boundary conditions */
    INPUT_DATA *str_value;
    for (i = 0; i < mod->grid->nnodes; i++) {
        if (mod->grid->node[i].string > NORMAL) {
            istr = mod->grid->node[i].string;
            
            if (mod->is_sediment_running) {
#ifdef _SEDIMENT
                str_value = &(mod->str_values[istr].sed[index]);
#endif
            } else {
                str_value = &(mod->str_values[istr].trans[index]);
            }
            
            if (str_value->bc_flag == BCT_DIR) {
                mod->bc_mask[i] = YES;
                isers = str_value->iu_0;
                c[i] = sseries_get_value(isers, mod->series_head,0) * c_inv;
            }
        }
    }
    
#ifdef _DEBUG
    if (DEBUG) {tl_check_all_pickets(__FILE__, __LINE__);}
#endif
    
#ifdef _MESSG
    comm_update_int(mod->bc_mask, 1, mod->grid->smpi);
#endif
    
}
