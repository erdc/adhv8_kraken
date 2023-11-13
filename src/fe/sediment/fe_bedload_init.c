/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Initializes bed load transport solutions and boundary conditions for a Newton iterate.
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

void fe_bedload_init(SMODEL *model) {
    
    int DEBUG = OFF;

    double c_inv;                 /* 1./con[itrn].property[0] */
    int i, ie;                    /* loop counter */
    int istr;                     /* the string number for evaluating boundary conditions */
    int isers;                    /* the series number for evaluating boundary conditions */
  
    // aliases
    int ised = model->ised;
    SSED *sed = model->sed;
    int nnodes = model->grid->nnodes_bed; // = nnodes in 2d
    c_inv = 1. / sed->grain[ised].reference_c;
    
    /* copies the old solution to the current solution */
    int inode = UNSET_INT;
    for (i = 0; i < nnodes; i++) {
        
        inode = i;
        if (model->grid->ndim == 3) {
            inode = model->grid->nodeID_2d_to_3d_bed[i];
        }
        
        // really?
        if (solv_isnan(sed->susload[ised].c[inode]) || solv_isinf(sed->susload[ised].c[inode])) {
            sed->bedload[ised].c[i] = 0.0;
        }
        
        model->bc_mask[inode] = NO;
        if (model->grid->node[inode].string > NORMAL) {
            istr = model->grid->node[inode].string;
            if (model->str_values[istr].sed[ised].bc_flag == BCT_CEQ) {
                model->bc_mask[inode] = YES;
            }
        }

        if (model->bc_mask[inode] == NO) {
            sed->bedload[ised].c[i] = sed->bedload[ised].old_c[i];
            if (model->grid->ndim == 2 && model->sw->d2->head[i] <= NOT_QUITE_SMALL) {
                sed->bedload[ised].c[i] = 0.0;
            }
        }
    }

    /* enforce the Dirichlet boundary conditions */
    for (i = 0; i < nnodes; i++) {
        
        inode = i;
        if (model->grid->ndim == 3) {
            inode = model->grid->nodeID_2d_to_3d_bed[i];
        }
        
        if (model->grid->node[inode].string > NORMAL) {
            istr = model->grid->node[inode].string;
            if (model->str_values[istr].sed[ised].bc_flag == BCT_DIR) {
                model->bc_mask[inode] = YES;
                isers = model->str_values[istr].sed[ised].iu_0;
                sed->bedload[ised].c[i] = sseries_get_value(isers, model->series_head,0) * c_inv * 1.E-6;
            }
        }
    }

#ifdef _DEBUG
    if (DEBUG) {
        tl_check_all_pickets(__FILE__, __LINE__);
    }
#endif

}
