#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes transport material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat_trn (SMAT_TRN **)  double pointer to a surface water material
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smat_trn_alloc_init(SMAT_TRN **mat_trn) {
    
    (*mat_trn) = (SMAT_TRN *) tl_alloc(sizeof(SMAT_TRN), 1);
    SMAT_TRN *mat = *mat_trn; // alias
    
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

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees transport material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat (SMAT_TRN *)  pointer to a transport material
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smat_trn_free(SMAT_TRN *mat) {
    mat = (SMAT_TRN *) tl_free(sizeof(SMAT_TRN), 1, mat);
}
