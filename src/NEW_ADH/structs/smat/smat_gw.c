#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes ground water material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat_gw (SMAT_GW **)  double pointer to a ground water material
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smat_gw_alloc_init(SMAT_GW **mat_gw) {
    
    (*mat_gw) = (SMAT_GW *) tl_alloc(sizeof(SMAT_GW), 1);
    SMAT_GW *mat = *mat_gw; // alias
    
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

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees ground water material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat (SMAT_GW *)  pointer to a ground water material
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smat_gw_free(SMAT_GW *mat) {
    mat = (SMAT_GW *) tl_free(sizeof(SMAT_GW), 1, mat);
}
