#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes surface water material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat_sw (SMAT_SW **)  double pointer to a surface water material
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smat_sw_alloc_init(SMAT_SW **mat_sw) {
    
    (*mat_sw) = (SMAT_SW *) tl_alloc(sizeof(SMAT_SW), 1);
    SMAT_SW *mat = *mat_sw; // alias
    
    mat->EEVF = OFF;
    mat->EVSF = OFF;
    mat->ev.xx = 0.;
    mat->ev.yy = 0.;
    mat->ev.xy = 0.;
    mat->ev.xz = 0.;
    mat->ev.yz = 0.;
    mat->ev.zz = 0.;
    mat->windatt = 1.0;
    mat->eev_mode = UNSET_INT;
    mat->fraction = 0.0;
    mat->turbulence_model_xy = UNSET_INT;
    mat->turbulence_model_z = UNSET_INT;
    mat->coriolis = 0.0;
    bed_disp_flag = 0;
    vor_flag = 0;
    eev_coef = 0;
    d_flag = 0;
    wind_flag = 0;
    smag_coeff = 0.0;
    supression_func = UNSET_INT;
    wall_func = UNSET_INT;
    min_tke = 0.0;
    min_tds = 0.0;
    len_max = 0.0;
    hyd_conductivity = 0.0;
    psi = 0.0;
    rooting_depth = 0.0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees surface water material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat (SMAT_SW *)  pointer to a surface water material
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smat_sw_free(SMAT_SW *mat) {
    mat = (SMAT_SW *) tl_free(sizeof(SMAT_SW), 1, mat);
}
