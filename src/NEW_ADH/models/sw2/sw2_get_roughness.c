/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes the eddy viscosity coefficients for the SW system.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] ev_st (double *) the streamline eddy viscosity coefficient
 * @param[out] ev_tr (double *) the transverse eddy viscosity coefficient
 * @param[in]  mod (SMODEL *) a pointer to a model struct
 * @param[in]  depth (double) the local elemental depth
 * @param[in]  velocity (double) the local elemental velocity tangent to the face
 * @param[in]  id (int) the elemental material string id
 * @param[in]  which (int) and ice friction function flag
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double fe_sw2_get_roughness(SMODEL_SUPER *mod, double depth, double velocity, int id, int which) {
    
    double zero_drag = 0.0;
    if(mod->str_values[id].fterms.sav_flag == YES) {
        return fr_sav_drag_coef(depth, mod->str_values[id].fterms.eqrheight,
                                mod->str_values[id].fterms.hghtstem);
    }
    else if(mod->str_values[id].fterms.urv_flag == YES) {
        return fr_urv_drag_coef(depth, mod->str_values[id].fterms.eqrheight,
                                mod->str_values[id].fterms.diamstem,
                                mod->str_values[id].fterms.densstem);
    }
    else if(mod->str_values[id].fterms.icethick != UNSET_FLT) {
        return fr_stationary_ice_coef(depth, velocity, mod->str_values[id].fterms.bedrhght,
                                      mod->str_values[id].fterms.icerhght, mod->density, which);
    }
    else {
        if(mod->str_values[id].fterms.eqrheight > SMALL) {
            return fr_bedshstr_drag_coef(depth, mod->str_values[id].fterms.eqrheight);
        }
        else {
            return zero_drag;
        }
    }
}
