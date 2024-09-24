/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_DB.c This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for reading transport constituents
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] *sm (SUPER_MODEL *)  a double pointer to an array of AdH supermodels
 * @param[in] **token (CHAR) a BC file line string token
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_bc_DB(SMODEL_SUPER *mod, char **token) {
    int ibc, itrn, iseries;
    
    if (strcmp(*token, "VEL") == 0) {
        ibc = get_id(token,mod->nstring,"Velocity String Error\n");
        mod->str_values[ibc].flow.bc_flag = BCT_VEL_DIR;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].flow.ivx = iseries;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].flow.ivy = iseries;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].flow.ivz = iseries;
    } else if (strcmp(*token, "PRS") == 0) {
        ibc = get_id(token,mod->nstring,"Pressure String Error\n");
        mod->str_values[ibc].pressure.bc_flag = BCT_PRS_DIR;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].pressure.iu_0 = iseries;
    } else if (strcmp(*token, "FRS") == 0) {
        ibc = get_id(token,mod->nstring,"Friction String Error\n");
        mod->str_values[ibc].displacement.bc_flag = BCT_FREE_DIR;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].displacement.iu_0 = iseries;
    } else if (strcmp(*token, "OVL") == 0) {
        ibc = get_id(token,mod->nstring,"Velocity String Error\n");
        mod->str_values[ibc].ol_flow.bc_flag = BCT_VEL_DIR;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].ol_flow.ivx = iseries;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].ol_flow.ivy = iseries;
    } else if (strcmp(*token, "OVH") == 0) {
        ibc = get_id(token,mod->nstring,"Pressure String Error\n");
        mod->str_values[ibc].ol_flow.bc_flag = BCT_VEL_PRS_DIR;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].ol_flow.ivx = iseries;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].ol_flow.ivy = iseries;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].ol_flow.iu_0 = iseries;
    } else if (strcmp(*token, "TRN") == 0) {
        ibc = get_id(token,mod->nstring,"Transport String Error\n");
        itrn = get_id(token,mod->nstring,"Transport String Error\n");
        mod->str_values[ibc].trans[itrn].bc_flag = BCT_DIR;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].trans[itrn].iu_0 = iseries;
    }
#ifdef _ADH_GROUNDWATER
    else if (strcmp(*token, "THD") == 0) {
        ibc = get_id(token,mod->nstring,"THD String Error\n");
        mod->str_values[ibc].flow.bc_flag = BCT_DIR;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].flow.iu_0 = iseries;
    } else if (strcmp(*token, "FLW") == 0) {
        ibc = get_id(token,mod->nstring,"FLW String Error\n");
        mod->str_values[ibc].flow.bc_flag = BCT_DIR;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].flow.iu_0 = iseries;
    } else if (strcmp(*token, "PSI") == 0) {
        ibc = get_id(token,mod->nstring,"PSE String Error\n");
        mod->str_values[ibc].flow.bc_flag = BCT_PRS_DIR;
        iseries = sseries_set_type(mod, token, TIME_SERIES);
        mod->str_values[ibc].flow.iu_0 = iseries;
    }
#endif

}
