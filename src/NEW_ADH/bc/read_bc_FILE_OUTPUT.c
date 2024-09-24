/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_FILE_OUTPUT.cThis file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for reading file outputs
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

void read_bc_FILE_OUTPUT(SMODEL_SUPER *mod, char **token) {
    
    if (strcmp(*token, "BEDVEL") == 0) {mod->file_output.bed_velocity = ON;}
    else if(strcmp(*token, "SURVEL") == 0) {mod->file_output.surface_velocity = ON;}
    else if(strcmp(*token, "AVGVEL") == 0) {mod->file_output.depth_avg_velocity = ON;}
    else if(strcmp(*token, "PRS") == 0) {mod->file_output.pressure = ON;}
    else if(strcmp(*token, "GSPEED") == 0) {mod->file_output.grid_speed = ON;}
    else if(strcmp(*token, "WIND") == 0) {mod->file_output.wind = ON;}
    else if(strcmp(*token, "WAVE") == 0) {mod->file_output.wave = ON;}
    else if(strcmp(*token, "VIS") == 0) {mod->file_output.hyd_vis = ON;}
    else if(strcmp(*token, "DIF") == 0) {mod->file_output.trn_dif = ON;}
    else if(strcmp(*token, "CHOP") == 0) {mod->file_output.chop= ON;}
    else if(strcmp(*token, "GRID") == 0) {mod->file_output.grid2dm = ON;}
    else if(strcmp(*token, "ADAPT") == 0) {
        get_next_token(token);
        if (strcmp(*token, "GRID") == 0) {
            mod->file_output.adaption = ON;
            mod->file_output.adapt_grid = ON;
        } else if (strcmp(*token, "SW") == 0) {
            mod->file_output.adaption = ON;
            mod->file_output.adapt_grid = ON;
            mod->file_output.adapt_sw = ON;
        } else if (strcmp(*token, "NS") == 0) {
            mod->file_output.adaption = ON;
            mod->file_output.adapt_grid = ON;
            mod->file_output.adapt_ns = ON;
        } else if (strcmp(*token, "CON") == 0) {
            mod->file_output.adaption = ON;
            mod->file_output.adapt_grid = ON;
            mod->file_output.adapt_con = ON;
        } else if (strcmp(*token, "SED") == 0) {
            mod->file_output.adaption = ON;
            mod->file_output.adapt_grid = ON;
            mod->file_output.adapt_sed = ON;
        } else if (strcmp(*token, "GW") == 0) {
            mod->file_output.adaption = ON;
            mod->file_output.adapt_grid = ON;
            mod->file_output.adapt_gw = ON;
        } else {
            mod->file_output.adaption = ON;
            mod->file_output.adapt_grid = ON;
            mod->file_output.adapt_sw = ON;
            mod->file_output.adapt_con = ON;
            mod->file_output.adapt_sed = ON;
        }
    }
    else if(strcmp(*token, "ALL") == 0) {
        mod->file_output.grid2dm = ON;
        mod->file_output.bed_velocity = ON;
        mod->file_output.surface_velocity = ON;
        mod->file_output.depth_avg_velocity = ON;
        mod->file_output.pressure = ON;
        mod->file_output.grid_speed = ON;
        mod->file_output.wind = ON;
        mod->file_output.wave = ON;
        mod->file_output.hyd_vis = ON;
        mod->file_output.trn_dif = ON;
    }
}
