/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_CN.c This file collects methods to read an AdH SuperModel boundary condition input file       */
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

void read_bc_CN(SMODEL_SUPER *mod, char **token) {
    
    int itrn = 0;
    
    if (strcmp(*token, "CON") == 0) {
        itrn = get_next_token_int(token) - 1;
        if (itrn > mod->ntransport - 1 || itrn < 0) {
            io_read_error(io, "The transport consistituent ID is either too big or too small.", TRUE);
        }
        mod->con[itrn].type = CON;
        mod->con[itrn].property[0] = get_next_token_dbl(token);
    } else if (strcmp(*token, "SAL") == 0) {
        itrn = get_next_token_int(token) - 1;
        if (itrn > mod->ntransport - 1 || itrn < 0) {
            io_read_error(io, "The transport consistituent ID is either too big or too small.", TRUE);
        }
        mod->con[itrn].type = SAL;
        mod->salinity_id = itrn;
        mod->flag.BAROCLINIC += 1;
        mod->con[itrn].property[0] = get_next_token_dbl(token);
    } else if (strcmp(*token, "TMP") == 0) {
        itrn = get_next_token_int(token) - 1;
        if (itrn > mod->ntransport - 1 || itrn < 0) {
            io_read_error(io, "The transport consistituent ID is either too big or too small.", TRUE);
        }
        mod->con[itrn].type = TMP;
        mod->temperature_id = itrn;
        mod->flag.BAROCLINIC += 10;
        mod->con[itrn].property[0] = get_next_token_dbl(token);
    } else if (strcmp(*token, "VOR") == 0) {
        itrn = get_next_token_int(token) - 1;
        if (itrn > mod->ntransport - 1 || itrn < 0) {
            io_read_error(io, "The transport consistituent ID is either too big or too small.", TRUE);
        }
        mod->flag.VORTICITY = ON;
        mod->con[itrn].type = VOR;
        mod->vorticity_id = itrn;
        mod->con[itrn].property[0] = get_next_token_dbl(token);
        mod->con[itrn].property[1] = get_next_token_dbl(token);
        mod->con[itrn].property[2] = get_next_token_dbl(token);
        if (mod->con[itrn].property[1] <= 0) {
            printf("Vorticity As term is specified as %lf. This is less than 0, defaulting to 5\n",mod->con[itrn].property[1]);
            mod->con[itrn].property[1] = 5.0;
        }
        if (mod->con[itrn].property[2] <= 0) {
            printf("Vorticity Ds term is specified as %lf. This is less than 0, defaulting to 0.5\n",mod->con[itrn].property[2]);
            mod->con[itrn].property[2] = 0.5;
        }
    }
    
}
