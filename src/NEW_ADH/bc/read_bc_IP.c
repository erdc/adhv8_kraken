/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_IP.c This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for reading solver parameters
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
void read_bc_IP(SMODEL_SUPER *mod, char **token) {
    
    if (strcmp(*token, "MIT") == 0) {
        mod->solver_info.max_lin_it = get_next_token_int(token);
        mod->solver_info.force_lin_it = NO;
    } else if (strcmp(*token, "NIT") == 0) {
        mod->max_nonlin_it = get_next_token_int(token);
        mod->solver_info.force_nonlin_it = NO;
    } else if (strcmp(*token, "ITL") == 0) {
        mod->inc_nonlin  = get_next_token_dbl(token);
    } else if (strcmp(*token, "NTL") == 0) {
        mod->tol_nonlin = get_next_token_dbl(token);
    }
    
}
