/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_TC.c  This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for reading model time controls
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] *sm (SUPER_MODEL *)  a double pointer to an array of AdH supermodels
 * @param[in] **token (CHAR) a BC file line string token
 *
 * \note CJT -- units on old AdH were optional
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_bc_TC(SMODEL_SUPER *mod, char **token) {

  int units;
    if (strcmp(*token, "T0") == 0) {
        mod->t_prev = mod->t_init = get_next_token_dbl(token);
        units = get_next_token_int(token);
        mod->t_init *= tc_conversion_factor(units, TO);
        
    } else if (strcmp(*token, "TF") == 0) {
        mod->t_final = get_next_token_dbl(token);
        units = get_next_token_int(token);
        mod->t_final *= tc_conversion_factor(units, TO);

    } else if (strcmp(*token, "NDP") == 0) {
        mod->t_adpt_flag = OFF;
    }
    
    
}
