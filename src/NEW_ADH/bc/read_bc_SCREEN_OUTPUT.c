/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_SCREEN_OUTPUT.c  This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for deciding screen output
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

void read_bc_SCREEN_OUTPUT(SMODEL_SUPER *mod, char **token) {

    if      (strcmp(*token, "RESID") == 0)  {screen_output.residuals = 1;}
    else if (strcmp(*token, "NLNODE") == 0) {screen_output.worse_node_nonlinear = 1;}
    else if (strcmp(*token, "LNODE") == 0)  {screen_output.worse_node_linear = 1;}
    else if (strcmp(*token, "MERROR") == 0) {screen_output.grid_mass_error  = 1;}
    else if (strcmp(*token, "ALL") == 0)    {
        screen_output.residuals = 1;
        screen_output.worse_node_nonlinear = 1;
        screen_output.worse_node_linear = 1;
    }
}
