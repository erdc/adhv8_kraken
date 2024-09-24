/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_DEBUG,c This file collects methods to read an AdH SuperModel boundary condition input file  for debugging      */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for reading debugging cards
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

void read_bc_DEBUG(SMODEL *mod, char **token) {
    if      (strcmp(*token, "DEBUG")  == 0) {debug.print = 1;}
    else if (strcmp(*token, "READBC") == 0) {debug.readBC = 1;}
    else if (strcmp(*token, "RHS")    == 0) {debug.rhs = 1;}
    else if (strcmp(*token, "RHSHV")  == 0) {debug.rhs_hvel = 1;}
    else if (strcmp(*token, "RHSWV")  == 0) {debug.rhs_wvel = 1;}
    else if (strcmp(*token, "MAT")    == 0) {debug.matrix = 1;}
    else if (strcmp(*token, "MATHV")  == 0) {debug.matrix_hvel = 1;}
    else if (strcmp(*token, "MATWV")  == 0) {debug.matrix_wvel = 1;}
    else if (strcmp(*token, "NEWTON") == 0) {debug.newton = 1;}
    else if (strcmp(*token, "RES")    == 0) {debug.residual = 1;}
    else if (strcmp(*token, "RESHV")  == 0) {debug.residual_hvel = 1;}
    else if (strcmp(*token, "RESWV")  == 0) {debug.residual_wvel = 1;}
    else if (strcmp(*token, "DIFF")   == 0) {debug.diffusion = 1;}
    else if (strcmp(*token, "CONV")   == 0) {debug.advection = 1;}
    else if (strcmp(*token, "SUPG")   == 0) {debug.supg = 1;}
    else if (strcmp(*token, "FRIC")   == 0) {debug.friction = 1;}
    else if (strcmp(*token, "TIME")   == 0) {debug.temporal = 1;}
    else if (strcmp(*token, "PRESS")  == 0) {debug.pressure = 1;}
    else if (strcmp(*token, "LOAD")   == 0) {debug.load = 1;}
    else if (strcmp(*token, "LOADHV") == 0) {debug.load_hvel = 1;}
    else if (strcmp(*token, "LOADWV") == 0) {debug.load_wvel = 1;}
    else if (strcmp(*token, "WAVE")   == 0) {debug.waves = 1;}
    else if (strcmp(*token, "WIND")   == 0) {debug.winds = 1;}
}
