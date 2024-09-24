
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_NOTERM.c  This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file to prevent certain terms from being included in the rhs in 3d, this is for the hvel (momentum equations)
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
void read_bc_NOTERM(SMODEL_SUPER *mod, char **token) {
    
    if (strcmp(*token, "HYDRO") == 0) {
        debug.no_hydro = ON;
        debug.u_vel = read_dbl_field(info, &subdata);
        debug.v_vel = read_dbl_field(info, &subdata);
        debug.w_vel = read_dbl_field(info, &subdata);
        debug.hard_disp = read_dbl_field(info, &subdata);
        mod->max_nsys = 1;
        mod->max_nsys_sq = 1;
    }
    else if (strcmp(*token, "DIFF") == 0)  {debug.no_diffusion = ON;}
    else if (strcmp(*token, "CONV") == 0)  {debug.no_advection = ON;}
    else if (strcmp(*token, "SUPG") == 0)  {debug.no_supg = ON;}
    else if (strcmp(*token, "FRIC") == 0)  {debug.no_friction = ON;}
    else if (strcmp(*token, "TIME") == 0)  {debug.no_temporal = ON;}
    else if (strcmp(*token, "PRESS") == 0) {debug.no_pressure = ON;}
    else if (strcmp(*token, "COR") == 0)   {debug.no_coriolis = ON;}
    
}
