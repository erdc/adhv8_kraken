/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_FR.c  This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for reading friction inputs
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

void read_bc_FR(SMODEL_SUPER *mod, char **token) {
    
    if (strcmp(*token, "MNG") == 0) {
        int ibc = get_id(token,mod->nstring,"MNG String Error\n");
        mod->str_values[ibc].roughness = get_next_token_dbl(token);
        mod->str_values[ibc].fterms.manningsn = mod->str_values[ibc].roughness;
        mod->str_values[ibc].fterms.mng_flag = YES;
        
    } else if (strcmp(*token, "ERH") == 0) {
        int ibc = get_id(token,mod->nstring,"MNG String Error\n");
        if (mod->str_values[ibc].fterms.eqrheight == UNSET_FLT) {
            mod->str_values[ibc].fterms.eqrheight = get_next_token_dbl(token);
        }
        mod->str_values[ibc].fterms.erh_flag = YES;
    }
    
}
