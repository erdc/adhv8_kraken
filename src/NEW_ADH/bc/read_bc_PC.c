/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_PC.c  This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for model print controls
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

void read_bc_PC(SMODEL_SUPER *mod, char **token) {
    
    if (strcmp(*token, "LVL") == 0) {
        mod->out_level = get_next_token_int(token);
    } else if (strcmp(*token, "XDF") == 0) {
#ifdef _ADH_HDF5
        mod->flag.PC_FILE_XDMF = ON;
#ifdef _MESSG
        if (mod->grid->smpi->myid==0)
#endif
            printf("\nSwitching output type to XDMF ...");
        
#else
        printf("\nWarning: Attempted to turn on XDMF output, but code was not compiled with HDF5 support.");
        printf("\nIgnoring XDMF card; Default AdH output will be used.");
        //tl_error("Attempted to turn on XDMF output, but code was not compiled with HDF5 support.\n");
#endif
    }
    
}
