/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_OP.c  This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for reading model operational inputs
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

void read_bc_OP(SMODEL_SUPER *mod, char **token) {
    
    if (strcmp(*token, "INC") == 0) {
        mod->nalloc_inc = get_next_token_int(token);
    } else if (strcmp(*token, "SW2") == 0) {
        mod->flag.SW_FLOW = ON;
        mod->flag.SW2_FLOW = ON;
        if (mod->ntransport > 0) mod->flag.SW2_TRANSPORT = ON;
    } else if (strcmp(*token, "SW3") == 0) {
        mod->flag.SW_FLOW = ON;
        mod->flag.SW3_FLOW = ON;
        mod->flag.MG = ON;
        if (mod->ntransport > 0) mod->flag.SW3_TRANSPORT = ON;
    } else if (strcmp(*token, "DIF") == 0) {
        mod->flag.SW_FLOW = ON;
        mod->flag.SW2_FLOW = ON;
        mod->flag.DIFFUSIVE_WAVE = ON;
        printf("Diffusive Wave Equation Set Is Active\n");
        printf("Velocities will be computed based on the water surface slope. Equation set is 1 DOF\n");
    } else if (strcmp(*token, "NS2") == 0) {
        mod->flag.NS_FLOW = ON;
        mod->flag.NS2_FLOW = ON;
        if (mod->ntransport > 0) mod->flag.NS2_TRANSPORT = ON;
    } else if (strcmp(*token, "NS3") == 0) {
        mod->flag.NS_FLOW = ON;
        mod->flag.NS3_FLOW = ON;
        if (mod->ntransport > 0) mod->flag.NS3_TRANSPORT = ON;
    } else if (strcmp(*token, "WAV") == 0) {
        printf("\n Short waves are turned on during this simulation. \n");
        printf("-- Wave data can be initialized through the hotstart file or given via CSTORM.\n");
        mod->flag.WAVE = ON;
    } else if (strcmp(*token, "WND") == 0) {
        printf("\n Wind stressing is turned on during this simulation. \n");
        mod->flag.WIND = ON;
    } else if (strcmp(*token, "TEM") == 0) {
        mod->tau_temporal = get_next_token_dbl(token);
        if (mod->tau_temporal < 0 || mod->tau_temporal > 1) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error("Tau_temporal must be between 0 and 1.");
        }
    } else if (strcmp(*token, "TPG") == 0) {
        mod->tau_pg = get_next_token_dbl(token);
        if (mod->tau_pg < 0 || mod->tau_pg > 0.5) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error("Tau_pg must be between 0 and .5");
        }
    } else if (strcmp(*token, "PRE") == 0) {
        mod->solver_info.prec_value = get_next_token_int(token);
    } else if (strcmp(*token, "MG") == 0) {
        mod->flag.MG = ON;
    } else if (strcmp(*token, "SEDLIB") == 0) {
#ifdef _SEDLIB
        read_bc_SEDLIB(mod, token);
#endif
    } else if (strcmp(*token, "EOS") == 0) {
        mod->flag.EOS = get_next_token_int(token);
        printf("AdH will use Equation of State %d (0: Linearlized 1: Full Equation)\n", mod->flag.EOS);
    }
    
}
