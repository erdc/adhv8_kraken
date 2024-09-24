/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_NB.c This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for reading natural boundary conditions
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

void read_bc_NB(SMODEL_SUPER *mod, char **token) {
    
    int itrn = UNSET_INT;
    int ibc = UNSET_INT, iseries = UNSET_INT;
    int nwer, ibc1;
    int nflp;
    int nslu;
    SIO info = *(mod->io);
    int nstring = mod->nstring;   // alias
    int nseries = mod->nseries;   // alias
    nwer = 0;
    nflp = 0;
    
    if (strcmp(*token, "OF") == 0) {
        ibc = get_id(token,mod->nstring,"OF String Error\n");
        if (mod->flag.SW2_FLOW || mod->flag.NS2_FLOW) {
            mod->str_values[ibc].ol_flow.bc_flag = BCT_OUTFLOW;
        } else if (mod->flag.SW3_FLOW || mod->flag.NS3_FLOW) {
            mod->str_values[ibc].flow.bc_flag = BCT_OUTFLOW;
        }
    } else if (strcmp(*token, "DIS") == 0) {
        ibc = get_id(token,mod->nstring,"DIS String Error\n");
        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
        if (mod->flag.SW2_FLOW || mod->flag.NS2_FLOW) {
            mod->str_values[ibc].ol_flow.bc_flag = BCT_DIS_NEU;
            mod->str_values[ibc].ol_flow.isigma = iseries;
            mod->flag.CONVEYANCE = YES;
        }  else if (mod->flag.SW3_FLOW || mod->flag.NS3_FLOW) {
            mod->str_values[ibc].flow.bc_flag = BCT_DIS_NEU;
            mod->str_values[ibc].flow.isigma = iseries;
        }
    } else if (strcmp(*token, "CPL") == 0) {
        ibc = get_id(token,mod->nstring,"CPL String Error\n");
        if (mod->flag.SW2_FLOW) {
            mod->str_values[ibc].ol_flow.bc_flag = BCT_FLUX_COUPLE;
        } else {
            mod->str_values[ibc].flow.bc_flag = BCT_FLUX_COUPLE;
        }
    } else if (strcmp(*token, "PRS") == 0) {
        ibc = get_id(token,mod->nstring,"PRS String Error\n");
        mod->str_values[ibc].pressure.bc_flag = BCT_PRS_NEU;
        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
        mod->str_values[ibc].pressure.isigma = iseries;
    } else if (strcmp(*token, "FRS") == 0) {
        ibc = get_id(token,mod->nstring,"FRS String Error\n");
        mod->str_values[ibc].flow.bc_flag = BCT_FRS;
        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
        mod->str_values[ibc].flow.isigma = iseries;
    } else if (strcmp(*token, "OTW") == 0) {
        ibc = get_id(token,mod->nstring,"OTW String Error\n");
        mod->str_values[ibc].ol_flow.bc_flag = BCT_PRS_NEU;
        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
        mod->str_values[ibc].ol_flow.isigma = iseries;
    } else if (strcmp(*token, "VEL") == 0) {
        ibc = get_id(token,mod->nstring,"VEL String Error\n");
        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
        if (mod->flag.SW2_FLOW || mod->flag.NS2_FLOW) {
            mod->str_values[ibc].ol_flow.bc_flag = BCT_VEL_NEU;
            mod->str_values[ibc].ol_flow.isigma = iseries;
        } else if (mod->flag.SW3_FLOW || mod->flag.NS3_FLOW
#ifdef _ADH_GROUNDWATER
                   || mod->flag.GW_FLOW
#endif
                   ) {
            mod->str_values[ibc].flow.bc_flag = BCT_VEL_NEU;
            mod->str_values[ibc].flow.isigma = iseries;
        }
    } else if (strcmp(*token, "SOURCE") == 0) {
        ibc = get_id(token,mod->nstring,"SOURCE String Error\n");
        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
        if (mod->flag.SW2_FLOW) {
            mod->str_values[ibc].ol_flow.bc_flag = BCT_WATER_SOURCE;
            mod->str_values[ibc].ol_flow.isigma = iseries;
        }  else if (mod->flag.SW3_FLOW) {
            mod->str_values[ibc].flow.bc_flag = BCT_WATER_SOURCE;
            mod->str_values[ibc].flow.isigma = iseries;
        }
    } else if (strcmp(*token, "OVH") == 0) {
        ibc = get_id(token,mod->nstring,"OVH String Error\n");
        mod->str_values[ibc].ol_flow.bc_flag = BCT_OVH_NEU;
        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
        mod->str_values[ibc].ol_flow.ivx = iseries;
        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
        mod->str_values[ibc].ol_flow.ivy = iseries;
        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
        mod->str_values[ibc].ol_flow.isigma = iseries;
    } else if (strcmp(*token, "ZDG") == 0) {
        ibc = get_id(token,mod->nstring,"ZDG String Error\n");
        mod->str_values[ibc].ol_flow.bc_flag = BCT_ZERO_DEPTH_GRAD;
    } else if (strcmp(*token, "CD") == 0) {
        ibc = get_id(token,mod->nstring,"CD String Error\n");
        mod->str_values[ibc].ol_flow.bc_flag = BCT_CRITICAL_DEPTH;
    } else if (strcmp(*token, "BED") == 0) {
        ibc = get_id(token,mod->nstring,"BED String Error\n");
        mod->str_values[ibc].bed.bc_flag = BCT_BED;
        mod->str_values[ibc].flow.bc_flag = BCT_BED;
        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
        mod->str_values[ibc].flow.isigma = iseries;
    } else if (strcmp(*token, "TRN") == 0) {
        ibc = get_id(token,mod->nstring,"TRN String Error\n");
        itrn = get_transport_id(info, &subdata, mod->ntransport);
        mod->str_values[ibc].trans[itrn].bc_flag = BCT_NEU;
        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
        mod->str_values[ibc].trans[itrn].isigma = iseries;
    } else if (strcmp(*token, "WRS") == 0) {
#ifdef _ADH_STRUCTURES
        nwer = read_int_field (info, &subdata);
        nwer --;
        ibc = get_id(token,mod->nstring,"WRS String Error\n"); /* Reference string upstream of upstream weir face */
        ibc1 = ibc;
        mod->str_values[ibc].weir_num = nwer;
        mod->str_values[ibc].ol_flow.bc_flag = BCT_WRSU;
        mod->weir[nwer].egsu = ibc;
        ibc = get_id(token,mod->nstring,"WRS String Error\n"); /* THe reference string downstream of the downstream weir face */
        if (ibc != ibc1)
        {
            mod->str_values[ibc].ol_flow.bc_flag = BCT_WRSD;
            mod->str_values[ibc].weir_num = nwer;
        }
        mod->weir[nwer].egsd = ibc;
        ibc = get_id(token,mod->nstring,"WRS String Error\n"); /* Weir face on the upstream */
        ibc1 = ibc;
        mod->str_values[ibc].ol_flow.bc_flag = BCT_WEIRU;
        mod->str_values[ibc].weir_num = nwer;
        ibc = get_id(token,mod->nstring,"WRS String Error\n"); /* Weir face on the downstream */
        if (ibc1 != ibc)
        {
            mod->str_values[ibc].ol_flow.bc_flag = BCT_WEIRD;
            mod->str_values[ibc].weir_num = nwer;
        }
        mod->weir[nwer].w1 = get_next_token_dbl(token);
        mod->weir[nwer].z1 = get_next_token_dbl(token);
        mod->weir[nwer].hw = get_next_token_dbl(token);
#endif
    } else if (strcmp(*token, "FGT") == 0) {
#ifdef _ADH_STRUCTURES
        nflp = read_int_field (info, &subdata); /* Flap number */
        nflp--;
        printf("Flap Number is %d\n", nflp + 1);
        ibc = read_int_field (info, &subdata);
        mod->flap[nflp].calc = ibc;
        ibc =  get_id(token,mod->nstring,"FGT String Error\n");    /* Ref String upstream of Flap */
        ibc1 = ibc;
        mod->str_values[ibc].flap_num = nflp;
        mod->str_values[ibc].ol_flow.bc_flag = BCT_FLPU;
        /*str_values[ibc].phys_flag = OFF; */
        mod->flap[nflp].egsu = ibc;
        
        ibc =  get_id(token,mod->nstring,"FGT String Error\n");  /* Ref string Downstrean of flap */
        if (ibc != ibc1) {
            mod->str_values[ibc].flap_num = nflp;   /* String Downstream of flap */
            mod->str_values[ibc].ol_flow.bc_flag = BCT_FLPD;
            /*str_values[ibc].phys_flag = OFF; */
        }
        mod->flap[nflp].egsd = ibc;
        ibc =  get_id(token,mod->nstring,"FGT String Error\n");    /* The Flap String on the u/s */
        ibc1 = ibc;
        mod->str_values[ibc].ol_flow.bc_flag = BCT_FLAPU;
        mod->str_values[ibc].flap_num = nflp;
        /*str_values[ibc].phys_flag = OFF; */
        ibc = get_id(token,mod->nstring,"FGT String Error\n");   /* The Flap String on the d/s */
        if (ibc != ibc1) {
            mod->str_values[ibc].ol_flow.bc_flag = BCT_FLAPD;
            mod->str_values[ibc].flap_num = nflp;
            
            /*str_values[ibc].phys_flag = OFF; */ }
        if (mod->flap[nflp].calc == 1) {
            printf("The Flap flow-stage difference relationship must be in terms of total flow\n");
            mod->flap[nflp].A = get_next_token_dbl(token);
            mod->flap[nflp].B = get_next_token_dbl(token);
            mod->flap[nflp].C = get_next_token_dbl(token);
            mod->flap[nflp].D = get_next_token_dbl(token);
            mod->flap[nflp].E = get_next_token_dbl(token);
            mod->flap[nflp].F = get_next_token_dbl(token);
            mod->flap[nflp].GG = get_next_token_dbl(token);
        }
        else if (mod->flap[nflp].calc == 2) {
            mod->flap[nflp].density = get_next_token_dbl(token);
            mod->flap[nflp].length = get_next_token_dbl(token);
            mod->flap[nflp].width = get_next_token_dbl(token);
        }
#endif
    } else if (strcmp(*token, "SLS") == 0) {
#ifdef _ADH_STRUCTURES
        nslu = read_int_field (info, &subdata);
        nslu --;
        ibc = get_id(token,mod->nstring,"SLS String Error\n"); /* Reference string upstream of upstream Sluice face */
        ibc1 = ibc;
        mod->str_values[ibc].sluice_num = nslu;
        mod->str_values[ibc].ol_flow.bc_flag = BCT_SLSU;
        mod->sluice[nslu].egsu = ibc;
        ibc = get_id(token,mod->nstring,"SLS String Error\n"); /* THe reference string downstream of the downstream Sluice face */
        if (ibc != ibc1)
        {
            mod->str_values[ibc].ol_flow.bc_flag = BCT_SLSD;
            mod->str_values[ibc].sluice_num = nslu;
        }
        mod->sluice[nslu].egsd = ibc;
        ibc = get_id(token,mod->nstring,"SLS String Error\n"); /* Sluice face on the upstream */
        ibc1 = ibc;
        mod->str_values[ibc].ol_flow.bc_flag = BCT_SLUICEU;
        mod->str_values[ibc].sluice_num = nslu;
        ibc = get_id(token,mod->nstring,"SLS String Error\n"); /* Sluice face on the downstream */
        if (ibc1 != ibc)
        {
            mod->str_values[ibc].ol_flow.bc_flag = BCT_SLUICED;
            mod->str_values[ibc].sluice_num = nslu;
        }
        mod->sluice[nslu].a = get_next_token_dbl(token);
        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
        mod->sluice[nslu].opening = iseries;
#endif
    } else if (strcmp(*token, "FLXNML") == 0) {
        if (mod->flag.SW2_FLOW) {
            printf("\nWarning: Ignoring card NB FLXNML since this is a 2D SW model.");
        } else if (mod->flag.SW3_FLOW) {
            printf("\nUsing flux weighted normals for Neumann boundaries.");
            mod->flag.FLUX_WEIGHTED_NORMALS = ON;
        }
    }
}
