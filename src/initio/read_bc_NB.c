#include "global_header.h"

void read_bc_NB(SMODEL *mod, char *data) {
    
    char line[MAXLINE];           /* the input line */
    char *subdata = NULL;         /* the data after the second card   is read */
    char *subsubdata = NULL;
    
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
    switch (parse_card(data, &subdata)) {
            
        case CARD_OF:
            ibc = get_string_id(info, &subdata, nstring);
            if (mod->flag.SW2_FLOW || mod->flag.NS2_FLOW) {
                mod->str_values[ibc].ol_flow.bc_flag = BCT_OUTFLOW;
            } else if (mod->flag.SW3_FLOW || mod->flag.NS3_FLOW) {
                mod->str_values[ibc].flow.bc_flag = BCT_OUTFLOW;
            }
            break;
        case CARD_DIS:
            ibc = get_string_id(info, &subdata, nstring);
            iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
            if (mod->flag.SW2_FLOW || mod->flag.NS2_FLOW) {
                mod->str_values[ibc].ol_flow.bc_flag = BCT_DIS_NEU;
                mod->str_values[ibc].ol_flow.isigma = iseries;
                mod->flag.CONVEYANCE = YES;
            }  else if (mod->flag.SW3_FLOW || mod->flag.NS3_FLOW) {
                mod->str_values[ibc].flow.bc_flag = BCT_DIS_NEU;
                mod->str_values[ibc].flow.isigma = iseries;
            }
            break;
            // cjt || FLUX_TAG
        case CARD_CPL:
            ibc = get_string_id(info, &subdata, nstring);
            if (mod->flag.SW2_FLOW) {
                mod->str_values[ibc].ol_flow.bc_flag = BCT_FLUX_COUPLE;
            } else {
                mod->str_values[ibc].flow.bc_flag = BCT_FLUX_COUPLE;
            }
            break;
        case CARD_PRS:
            ibc = get_string_id(info, &subdata, nstring);
            mod->str_values[ibc].pressure.bc_flag = BCT_PRS_NEU;
            iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
            mod->str_values[ibc].pressure.isigma = iseries;
            break;
        case CARD_FRS:
            ibc = get_string_id(info, &subdata, nstring);
            mod->str_values[ibc].flow.bc_flag = BCT_FRS;
            iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
            mod->str_values[ibc].flow.isigma = iseries;
            break;
        case CARD_OTW:
            ibc = get_string_id(info, &subdata, nstring);
            mod->str_values[ibc].ol_flow.bc_flag = BCT_PRS_NEU;
            iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
            mod->str_values[ibc].ol_flow.isigma = iseries;
            break;
        case CARD_VEL:
            ibc = get_string_id(info, &subdata, nstring);
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
            break;
        case CARD_SOURCE:
            ibc = get_string_id(info, &subdata, nstring);
            iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
            if (mod->flag.SW2_FLOW) {
                mod->str_values[ibc].ol_flow.bc_flag = BCT_WATER_SOURCE;
                mod->str_values[ibc].ol_flow.isigma = iseries;
            }  else if (mod->flag.SW3_FLOW) {
                mod->str_values[ibc].flow.bc_flag = BCT_WATER_SOURCE;
                mod->str_values[ibc].flow.isigma = iseries;
            }
            break;
        case CARD_OVH:
            ibc = get_string_id(info, &subdata, nstring);
            mod->str_values[ibc].ol_flow.bc_flag = BCT_OVH_NEU;
            iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
            mod->str_values[ibc].ol_flow.ivx = iseries;
            iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
            mod->str_values[ibc].ol_flow.ivy = iseries;
            iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
            mod->str_values[ibc].ol_flow.isigma = iseries;
            break;
        case CARD_BED:
            ibc = get_string_id(info, &subdata, nstring);
            mod->str_values[ibc].bed.bc_flag = BCT_BED;
            mod->str_values[ibc].flow.bc_flag = BCT_BED;
            iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
            mod->str_values[ibc].flow.isigma = iseries;
            break;
        case CARD_TRN:
            ibc = get_string_id(info, &subdata, nstring);
            itrn = get_transport_id(info, &subdata, mod->ntransport);
            mod->str_values[ibc].trans[itrn].bc_flag = BCT_NEU;
            iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
            mod->str_values[ibc].trans[itrn].isigma = iseries;
            break;
        case CARD_WRS:
#ifdef _ADH_STRUCTURES
            nwer = read_int_field (info, &subdata);
            nwer --;
            ibc = get_string_id(info, &subdata, nstring); /* Reference string upstream of upstream weir face */
            ibc1 = ibc;
            mod->str_values[ibc].weir_num = nwer;
            mod->str_values[ibc].ol_flow.bc_flag = BCT_WRSU;
            mod->weir[nwer].egsu = ibc;
            ibc = get_string_id(info, &subdata, nstring); /* THe reference string downstream of the downstream weir face */
            if (ibc != ibc1)
            {
                mod->str_values[ibc].ol_flow.bc_flag = BCT_WRSD;
                mod->str_values[ibc].weir_num = nwer;
            }
            mod->weir[nwer].egsd = ibc;
            ibc = get_string_id(info, &subdata, nstring); /* Weir face on the upstream */
            ibc1 = ibc;
            mod->str_values[ibc].ol_flow.bc_flag = BCT_WEIRU;
            mod->str_values[ibc].weir_num = nwer;
            ibc = get_string_id(info, &subdata, nstring); /* Weir face on the downstream */
            if (ibc1 != ibc)
            {
                mod->str_values[ibc].ol_flow.bc_flag = BCT_WEIRD;
                mod->str_values[ibc].weir_num = nwer;
            }
            mod->weir[nwer].w1 = read_dbl_field(info, &subdata);
            mod->weir[nwer].z1 = read_dbl_field(info, &subdata);
            mod->weir[nwer].hw = read_dbl_field(info, &subdata);
#endif
            break;
        case CARD_FGT:
#ifdef _ADH_STRUCTURES
            nflp = read_int_field (info, &subdata); /* Flap number */
            nflp--;
            printf("Flap Number is %d\n", nflp + 1);
            ibc = read_int_field (info, &subdata);
            mod->flap[nflp].calc = ibc;
            ibc =  get_string_id(info, &subdata, nstring);    /* Ref String upstream of Flap */
            ibc1 = ibc;
            mod->str_values[ibc].flap_num = nflp;
            mod->str_values[ibc].ol_flow.bc_flag = BCT_FLPU;
            /*str_values[ibc].phys_flag = OFF; */
            mod->flap[nflp].egsu = ibc;
            
            ibc =  get_string_id(info, &subdata, nstring);  /* Ref string Downstrean of flap */
            if (ibc != ibc1) {
                mod->str_values[ibc].flap_num = nflp;   /* String Downstream of flap */
                mod->str_values[ibc].ol_flow.bc_flag = BCT_FLPD;
                /*str_values[ibc].phys_flag = OFF; */
            }
            mod->flap[nflp].egsd = ibc;
            ibc =  get_string_id(info, &subdata, nstring);    /* The Flap String on the u/s */
            ibc1 = ibc;
            mod->str_values[ibc].ol_flow.bc_flag = BCT_FLAPU;
            mod->str_values[ibc].flap_num = nflp;
            /*str_values[ibc].phys_flag = OFF; */
            ibc = get_string_id(info, &subdata, nstring);   /* The Flap String on the d/s */
            if (ibc != ibc1) {
                mod->str_values[ibc].ol_flow.bc_flag = BCT_FLAPD;
                mod->str_values[ibc].flap_num = nflp;
                
                /*str_values[ibc].phys_flag = OFF; */ }
            if (mod->flap[nflp].calc == 1) {
                printf("The Flap flow-stage difference relationship must be in terms of total flow\n");
                mod->flap[nflp].A = read_dbl_field(info, &subdata);
                mod->flap[nflp].B = read_dbl_field(info, &subdata);
                mod->flap[nflp].C = read_dbl_field(info, &subdata);
                mod->flap[nflp].D = read_dbl_field(info, &subdata);
                mod->flap[nflp].E = read_dbl_field(info, &subdata);
                mod->flap[nflp].F = read_dbl_field(info, &subdata);
                mod->flap[nflp].GG = read_dbl_field(info, &subdata);
            }
            else if (mod->flap[nflp].calc == 2) {
                mod->flap[nflp].density = read_dbl_field(info, &subdata);
                mod->flap[nflp].length = read_dbl_field(info, &subdata);
                mod->flap[nflp].width = read_dbl_field(info, &subdata);
            }
#endif
            break;
            
        case CARD_SLS:
#ifdef _ADH_STRUCTURES
            nslu = read_int_field (info, &subdata);
            nslu --;
            ibc = get_string_id(info, &subdata, nstring); /* Reference string upstream of upstream Sluice face */
            ibc1 = ibc;
            mod->str_values[ibc].sluice_num = nslu;
            mod->str_values[ibc].ol_flow.bc_flag = BCT_SLSU;
            mod->sluice[nslu].egsu = ibc;
            ibc = get_string_id(info, &subdata, nstring); /* THe reference string downstream of the downstream Sluice face */
            if (ibc != ibc1)
            {
                mod->str_values[ibc].ol_flow.bc_flag = BCT_SLSD;
                mod->str_values[ibc].sluice_num = nslu;
            }
            mod->sluice[nslu].egsd = ibc;
            ibc = get_string_id(info, &subdata, nstring); /* Sluice face on the upstream */
            ibc1 = ibc;
            mod->str_values[ibc].ol_flow.bc_flag = BCT_SLUICEU;
            mod->str_values[ibc].sluice_num = nslu;
            ibc = get_string_id(info, &subdata, nstring); /* Sluice face on the downstream */
            if (ibc1 != ibc)
            {
                mod->str_values[ibc].ol_flow.bc_flag = BCT_SLUICED;
                mod->str_values[ibc].sluice_num = nslu;
            }
            mod->sluice[nslu].a = read_dbl_field(info, &subdata);
            iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
            mod->sluice[nslu].opening = iseries;
#endif
            break;
        case CARD_FLXNML:
            if (mod->flag.SW2_FLOW) {
                printf("\nWarning: Ignoring card NB FLXNML since this is a 2D SW model.");
            } else if (mod->flag.SW3_FLOW) {
                printf("\nUsing flux weighted normals for Neumann boundaries.");
                mod->flag.FLUX_WEIGHTED_NORMALS = ON;
            }
            break;
        default:
            break;
            
    }
}
