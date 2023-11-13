
#include "global_header.h"

void read_bc_OP(SMODEL *mod, char *data) {
    
    char line[MAXLINE];           /* the input line */
    char *subdata = NULL;         /* the data after the second card is read */
    char *subsubdata = NULL;
    
    SIO info = *(mod->io);
    
    switch (parse_card(data, &subdata)) {
            
        case CARD_TRN:
            break;
            
#ifdef _SEDLIB
        case CARD_SEDLIB:
            read_bc_SEDLIB(mod, data);
            break;
#endif
            
        case CARD_WQNSM:
            mod->flag.NSM = TRUE;
            printf("WQ Simulations with NSM active\n");
            read_bc_NSM(mod, data);
            break;
            
        case CARD_INC:
            mod->nalloc_inc = read_int_field(info, &subdata);
            break;
            
        case CARD_BLK:
            mod->nblock = read_int_field(info, &subdata);
            if (mod->nblock < 1) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error("BLK card value has to be positive.\n");
            }
            break;
        case CARD_EOS:
            mod->flag.EOS = read_int_field(info, &subdata); /* This reads the type of EOS to use in AdH Computations */
            printf("AdH will use Equation of State %d (0: Linearlized 1: Full Equation)\n", mod->flag.EOS);
            break;
        case CARD_SW2:
            mod->flag.SW_FLOW = ON;
            mod->flag.SW2_FLOW = ON;
            if (mod->ntransport > 0)
                mod->flag.SW2_TRANSPORT = ON;
            break;
        case CARD_DIF:
            mod->flag.SW_FLOW = ON;
            mod->flag.SW2_FLOW = ON;
            mod->flag.DIFFUSIVE_WAVE = ON;
            printf("Diffusive Wave Equation Set Is Active\n");
            printf("Velocities will be computed based on the water surface slope. Equation set is 1 DOF\n");
            break;
        case CARD_SW3:
            mod->flag.SW_FLOW = ON;
            mod->flag.SW3_FLOW = ON;
            mod->flag.MG = ON;
            if (mod->ntransport > 0)
                mod->flag.SW3_TRANSPORT = ON;
            break;
        case CARD_NS3:
			mod->flag.NS_FLOW = ON;
            mod->flag.NS3_FLOW = ON;
            if (mod->ntransport > 0)
                mod->flag.NS3_TRANSPORT = ON;
            break; 
        case CARD_NS2:
            mod->flag.NS_FLOW = ON;
			mod->flag.NS2_FLOW = ON;
            if (mod->ntransport > 0)
                mod->flag.NS2_TRANSPORT = ON;
            break;

        case CARD_WAV:
            printf("\n Short waves are turned on during this simulation. \n");
            printf("-- Wave data can be initialized through the hotstart file or given via CSTORM.\n");
            mod->flag.WAVE = ON;
            break;
            
        case CARD_WND:
            printf("\n Wind stressing is turned on during this simulation. \n");
            mod->flag.WIND = ON;
            break;
            
        case CARD_TEM:
            mod->tau_temporal = read_dbl_field(info, &subdata);
            if (mod->tau_temporal < 0 || mod->tau_temporal > 1) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error("Tau_temporal must be between 0 and 1.");
            }
            break;
            
        case CARD_TPG:
            mod->tau_pg = read_dbl_field(info, &subdata);
            if (mod->tau_pg < 0 || mod->tau_pg > 0.5) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error("Tau_pg must be between 0 and .5");
            }
            break;
            
        case CARD_PRE:
            mod->solver_info.prec_value = read_int_field(info, &subdata);
            break;
            
        case CARD_MG:
            mod->flag.MG = ON;
            break;
            
            
            //////////////////////////////////////////////////////////////////////////////////////////////////
            // FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //
        case CARD_WNDLIB:
#ifdef WINDLIB
        {
            if (mod->flag.WIND == OFF) {
                printf("\n Wind stressing must be turned on using OP WND in *.bc file BEFORE using card OP WNDLIB. Exiting program. \n");
                exit(0);
            }
            mod->flag.WIND_LIBRARY = ON;
            printf("\n Wind library will be used in simulation. Reading wind library input. \n");
            read_bc_WINDLIB(mod, data);
            swindlib_printScreen(mod->windlib);
            printf("Finished reading user wind library input\n");
        }
#else
        {
            printf("\n Attempted to access wind library. No wind library installed! Exiting program. \n");
            exit(0);
        }
#endif
            break;
            // ABOVE LINES ADDED BY GAJANAN                                                                 //
            //////////////////////////////////////////////////////////////////////////////////////////////////
            
            
        default:
            break;
            
    }
}
