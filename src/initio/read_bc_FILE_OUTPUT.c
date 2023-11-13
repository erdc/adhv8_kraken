#include "global_header.h"

// determines which optional files are output

void read_bc_FILE_OUTPUT(SMODEL * mod, char *data) {
    
    char line[MAXLINE];         /* the input line */
    char *subdata = NULL;       /* the data after the second card is read */
    char *subsubdata = NULL;
    
    switch (parse_card(data, &subdata)) {
        case CARD_BEDVEL:
            mod->file_output.bed_velocity = ON;
            break;
            
        case CARD_SURVEL:
            mod->file_output.surface_velocity = ON;
            break;
            
        case CARD_AVGVEL:
            mod->file_output.depth_avg_velocity = ON;
            break;
            
        case CARD_PRS:
            mod->file_output.pressure = ON;
            break;
            
        case CARD_GSPEED:
            mod->file_output.grid_speed = ON;
            break;
            
        case CARD_WIND:
            mod->file_output.wind = ON;
            break;
            
        case CARD_WAVE:
            mod->file_output.wave = ON;
            break;
            
        case CARD_VIS:
            mod->file_output.hyd_vis = ON;
            break;
            
        case CARD_DIF:
            mod->file_output.trn_dif = ON;
            break;
            
        case CARD_CHOP:
            mod->file_output.chop = ON;
            break;
            
        case CARD_ADAPT:
            switch (parse_card(subdata, &subsubdata)) {
                case CARD_GRID:
                    mod->file_output.adaption = ON;
                    mod->file_output.adapt_grid = ON;
                    break;
                case CARD_SW:
                    mod->file_output.adaption = ON;
                    mod->file_output.adapt_grid = ON;
                    mod->file_output.adapt_sw = ON;
                    break;
                case CARD_NS:
                    mod->file_output.adaption = ON;
                    mod->file_output.adapt_grid = ON;
                    mod->file_output.adapt_ns = ON;
                    break;
                case CARD_CON:
                    mod->file_output.adaption = ON;
                    mod->file_output.adapt_grid = ON;
                    mod->file_output.adapt_con = ON;
                    break;
                case CARD_SED:
                    mod->file_output.adaption = ON;
                    mod->file_output.adapt_grid = ON;
                    mod->file_output.adapt_sed = ON;
                    break;
#ifdef _ADH_GROUNDWATER
                case CARD_GW:
                    mod->file_output.adaption = ON;
                    mod->file_output.adapt_grid = ON;
                    mod->file_output.adapt_gw = ON;
                    break;
#endif
                default:
                    mod->file_output.adaption = ON;
                    mod->file_output.adapt_grid = ON;
                    mod->file_output.adapt_sw = ON;
                    mod->file_output.adapt_con = ON;
                    mod->file_output.adapt_sed = ON;
                    break;
            }
            break;
            
        case CARD_GRID:
            mod->file_output.grid2dm = ON;
            break;
            
        case CARD_ALL:
            mod->file_output.grid2dm = ON;
            mod->file_output.bed_velocity = ON;
            mod->file_output.surface_velocity = ON;
            mod->file_output.depth_avg_velocity = ON;
            mod->file_output.pressure = ON;
            mod->file_output.grid_speed = ON;
            mod->file_output.wind = ON;
            mod->file_output.wave = ON;
            mod->file_output.hyd_vis = ON;
            mod->file_output.trn_dif = ON;
            break;
            
        default:
            break;
    }
}
