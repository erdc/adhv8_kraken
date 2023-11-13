#include "global_header.h"

void comm_update_swaves(SGRID *grid, SWAVE *waves, int cstorm_flag) {

#ifdef _MESSG
    int i,j,k;
    
    // determine number of doubles to pack
    int num_doubles = -1;
    switch (cstorm_flag) {
        case 7:
            num_doubles = 3;
            break;
        case 12:
            num_doubles = 8;
            tl_error("c2shore CStorm option not supported yet in comm_update_swaves");
            break;
        case 14:
            num_doubles = 2;
            break;
        default:
            tl_error("using a bad cstorm wsid flag in comm_update_swaves");
            break;
    }
    
    // pack it
    double *packed = (double *)tl_alloc(sizeof(double), grid->nnodes * num_doubles);
    for(j=0, i=0; i<grid->nnodes; i++) {
        switch (cstorm_flag) {
            case 7:
                packed[j++] = waves[i].rads.xx;
                packed[j++] = waves[i].rads.xy;
                packed[j++] = waves[i].rads.yy;
                break;
            case 12:
                tl_error("C2shore CStorm option not supported yet in comm_update_swaves");
                break;
            case 14:
                packed[j++] = waves[i].stress.x;
                packed[j++] = waves[i].stress.y;
                break;
            default:
                tl_error("using a bad cstorm wsid flag in comm_update_swaves");
                break;
        }
    }
    
    // update core arrays
    comm_update_double(packed, num_doubles, grid->smpi);
  
    // unpack it
    for(j=0, i=0; i<grid->nnodes; i++) {
        switch (cstorm_flag) {
            case 7:
                waves[i].rads.xx = packed[j++];
                waves[i].rads.xy = packed[j++];
                waves[i].rads.yy = packed[j++];
                break;
            case 12:
                tl_error("C2shore CStorm option not supported yet in comm_update_swaves");
                break;
            case 14:
                waves[i].stress.x = packed[j++];
                waves[i].stress.y = packed[j++];
                break;
            default:
                tl_error("using a bad cstorm wsid flag in comm_update_swaves");
                break;
        }
    }
    
    packed = (double *)tl_free(sizeof(double), grid->nnodes * num_doubles, packed);
#endif
}
