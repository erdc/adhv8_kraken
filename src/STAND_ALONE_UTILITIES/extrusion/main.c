#include "extrusion.h"

int main(int argc, char *argv[]) {
    
    int i;
    int mesh_type = TET_MESH;
	int flag_bins = OFF;
	int flag_min_layer = OFF;
    int minimum_layers = 1;
    
    //------------------------------------------------------------------
    //--- COMMAND LINE INFO --------------------------------------------
    
    for (i=0; i<argc; i++) {
        printf("%s: ",argv[i]);
        if (strcmp(argv[i], "-mixed") == AGREE) {mesh_type = MIXED_ELEMENT_MESH;}
        if (strcmp(argv[i], "-bin") == AGREE) {
            printf("\n LAYER BIN OPTION USED.\n");
            flag_bins = ON;
        }
        if (strcmp(argv[i], "-minLayer") == AGREE) {
            printf("\n MINIMUM LAYER OPTION USED.\n");
            flag_min_layer = ON;
            minimum_layers = atoi(argv[i+1]);
            i++;
            if (minimum_layers <= 0) {
                fprintf(stderr, "Error :: minimum_#_bed_layers must be >= 0\n");
                exit(-1);
            }
        }
        
    }
    
    if (flag_min_layer == ON && flag_bins == ON) {
        fprintf(stderr, "Error :: Cannot use the minimum layer and bin options simultaneously.\n");
        exit(-1);
    }
    if (flag_min_layer == OFF && flag_bins == OFF) {
        fprintf(stderr, "Error :: Must use either -bin or -minLayer options.\n");
        exit(-1);
    }

    char adh_root[MAXLINE];
    strcpy(adh_root, argv[1]);
 
    int ERROR = extrudeAdH(flag_bins,flag_min_layer,minimum_layers,adh_root,mesh_type);
    return ERROR;
}
