#include "global_header.h"
#include "xdmf_generator.h"
#ifdef _ADH_HDF5
#include "adh_hdf5.h"
#endif

int main(int argc, char *argv[]) {
 
    if (argc < 2 || argc > 3) {
        fprintf(stderr, "\n Missing argument.\n Usage:\n" "   To convert 3dm style files to XDMF\n     xdmf_convert file_base \n");
        exit(0);
    }

    char adh_root[MAXLINE];
    strcpy(adh_root, argv[1]);
    
    geo_xdmf_write(adh_root);

    return 0;
}

