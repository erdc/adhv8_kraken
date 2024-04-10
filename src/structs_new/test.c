#include "local_header.h"

int main(int argc, char *argv[]) {
    
    DESIGN_MODEL *dmod;
    int nSuperModels = 2;
    int nSubModels[2];
    nSubModels[0] = 1;
    nSubModels[1] = 3;
    
    // create a design model with 2 supermodels
    design_model_alloc_init(&(dmod),nSuperModels);
    design_model_printScreen(&(dmod[0]));
    
    sgrid_read(&(dmod->grid), NULL, NULL);
    sgrid_printScreen(dmod[0].grid);
    
    return 0;
}
    
