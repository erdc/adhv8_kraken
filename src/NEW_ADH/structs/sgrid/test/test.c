#include "global_header.h"

int main(int argc, char *argv[]) {
    int ierror = 0;
    
#ifdef _MESSG
    MPI_Comm ADH_COMM;
    ADH_COMM = MPI_COMM_WORLD;
#endif
    
#ifndef _MESSG
    debug_initialize();
#else
    ierror = MPI_Init(&argc, &argv);
    debug_initialize(ADH_COMM); // uses COMM_DUPLICATE
#endif
    
    DESIGN_MODEL *dmod;
    int nSuperModels = 2;
    int nSubModels[2];
    nSubModels[0] = 1;
    nSubModels[1] = 3;
    
    // create a design model with 2 supermodels
    design_model_alloc_init(&(dmod),nSuperModels);
    design_model_printScreen(&(dmod[0]));
    
#ifndef _MESSG
    sgrid_read(&(dmod->grid), "test");
#else
    sgrid_read(&(dmod->grid), "test", ADH_COMM);
#endif
    sgrid_printScreen(dmod[0].grid);

    //call grid writer
    init_hdf5_file(dmod->grid);
    sgrid_write_hdf5(dmod->grid);
    sgrid_write_xdmf(dmod->grid);
    sgrid_write_nodal_pe(dmod->grid);
    sgrid_write_elemental_pe(dmod->grid);
    sgrid_write_xdmf_nodal_pe(dmod->grid);
    sgrid_write_xdmf_elemental_pe(dmod->grid);
    sgrid_read_nodal_attribute(dmod->grid);
    
#ifdef _MESSG
    MPI_Finalize();
#endif
    
    return 0;
}
    
