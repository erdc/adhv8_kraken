
#include "test.h"
// need to include SGRID in this to test - doesn't work
int main(int argc, char *argv[]) {
    
    SMAT_PHYSICS *mat;
    int elem3d_mat_ids[10];
    int elem2d_mat_ids[5];
    int elem1d_mat_ids[2];
    
    smat_physics_allocate_read(&mat,elem3d_mat_ids,elem2d_mat_ids,elem1d_mat_ids);
    
    return 0;
}
