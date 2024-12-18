#ifndef H_SMAT_PHYSICS_
#define H_SMAT_PHYSICS_

typedef struct {
    int ntrns; // the number of transport constituents on this material
    char elemVarCode[4]; //codes read in by mesh, point to physics routines not var codes per say

    bool SW_FLOW;   // 1, 2, 3
    bool SW1_FLOW;  // 1
    bool SW2_FLOW;  // 2
    bool SW3_FLOW;  // 3
    bool NS_FLOW;   // 4
    bool NS3_FLOW;  // 4
    bool NS3_SPLIT; // 5
    bool DW_FLOW;   // 6
    bool WVEL_SPLIT;// 7
    bool PRESSURE;  // 8
    bool GW_FLOW;   // 9
    bool *TRANSPORT; // ntrans
    
    //convenitent to have total nvar to solve
    int nvar;
    //needs var codes
    int *vars;//int elem_vars[nvar];
    //similarly, how many physics routines to call
    int nSubmodels;
    //array of model structs, all they contain are pointers to physics routines
    SMODEL *model; // [nSubModels] length array
//    bool VORTICITY;
//    bool SEDIMENT;
//    bool SEDLIB;
//    bool ICM;
//    bool NSM;
//    bool WAVE;
//    bool WIND;

} SMAT_PHYSICS;

// methods
//void smat_physics_alloc_init_array(SMAT_PHYSICS **mat_physics, int nmat, int *ntrns);
void smat_physics_alloc_init_array(SMAT_PHYSICS **mat_physics, int nmat, int *ntrns, int *nvars, int *nSubMods, int **subMod_nvars);
void smat_physics_alloc_init(SMAT_PHYSICS *mat, int ntrns, int nvar, int nSubMods, int *nSubMod_nvar);
void smat_physics_free_array(SMAT_PHYSICS *mat, int nmat);
//void smat_physics_allocate_read(SMAT_PHYSICS **mat, SGRID *grid);
#endif
