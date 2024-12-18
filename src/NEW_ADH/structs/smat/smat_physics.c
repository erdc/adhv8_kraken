//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <stdbool.h>
//#include <assert.h>
//#include "smat_physics.h"
//#include "debug.h"
//#include "header_tl_alloc.h"
//#include "define.h"
//#include "tokens.h"
//#include "sgrid.h"

#include "adh.h"

static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes physics material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat_physics (SMAT_PHYSICS**)  double pointer to a physics material
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smat_physics_alloc_init_array(SMAT_PHYSICS **mat_physics, int nmat, int *ntrns, int *nvars, int *nSubMods, int **subMod_nvars) {
    
    assert(nmat>=0);
    if (nmat == 0){
        mat_physics = NULL;
        return;
    }
    else{
        //allocate array of mat physics
        (*mat_physics) = (SMAT_PHYSICS *) tl_alloc(sizeof(SMAT_PHYSICS), nmat);
        int imat;
        SMAT_PHYSICS *mat; // for alias
    
        for (imat=0; imat<nmat; imat++) {
            mat = &(*mat_physics)[imat]; // alias
            smat_physics_alloc_init(mat, ntrns[imat], nvars[imat], nSubMods[imat], subMod_nvars[imat]);
        }
        return;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes physics material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat_physics (SMAT_PHYSICS**)  double pointer to a physics material
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smat_physics_alloc_init(SMAT_PHYSICS *mat, int ntrns, int nvar, int nSubMods, int *nSubMod_nvar){
    int i;
    assert(nvar>0);
    //allocate and initialize each mat physics
    mat->ntrns = ntrns;
    mat->nvar = nvar;
    mat->nSubmodels = nSubMods;

    if (ntrns > 0){
        (mat->TRANSPORT) = (bool *) tl_alloc(sizeof(bool), ntrns);
        for(i=0;i<ntrns;i++){
            mat->TRANSPORT[i] = false;
        }
    }else{
        mat->TRANSPORT = NULL;
    }
    mat->SW_FLOW = OFF;
    mat->SW1_FLOW = OFF;
    mat->SW2_FLOW = OFF;
    mat->SW3_FLOW = OFF;
    mat->NS_FLOW = OFF;
    mat->NS3_FLOW = OFF;
    mat->NS3_SPLIT = OFF;
    mat->GW_FLOW = OFF;
    mat->DW_FLOW = OFF;
    mat->WVEL_SPLIT = OFF;
    mat->PRESSURE = OFF;
    (mat->vars) = (int *) tl_alloc(sizeof(int), nvar);
    //initalize to 0, or different value
    for(i=0;i<nvar;i++){
        mat->vars[i] = 0;
    }
    //allocate and initialize array of elem_phyiscs
    if(nSubMods>0){
        smodel_alloc_init_array(&(mat->model), nSubMods, nSubMod_nvar); 
    }else{
        mat->model = NULL;
    }

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees physics material properties
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] mat (SMAT_PHYSICS *)  pointer to a transport material
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smat_physics_free_array(SMAT_PHYSICS *mat, int nmat) {
    
    int imat;
    for (imat=0; imat<nmat; imat++) {
        if(mat[imat].TRANSPORT!=NULL){
            mat[imat].TRANSPORT = (bool *) tl_free(sizeof(bool), mat[imat].ntrns, mat[imat].TRANSPORT);
        }
        if(mat[imat].vars!=NULL){
            mat[imat].vars = (int *) tl_free(sizeof(int), mat[imat].nvar, mat[imat].vars);
        }
        if(mat[imat].model!=NULL){
            smodel_free_array(mat[imat].model, mat[imat].nSubmodels); 
        }

    }
    mat = (SMAT_PHYSICS *) tl_free(sizeof(SMAT_PHYSICS), nmat, mat);
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and reads physics material property files
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout]
 *
 * \note Requires allocated SGRID
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*
void smat_physics_allocate_read(SMAT_PHYSICS **mat, SGRID *grid) {
    
    int imat, itrns, ielem, nmat, flag = 0;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *token;
    
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Open and read physics material file
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    FILE *fp = fopen("physics.mat", "r");
    if (fp == 0) {
        printf("WARNING: no physics.mat input file found.\n");
        return;
    }
    
    //Mark, weird build error

    // read through the element list first for counting materials
    nmat = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        get_compare_token_increment(line,&token,"MAT",&nmat);
    }
    printf("physics.mat read: total # of materials found: %d\n",nmat);
    rewind(fp);

    // now count transport materials
    int ntrns[nmat]; for (imat=0; imat<nmat; imat++) {ntrns[imat] = 0;}
    while ((read = getline(&line, &len, fp)) != -1) {
        if (get_compare_token(line,&token,"MAT") == 1) {
            imat = get_next_token_int(&token);
            while (token != NULL) {
                get_compare_token_increment(line,&token,"TRN",&ntrns[imat-1]);
            }
        }
    }
    for (imat=0; imat<nmat; imat++) {
        printf("physics.mat read: total # of transport materials found for mat[%d]: %d\n",imat,ntrns[imat]);
    }
    rewind(fp);

    // allocate physics materials
    smat_physics_alloc_init(mat, nmat, ntrns);
    SMAT_PHYSICS *mat_phys = *mat; // alias

    // read and store hydrodynamic material physics headers
    while ((read = getline(&line, &len, fp)) != -1) {
        if (get_compare_token(line,&token,"MAT") == 1) {
            imat = get_next_token_int(&token) - 1; // the element ID
            flag = 0;
            while (1) {
                get_next_token(&token);
                //printf("token: %s\n",token);
                if (token == NULL) break;
                //printf("%d\n",strcmp(token, "DW"));
                if (strcmp(token, "SW1") == 0) {
                    mat_phys[imat].SW_FLOW = TRUE;
                    mat_phys[imat].SW1_FLOW = TRUE;
                    strcpy(&mat_phys[imat].elemVarCode[0],"1");
                    flag = 1;
                    
                } else if (strcmp(token, "SW2") == 0) {
                    mat_phys[imat].SW_FLOW = TRUE;
                    mat_phys[imat].SW2_FLOW = TRUE;
                    strcpy(&mat_phys[imat].elemVarCode[0],"2");
                    flag = 1;
                    
                } else if (strcmp(token, "SW3") == 0) {
                    mat_phys[imat].SW_FLOW = TRUE;
                    mat_phys[imat].SW3_FLOW = TRUE;
                    strcpy(&mat_phys[imat].elemVarCode[0],"3");
                    flag = 1;
                    
                } else if (strcmp(token, "NS3") == 0) {
                    mat_phys[imat].NS_FLOW = TRUE;
                    mat_phys[imat].NS3_FLOW = TRUE;
                    strcpy(&mat_phys[imat].elemVarCode[0],"4");
                    flag = 1;
                    
                } else if (strcmp(token, "NS3_SPLIT") == 0) {
                    mat_phys[imat].NS_FLOW = TRUE;
                    mat_phys[imat].NS3_SPLIT = TRUE;
                    strcpy(&mat_phys[imat].elemVarCode[0],"5");
                    flag = 1;
                    
                } else if (strcmp(token, "DW\n") == 0) {
                    mat_phys[imat].DW_FLOW = TRUE;
                    strcpy(&mat_phys[imat].elemVarCode[0],"6");
                    flag = 1;
                    
                } else if (strcmp(token, "WVEL") == 0) {
                    mat_phys[imat].WVEL_SPLIT = TRUE;
                    strcpy(&mat_phys[imat].elemVarCode[0],"7");
                    flag = 1;
                    
                } else if (strcmp(token, "PRESSURE") == 0) {
                    mat_phys[imat].PRESSURE = TRUE;
                    strcpy(&mat_phys[imat].elemVarCode[0],"8");
                    flag = 1;
                    
                }
            }
            if (flag == 0) strcpy(&mat_phys[imat].elemVarCode[0],"0");
        }
    }
    rewind(fp);

//    for (imat=0; imat<nmat; imat++) {
//        printf("mat: %d varcode: %s\n",imat,mat_phys[imat].elemVarCode);
//    }

    // read and store groundwater material physics headers
    while ((read = getline(&line, &len, fp)) != -1) {
        if (get_compare_token(line,&token,"MAT") == 1) {
            imat = get_next_token_int(&token) - 1; // the element ID
            flag = 0;
            while (1) {
                get_next_token(&token);
                if (token == NULL) break;
                if (strcmp(token, "GW") == 0) {
                    mat_phys[imat].GW_FLOW = TRUE;
                    strcpy(&mat_phys[imat].elemVarCode[1],"1");
                    flag = 1;
                }
            }
            if (flag == 0) strcpy(&mat_phys[imat].elemVarCode[1],"0");
        }
    }
    rewind(fp);
    
//    for (imat=0; imat<nmat; imat++) {
//        printf("mat: %d varcode: %s\n",imat,mat_phys[imat].elemVarCode);
//    }

    // read and store transport material physics headers
    while ((read = getline(&line, &len, fp)) != -1) {
        if (get_compare_token(line,&token,"MAT") == 1) {
            imat = get_next_token_int(&token) - 1; // the element ID
            flag = 0;
            while (1) {
                get_next_token(&token);
                if (token == NULL) break;
                if (strcmp(token, "TRN") == 0) {
                    itrns = get_next_token_int(&token) - 1;
                    mat_phys[imat].TRANSPORT[itrns] = TRUE;
                    flag = 1;
                }
            }
            if (flag == 0) {
                //strcat(varCode[imat], "0");
                strcpy(&mat_phys[imat].elemVarCode[2],"0");
            } else {
                char tmp[2];
                sprintf(tmp, "%d", ntrns[imat]);
                //strcat(varCode[imat], tmp);
                strcpy(&mat_phys[imat].elemVarCode[2],tmp);
            }
        }
    }
    fclose(fp);
    
    for (imat=0; imat<nmat; imat++) {
        printf("mat: %d varcode: %s\n",imat,mat_phys[imat].elemVarCode);
    }
    

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    // read and store element material IDs
    int ie3d_count = 0, ie2d_count = 0, ie1d_count = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token); // the element type
        if (strcmp(token, "TET") == 0 || strcmp(token, "PRISM") == 0) {
            ie3d_count++;
            ielem = get_next_token_int(&token)-1; // the element ID
            imat  = get_next_token_int(&token)-1; // the element ID
            assert(ie3d_count < grid->nelems3d);
            assert(ielem > 0);
            assert(ielem < grid->nelems3d);
            assert(imat > 0);
            assert(imat < nmat);
            grid->elem3d[ielem].mat = imat;
        } else if (strcmp(token, "TRI") == 0 || strcmp(token, "QUAD") == 0) {
            ie2d_count++;
            ielem = get_next_token_int(&token)-1; // the element ID
            imat  = get_next_token_int(&token)-1; // the element ID
            assert(ie2d_count < grid->nelems2d);
            assert(ielem > 0);
            assert(ielem < grid->nelems2d);
            assert(imat > 0);
            assert(imat < nmat);
            grid->elem2d[ielem].mat = imat;
        } else if (strcmp(token, "SEG") == 0) {
            ie1d_count++;
            ielem = get_next_token_int(&token)-1; // the element ID
            imat  = get_next_token_int(&token)-1; // the element ID
            assert(ie1d_count < grid->nelems1d);
            assert(ielem > 0);
            assert(ielem < grid->nelems1d);
            assert(imat > 0);
            assert(imat < nmat);
            grid->elem1d[ielem].mat = imat;
        }
    }
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++

}
*/

