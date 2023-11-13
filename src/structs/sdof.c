#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "sdof.h"

static int printFieldWidth = 30;
static int printPrecision  = 20;

/***********************************************************/
/***********************************************************/
/***********************************************************/

// initialization routines

void dof2_init(DOF_2 *dof) { dof->x_eq = 0.; dof->y_eq = 0.;}
void dof3_init(DOF_3 *dof) { dof->x_eq = 0.; dof->y_eq = 0.; dof->c_eq = 0.;}
void dof4_init(DOF_4 *dof) { dof->x_eq = 0.; dof->y_eq = 0.; dof->z_eq = 0.; dof->c_eq = 0.;}

void dof2_init_array(DOF_2 *dof, int array_size) {
    int i;
    for (i=0; i<array_size; i++) {
        dof2_init(&dof[i]);
    }
}

void dof3_init_array(DOF_3 *dof, int array_size) {
    int i;
    for (i=0; i<array_size; i++) {
        dof3_init(&dof[i]);
    }
}

void dof4_init_array(DOF_4 *dof, int array_size) {
    int i;
    for (i=0; i<array_size; i++) {
        dof4_init(&dof[i]);
    }
}

// debug
void dof3_debug(DOF_3 *dof, int n, char *filename, int linenumber) {
    int i = 0;
    for(i = 0; i < n; i++) {
        Is_Double_Inf_or_NaN(dof[i].c_eq, filename, linenumber);
        Is_Double_Inf_or_NaN(dof[i].x_eq, filename, linenumber);
        Is_Double_Inf_or_NaN(dof[i].y_eq, filename, linenumber);
    }
    return;
}

void dof4_debug(DOF_4 *dof, int n, char *filename, int linenumber) {
    int i = 0;
    for(i = 0; i < n; i++) {
        Is_Double_Inf_or_NaN(dof[i].c_eq, filename, linenumber);
        Is_Double_Inf_or_NaN(dof[i].x_eq, filename, linenumber);
        Is_Double_Inf_or_NaN(dof[i].y_eq, filename, linenumber);
		Is_Double_Inf_or_NaN(dof[i].z_eq, filename, linenumber);
    }
    return;
}

int dof3_debug_continue(DOF_3 *dof, int n) {
    int i = 0, FLAG = 0;
    for(i = 0; i < n; i++) {
        if (solv_isinf(dof[i].c_eq) != 0 || solv_isnan(dof[i].c_eq) != 0) FLAG = 1;
        if (solv_isinf(dof[i].x_eq) != 0 || solv_isnan(dof[i].x_eq) != 0) FLAG = 1;
        if (solv_isinf(dof[i].y_eq) != 0 || solv_isnan(dof[i].y_eq) != 0) FLAG = 1;
    }
    return FLAG;
}

int dof4_debug_continue(DOF_4 *dof, int n) {
    int i = 0, FLAG = 0;
    for(i = 0; i < n; i++) {
        if (solv_isinf(dof[i].c_eq) != 0 || solv_isnan(dof[i].c_eq) != 0) FLAG = 1;
        if (solv_isinf(dof[i].x_eq) != 0 || solv_isnan(dof[i].x_eq) != 0) FLAG = 1;
        if (solv_isinf(dof[i].y_eq) != 0 || solv_isnan(dof[i].y_eq) != 0) FLAG = 1;
		if (solv_isinf(dof[i].z_eq) != 0 || solv_isnan(dof[i].z_eq) != 0) FLAG = 1;
    }
    return FLAG;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

// copy routines

void dof2_copy(DOF_2 *dof_to, DOF_2 dof_from) { 
    dof_to->x_eq = dof_from.x_eq; dof_to->y_eq = dof_from.y_eq; 
}
void dof3_copy(DOF_3 *dof_to, DOF_3 dof_from) { 
    dof_to->x_eq = dof_from.x_eq; dof_to->y_eq = dof_from.y_eq; dof_to->c_eq = dof_from.c_eq; 
}
void dof4_copy(DOF_4 *dof_to, DOF_4 dof_from) { 
    dof_to->x_eq = dof_from.x_eq; dof_to->y_eq = dof_from.y_eq; dof_to->z_eq = dof_from.z_eq; dof_to->c_eq = dof_from.c_eq; 
}

void dof2_copy_array(DOF_2 *dof_to, DOF_2 *dof_from, int array_size) {
    int i;
    for (i=0; i<array_size; i++) {
        dof2_copy(&dof_to[i], dof_from[i]);
    }
}

void dof3_copy_array(DOF_3 *dof_to, DOF_3 *dof_from, int array_size) {
    int i;
    for (i=0; i<array_size; i++) {
        dof3_copy(&dof_to[i], dof_from[i]);
    }
}

void dof4_copy_array(DOF_4 *dof_to, DOF_4 *dof_from, int array_size) {
    int i;
    for (i=0; i<array_size; i++) {
        dof4_copy(&dof_to[i], dof_from[i]);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

// print routines

void dof1_printScreen(double dof) {
    printf("x_eq: %*.*e \n",printFieldWidth,printPrecision,dof);
}

void dof1_printScreen_array(char * descript, double *dof, int array_size, int linenumber, char *filename) {
    int i;
    printf("\n");
    printf("printing dof: %s @ line %s:%d :: size: %d \n",descript, filename, linenumber, array_size);
    for (i=0; i<array_size; i++) {
        printf("[%d] \t",i);
        dof1_printScreen(dof[i]);
    }
}

void dof2_printScreen(DOF_2 dof) {
    printf("x_eq: %*.*e \t y_eq: %*.*e\n",
            printFieldWidth,printPrecision,dof.x_eq,
            printFieldWidth,printPrecision,dof.y_eq);
}

void dof2_printScreen_array(char * descript, DOF_2 *dof, int array_size, int linenumber, char *filename) {
    int i;
    printf("\n");
    printf("printing dof2: %s @ line %s:%d :: size: %d \n",descript, filename, linenumber, array_size);
    for (i=0; i<array_size; i++) {
        printf("[%d] \t",i);
        dof2_printScreen(dof[i]);
    }
}

void dof3_printScreen(DOF_3 dof) {
    printf("c_eq: %*.*e \t x_eq: %*.*e \t y_eq: %*.*e\n",
            printFieldWidth,printPrecision,dof.c_eq,
            printFieldWidth,printPrecision,dof.x_eq,
            printFieldWidth,printPrecision,dof.y_eq);
}

void dof3_printScreen_array(char * descript, DOF_3 *dof, int array_size, int linenumber, char *filename) {
    int i;
    printf("\n");
    printf("printing dof3: %s @ line %s:%d :: size: %d \n",descript, filename, linenumber, array_size);
    for (i=0; i<array_size; i++) {
        printf("[%d] \t",i);
        dof3_printScreen(dof[i]);
    }
}

void dof4_printScreen(DOF_4 dof, int printFieldWidth, int printPrecision) {
     printf("c_eq: %+-*.*e  x_eq: %+-*.*e  y_eq: %+-*.*e  z_eq: %+-*.*e\n",
             printFieldWidth,printPrecision,dof.c_eq, 
             printFieldWidth,printPrecision,dof.x_eq, 
             printFieldWidth,printPrecision,dof.y_eq, 
             printFieldWidth,printPrecision,dof.z_eq);
}

void dof4_printScreen_array(char * descript, DOF_4 *dof, int array_size, int linenumber, char *filename, int printFieldWidth, int printPrecision) {
    int i;
    printf("\n");
    printf("printing dof4: %s @ line %s:%d :: size: %d \n",descript, filename, linenumber, array_size);
    for (i=0; i<array_size; i++) {
        printf("[%d] \t",i);
        dof4_printScreen(dof[i],printFieldWidth,printPrecision);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

// subtract

DOF_2 dof2_subtract(DOF_2 dof1, DOF_2 dof2) {
    DOF_2 result;
    result.x_eq = dof1.x_eq - dof2.x_eq; 
    result.y_eq = dof1.y_eq - dof2.y_eq;
    return result;
}

DOF_3 dof3_subtract(DOF_3 dof1, DOF_3 dof2) {
    DOF_3 result;
    result.c_eq = dof1.c_eq - dof2.c_eq;
    result.x_eq = dof1.x_eq - dof2.x_eq;
    result.y_eq = dof1.y_eq - dof2.y_eq;
    return result;
}

DOF_4 dof4_subtract(DOF_4 dof1, DOF_4 dof2) {
    DOF_4 result;
    result.c_eq = dof1.c_eq - dof2.c_eq;
    result.x_eq = dof1.x_eq - dof2.x_eq;
    result.y_eq = dof1.y_eq - dof2.y_eq;
    result.z_eq = dof1.z_eq - dof2.z_eq;
    return result;
}

void dof2_subtract_array(DOF_2 *diff, DOF_2 *dof1, DOF_2 *dof2, int array_size) {
    int i;
    for (i=0; i<array_size; i++) {
        diff[i] = dof2_subtract(dof1[i], dof2[i]);
    }
}

void dof3_subtract_array(DOF_3 *diff, DOF_3 *dof1, DOF_3 *dof2, int array_size) {
    int i;
    for (i=0; i<array_size; i++) {
        diff[i] = dof3_subtract(dof1[i], dof2[i]);
    }
}

void dof4_subtract_array(DOF_4 *diff, DOF_4 *dof1, DOF_4 *dof2, int array_size) {
    int i;
    for (i=0; i<array_size; i++) {
        diff[i] = dof4_subtract(dof1[i], dof2[i]);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

// add

DOF_2 dof2_add(DOF_2 dof1, DOF_2 dof2) {
    DOF_2 result;
    result.x_eq = dof1.x_eq + dof2.x_eq;
    result.y_eq = dof1.y_eq + dof2.y_eq;
    return result;
}

DOF_3 dof3_add(DOF_3 dof1, DOF_3 dof2) {
    DOF_3 result;
    result.c_eq = dof1.c_eq + dof2.c_eq;
    result.x_eq = dof1.x_eq + dof2.x_eq;
    result.y_eq = dof1.y_eq + dof2.y_eq;
    return result;
}

DOF_4 dof4_add(DOF_4 dof1, DOF_4 dof2) {
    DOF_4 result;
    result.c_eq = dof1.c_eq + dof2.c_eq;
    result.x_eq = dof1.x_eq + dof2.x_eq;
    result.y_eq = dof1.y_eq + dof2.y_eq;
    result.z_eq = dof1.z_eq + dof2.z_eq;
    return result;
}

void dof2_add_array(DOF_2 *diff, DOF_2 *dof1, DOF_2 *dof2, int array_size) {
    int i;
    for (i=0; i<array_size; i++) {
        diff[i] = dof2_add(dof1[i], dof2[i]);
    }
}

void dof3_add_array(DOF_3 *diff, DOF_3 *dof1, DOF_3 *dof2, int array_size) {
    int i;
    for (i=0; i<array_size; i++) {
        diff[i] = dof3_add(dof1[i], dof2[i]);
    }
}

void dof3_add_replace_array(DOF_3 *dof1, DOF_3 *dof2, int array_size) {
    int i;
    for (i=0; i<array_size; i++) {
        dof1[i] = dof3_add(dof1[i], dof2[i]);
    }
}

void dof4_add_array(DOF_4 *diff, DOF_4 *dof1, DOF_4 *dof2, int array_size) {
    int i;
    for (i=0; i<array_size; i++) {
        diff[i] = dof4_add(dof1[i], dof2[i]);
    }
}

        

