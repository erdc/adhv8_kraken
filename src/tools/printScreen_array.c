#include <stdio.h>
#include <string.h>

void printScreen_dble_array(char * descript, double *array, int size,  int linenumber, char *filename) {
    int i;
    printf("\n");
    printf("printing array: %s @ line %s:%d :: size: %d \n",descript,filename, linenumber, size);
    for (i=0; i<size; i++) {
        printf("[%d] \t %30.20f \n",i,array[i]);
    }
}

void printScreen_int_array(char * descript, int *array, int size, int linenumber, char *filename) {
    int i;
    printf("\n");
    printf("printing array: %s @ line %s:%d :: size: %d \n",descript,filename, linenumber, size);
    for (i=0; i<size; i++) {
        printf("[%d] \t %20d \n",i,array[i]);
    }
}
