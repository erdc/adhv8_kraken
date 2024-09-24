#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int intFromBinary(const char *s) {
  return (int) strtol(s, NULL, 2);
}


void encode(char *num) {
    int i;
    int len = strlen(num);
    int result = 0;
    for(i=0; i<len; i++){
        result = result * 10 + ( num[i] - '0' );
    }
    printf("result: %d\n",result);
}

int count_vars(int score, char *vars) {
    int i, count = 0, digit, div, hflag = 0;;
    for (div = 1; div <= score; div *= 10);


    // surface water
    div /= 10;
    digit = score/div;
    score %= div;
    printf("surface_water: %d\n",digit);
    if (digit == 1) {
        count  += 2; // 1D SW
        strcat(vars,"hu");
        hflag = 1;
    } else if (digit == 2) {
        count  += 3; // 2D SW
        strcat(vars,"huv");
        hflag = 1;
    } else if (digit == 3) {
        count  += 3; // 3D SW - SPLIT - HVEL -  {h,u,v}
        strcat(vars,"huv");
        hflag = 1;
    } else if (digit == 4) {
        count  += 4; // 3D NS {u,v,w,p}
        strcat(vars,"uvwp");
    } else if (digit == 5) {
        count  += 3; // 3D NS - SPLIT - {u,v,w}
        strcat(vars,"uvw");
    } else if (digit == 6) {
        count  += 1; // 2D OF
        strcat(vars,"h");
        hflag = 1;
    } else if (digit == 7) {
        count  += 1; // SW 3D - SPLIT - w
        strcat(vars,"w");
    } else if (digit == 8) {
        count  += 1; // 3D NS - SPLIT - p
        strcat(vars,"p");
    } else {
        //tl_error("physics material error\n");
    }
    
    // ground water
    if (hflag == 0 ) { // not already included in surface water
        div /= 10;
        digit = score/div;
        score %= div;
        printf("ground_water: %d\n",digit);
        strcat(vars,"h");
        count += digit; // 1D/2D/3D GROUNDWATER
    } else {
        div /= 10;
        digit = score/div;
        score %= div;
    }
    
    // transport
    div /= 10;
    digit = score/div;
    score %= div;
    printf("transport: %d\n",digit);
    count += digit; // 1D/2D/3D TRANSPORT
    for (i=0; i<digit; i++) {
        char buffer[10];
        sprintf(buffer,"%d",i);
        strcat(vars,buffer);
    }
    
    return count;
}


int main() {
    int i;
    
    int score = 215;
    char vars[10] = "";
    int nvars = count_vars(score,vars);
    
    printf("nvars: %d \t vars: %s\n",nvars,vars);
    for (i = 0; i < strlen(vars); i++) {
        printf("var[%d]: %c\n",i,vars[i]);
    }
    
    
    
    
//    encode("huvw123456789");
    
//    int val;
//    char *str;
//    str = "1509.10E";
//    val = atoi(str);
//    printf("integral number = %d", val);
//
//
//    char *vars;
//    str = "huvwpd"; //123456789"; //"huvwd"; // hardwired in ADH
//    val = atoi(vars);
//    printf("varCode Integer: %d\n",val);
    
    
    
//    char *varCodeStr = "110101010101"; // found by element physics
//    //printf("varCode Integer: %d\n",intFromBinary(varCodeStr));
//
//    // which DOFs
//    int count = 0;
//    for (i = 0; i < strlen(varCodeStr); i++) {
//        char tmp = varCodeStr[i];
//        int varCode = atoi(&tmp);
//        char dof[20] = "";
//        if (i<5) {
//            strncpy(dof, vars + i, 1);
//        } else {
//            char num[3]="";
//            sprintf(num, "%d", i - 5);
//            strcpy(dof, "c");
//            strcat(dof, num);
//        }
//        printf("dof: %s || switch: %d\n",dof,varCode);
//    }
    return 0;
}


// "huvwpd123456789"
