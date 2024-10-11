//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <stdbool.h>

#include "adh.h"

static char delim[] = " ";

void get_token(char *line, char **token) {
    (*token) = strtok(line,delim);
    if (*token == NULL) return;
    token[strcspn((*token), "\r\n")] = 0;
}

void get_next_token(char **token) {
    (*token) = strtok(NULL,delim);
    if (*token == NULL) return;
    token[strcspn((*token), "\r\n")] = 0;
}

int get_next_token_int(char **token) {
    (*token) = strtok(NULL,delim);
    if (*token == NULL) return -999;
    token[strcspn((*token), "\r\n")] = 0;
    int itmp;
    sscanf((*token), "%d", &itmp);
    return itmp;
}

double get_next_token_dbl(char **token) {
    (*token) = strtok(NULL,delim);
    if (*token == NULL) return -999;
    token[strcspn((*token), "\r\n")] = 0;
    double dtmp;
    sscanf((*token), "%lf", &dtmp);
    return dtmp;
}


int get_compare_token(char *line, char **token, char *cmp) {
    get_token(line,token);
    if (*token == NULL) return -1; //empty line
    if (strcmp(*token, cmp) == 0) {
        return 1;
    } else {
        return 0;
    }
}

void get_compare_token_increment(char *line, char **token, char *cmp, int *idum) {
    get_token(line,token);
    if (*token == NULL) return; //empty line
    if (strcmp(*token, cmp) == 0) {
        (*idum)++;
        return;
    }
}


int get_compare_next_token(char **token, char *cmp) {
    get_next_token(token);
    if (*token == NULL) return -1; //empty line
    if (strcmp(*token, cmp) == 0) {
        return 1;
    } else {
        return 0;
    }
}
