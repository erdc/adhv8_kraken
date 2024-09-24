#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char delim[] = " ";

void get_token(char *line, char **token) {
    (*token) = strtok(line,delim);
    token[strcspn((*token), "\n")] = 0;
}

void get_next_token(char **token) {
    (*token) = strtok(NULL,delim);
    token[strcspn((*token), "\n")] = 0;
}

int get_next_token_int(char **token) {
    (*token) = strtok(NULL,delim);
    token[strcspn((*token), "\n")] = 0;
    int itmp;
    sscanf((*token), "%d", &itmp);
    return itmp;
}

double get_next_token_dbl(char **token) {
    (*token) = strtok(NULL,delim);
    token[strcspn((*token), "\n")] = 0;
    double dtmp;
    sscanf((*token), "%lf", &dtmp);
    return dtmp;
}
