#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void get_token(char *line, char **token);
void get_next_token(char **token);
int  get_next_token_int(char **token);
double get_next_token_dbl(char **token);

void read_next_level(char **token) {
    int idum = -3;
    if (strcmp((*token),"SW2") == 0) {
        idum = get_next_token_int(token);
        printf("idum: %d \n",idum);
    }
    
}
