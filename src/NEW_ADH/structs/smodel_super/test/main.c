
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void read_next_level(char **token);
void get_token(char *line, char **token);
void get_next_token(char **token);
int  get_next_token_int(char **token);
double get_next_token_dbl(char **token);

int main() {
    
    char *line = NULL, cdum[20], line_save[100];
    size_t len = 0;
    ssize_t read;
    char *token, *subtoken, *subsubtoken;
    char delim[] = " ";
    int idum = -3;
    double fdum;
    
    
    FILE *fp = fopen("test.bc", "r");
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token);
        printf("token: %s \n",token);
        if (strcmp(token, "OP") == 0) {
            
            get_next_token(&token);
            printf("subtoken: %s\n",token);
            read_next_level(&token);
            //if (strcmp(token,"SW2") == 0) {
            //    idum = get_next_token_int(&token);
            //    printf("idum: %d \n",idum);
            //}
        }
        if (strcmp(token, "ND") == 0) {
            double x=0,y=0,z=0;
            x = get_next_token_dbl(&token);
            y = get_next_token_dbl(&token);
            z = get_next_token_dbl(&token);
            printf("ND: %lf, %lf %lf \n",x,y,z);
        }
    }

    return 0;
}




