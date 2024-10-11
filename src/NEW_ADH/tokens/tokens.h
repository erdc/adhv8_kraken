#ifndef H_TOKENS_
#define H_TOKENS_
// token
void get_token(char *line, char **token);
void get_next_token(char **token);
int get_next_token_int(char **token);
double get_next_token_dbl(char **token);
int get_compare_token(char *line, char **token, char *cmp);
void get_compare_token_increment(char *line, char **token, char *cmp, int *idum);
int get_compare_next_token(char **token, char *cmp);

#endif
