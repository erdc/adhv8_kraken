#include <global_header.h>
#include <ctype.h>

size_t strlen(const char *str) {
    const char *s;
    for (s = str; *s; ++s);
    return(s - str);
}

int main ( int argc, char *argv[] ) {
    
    int i, ia;
    char *name;
    char **data;
    
    if ( argc != 2 ) {
        /* We print argv[0] assuming it is the program name */
        printf( "usage: %s card", argv[0] );
        return -1;
    }
    
    name = argv[1];
    int length = strlen(name);
    printf("card length %d\n",length);
    

    ia = parse_card(name,data);

    printf("CARD_%s=%d\n",name,ia);
    return 0;
}
