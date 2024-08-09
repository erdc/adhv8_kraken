#ifndef H_TL_ALLOC_
#define H_TL_ALLOC_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#ifdef _MESSG
#include "mpi.h"
MPI_Comm debug_comm; // cjt :: can't use COMM_WORLD for CSTORM :: only works for 1 grid in CSTORM
#endif

#ifdef _MESSG
void debug_initialize(MPI_Comm);
#else
void debug_initialize(void);
#endif
void debug_finalize(void);

#ifndef _DEBUG
void tag(void);
void tl_error(char *);
void *tl_alloc(size_t,int);
void *tl_realloc(size_t, int, int, void *);
void *tl_free(size_t, int, void *);
void *tl_free2(size_t, int, size_t, int, void **);
void *tl_free3(size_t, int, size_t, int, size_t, int, void ***);

#else

#ifdef _MESSG
#define tag(A) tag_debug(__LINE__,__FILE__,A)
#else
#define tag() tag_debug(__LINE__,__FILE__)
#endif

#define tl_error(A) tl_error_debug(A,__LINE__,__FILE__)
#define tl_alloc(A,B) tl_alloc_debug(A,B,__LINE__,__FILE__)
#define tl_realloc(A,B,C,D) tl_realloc_debug(A,B,C,__LINE__,__FILE__,D)
#define tl_free(A,B,C) tl_free_debug(A,B,__LINE__,__FILE__,C)
#define tl_free2(A,B,C,D,E) tl_free2_debug(A,B,C,D,__LINE__,__FILE__,E)
#define tl_free3(A,B,C,D,E,F,G) tl_free3_debug(A,B,C,D,E,F,__LINE__,__FILE__,G)

#ifdef _MESSG
void tag_debug(int ,char *,MPI_Comm); // cjt
#else
void tag_debug(int ,char *);
#endif
void tl_error_debug(char *, int, char *); // cjt
void *tl_alloc_debug(size_t, int, int, char *);
void *tl_realloc_debug(size_t, int, int, int, char *, void *);
void *tl_free_debug(size_t, int, int, char *, void *);
void *tl_free2_debug(size_t, int, size_t, int, int, char *, void **);
void *tl_free3_debug(size_t, int, size_t, int, size_t, int, int, char *, void ***);
void tl_check_all_pickets(char *, int);
#endif

#endif // ifndef H_TL_ALLOC
