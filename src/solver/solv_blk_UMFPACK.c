#include <stdio.h>
#ifdef _UMFPACK
#include "umfpack.h"
#include "global_header.h"

/* this file interfaces to UMFPACK */
/* The distribution site is http://www.cise.ufl.edu/research/sparse/umfpack/ */
/* ( A user guide is included in the distribution) */

static double Info[UMFPACK_INFO], Control[UMFPACK_CONTROL];
static void *Symbolic, *Numeric;
static long status;
static int solve_number;
static int debug_level = 0;
static int myid_u, npes_u;
static double sum, avg, max, bal;

#ifdef _MESSG
int solv_blk_UMFPACK_fact(int, int, long *, long *, double *, char *, MPI_Comm);
int solv_blk_UMFPACK_solve(long *, long *, double *, int, double *, double *, MPI_Comm);
#else
int solv_blk_UMFPACK_fact(int, int, long *, long *, double *, char *);
int solv_blk_UMFPACK_solve(long *, long *, double *, int, double *, double *);
#endif

/*performs LU factorization of sparse matrix */
#ifdef _MESSG
int solv_blk_UMFPACK_fact(int n, int nnz, long *Ap, long *Ai, double *Ax, char *form, MPI_Comm ADH_COMM)
#else
int solv_blk_UMFPACK_fact(int n, int nnz, long *Ap, long *Ai, double *Ax, char *form)
#endif
  /* on input, {Ap,Ai,Ax} are in either "triplet" form or compressed column form. */
  /* on output, {Ap,Ai,Ax} are in std compressed column form */
{
    long *Map = NULL;
    long *Ti = NULL;
    long *Tj = NULL;
    double *Tx = NULL;
    int i = 0;
    int UMFail = NO;

    solve_number = 0;

    umfpack_dl_defaults(Control);    /* default control parameters */

    Control[UMFPACK_PRL] = 1;   /* print level */

    if (strncmp(form, "triplet", 7) == 0) {

        /* temporary storage for triplet form */

        Ti = (long *) tl_alloc(sizeof(long), nnz);
        Tj = (long *) tl_alloc(sizeof(long), nnz);
        Tx = (double *) tl_alloc(sizeof(double), nnz);

        for (i = 0; i < nnz; i++) {
            Tj[i] = Ai[i];
            Ti[i] = Ap[i];
            Tx[i] = Ax[i];
        }
        /* convert triplet form to column form */

        (void) umfpack_dl_report_triplet((long)(n), (long) n, (long) nnz, Ti, Tj, Tx, Control);
        status = umfpack_dl_triplet_to_col((long)(n), (long) n, (long) nnz, Ti, Tj, Tx, Ap, Ai, Ax, Map);

        if (status != UMFPACK_OK) {
            umfpack_dl_report_status(Control, status);
            /*tl_error
               ("Fatal error in umfpack_l_triplet_to_col, called from solv_blk_UMFPACK_fact."); */
            UMFail = YES;
        }
        Tx = (double *) tl_free(sizeof(double), nnz, Tx);
        Ti = (long *) tl_free(sizeof(long), nnz, Ti);
        Tj = (long *) tl_free(sizeof(long), nnz, Tj);
    }

    /* print the column-form of A */

    (void) umfpack_dl_report_matrix((long)(n), (long) (n), Ap, Ai, Ax, (long)(1), Control);

    status = umfpack_dl_symbolic((long)(n), (long)(n), Ap, Ai, Ax, &Symbolic, Control, Info);
    if (status != UMFPACK_OK) {
        umfpack_dl_report_info(Control, Info);
        umfpack_dl_report_status(Control, status);
        /*tl_error
           ("Fatal error in umfpack_l_symbolic, called from solv_blk_UMFPACK_fact."); */
        UMFail = YES;
    }
    status = umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
    if (status != UMFPACK_OK) {
        umfpack_dl_report_info(Control, Info);
        umfpack_dl_report_status(Control, status);
        /*tl_error
           ("Fatal error in umfpack_l_numeric, called from solv_blk_UMFPACK_fact."); */
        UMFail = YES;
    }
#ifdef _MESSG
    if (debug_level > 0) {
        npes_u = messg_comm_size(ADH_COMM);
        myid_u = messg_comm_rank(ADH_COMM);

        sum = messg_dsum(Info[UMFPACK_NROW], ADH_COMM);
        max = messg_dmax(Info[UMFPACK_NROW], ADH_COMM);

        if (myid_u == 0) {
            avg = sum / npes_u;
            bal = max / avg;
            printf("\nFine-Grid preconditioner factorizatiuon statistics (solv_blk_UMFPACK_fact)\n");
            printf("matrix size :\t max %8d \t avg %8d \t max/avg %8.2f\n", (int) max, (int) avg, bal);
        }
        sum = messg_dsum(Info[UMFPACK_NZ], ADH_COMM);
        max = messg_dmax(Info[UMFPACK_NZ], ADH_COMM);

        if (myid_u == 0) {
            avg = sum / npes_u;
            bal = max / avg;
            printf("nonzeros   :\t max %8d \t avg %8d \t max/avg %8.2f\n", (int) max, (int) avg, bal);
        }
        sum = messg_dsum(Info[UMFPACK_LU_ENTRIES], ADH_COMM);
        max = messg_dmax(Info[UMFPACK_LU_ENTRIES], ADH_COMM);

        if (myid_u == 0) {
            avg = sum / npes_u;
            bal = max / avg;
            printf("LU nonzeros:\t max %8d \t avg %8d \t max/avg %8.2f\n", (int) max, (int) avg, bal);
        }
        sum = messg_dsum(Info[UMFPACK_NUMERIC_TIME], ADH_COMM);
        max = messg_dmax(Info[UMFPACK_NUMERIC_TIME], ADH_COMM);

        if (myid_u == 0) {
            avg = sum / npes_u;
            if (sum > 0.0)
                bal = max / avg;
            else
                bal = 0.0;
            printf("numeric sec:\t max %8.2e \t avg %8.2e \t max/avg %8.2f\n", max, avg, bal);
        }
    }
#endif
    return UMFail;

}

#ifdef _MESSG
int solv_blk_UMFPACK_solve     /* solve linear system with LU factorization */
    (long *Ap, long *Ai, double *Ax, int n, double *b, double *x, MPI_Comm ADH_COMM)
#else
int solv_blk_UMFPACK_solve     /* solve linear system with LU factorization */
    (long *Ap, long *Ai, double *Ax, int n, double *b, double *x)
#endif
  /* this routine overwrits the rhs with the solution */
{
    int i;
    int UMFail = NO;

    status = umfpack_dl_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, Control, Info);
    if (status != UMFPACK_OK) {
        umfpack_dl_report_info(Control, Info);
        umfpack_dl_report_status(Control, status);
        /*tl_error
           ("Fatal error in umfpack_l_solve, called from solv_blk_UMFPACK_solve."); */
        UMFail = YES;
    }
    for (i = 0; i < n; i++)
        b[i] = x[i];
#ifdef _MESSG
    if ((debug_level > 0) && (++solve_number == 1)) {   /* cjt edited */
        npes_u = messg_comm_size(ADH_COMM);
        myid_u = messg_comm_rank(ADH_COMM);
        sum = messg_dsum(Info[UMFPACK_SOLVE_TIME], ADH_COMM);
        max = messg_dmax(Info[UMFPACK_SOLVE_TIME], ADH_COMM);

        if (myid_u == 0) {
            avg = sum / npes_u;
            if (sum > 0.0)
                bal = max / avg;
            else
                bal = 0.0;
      /**printf ("\nFine Grid preconditioner statistics for first solve\n");
	  printf ("solve sec:\t max %8e \t avg %8e \t max/avg %8.2f\n", max,
		  avg, bal);**/
        }
    }
#endif
    return UMFail;
}

void solv_blk_UMFPACK_free()
{  /* deallocate data structures */
    umfpack_dl_free_symbolic(&Symbolic);
    umfpack_dl_free_numeric(&Numeric);
}
#else
void solv_blk_UMFPACK_free();
#endif /* ifdef UMFPACK */
