/*!
   \file solv_type.h
   \brief Data Structures Needed by Linear Solver Routine
 */

#ifdef SOLV_INIT
#define EXTERN_SOLV
#else
#define EXTERN_SOLV extern
#endif

#ifndef _UMFPACK_INT
#define _UMFPACK_INT 0
#endif

typedef struct {
  double *value;        /* the values of the sparse vector */
  int *index;           /* the indices of the sparse vector */
  int size;         /* the number of entries in the sparse vector */
  int max_size;         /* the number of allocated entries in the sparse vector */
} SPARSE_VECT;          /* a sparse vector */

typedef struct
{
  int begin;            /* the beginning of the band */
  int end;          /* the end of the band */
  int size;         /* the allocated size of the band */
  double *value;        /* the values */
} BAND_VECT;            /* a band storage sparse vector */

/*!
   \struct SPARSE_CSR 
   \brief Compressed Sparse Row, matrix storage format 

   \param num_row_alloc Number of Rows Allocated
   \param nnz_alloc Number of Non-Zero Entries Allocated
   \param num_row Number of Rows in Matrix 
   \param num_col Number of Columns in Matrix 
   \param nnz Number of Non-Zero Entries
   \param A Pointer to Entries of Matrix
   \param C_ind Pointer to Column Indices
   \param R_ptr Pointer to Row Pointers
 */
typedef struct
{
  int row_alloc;		/* Number of Rows Allocated */
  int nnz_alloc;		/* Number of Non-zero Entries Allocated */
  int num_row;			/* Number of Rows of Matrix */
  int num_col;			/* Number of Columns of Matrix */
  int nnz;			/* Number of Non-zero Entries of Matrix */
  double *A;			/* Matrix Entries */
  int *C_ind;			/* Column Indices */
  int *R_ptr;			/* Row pointer */
} SPARSE_CSR;

/*!
   \struct SPARSE_CSC 
   \brief Compressed Sparse Column, matrix storage format 

   \param num_row Number of Rows
   \param num_col Number of Columns
   \param nnz Number of Non-Zero Entries
   \param A Pointer to Entries of Matrix
   \param R_ind Pointer to Row Indices
   \param C_ptr Pointer to Column Pointers
 */
typedef struct
{
  int col_alloc;		/* Number of Rows Allocated */
  int nnz_alloc;		/* Number of Non-zero Entries Allocated */
  int num_row;			/* Number of Rows */
  int num_col;			/* Number of Columns */
  int nnz;			/* Number of Non-zero Entries */
  double *A;			/* Matrix Entries */
  int *R_ind;			/* Row Indices */
  int *C_ptr;			/* Column pointer */
} SPARSE_CSC;

/*!
   \struct SPARSE_AIJ 
   \brief Coordinate (or AIJ) Format, matrix storage format 

   \param num_row Number of Rows
   \param num_col Number of Columns
   \param nnz Number of Non-Zero Entries
   \param A Pointer to Entries of Matrix
   \param IR Pointer to Row Indices
   \param JC Pointer to Column Indices
 */
typedef struct
{
  int nnz_alloc;		/* Number of Non-zero Entries Allocated */
  int num_row;			/* Number of Rows */
  int num_col;			/* Number of Columns */
#if _UMFPACK_INT == 32
  int nnz;			/* Number of Non-zero Entries */
  int *IR;			/* Row Indices */
  int *JC;			/* Column Indices */
#elif _UMFPACK_INT == 64
  long nnz;			/* Number of Non-zero Entries */
  long *IR;			/* Row Indices */
  long *JC;			/* Column Indices */
#endif
  double *A;			/* Matrix Entries */
} SPARSE_AIJ;


typedef struct {
    int size; //= 0;               /* the size of the profile matrix */
    BAND_VECT *rows; // = NULL;    /* the rows of the profile matrix */
    BAND_VECT *cols; // = NULL;    /* the columns of the profile matrix */
} Profile_Info;

/*!
   \struct Solver_Matrices
   \brief Pointers to Matrix, Left and Right Preconditioner abstracts
 */
typedef struct
{
  int S_MAT_KIND;
  int A_MAT_KIND;
  int M_MAT_KIND;
  SPARSE_AIJ S_AIJ;
  SPARSE_CSR S_CSR;
  SPARSE_AIJ A_AIJ;
  SPARSE_CSR A_CSR;
  SPARSE_AIJ M_AIJ;
  SPARSE_CSR M_CSR;
} Solver_Matrices;

/*!
   \struct Solver_Info
   \brief General Parameters Passed to Iterative Solver
 */
typedef struct
{
  double tol_lin;                   /* Linear Tolerance */
  double tol_nonlin;                /* Non-linear Tolerance */
  int prec_value;                   /* Which Preconditioner */
  int it_count_lin;                 /* Total number of linear steps taken */
  int it_count_nonlin;              /* Total number of nonlinear steps taken */
  int it_count_lin_failed;          /* Linear steps lost to failure */
  int it_count_nonlin_failed;       /* Non-Linear steps lost to failure */
  int lin_fail;                     /* Indication of linear failure YES/NO */
  int nonlin_fail;                  /* Indication of non-linear failure YES/NO */
  int max_lin_it;                   /* Max Num of linear iter. per step */
  int max_nonlin_it;                /* Max Num of non-linear iter. per step */
  int blocksize;                    /* Block Size */
  int refresh;                      /* Refresh preconditioner */
  int force_lin_it;                 /* Force Linear solver to Accept End result */
  int force_nonlin_it;              /* Force Non-Linear solver to Accept End result */
  int max_nonlin_linesearch_cuts;	/* Max Num of non-linear line search cuts */
  int LINEAR_PROBLEM;               /* flag whether a linear problem or not */
  int PRN_NEWTON;                   /* if ON, print the resid info in fe_newton */
  int UMFail;

  int *node_block;        /* which block a node is assigned */
  int *solv_node_block;
  Solver_Matrices solver_matrices;
  Profile_Info profile;

} Solver_Info;

/*********************************************************/
/* struct methods -------------------------------------- */


void spv_init(SPARSE_VECT, int);
void spv_alloc(SPARSE_VECT *, int);
void spv_dot(SPARSE_VECT, double *, double *, int);
void spv_load(SPARSE_VECT *, int, double *, int, int, int);
void spv_norm(SPARSE_VECT, double *, int);
void spv_reset(SPARSE_VECT *);
void spv_vscale(SPARSE_VECT *, double *, double *, int);

void bv_alloc(BAND_VECT *);
double bv_dot(BAND_VECT *, BAND_VECT *, int, int);
void bv_init(BAND_VECT *);
void bv_load(BAND_VECT *, double, int);

void solv_blk_load_profile(int *, SPARSE_VECT *, double *, BAND_VECT *, BAND_VECT *, int, int, int);
void solv_blk_load_sparse(int *, SPARSE_VECT *, double *, int, int, int, int, int, long *, long *, double *);
void solv_blk_factor_profile(BAND_VECT *, BAND_VECT *, int);
void solv_blk_set_profile(Solver_Info *, Profile_Info *, SPARSE_VECT *, double *, int, int, int);
void solv_blk_set_clean(Profile_Info *);
void solv_blk_free_sparse(void);


void solv_bcgstab_clean(void);
#ifdef _MESSG
int  solv_bcgstab(Solver_Info *, int *, SPARSE_VECT *, double *, double *, double *, double *, int, int, int, SMPI *);
void solv_prec(Solver_Info *, int *, SPARSE_VECT *, double *, double *, double *, double *, int, int, int, int, SMPI *);
void solv_amult(Solver_Info *, double *, double *, double *, SPARSE_VECT *, int, int, SMPI *);
void solv_return_old(Solver_Info *, double *, double *, double *, int, int, SMPI *);
int solv_blk_solve(Profile_Info, double *, int, MPI_Comm);
int solv_blk_set_sparse(int *, SPARSE_VECT *, double *, int, int, int, MPI_Comm);
void solv_blk_set(Solver_Info *, SPARSE_VECT *, double *, int, int, int, MPI_Comm);
void spv_jacobi(SPARSE_VECT *, double *, double *, double *, double *, int, int, int, SMPI *);
#else
int  solv_bcgstab(Solver_Info *, int *, SPARSE_VECT *, double *, double *, double *, double *, int, int, int);
void solv_prec(Solver_Info *, int *, SPARSE_VECT *, double *, double *, double *, double *, int, int, int, int);
void solv_amult(Solver_Info *, double *, double *, double *, SPARSE_VECT *, int, int);
void solv_return_old(Solver_Info *, double *, double *, double *, int, int);
int solv_blk_solve(Profile_Info, double *, int);
int solv_blk_set_sparse(int *, SPARSE_VECT *, double *, int, int, int);
void solv_blk_set(Solver_Info *, SPARSE_VECT *, double *, int, int, int);
void spv_jacobi(SPARSE_VECT *, double *, double *, double *, double *, int, int, int);
#endif

void solv_coarse_solve(int *, double *, int, int);

void quickSort(int *, double *, int);
void q_sort(int *, double *, int, int);

int  solv_bicgstab(Solver_Info *, int *, SPARSE_CSR *, double *, double *, int, int, int, double *, int *);
void solv_bicgstab_clean(void);
void solv_csr_initCSR(SPARSE_CSR *);
void solv_csr_freeCSR(SPARSE_CSR *);
void solv_csr_MatVecMult(SPARSE_CSR *, double *, double *, int);
void solv_csr_print(SPARSE_CSR *);
void solv_initialize(Solver_Info *);    // cjt
void solv_finalize(void);
int  solv_lapack_direct(SPARSE_VECT *, double *, double *, double *, int, int);
void solv_left_prec(int *, SPARSE_CSR *, double *, double *, double *, int, int);
#ifdef _MESSG
void solv_linear_sys_setup(Solver_Info *, int *, SPARSE_VECT *, double *, double *, double *, double *, int, int, int, SMPI *);
int  solv_linear_sys_solve(Solver_Info *, int *, SPARSE_VECT *, double *, double *, double *, double *, int, int, int, SMPI *);
#else
void solv_linear_sys_setup(Solver_Info *, int *, SPARSE_VECT *, double *, double *, double *, double *, int, int, int);
int  solv_linear_sys_solve(Solver_Info *, int *, SPARSE_VECT *, double *, double *, double *, double *, int, int, int);
#endif
void solv_point_jacobi(SPARSE_CSR *, double *, double *, double *, int);
void solv_return(double *, double *, int, int);
void solv_right_prec(int *, SPARSE_CSR *, double *, double *, double *, int, int);
int  solv_spv_2_aij(SPARSE_VECT *, double *, int, int, int, SPARSE_AIJ *);
int  solv_spv_2_csr(SPARSE_VECT *, double *, int, int, int, SPARSE_CSR *);

#ifdef _MESSG
int solv_blk_solve_sparse(double *, int, MPI_Comm);
int solv_blk_UMFPACK_fact(int, int, long *, long *, double *, char *, MPI_Comm);
int solv_blk_UMFPACK_solve(long *, long *, double *, int, double *, double *, MPI_Comm);
#else
int solv_blk_solve_sparse(double *, int);
int solv_blk_UMFPACK_fact(int, int, long *, long *, double *, char *);
int solv_blk_UMFPACK_solve(long *, long *, double *, int, double *, double *);
#endif
void solv_blk_load_sparse(int *, SPARSE_VECT *, double *, int, int, int, int, int, long *, long *, double *);
void solv_blk_UMFPACK_free();


