/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  tl_matrix.c This file collects functions that operate on the AdH matrix          */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Initializes the AdH Newton Jacobi matrix.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] matrix the sparse AdH matrix
 *  @param[in,out] diagonal the sparse AdH matrix diagonal blocks
 *  @param[in,out] nnodes the total number of grid nodes
 *  @param[in,out] max_nsys_sq the maximum number of equations in system squared
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void init_adh_matrix(int nnodes, int max_nsys_sq, SPARSE_VECT *matrix, double *diagonal) {
    
    int i, iend;

    for (i = 0, iend = nnodes * max_nsys_sq; i < iend; i++) {
        diagonal[i] = 0.0;
    }
    
    for (i = 0; i < nnodes; i++) {
        spv_init(matrix[i], max_nsys_sq);
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Checks the diagonal of the AdH Jacobi Matrix for "non-existant" entries and "zeros" them.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] diagonal the sparse AdH matrix diagonal blocks
 *  @param[in] nnodes the total number of grid nodes
 *  @param[in] max_nsys_sq the maximum number of equations in system squared
 *  @param[in] nsys the number of equations in system
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void check_matrix_diagonal_for_nonexistant_entries(int nnodes, int max_nsys_sq, int nsys, double *diagonal) {
    int i,j, isys;
    for (i=0; i<nnodes; i++) {
        j = i * max_nsys_sq;
        for (isys=0; isys<nsys; isys++) {
            if (fabs(diagonal[i]) < NOT_QUITE_SMALL) {diagonal[i] = NOT_QUITE_SMALL;}
            j += (nsys + 1);
        }
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Prints global matrix to screen
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in] descript (char *) the matrix print desciption
 *  @param[in] diagonal (double *) the block diagonal array of the matrix
 *  @param[in] matrix (SPARSE_VECT *) the off-diagonal blocks of the matrix
 *  @param[in] nnode (int) the number of global nodes
 *  @param[in] nsys_sq (int) the number of equations in system squared
 *  @param[in] linenumber (int) the line number of the calling routine
 *  @param[in] filename (int) the file name of the calling routine
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void printScreen_matrix(char * descript, double *diagonal, SPARSE_VECT * matrix, int nnode, int nsys_sq, int linenumber, char *filename)
{
    int i, j, k;
    int ierr_code, myid=0;
#ifdef _MESSG
    ierr_code = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
    char fn[30+1];
    snprintf(fn, 30, "matrix_world_pe%d.txt", myid);
    FILE *fp = fopen(fn, "w");
    fflush(fp);

    fprintf(fp,"printing matrix: %s @ line %s:%d \n",descript,filename, linenumber);
    for (i = 0; i < nnode; i++) {
        for (j = 0; j < nsys_sq; j++) {
            fprintf(fp," i=%d, diagonal %12.4e \n", i + 1, diagonal[i * nsys_sq + j]);
        }
        fprintf(fp,"\n\n----------row %d--------------\n", (i + 1));
        for (j = 0; j < matrix[i].size; j++) {
            fprintf(fp,"**col %d**\n", (matrix[i].index[j] + 1));
            for (k = j * nsys_sq; k < (j + 1) * nsys_sq; k++)
                fprintf(fp,"%3d %12.4e\n", k - j * nsys_sq, matrix[i].value[k]);
        }
    }

    fflush(fp);
    fclose(fp);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Prints global matrix to FILE
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in] fp (FILE *) a file pointer
 *  @param[in] descript (char *) the matrix print desciption
 *  @param[in] diagonal (double *) the block diagonal array of the matrix
 *  @param[in] matrix (SPARSE_VECT *) the off-diagonal blocks of the matrix
 *  @param[in] nnode (int) the number of global nodes
 *  @param[in] nsys_sq (int) the number of equations in system squared
 *  @param[in] linenumber (int) the line number of the calling routine
 *  @param[in] filename (int) the file name of the calling routine
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void printFile_matrix(FILE *fp, char *descript, double *diagonal, SPARSE_VECT * matrix, int nnode, int nsys_sq, int linenumber, char *filename) {
    int i, j, k;
    fprintf(fp,"\n");
    fprintf(fp,"printing matrix: %s @ line %s:%d \n",descript,filename, linenumber);
    for (i = 0; i < nnode; i++) {
        for (j = 0; j < nsys_sq; j++) {
            fprintf(fp," i=%d, diagonal %30.20e \n", i + 1, diagonal[i * nsys_sq + j]);
        }
        fprintf(fp,"\n\n----------row %d--------------\n", (i + 1));
        for (j = 0; j < matrix[i].size; j++) {
            fprintf(fp,"**col %d**\n", (matrix[i].index[j] + 1));
            for (k = j * nsys_sq; k < (j + 1) * nsys_sq; k++)
                fprintf(fp,"%3d %30.20e\n", k - j * nsys_sq, matrix[i].value[k]);
        }
    }
    
}




