#include "adh.h"

//*****************************************************//
//*****************************************************//
void sarray_init_int(int *array, int size) {
    int i;
    for (i=0; i<size; i++) {
        array[i] = 0;
    }
}
//----------------------------------------------------//
void sarray_init_dbl(double *array, int size) {
    int i;
    for (i=0; i<size; i++) {
        array[i] = 0.;
    }
}
//*****************************************************//
//*****************************************************//
void sarray_init_value_int(int *array, int size, int value) {
    int i;
    for (i=0; i<size; i++) {
        array[i] = value;
    }
}
//----------------------------------------------------//
void sarray_init_value_dbl(double *array, int size, double value) {
    int i;
    for (i=0; i<size; i++) {
        array[i] = value;
    }
}
//*****************************************************//
//*****************************************************//
void sarray_init_value_range_int(int *array, int value, int start, int end) {
    int i;
    for (i=start; i<end; i++) {
        array[i] = value;
    }
}
//----------------------------------------------------//
void sarray_init_value_range_dbl(double *array, double value, int start, int end) {
    int i;
    for (i=start; i<end; i++) {
        array[i] = value;
    }
}
//*****************************************************//
//*****************************************************//
void sarray_init2one_int(int *array, int size) {
    int i;
    for (i=0; i<size; i++) {
        array[i] = 1;
    }
}
//----------------------------------------------------//
void sarray_init2one_dbl(double *array, int size) {
    int i;
    for (i=0; i<size; i++) {
        array[i] = 1.;
    }
}
//*****************************************************//
//*****************************************************//
void sarray_copy_int(int *to, int *from, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        to[i] = from[i];
    }
}
//----------------------------------------------------//
//*****************************************************//
void sarray_init_int_2d(int **to, int nrows, int ncols) {
    int i,j;
    for (i=0; i<nrows; i++) {
        for (j=0; j<ncols; j++) {
            to[i][j] = 0;
        }
    }
}
//----------------------------------------------------//
//*****************************************************//
void sarray_init_double_2d(double **to, int nrows, int ncols) {
    int i,j;
    for (i=0; i<nrows; i++) {
        for (j=0; j<ncols; j++) {
            to[i][j] = 0.0;
        }
    }
}
//----------------------------------------------------//
//Mark, is memcpy faster??
//void sarray_copy_dbl(double *to, double *from, int nsize) {
//    int i;
//    for (i=0; i<nsize; i++) {
//        to[i] = from[i];
//    }
//}
void sarray_copy_dbl(double *to, double *from, int nsize){
#ifdef _BLAS
    int IONE = 1;
    dcopy_(&nsize, from, &IONE, to, &IONE);
#else
    memcpy(to, from, ((size_t) (nsize)) * sizeof(double));
#endif
    return;
}
//*****************************************************//
//*****************************************************//
int sarray_sum_int(int *array, int nsize) {
    int i=0;
    int sum = 0.;
    for (i=0; i<nsize; i++) {
        sum += array[i];
    }
    return sum;
}
//----------------------------------------------------//
double sarray_sum_dbl(double *array, int nsize) {
    int i=0;
    double sum = 0.;
    for (i=0; i<nsize; i++) {
        sum += array[i];
    }
    return sum;
}
/*!
 \brief Returns the inner product, \f$(x \cdot y)\f$, of two (double) vectors
 
 \param n the length of the vectors (int)
 \param x pointer to the first vector (double, length \a n)
 \param y pointer to the second vector (double, length \a n)
 */
double sarray_dot_dbl(double *x, double *y,int n)
{
    double value = DZERO;       /* the partial sum */
#ifdef _BLAS
    int IONE = 1;
    
    value = ddot_(&n, x, &IONE, y, &IONE);
#else
    int i = 0;          /* loop counter */
    
    for(i = 0; i < n; i++)
        value += x[i] * y[i];
#endif

    /* returns the sum */
    return (value);
}

double sarray_l_infty_norm(double *v1, int n)
{
    double value = DZERO;   /* the partial sum */
    int ii = 0;     /* loop counter */
    
    /* computes the value for this processor */
    for(ii = 0; ii < n; ii++)
    {
        value = max_dbl(value, fabs(v1[ii]));
    }
    
    /* returns the maximum */
    return value;
}

//*****************************************************//
//*****************************************************//
void sarray_add_replace_int(int *array1, int *array2, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        array1[i] += array2[i];
    }
}
//----------------------------------------------------//
void sarray_add_replace_dbl(double *array1, double *array2, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        array1[i] += array2[i];
    }
}
//*****************************************************//
//*****************************************************//
void sarray_add_int(int *result, int *array1, int *array2, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        result[i] = array1[i] + array2[i];
    }
}
//----------------------------------------------------//
void sarray_add_dbl(double *result, double *array1, double *array2, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        result[i] = array1[i] + array2[i];
    }
}
//*****************************************************//
//*****************************************************//
void sarray_subtract_replace_int(int *array1, int *array2, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        array1[i] -= array2[i];
    }
}
//----------------------------------------------------//
void sarray_subtract_replace_dbl(double *array1, double *array2, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        array1[i] -= array2[i];
    }
}
//*****************************************************//
//*****************************************************//
void sarray_subtract_int(int *result, int *array1, int *array2, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        result[i] = array1[i] - array2[i];
    }
}
//----------------------------------------------------//
void sarray_subtract_dbl(double *result, double *array1, double *array2, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        result[i] = array1[i] - array2[i];
    }
}
//*****************************************************//
//*****************************************************//
void sarray_add_scalar_int(int *array, int scalar, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        array[i] += scalar;
    }
}
//----------------------------------------------------//
void sarray_add_scalar_dbl(double *array, double scalar, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        array[i] += scalar;
    }
}
//*****************************************************//
//*****************************************************//
void sarray_subtract_scalar_int(int *array, int scalar, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        array[i] -= scalar;
    }
}
//----------------------------------------------------//
void sarray_subtract_scalar_dbl(double *array, double scalar, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        array[i] -= scalar;
    }
}
//*****************************************************//
//*****************************************************//
double sarray_avg_int(int *array, int nsize) {
    assert(nsize > 0);
    int i;
    double avg = 0.;
    for (i=0; i<nsize; i++) {
        avg += array[i];
    }
    return avg/(double) nsize;
}
//----------------------------------------------------//
double sarray_avg_dbl(double *array, int nsize) {
    assert(nsize > 0);
    int i;
    double avg = 0.;
    for (i=0; i<nsize; i++) {
        avg += array[i];
    }
    return avg/(double) nsize;
}
//*****************************************************//
//*****************************************************//
void sarray_scale_int(int *scaled_array, int *array, int scale, int nsize) {
    int i=0;
    for (i=0; i<nsize; i++) {
        scaled_array[i] = scale * array[i];
    }
}
//*****************************************************//
//returns 0 if not in array, 1 if found
int sarray_is_in_int(int *arr, int size, int val){
    for(int i=0;i<size;i++){
        if (arr[i] == val){
            return 1;
        }
    }
    return 0;
}
//*****************************************************//
//----------------------------------------------------//
void sarray_scale_dbl(double *scaled_array, double *array, double scale, int nsize) {
    int i=0;
    for (i=0; i<nsize; i++) {
        scaled_array[i] = scale * array[i];
    }
}
//*****************************************************//
//*****************************************************//
void sarray_scale_replace_int(int *array, int scale, int nsize) {
    int i=0;
    for (i=0; i<nsize; i++) {
        array[i] = scale * array[i];
    }
}
//----------------------------------------------------//
//void sarray_scale_replace_dbl(double *array, double scale, int nsize) {
//    int i=0;
//    for (i=0; i<nsize; i++) {
//        array[i] = scale * array[i];
//    }
//}
/*!
 \brief Multiplies all entries of a (double) vector by a scalar
 
 \param n the length of the vectors (int)
 \param scale scalar multiplier (double)
 \param array pointer to the vector (double, length n)
 */
void sarray_scale_replace_dbl(double *array, double scale, int nsize)
{
#ifdef _BLAS
    int IONE = 1;
    
    dscal_(&nsize, &scale, array, &IONE);
#else
    int i = 0;          /* loop counter */
    
    /* multiplies by a scalar */
    for(i = 0; i < nsize; i++)
        array[i] *= scale;
#endif
    return;
}

/*!
 \brief Performs daxpy operation for (double) vectors, \f$y = \alpha x  + y\f$
 
 \param n the length of the vectors (int)
 \param alpha scalar multiplier (double)
 \param x pointer to the first vector (double, length \a n)
 \param y pointer to the second vector (double, length \a n), changed on output
 */
void sarray_y_plus_ax_dbl(double *y, double alpha, double *x,  int n)
{
#ifdef _BLAS
    int IONE = 1;
    
    daxpy_(&n, &alpha, x, &IONE, y, &IONE);
#else
    int i = 0;          /* loop counter */
    
    /* adds the arrays */
    for(i = 0; i < n; i++)
        y[i] += alpha * x[i];
#endif
    return;
}


int sarray_unique_int(int *arr, int size){
    int unique_size = 1; // We'll start with the assumption that the first element is unique
    for (int i = 1; i < size; i++) {
        int is_unique = 1;
        for (int j = 0; j < unique_size; j++) {
            if (arr[i] == arr[j]) {
                is_unique = 0; // We found a duplicate
                break;
            }
        }
        if (is_unique) {
            arr[unique_size] = arr[i];
            unique_size++;
        }
    }
    return unique_size;
}

int sarray_argsort_int(int *my_int_arr, int *my_index_arr, int size){

    bool switched;
    int temp1;
    //fill my_index array (assumes it is already allocated)
    for(int i=0;i<size;i++){my_index_arr[i]=i;}
do
{
    switched = false;
    for(int i = 1; i < size; i++)
    {
        if(my_int_arr[my_index_arr[i-1]] > my_int_arr[my_index_arr[i]])
        {
            temp1 = my_index_arr[i];
            my_index_arr[i] = my_index_arr[i - 1];
            my_index_arr[i - 1] = temp1;
            
            switched = true;
        }
    }
}
while(switched);
//Also switch actual array
int temp[size];
for(int i=0;i<size;i++){temp[i] = my_int_arr[i];}
for(int i=0;i<size;i++){my_int_arr[i] =temp[my_index_arr[i]];}


    return 0;
}

int sarray_reverse_argsort_int(int *my_int_arr, int *my_index_arr, int size){

    bool switched;
    int temp1;
    //fill my_index array (assumes it is already allocated)
    for(int i=0;i<size;i++){my_index_arr[i]=i;}
do
{
    switched = false;
    for(int i = 1; i < size; i++)
    {
        if(my_int_arr[my_index_arr[i-1]] < my_int_arr[my_index_arr[i]])
        {
            temp1 = my_index_arr[i];
            my_index_arr[i] = my_index_arr[i - 1];
            my_index_arr[i - 1] = temp1;
            
            switched = true;
        }
    }
}
while(switched);
//Also switch actual array
int temp[size];
for(int i=0;i<size;i++){temp[i] = my_int_arr[i];}
for(int i=0;i<size;i++){my_int_arr[i] =temp[my_index_arr[i]];}


    return 0;
}

int sarray_shuffle_int(int *my_int_arr, int *my_index_arr, int size){
    //takes an array and shuffles the array according to index_arr
    int temp[size];
    for (int i = 0;i<size; i++){temp[i]=my_int_arr[i];}
    for(int i = 0; i<size; i++){
        my_int_arr[i] = temp[my_index_arr[i]];
    }

    return 0;
}

//*****************************************************//
//*****************************************************//
void sarray_integrity_check_int(int *array, int nsize, int linenumber, char *filename) {
    // nothing for now
}
//----------------------------------------------------//
void sarray_integrity_check_dbl(double *array, int nsize, int linenumber, char *filename) {
    
    //Is_DoubleArray_Inf_or_NaN(array, nsize, filename, linenumber);
}
//*****************************************************//
//*****************************************************//
void sarray_printScreen_int(int *array, int nsize, char *name) {
    int i=0;
    for (i=0; i<nsize; i++) {
        printf("%s[%d]: %d\n",name,i,array[i]);
    }
}
//----------------------------------------------------//
void sarray_printScreen_dbl(double *array, int nsize, char *name) {
    int i=0;
    for (i=0; i<nsize; i++) {
        printf("%s[%d]: %20.10e\n",name,i,array[i]);
    }
}

//*****************************************************//
//*****************************************************//

void sarray_vector_matrix_multiply(double *v, double **matrix, double *vout, int ndim) {
    int i,j,k;
    for(i = 0; i < ndim; i++) {
        vout[i] = 0;
        for(j = 0; j < ndim; j++) {
            vout[i] += matrix[i][j] * v[j];
        }
    }
}

//*****************************************************//
//*****************************************************//

void sarray_matrix_matrix_multiply(double **matrix1, double **matrix2, double **matrix_out, int ndim1, int ndim2) {
    int i,j,k;
    
    for(i = 0; i < ndim1; i++ ){
        for(j = 0; j < ndim2; j++){
            matrix_out[i][j] = 0;
            for(k = 0; k < ndim1; k++){
                matrix_out[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
}

//*****************************************************//
//*****************************************************//
void sarray_print_int(SGRID *grid, FILE * fp_out1, FILE *fp_out2, int *data, int **ndata, int my_nnode_max, int *my_nnode_ext, int flag){
    int i;                        /* loop counter */
    int ip, max_nnode;
    int *gdata, *edata;
    int ierr;
    FILE *fp_out;
#ifdef _MESSG
    MPI_Status *msg_status= grid->smpi->msg_status;
#endif
    
    if (grid->smpi->myid <= 0)
    {
        
        gdata = (int *) tl_alloc(sizeof(int), grid->macro_nnodes);
        edata = (int *) tl_alloc(sizeof(int), my_nnode_max);
        for (i = 0; i < my_nnode_ext[0]; i++){
            gdata[ndata[0][i]] = data[i];
        }
#ifdef _MESSG
        
        for (ip = 1; ip < grid->smpi->npes; ip++)
        {
            ierr = MPI_Recv(edata, my_nnode_ext[ip], MPI_INT, ip, 999, grid->smpi->ADH_COMM, msg_status);
            if (ierr != MPI_SUCCESS)
                messg_err(ierr);
            
            for (i = 0; i < my_nnode_ext[ip]; i++){
                gdata[ndata[ip][i]] = edata[i];
            }
        }
#endif
        for (ip=0;ip<=flag;ip++){
            if(ip==0){
                fp_out=fp_out1;
                max_nnode=grid->orig_macro_nnodes;
            }else{
                fp_out=fp_out2;
                max_nnode=grid->macro_nnodes;
            }
            for (i = 0; i < max_nnode; i++){
                fprintf(fp_out, "%d\n", gdata[i]);
            }
        }
        
        gdata = (int *) tl_free(sizeof(int), grid->macro_nnodes, gdata);
        edata = (int *) tl_free(sizeof(int), my_nnode_max, edata);
    }
    else
    {
#ifdef _MESSG
        ierr = MPI_Send(data, my_nnode_ext[grid->smpi->myid], MPI_INT, 0, 999, grid->smpi->ADH_COMM);
        if (ierr != MPI_SUCCESS)
            messg_err(ierr);
#endif
    }
}
//----------------------------------------------------//
void sarray_print_dbl(SGRID *grid, FILE * fp_out1, FILE *fp_out2, double *data, int **ndata, int my_nnode_max, int *my_nnode_ext, int flag){
    int i;                        /* loop counter */
    int ip, max_nnode;
    double *gdata, *edata;
    int ierr;
    FILE *fp_out;
#ifdef _MESSG
    MPI_Status *msg_status= grid->smpi->msg_status;
#endif
    
    if (grid->smpi->myid <= 0)
    {
        
        gdata = (double *) tl_alloc(sizeof(double), grid->macro_nnodes);
        edata = (double *) tl_alloc(sizeof(double), my_nnode_max);
        for (i = 0; i < my_nnode_ext[0]; i++){
            gdata[ndata[0][i]] = data[i];
        }
#ifdef _MESSG
        
        for (ip = 1; ip < grid->smpi->npes; ip++)
        {
            ierr = MPI_Recv(edata, my_nnode_ext[ip], MPI_DOUBLE, ip, 999, grid->smpi->ADH_COMM, msg_status);
            if (ierr != MPI_SUCCESS)
                messg_err(ierr);
            
            for (i = 0; i < my_nnode_ext[ip]; i++){
                gdata[ndata[ip][i]] = edata[i];
            }
        }
#endif
        for (ip=0;ip<=flag;ip++){
            if(ip==0){
                fp_out=fp_out1;
                max_nnode=grid->orig_macro_nnodes;
            }else{
                fp_out=fp_out2;
                max_nnode=grid->macro_nnodes;
            }
            for (i = 0; i < max_nnode; i++){
                fprintf(fp_out, "%16.8e\n", gdata[i]);
            }
        }
        
        gdata = (double *) tl_free(sizeof(double), grid->macro_nnodes, gdata);
        edata = (double *) tl_free(sizeof(double), my_nnode_max, edata);
    }
    else
    {
#ifdef _MESSG
        
        ierr = MPI_Send(data, my_nnode_ext[grid->smpi->myid], MPI_DOUBLE, 0, 999, grid->smpi->ADH_COMM);
        if (ierr != MPI_SUCCESS)
            messg_err(ierr);
#endif
    }
}
//*****************************************************//
//*****************************************************//
void sarray_2column_print_int(FILE *fout, int *array1, int *array2, int size){
    int i=0;
    for (i=0; i<size; i++) {
        fprintf(fout, "%d %d\n", array1[i],array2[i]);
    }
}
//----------------------------------------------------//
void sarray_2column_print_dbl(FILE *fout, double *array1, double *array2, int size){
    int i=0;
    for (i=0; i<size; i++) {
        fprintf(fout, "%16.8e %16.8e\n", array1[i],array2[i]);
    }
}
//*****************************************************//
//*****************************************************//
void sarray_3column_print_int(FILE *fout, int *array1, int *array2, int *array3, int size){
    int i=0;
    for (i=0; i<size; i++) {
        fprintf(fout, "%d %d %d\n", array1[i],array2[i],array3[i]);
    }
}
//----------------------------------------------------//
void sarray_3column_print_dbl(FILE *fout, double *array1, double *array2, double *array3, int size){
    int i=0;
    for (i=0; i<size; i++) {
        fprintf(fout, "%16.8e %16.8e %16.8e\n", array1[i],array2[i],array3[i]);
    }
}

//*****************************************************//
//*****************************************************//

int sarray_max_int(int a[], int num_elements) {
    int i, max = a[0];
    for (i=1; i<num_elements; i++) {
        if (a[i] > max) max = a[i];
    }
    return (max);
}

double sarray_max_dbl(double a[], int num_elements) {
    int i;
	double max = a[0];
    for (i=1; i<num_elements; i++) {
        if (a[i] > max) max = a[i];
    }
    return (max);
}
// Function to find the index of the minimum element in an integer array.
// Returns -1 if the array is NULL or has zero elements.
int sarray_argmin_int(int *arr, int n) {
    if (arr == NULL || n <= 0) {
        return -1; // Indicate an invalid input.
    }

    int min_index = 0;
    int min_value = arr[0];

    for (int i = 1; i < n; i++) {
        if (arr[i] < min_value) {
            min_value = arr[i];
            min_index = i;
        }
    }

    return min_index;
}
