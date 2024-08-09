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
void sarray_copy_dbl(double *to, double *from, int nsize) {
    int i;
    for (i=0; i<nsize; i++) {
        to[i] = from[i];
    }
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
void sarray_scale_replace_dbl(double *array, double scale, int nsize) {
    int i=0;
    for (i=0; i<nsize; i++) {
        array[i] = scale * array[i];
    }
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
