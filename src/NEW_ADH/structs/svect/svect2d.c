#include "adh.h"

/**********************************************/
/**********************************************/

void svect2d_init(SVECT2D *v) {
    v->x = 0.;
    v->y = 0.;
}

void svect2d_init2one(SVECT2D *v) {
    v->x = 1.;
    v->y = 1.;
}

void svect2d_init_array(SVECT2D *v, int size) {
    int i=0;
    for (i=0; i<size; i++) {
        v[i].x = 0.;
        v[i].y = 0.;
    }
}

void svect2d_init_array_value_range(SVECT2D *v, double xValue, double yValue, int start, int end) {
    int i=0;
    for (i=start; i<end; i++) {
        v[i].x = xValue;
        v[i].y = yValue;
    }
}


void svect2d_init2one_array(SVECT2D *v, int size) {
    int i=0;
    for (i=0; i<size; i++) {
        v[i].x = 1.;
        v[i].y = 1.;
    }
}

/**********************************************/
/**********************************************/

void svect2d_print(FILE *fout, SVECT2D v) {
    
    if (fout == NULL) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> trying to print to bad file!");
    }
    else {
        fprintf(fout, "%16.8e %16.8e %16.8e\n", v.x, v.y, 0.);
    }
}

void svect2d_print_array(FILE *fout, SVECT2D *v, int nnodes) {
    int inode;
    if (v == NULL) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> trying to print non-allocated array!");
    }
    
    for(inode=0; inode < nnodes; inode++) {
        svect2d_print(fout, v[inode]);
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes a 2D AdH vector array to screen with global IDS
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in] descript variable name string
 *  @param[in] f the vector array
 *  @param[in] n the number of array elements
 *  @param[in] global_nd_ids the global IDs of the array nodes
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void printScreen_debug_svec2d(char *descript, SVECT2D *v, int n, int *global_nd_ids) {
    int i;
    for (i=0; i<n; i++) printf("%s[%d] {x,y} gnode id: %d || {%30.20e,%30.20e}\n",descript,i,global_nd_ids[i],v[i].x,v[i].y);
}


void svect2d_print_array_MPI(SGRID *grid, FILE * fp_out1, FILE *fp_out2, SVECT2D *data, int **ndata, int my_nnode_max, int *my_nnode_ext, int flag)
{
//    int i;                        /* loop counter */
//    int ip, flag2=0;
//    SVECT2D *gdata, *edata;
//    int ierr;
//    int myid = grid->smpi->myid;
//    int npes = grid->smpi->npes;
//    int global_nnode, max_nnode;
//    FILE *fp_out;
//#ifdef _MESSG
//    MPI_Status * msg_status = grid->smpi->msg_status;
//
//    MPI_Datatype MPI_VECT2D;
//    MPI_Datatype Dt_type[2] = { MPI_DOUBLE, MPI_DOUBLE };
//    int Dt_block[2] = { 1, 1 };
//    MPI_Aint Dt_disp[2];
//
//    ierr = MPI_Get_address(data, &Dt_disp[0]);
//    if (ierr != MPI_SUCCESS)
//        messg_err(ierr);
//    ierr = MPI_Get_address(&data[0].y, &Dt_disp[1]);
//    if (ierr != MPI_SUCCESS)
//        messg_err(ierr);
//    for (i = 1; i >= 0; i--)
//        Dt_disp[i] -= Dt_disp[0];
//    ierr = MPI_Type_create_struct(2, Dt_block, Dt_disp, Dt_type, &MPI_VECT2D);
//    if (ierr != MPI_SUCCESS)
//        messg_err(ierr);
//    ierr = MPI_Type_commit(&MPI_VECT2D);
//    if (ierr != MPI_SUCCESS)
//        messg_err(ierr);
//#endif
//    if(grid->ndim==3){
//        global_nnode=grid->macro_nnodes_sur;
//    }else{
//        global_nnode=grid->macro_nnodes;
//    }
//    if (myid <= 0)
//    {
//        max_nnode=0;
//      	for(i=0;i<grid->smpi->npes;i++)max_nnode+=my_nnode_ext[i];
//      	if((grid->ndim==3) && (max_nnode<grid->macro_nnodes)) flag2=1;
//        gdata = (SVECT2D *) tl_alloc(sizeof(SVECT2D), global_nnode);
//        edata = (SVECT2D *) tl_alloc(sizeof(SVECT2D), my_nnode_max);
//        for (i = 0; i < my_nnode_ext[0]; i++)
//        {
//            gdata[ndata[0][i]].x = data[i].x;
//            gdata[ndata[0][i]].y = data[i].y;
//
//        }
//#ifdef _MESSG
//        for (ip = 1; ip < npes; ip++)
//        {
//            ierr = MPI_Recv(edata, my_nnode_ext[ip], MPI_VECT2D, ip, 999, grid->smpi->ADH_COMM, msg_status);
//            if (ierr != MPI_SUCCESS)
//                messg_err(ierr);
//
//            for (i = 0; i < my_nnode_ext[ip]; i++){
//                gdata[ndata[ip][i]].x = edata[i].x;
//                gdata[ndata[ip][i]].y = edata[i].y;
//            }
//        }
//#endif
//        for (ip=0;ip<=flag;ip++){
//            if(ip==0){
//                fp_out=fp_out1;
//                if(flag2==0){max_nnode=grid->orig_macro_nnodes;}
//          			else{max_nnode=grid->orig_macro_nnodes_sur;}
//            }else{
//                fp_out=fp_out2;
//								if(flag2==0){max_nnode=grid->macro_nnodes;}
//          			else{max_nnode=grid->macro_nnodes_sur;}
//            }
//            svect2d_print_array(fp_out, gdata, max_nnode);
//        }
//
//        gdata = (SVECT2D *) tl_free(sizeof(SVECT2D), global_nnode, gdata);
//        edata = (SVECT2D *) tl_free(sizeof(SVECT2D), my_nnode_max, edata);
//    }
//    else
//    {
//#ifdef _MESSG
//        ierr = MPI_Send(data, my_nnode_ext[myid], MPI_VECT2D, 0, 999, grid->smpi->ADH_COMM);
//        if (ierr != MPI_SUCCESS)
//            messg_err(ierr);
//
//#endif
//    }
//#ifdef _MESSG
//    ierr = MPI_Type_free(&MPI_VECT2D);
//    if (ierr != MPI_SUCCESS)
//        messg_err(ierr);
//#endif
}

/**********************************************/
/**********************************************/

void svect2d_printScreen(SVECT2D v) {
    printf("x: %30.20e \t y: %30.20e\n", v.x, v.y);
}

void svect2d_printScreen_array(char * descript, SVECT2D *v, int size, int linenumber, char *filename) {
    int i;
    printf("\n");
    printf("printing vector: %s @ line %s:%d :: size: %d \n",descript, filename, linenumber,  size);
    for (i=0; i<size; i++) {
        printf("[%d] ",i);
        svect2d_printScreen(v[i]);
    }
}

/**********************************************/
/**********************************************/

void svect2d_copy_array(SVECT2D *vdest, SVECT2D *vorig, int nnodes) {
    int inode;
    if (vdest == NULL) {
        vorig = NULL;
        return;
    }
    for(inode=0; inode < nnodes; inode++) {
        vdest[inode].x = vorig[inode].x;
        vdest[inode].y = vorig[inode].y;
    }
}
/**********************************************/
/**********************************************/

double svect2d_mag(SVECT2D v) {
    return sqrt(v.x * v.x + v.y * v.y);
}
double svect2d_mag_safe(SVECT2D v) {
        return sqrt(v.x * v.x + v.y * v.y  + SMALL);
}

SVECT2D svect2d_subtract(SVECT2D v1, SVECT2D v2) {
    SVECT2D v3;
    v3.x = v1.x - v2.x;
    v3.y = v1.y - v2.y;
    return v3;
}

SVECT2D svect2d_add(SVECT2D v1, SVECT2D v2) {
    SVECT2D v3;
    v3.x = v1.x + v2.x;
    v3.y = v1.y + v2.y;
    return v3;
}

void svect2d_add_array(SVECT2D *v, SVECT2D *v1, SVECT2D *v2, int size) {
    int i = 0;
    for (i = 0; i < size; i++) {
        v[i] = svect2d_add(v1[i], v2[i]);
    }
}

SVECT2D svect2d_scale(SVECT2D v1, double scale) {
    SVECT2D v2;
    v2.x = scale * v1.x;
    v2.y = scale * v1.y;
    return v2;
}

void svect2d_scale_array(SVECT2D *v, int size, double scale) {
    int i = 0;
    for (i = 0; i < size; i++) {
        v[i].x *= scale;
        v[i].y *= scale;
    }
}

void svect2d_scale_replace_array(SVECT2D *v, double scale, int size) {
    int i;
    for (i=0; i<size; i++) {
        v[i].x *= scale;
        v[i].y *= scale;
    }
}

void svect2d_nscale_array(SVECT2D *v, int size, double *scale) {
    int i = 0;
    for (i = 0; i < size; i++) {
        v[i].x *= scale[i];
        v[i].y *= scale[i];
    }
}

double svect2d_dotp(SVECT2D v1, SVECT2D v2) {
    return (v1.x * v2.x + v1.y * v2.y);
}

SVECT2D svect2d_avg(SVECT2D v1, SVECT2D v2) {
    SVECT2D v3;
    v3.x = 0.5 * (v1.x + v2.x);
    v3.y = 0.5 * (v1.y + v2.y);
    return v3;
}

SVECT2D svect2d_average_array(SVECT2D *vect, int size) {
    SVECT2D avg; svect2d_init(&avg);
    int i = 0;
    for (i=0; i<size; i++) {
        avg.x += vect[i].x; avg.y += vect[i].y;
    }
    avg.x *= (1./(double) size);  // apparently, doing this instead of dividing
    avg.y *= (1./(double) size);  // matters a whole lot, and stabilized source test
    return avg;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     return two double arrays from a 2D AdH vector array
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] u a double array of x-velocities
 *  @param[out] v a double array of y-velocities
 *  @param[in] vel the 2D vector array
 *  @param[in] nnodes the total number of elements in the 2D vector array
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void dumpVector2D(SVECT2D *vel, int nnodes, double *u, double *v) {
    int i;
    for (i=0; i<nnodes; i++) {
        u[i] = vel[i].x;
        v[i] = vel[i].y;
    }
}

/**********************************************/
/**********************************************/

void svect2d_integrity_check(SVECT2D vect, int linenumber, char *filename) {
//    Is_Double_Inf_or_NaN(vect.x, filename, linenumber);
//    Is_Double_Inf_or_NaN(vect.y, filename, linenumber);
}

void svect2d_integrity_check_array(SVECT2D *vect, int nsize, int linenumber, char *filename) {
    int i=0;
    for (i=0; i<nsize; i++) {
        svect2d_integrity_check(vect[i], linenumber, filename);
    }
}
