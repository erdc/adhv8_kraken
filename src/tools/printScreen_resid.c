/* This routine prints out the residuals by node and equation */
/* called from resid routines in a DEBUG section */

#include "global_header.h"

void printScreen_resid(char * descript, double *resid, int nnode, int nsys, int linenumber, char *filename)
{
    int i, j, k;
    FILE *fp = stdout;
#ifdef _MESSG
    int ierr_code, myid=0;
    ierr_code = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    char fn[30+1];
    snprintf(fn, 30, "residual_pe%d.txt", myid);
    fp = fopen(fn, "w");
#endif
    fflush(fp);

    fprintf(fp,"\n");
    fprintf(fp,"printing residual: %s @ line %s:%d \n",descript,filename, linenumber);
    for (i = 0; i < nnode; i++) {
        j = i * nsys;
        fprintf(fp,"node: %-8d \t ", (i + 1));
        if (nsys == 1) {
            fprintf(fp,"  Rc_eq: %+12.10e ", resid[j + 0]);
        } else if (nsys == 2) {
            fprintf(fp," Rmx_eq: %+-12.10e ", resid[j + 0]);
            fprintf(fp,"  Rc_eq: %+-12.10e ", resid[j + 1]);
        } else if (nsys == 3) {
            fprintf(fp," Rmx_eq: %+-12.10e ", resid[j + 0]);
            fprintf(fp," Rmy_eq: %+-12.10e ", resid[j + 1]);
            fprintf(fp,"  Rc_eq: %+-12.10e ", resid[j + 2]);
        } else if (nsys == 4)  {
            fprintf(fp," Rmx_eq: %+-12.10e ", resid[j + 0]);
            fprintf(fp," Rmy_eq: %+-12.10e ", resid[j + 1]);
            fprintf(fp," Rmz_eq: %+-12.10e ", resid[j + 2]);
            fprintf(fp,"  Rc_eq: %+-12.10e ", resid[j + 3]);
        }
        
        fprintf(fp,"\n");
    }
}
