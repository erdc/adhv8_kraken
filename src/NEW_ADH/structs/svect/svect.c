#include "adh.h"

/**********************************************/
/**********************************************/

void svect_init(SVECT *v) {
    v->x = 0.;
    v->y = 0.;
    v->z = 0.;
}

void svect_init_value(SVECT *v, double x, double y, double z) {
    v->x = x;
    v->y = y;
    v->z = z;
}

void svect_init_array(SVECT *v, int size) {
    int i=0;
    for (i=0; i<size; i++) {
        v[i].x = 0.;
        v[i].y = 0.;
        v[i].z = 0.;
    }
}


SVECT svect_avg(SVECT v1, SVECT v2) {
    SVECT v3;
    v3.x = 0.5 * (v1.x + v2.x);
    v3.y = 0.5 * (v1.y + v2.y);
    v3.z = 0.5 * (v1.z + v2.z);
    return v3;
}

/**********************************************/
/**********************************************/

void svect_print(FILE *fout, SVECT v) {
  if (fout == NULL) {
       printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
       tl_error(">> trying to print to bad file!");
  }
  else {
      fprintf(fout, "%15.8e %15.8e %15.8e\n", v.x, v.y, v.z);
  }
}

void svect_print_array(FILE *fout, SVECT *v, int nnodes) {
      int inode;
      if (v == NULL) {
          printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
          tl_error(">> trying to print non-allocated array!");
      }

      for(inode=0; inode < nnodes; inode++) {
          svect_print(fout, v[inode]);
      }
}

/**********************************************/
/**********************************************/

void svect_printScreen(SVECT v, char *name) {
    printf("\t %s.x: %30.20e \t %s.y: %30.20e \t %s.z: %30.20e\n", name, v.x, name, v.y, name, v.z);
}

void svect_printScreen_array(char * descript, SVECT *v, char *name, int size, int linenumber, char *filename) {
    int i;
    printf("\n");
    printf("printing vector: %s @ line %s:%d :: size: %d \n",descript, filename, linenumber, size); 
    for (i=0; i<size; i++) {
        printf("[%d] ",i);
        svect_printScreen(v[i],name);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes a 3D AdH vector array to screen
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in] descript variable name string
 *  @param[in] f the vector array
 *  @param[in] n the number of array elements
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void printScreen_debug_vec(char *descript, SVECT *f, int n) {
    int i;
    printf("%s || x || ",descript);
    for (i=0; i<n; i++) {
        printf("%30.20f \t ",f[i].x);
    }
    printf("\n");

    printf("%s || y || ",descript);
    for (i=0; i<n; i++) {
        printf("%30.20f \t ",f[i].y);
    }
    printf("\n");

    printf("%s || z || ",descript);
    for (i=0; i<n; i++) {
        printf("%30.20f \t ",f[i].z);
    }
    printf("\n");
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes a 3D AdH vector array to screen with global IDS
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
inline void printScreen_debug_svect(char *descript, SVECT *v, int n, int *global_nd_ids) {
    int i;
    for (i=0; i<n; i++) printf("%s[%d] {x,y,z} gnode id: %d || {%20.10e,%20.10e,%20.10e}\n",descript,i,global_nd_ids[i],v[i].x,v[i].y,v[i].z);
}


void svect_print_array_MPI(SGRID *grid, FILE * fp_out1, FILE *fp_out2, SVECT *data, int **ndata, int my_nnode_max, int *my_nnode_ext, int flag)
{
//  int i;                        /* loop counter */
//  int ip, max_nnode, flag2=0;
//  SVECT *gdata, *edata;
//  int ierr;
//  int myid = grid->smpi->myid;
//  int npes = grid->smpi->npes;
//  FILE *fp_out;
//#ifdef _MESSG
//  MPI_Status * msg_status = grid->smpi->msg_status;
//
//  MPI_Datatype MPI_VECT;
//  MPI_Datatype Dt_type[3] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
//  int Dt_block[3] = { 1, 1, 1 };
//  MPI_Aint Dt_disp[3];
//
//  ierr = MPI_Get_address(data, &Dt_disp[0]);
//  if (ierr != MPI_SUCCESS)
//    messg_err(ierr);
//  ierr = MPI_Get_address(&data[0].y, &Dt_disp[1]);
//  if (ierr != MPI_SUCCESS)
//    messg_err(ierr);
//  ierr = MPI_Get_address(&data[0].z, &Dt_disp[2]);
//  if (ierr != MPI_SUCCESS)
//    messg_err(ierr);
//  for (i = 2; i >= 0; i--)
//    Dt_disp[i] -= Dt_disp[0];
//  ierr = MPI_Type_create_struct(3, Dt_block, Dt_disp, Dt_type, &MPI_VECT);
//  if (ierr != MPI_SUCCESS)
//    messg_err(ierr);
//  ierr = MPI_Type_commit(&MPI_VECT);
//  if (ierr != MPI_SUCCESS)
//    messg_err(ierr);
//#endif
//
//  if (myid <= 0)
//    {
//      max_nnode=0;
//      for(i=0;i<grid->smpi->npes;i++)max_nnode+=my_nnode_ext[i];
//      if((grid->ndim==3) && (max_nnode<grid->macro_nnodes)) flag2=1;
//      gdata = (SVECT *) tl_alloc(sizeof(SVECT), grid->macro_nnodes);
//      edata = (SVECT *) tl_alloc(sizeof(SVECT), my_nnode_max);
//      for (i = 0; i < my_nnode_ext[0]; i++)
//      {
//        gdata[ndata[0][i]].x = data[i].x;
//        gdata[ndata[0][i]].y = data[i].y;
//        gdata[ndata[0][i]].z = data[i].z;
//      }
//#ifdef _MESSG
//      for (ip = 1; ip < npes; ip++)
//        {
//          ierr = MPI_Recv(edata, my_nnode_ext[ip], MPI_VECT, ip, 999, grid->smpi->ADH_COMM, msg_status);
//          if (ierr != MPI_SUCCESS)
//            messg_err(ierr);
//          for (i = 0; i < my_nnode_ext[ip]; i++){
//            gdata[ndata[ip][i]].x = edata[i].x;
//            gdata[ndata[ip][i]].y = edata[i].y;
//            gdata[ndata[ip][i]].z = edata[i].z;
//          }
//        }
//#endif
//      for (ip=0;ip<=flag;ip++){
//        if(ip==0){
//          fp_out=fp_out1;
//          if(flag2==0){max_nnode=grid->orig_macro_nnodes;}
//          else{max_nnode=grid->orig_macro_nnodes_sur;}
//        }else{
//          fp_out=fp_out2;
//          if(flag2==0){max_nnode=grid->macro_nnodes;}
//          else{max_nnode=grid->macro_nnodes_sur;}
//        }
//        svect_print_array(fp_out, gdata, max_nnode);
//      }
//
//      gdata = (SVECT *) tl_free(sizeof(SVECT), grid->macro_nnodes, gdata);
//      edata = (SVECT *) tl_free(sizeof(SVECT), my_nnode_max, edata);
//    }
//  else
//    {
//#ifdef _MESSG
//      ierr = MPI_Send(data, my_nnode_ext[myid], MPI_VECT, 0, 999, grid->smpi->ADH_COMM);
//      if (ierr != MPI_SUCCESS)
//        messg_err(ierr);
//#endif
//    }
//#ifdef _MESSG
//  ierr = MPI_Type_free(&MPI_VECT);
//  if (ierr != MPI_SUCCESS)
//    messg_err(ierr);
//
//#endif
}
/**********************************************/
/**********************************************/

void svect_copy_array(SVECT *vdest, SVECT *vorig, int nnodes) {
      int inode;
      if (vdest == NULL) {
          vorig = NULL;
          return;
      }

      for(inode=0; inode < nnodes; inode++) {
          vdest[inode].x = vorig[inode].x;
          vdest[inode].y = vorig[inode].y;
          vdest[inode].z = vorig[inode].z;
      }
}

void svect_copy(SVECT *vdest, SVECT vorig) {
      vdest->x = vorig.x;
      vdest->y = vorig.y;
      vdest->z = vorig.z;
}

/**********************************************/
/**********************************************/
double svect_mag(SVECT v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

double svect_mag_safe(SVECT v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z + NOT_QUITE_SMALL);
}

SVECT svect_subtract(SVECT v1, SVECT v2) {
    SVECT v3;
    v3.x = v1.x - v2.x;
    v3.y = v1.y - v2.y;
    v3.z = v1.z - v2.z;
    return v3;
}


SVECT svect_add(SVECT v1, SVECT v2) {
    SVECT v3;
    v3.x = v1.x + v2.x;
    v3.y = v1.y + v2.y;
    v3.z = v1.z + v2.z;
    return v3;
}

SVECT svect_scale(SVECT v1, double scale) {
    SVECT v2;
    v2.x = scale * v1.x;
    v2.y = scale * v1.y;
    v2.z = scale * v1.z;
    return v2;
}

void svect_scale_replace(SVECT *v, double scale) {
    v->x *= scale;
    v->y *= scale;
    v->z *= scale;
}

void svect_scale_replace_array(SVECT *v, double scale, int size) {
    int i;
    for (i=0; i<size; i++) {
        v[i].x *= scale;
        v[i].y *= scale;
        v[i].z *= scale;
    }
}

SVECT svect_sum_array(SVECT *vect, int size) {
    SVECT sum; svect_init(&sum);
    int i = 0;
    for (i = 0; i < size; i++) {
        sum.x += vect[i].x;
        sum.y += vect[i].y;
        sum.z += vect[i].z;
    }
    return sum;
}

SVECT svect_subtract_array(SVECT *vect, int size) {
    SVECT sum; svect_init(&sum);
    int i = 0;
    for (i = 0; i < size; i++) {
        sum.x -= vect[i].x;
        sum.y -= vect[i].y;
        sum.z -= vect[i].z;
    }
    return sum;
}

SVECT svect_average_array(SVECT *vect, int size) {
    SVECT avg; svect_init(&avg);
    int i = 0;
    for (i = 0; i < size; i++) {
        avg.x += vect[i].x; avg.y += vect[i].y; avg.z += vect[i].z;
    }
    avg.x *= (1./((double) size));
    avg.y *= (1./((double) size));
    avg.z *= (1./((double) size));
    return avg;
}

void svect_add_array2(SVECT *v, SVECT *v1, SVECT *v2, int size) {
    svect_init(v);
    int i = 0;
    for (i = 0; i < size; i++) {
        v[i] = svect_add(v1[i], v2[i]);
    }
}

void svect_subtract_array2(SVECT *v, SVECT *v1, SVECT *v2, int size) {
    svect_init(v);
    int i = 0;
    for (i = 0; i < size; i++) {
        v[i] = svect_subtract(v1[i], v2[i]);
    }
}

double svect_dotp(SVECT v1, SVECT v2) {
    return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

void svect_dotp_array(SVECT *v1, SVECT *v2, int size, double *result) {
    int i = 0;
    for (i = 0; i < size; i++) {
        result[i] = v1[i].x * v2[i].x + v1[i].y * v2[i].y + v1[i].z * v2[i].z;
    }
}

SVECT svect_cross(SVECT v1, SVECT v2) {
    SVECT result;
    result.x = v1.y*v2.z - v1.z*v2.y;
    result.y = v1.z*v2.x - v1.x*v2.z;
    result.z = v1.x*v2.y - v1.y*v2.x;
    return result;
}

/**********************************************/
/**********************************************/

void svect_integrity_check(SVECT vect, int linenumber, char *filename) {
//    Is_Double_Inf_or_NaN(vect.x, filename, linenumber);
//    Is_Double_Inf_or_NaN(vect.y, filename, linenumber);
//    Is_Double_Inf_or_NaN(vect.z, filename, linenumber);
}

void svect_integrity_check_array(SVECT *vect, int nsize, int linenumber, char *filename) {
    int i=0;
    for (i=0; i<nsize; i++) {
        svect_integrity_check(vect[i], linenumber, filename);
    }
}
