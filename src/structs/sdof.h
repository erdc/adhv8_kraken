#ifndef H_SDOF_
#define H_SDOF_


/************************************************************/
/* structures ----------------------------------------------*/

typedef struct {
  double x_eq, y_eq, z_eq, c_eq;  /* the particular equation x,y,z, or continuity */
} DOF_4;

typedef struct {
  double x_eq, y_eq, c_eq;  /* the particular equation x,y or continuity */
} DOF_3;

typedef struct {
  double x_eq, y_eq;        /* the particular equation x or y */
} DOF_2;


/*********************************************************/
/* struct methods -------------------------------------- */

void dof2_init(DOF_2 *);
void dof3_init(DOF_3 *);
void dof4_init(DOF_4 *);
void dof2_init_array(DOF_2 *, int);
void dof3_init_array(DOF_3 *, int);
void dof4_init_array(DOF_4 *, int);
void dof2_copy(DOF_2 *, DOF_2);
void dof3_copy(DOF_3 *, DOF_3);
void dof4_copy(DOF_4 *, DOF_4);
void dof2_copy_array(DOF_2 *, DOF_2 *, int);
void dof3_copy_array(DOF_3 *, DOF_3 *, int);
void dof4_copy_array(DOF_4 *, DOF_4 *, int);
void dof2_printScreen(DOF_2);
void dof3_printScreen(DOF_3);
void dof4_printScreen(DOF_4, int, int);
void dof2_printScreen_array(char *, DOF_2 *, int, int, char *);
void dof3_printScreen_array(char *, DOF_3 *, int, int, char *);
void dof4_printScreen_array(char *, DOF_4 *, int, int, char *, int , int);
DOF_2 dof2_subtract(DOF_2, DOF_2);
DOF_3 dof3_subtract(DOF_3, DOF_3);
DOF_4 dof4_subtract(DOF_4, DOF_4);
void dof2_subtract_array(DOF_2 *, DOF_2 *, DOF_2 *, int);
void dof3_subtract_array(DOF_3 *, DOF_3 *, DOF_3 *, int);
void dof4_subtract_array(DOF_4 *, DOF_4 *, DOF_4 *, int);
void dof3_add_replace_array(DOF_3 *dof1, DOF_3 *dof2, int array_size);
void dof1_printScreen(double);
void dof1_printScreen_array(char *, double *, int, int, char *);
void dof3_debug(DOF_3 *dof, int n, char *filename, int linenumber);
void dof4_debug(DOF_4 *dof, int n, char *filename, int linenumber);
int dof3_debug_continue(DOF_3 *dof, int n);
int dof4_debug_continue(DOF_4 *dof, int n);

#endif

