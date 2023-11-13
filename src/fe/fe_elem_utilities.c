/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_elem_utilties.c This file collects functions that are useful for elemental operations */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes 3 degree of freedom elemental residual contribution to screen
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in] string name of the residual
 *  @param[in] nnodes the total number of nodes on the element
 *  @param[in] ie the global element ID
 *  @param[in] nodes an array of global node IDs on the element
 *  @param[in] elem_rhst the residual w/o contribution to print
 *  @param[in] elem_rhs the residual including the contribution to print
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void rhs_contrib_3dof(char *string, int nnodes, int ie, int *nodes, DOF_3 *elem_rhst, DOF_3 *elem_rhs) {
    int i;
    DOF_3 elem_rhs_diff[nnodes];
    dof3_subtract_array(elem_rhs_diff, elem_rhs, elem_rhst, nnodes);
    printf("\nRHS contribution from: %s ----------------------- \n",string);
    
    for (i=0; i<nnodes; i++) {
        printf("rhs contrib %s[%d] ie: %d node: %d :: ",string+1,i+1,ie+1,nodes[i]+1);
        dof3_printScreen(elem_rhs_diff[i]);
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes 3 degree of freedom elemental residual to screen
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in] string name of the residual
 *  @param[in] nnodes the total number of nodes on the element
 *  @param[in] ie the global element ID
 *  @param[in] nodes an array of global node IDs on the element
 *  @param[in] elem_rhs the residual to print
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void rhs_3dof(char *string, int nnodes, int ie, int *nodes, DOF_3 *elem_rhs) {
    int i;
    printf("\n***************************************\n%s\n",string);
    for (i=0; i<nnodes; i++) {
        printf("RHS contrib from %s[%d] ie: %d node: %-10d :: ",string,i+1,ie+1,nodes[i]+1);
        dof3_printScreen(elem_rhs[i]);
    }
    //printf("**************************************\n");
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes 4 degree of freedom elemental residual to screen
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in] string name of the residual
 *  @param[in] nnodes the total number of nodes on the element
 *  @param[in] ie the global element ID
 *  @param[in] nodes an array of global node IDs on the element
 *  @param[in] elem_rhs the residual to print
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void rhs_4dof(char *string, int nnodes, int ie, int *nodes, DOF_4 *elem_rhs, int printFieldWidth, int printPrecision) {
    int i;
    printf("\n***************************************\n%s\n",string);
    for (i=0; i<nnodes; i++) {
        printf("RHS contrib from %s[%d] ie: %d node: %-10d :: ",string,i+1,ie+1,nodes[i]+1);
        dof4_printScreen(elem_rhs[i],printFieldWidth,printPrecision);
    }
    //printf("**************************************\n");
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes 1 degree of freedom elemental residual to screen
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in] string name of the residual
 *  @param[in] nnodes the total number of nodes on the element
 *  @param[in] ie the global element ID
 *  @param[in] nodes an array of global node IDs on the element
 *  @param[in] elem_rhst the residual w/o contribution to print
 *  @param[in] elem_rhs the residual including the contribution to print
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void rhs_contrib_1dof(char *string, int nnodes, int ie, int *nodes, double *elem_rhst, double *elem_rhs) {
    int i;
    double elem_rhs_diff[nnodes];
    sarray_subtract_dbl(elem_rhs_diff, elem_rhs, elem_rhst, nnodes);
    printf("\nRHS contribution from: %s ----------------------- \n",string);
    
    for (i=0; i<nnodes; i++) {
        printf("rhs contrib %s[%d] ie: %d node: %d :: ",string,i+1,ie+1,nodes[i]);
        dof1_printScreen(elem_rhs_diff[i]);
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes 1 degree of freedom elemental residual to screen
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in] string name of the residual
 *  @param[in] nnodes the total number of nodes on the element
 *  @param[in] ie the global element ID
 *  @param[in] nodes an array of global node IDs on the element
 *  @param[in] elem_rhs the residual to print
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void rhs_1dof(char *string, int nnodes, int ie, int *nodes, double *elem_rhs) {
    int i;
    printf("\n%s ****************************\n",string);
    for (i=0; i<nnodes; i++) {
        printf("RHS contrib from %s[%d] ie: %d node: %d :: ",string,i+1,ie+1,nodes[i]);
        dof1_printScreen(elem_rhs[i]);
    }
    //printf("**************************************\n");
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes an integer array to screen
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in] descript variable name string
 *  @param[in] f the integer array
 *  @param[in] n the number of array elements
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void printScreen_debug_int(char *descript, int *f, int n) {
    int i;
    printf("%s || ",descript);
    for (i=0; i<n; i++) {
        printf("%d \t ",f[i]);
    }
    printf("\n");
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes a double array to screen
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in] descript variable name string
 *  @param[in] f the integer array
 *  @param[in] n the number of array elements
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void printScreen_debug_dbl(char *descript, double *f, int n) {
    int i;
    printf("%s || ",descript);
    for (i=0; i<n; i++) {
        printf("%30.20e \t ",f[i]);
    }
    printf("\n");
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes a double array to screen with global IDs
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in] descript variable name string
 *  @param[in] f the integer array
 *  @param[in] n the number of array elements
 *  @param[in] global_nd_ids the global IDs of the array nodes
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void printScreen_debug2_dbl(char *descript, double *f, int n, int *global_nd_ids) {
    int i;
    for (i=0; i<n; i++) printf("%s[%d] || gnode id: %d || %40.30e\n",descript,i,global_nd_ids[i],f[i]);
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


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes a 2D AdH vector to screen
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
inline void printScreen_debug_vec2d(char *descript, SVECT2D *f, int n) {
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
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates a scaled array of second-order temporal values at t + (3/2)*dt
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] local1 an array for output
 *  @param[in] local2 a double array for variable 1
 *  @param[in] local3 a double array for variable 2
 *  @param[in] scale a double scaling coefficient
 *  @param[in] time1 (not sure)
 *  @param[in] time2 (not sure)
 *  @param[in] n the number of elements in the array
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void elem_get_tposition(double *local1, double *local2, double *local3, double scale, double time1, double time2, int n) {
    int i=0;
    for (i=0; i<n; i++) {
        local1[i]=scale*(1.5*local2[i]-0.5*local3[i])+(1.-scale)*time1*local2[i]/time2;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates a scaled second-order temporal value at t + (3/2)*dt
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  \returns variable value at t+(3/2)*dt
 *
 *  @param[in] local2 a double array for variable 1
 *  @param[in] local3 a double array for variable 2
 *  @param[in] scale a double scaling coefficient
 *  @param[in] time1 (not sure)
 *  @param[in] time2 (not sure)
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double elem_get_tposition_scalar(double local2, double local3, double scale, double time1, double time2) {
    double var = 0.;
    var = scale*(1.5*local2-0.5*local3)+(1.-scale)*time1*local2/time2;
    return var;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates a scaled array of second-order 3D AdH vector temporal values at t + (3/2)*dt
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] local1 a vector array for output
 *  @param[in] local2 a vector array for variable 1
 *  @param[in] local3 a vector array for variable 2
 *  @param[in] scale a double scaling coefficient
 *  @param[in] time1 (not sure)
 *  @param[in] time2 (not sure)
 *  @param[in] n the number of elements in the vector array
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void elem_get_tposition_vect(SVECT *local1, SVECT *local2, SVECT *local3, double scale, double time1, double time2, int n) {
    int i=0;
    for (i=0; i<n; i++) {
        local1[i].x = scale*(1.5*local2[i].x - 0.5*local3[i].x) + (1.-scale)*time1*local2[i].x/time2;
        local1[i].y = scale*(1.5*local2[i].y - 0.5*local3[i].y) + (1.-scale)*time1*local2[i].y/time2;
        local1[i].z = scale*(1.5*local2[i].z - 0.5*local3[i].z) + (1.-scale)*time1*local2[i].z/time2;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates a second-order temporal value at t + (3/2)*dt
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  \returns variable value at t+(3/2)*dt
 *
 *  @param[in] local1 variable 1
 *  @param[in] local2 variable 2
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double get_second_order(double local1, double local2) {
    return (1.5*local1 - 0.5*local2);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates an array of second-order temporal values at t + (3/2)*dt
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  \returns variable value at t+(3/2)*dt
 *
 *  @param[in] local  the output array
 *  @param[in] local1 a double array of variable 1
 *  @param[in] local2 a double array variable 2
 *  @param[in] n the number of elements in the  array
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void get_second_order_array(double *local, double *local1, double *local2, int n) {
    int i=0;
    for (i=0; i<n; i++) {
        local[i] = get_second_order(local1[i], local2[i]);
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates an array of second-order 3D AdH vector temporal values at t + (3/2)*dt
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] local1 a vector array for output
 *  @param[in] local2 a vector array for variable 1
 *  @param[in] local3 a vector array for variable 2
 *  @param[in] n the number of elements in the vector array
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void elem_get_vect_second_order(SVECT *local1, SVECT *local2, SVECT *local3, int n) {
    int i=0;
    for (i=0; i<n; i++) {
        local1[i].x = get_second_order(local2[i].x, local3[i].x);
        local1[i].y = get_second_order(local2[i].y, local3[i].y);
        local1[i].z = get_second_order(local2[i].z, local3[i].z);
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates a 3 DOF first order finite difference of a Jacobi matrix elemental block
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] mat a 3 DOF matrix which stores the gradients
 *  @param[in] indx the block starting position
 *  @param[in] nnodes the total number of nodes on the element
 *  @param[in] local1 a 3 DOF residual with a (+) perturbation
 *  @param[in] local2 a 3 DOF residual with a (-) perturbation
 *  @param[in] diff_ep two times the perturbations size
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void elem_matrix_deriv(int indx, int nnodes, DOF_3 *local1, DOF_3 *local2, DOF_3 *mat, double diff_ep) {
    int inode;
    for (inode=0; inode<nnodes; inode++) {
        mat[ indx + inode*nnodes ].x_eq = (local1[inode].x_eq - local2[inode].x_eq)/diff_ep;
        mat[ indx + inode*nnodes ].y_eq = (local1[inode].y_eq - local2[inode].y_eq)/diff_ep;
        mat[ indx + inode*nnodes ].c_eq = (local1[inode].c_eq - local2[inode].c_eq)/diff_ep;
    }
}
inline void elem_matrix_deriv_4dof(int indx, int nnodes, DOF_4 *local1, DOF_4 *local2, DOF_4 *mat, double diff_ep) {
    int inode;
    for (inode=0; inode<nnodes; inode++) {
        mat[ indx + inode*nnodes ].x_eq = (local1[inode].x_eq - local2[inode].x_eq)/diff_ep;
        mat[ indx + inode*nnodes ].y_eq = (local1[inode].y_eq - local2[inode].y_eq)/diff_ep;
		mat[ indx + inode*nnodes ].z_eq = (local1[inode].z_eq - local2[inode].z_eq)/diff_ep;
        mat[ indx + inode*nnodes ].c_eq = (local1[inode].c_eq - local2[inode].c_eq)/diff_ep;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates a 3 DOF first order finite difference of a Jacobi matrix elemental block
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] mat a 3 DOF matrix which stores the gradients
 *  @param[in] indx the block starting position
 *  @param[in] nnodes the total number of nodes on the element
 *  @param[in] local1 a 3 DOF residual with a (+) perturbation
 *  @param[in] local2 a 3 DOF residual with a (-) perturbation
 *  @param[in] the inverse of diff_ep two times the perturbations size
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void elem_matrix_deriv2(int indx, int nnodes, DOF_3 *local1, DOF_3 *local2, DOF_3 *mat, double diff_ep_inv) {
    int inode;
    for (inode=0; inode<nnodes; inode++) {
        mat[ indx + inode*nnodes ].x_eq = (local1[inode].x_eq - local2[inode].x_eq) * diff_ep_inv;
        mat[ indx + inode*nnodes ].y_eq = (local1[inode].y_eq - local2[inode].y_eq) * diff_ep_inv;
        mat[ indx + inode*nnodes ].c_eq = (local1[inode].c_eq - local2[inode].c_eq) * diff_ep_inv;
    }
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates a 1 DOF first order finite difference of a Jacobi matrix elemental block
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] mat a 1 DOF matrix which stores the gradients
 *  @param[in] indx the block starting position
 *  @param[in] nnodes the total number of nodes on the element
 *  @param[in] local1 a 1 DOF residual with a (+) perturbation
 *  @param[in] local2 a 1 DOF residual with a (-) perturbation
 *  @param[in] diff_ep two times the perturbations size
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void elem_matrix_deriv_1dof(int indx, int nnodes, double *local1, double *local2, double *mat, double diff_ep) {
    int inode;
    for (inode=0; inode<nnodes; inode++) {
        mat[ indx + inode*nnodes ] = (local1[inode] - local2[inode])/diff_ep;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates grad(phi) times a 3D vector
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] grad_x a 3D vector array of values for grad(phi) * u
 *  @param[out] grad_y a 3D vector array of values for grad(phi) * v
 *  @param[out] grad_z a 3D vector array of values for grad(phi) * w
 *  @param[in] grad_phi an array of cartesian space shape function gradients on the element
 *  @param[in] v a 3D vector array
 *  @param[in] ndof the number of nodes on the element
 *
 * \f{eqnarray*}{
 *  \bbnabla u = \bbnabla(\sum\limits_{i}{\phidd{i} * u_i}) = \sum\limits_{i}{\bbnabla\phidd{i} * u_i}
 *  \bbnabla v = \bbnabla(\sum\limits_{i}{\phidd{i} * v_i}) = \sum\limits_{i}{\bbnabla\phidd{i} * v_i}
 *  \bbnabla w = \bbnabla(\sum\limits_{i}{\phidd{i} * w_i}) = \sum\limits_{i}{\bbnabla\phidd{i} * w_i}
 * \f}
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void grad_phi_dot_v(SVECT *grad_phi, SVECT *v, SVECT *grad_x, SVECT *grad_y, SVECT *grad_z, int ndof) {
    svect_init(grad_x);  svect_init(grad_y);  svect_init(grad_z);
    int idof=0;
    for (idof=0; idof<ndof; idof++) {
        grad_x->x += v[idof].x * grad_phi[idof].x;  grad_x->y += v[idof].x * grad_phi[idof].y;  grad_x->z += v[idof].x * grad_phi[idof].z;
        grad_y->x += v[idof].y * grad_phi[idof].x;  grad_y->y += v[idof].y * grad_phi[idof].y;  grad_y->z += v[idof].y * grad_phi[idof].z;
        grad_z->x += v[idof].z * grad_phi[idof].x;  grad_z->y += v[idof].z * grad_phi[idof].y;  grad_z->z += v[idof].z * grad_phi[idof].z;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates grad(phi) times a function
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] grad_x a 3D vector array of values for grad(phi) * f
 *  @param[in] grad_phi an array of cartesian space shape function gradients on the element
 *  @param[in] f a double array
 *  @param[in] ndof the number of nodes on the element
 *
 * \f{eqnarray*}{
 *  \bbnabla f = \bbnabla(\sum\limits_{i}{\phidd{i} * f_i}) = \sum\limits_{i}{\bbnabla\phidd{i} * f_i}
 * \f}
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void grad_phi_f(SVECT *grad_phi, double *f, SVECT *grad, int ndof) {
    svect_init(grad);
    int idof=0;
    for (idof=0; idof<ndof; idof++) {
        grad->x += f[idof] * grad_phi[idof].x;  grad->y += f[idof] * grad_phi[idof].y;  grad->z += f[idof] * grad_phi[idof].z;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates grad(phi) times a 2D vector
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] grad_x a 2D vector array of values for grad(phi) * u
 *  @param[out] grad_y a 2D vector array of values for grad(phi) * v
 *  @param[out] grad_z a 2D vector array of values for grad(phi) * w
 *  @param[in] grad_phi an array of cartesian space shape function gradients on the element
 *  @param[in] v a 2D vector array
 *  @param[in] ndof the number of nodes on the element
 *
 * \f{eqnarray*}{
 *  \bbnabla u = \bbnabla(\sum\limits_{i}{\phidd{i} * u_i}) = \sum\limits_{i}{\bbnabla\phidd{i} * u_i}
 *  \bbnabla v = \bbnabla(\sum\limits_{i}{\phidd{i} * v_i}) = \sum\limits_{i}{\bbnabla\phidd{i} * v_i}
 *  \bbnabla w = \bbnabla(\sum\limits_{i}{\phidd{i} * w_i}) = \sum\limits_{i}{\bbnabla\phidd{i} * w_i}
 * \f}
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void grad2d_phi_dot_v(SVECT2D *grad_phi, SVECT2D *v, SVECT2D *grad_x, SVECT2D *grad_y, int ndof) {
    svect2d_init(grad_x);  svect2d_init(grad_y);
    int idof=0;
    for (idof=0; idof<ndof; idof++) {
        grad_x->x += v[idof].x * grad_phi[idof].x;  grad_x->y += v[idof].x * grad_phi[idof].y;
        grad_y->x += v[idof].y * grad_phi[idof].x;  grad_y->y += v[idof].y * grad_phi[idof].y;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates 2D grad(phi) times a function
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] grad a 2D vector array of values for grad(phi) * f
 *  @param[in] grad_phi an array of cartesian space shape function gradients on the element
 *  @param[in] f a double array
 *  @param[in] ndof the number of nodes on the element
 *
 * \f{eqnarray*}{
 *  \bbnabla f = \bbnabla(\sum\limits_{i}{\phidd{i} * f_i}) = \sum\limits_{i}{\bbnabla\phidd{i} * f_i}
 * \f}
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void grad2d_phi_f(SVECT2D *grad_phi, double *f, SVECT2D *grad, int ndof) {
    svect2d_init(grad);
    int idof=0;
    for (idof=0; idof<ndof; idof++) {
        grad->x += f[idof] * grad_phi[idof].x;  grad->y += f[idof] * grad_phi[idof].y;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates a 3D tensor 3D vector product
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  \returns product
 *  @param[in] tens a 3D tensor structure
 *  @param[in] vect a 3D vector
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline SVECT tensor_vector_product(STENSOR tens, SVECT vect) {
    SVECT result; svect_init(&result);
    result.x = tens.xx*vect.x + tens.xy*vect.y + tens.xz*vect.z; 
    result.y = tens.xy*vect.x + tens.yy*vect.y + tens.yz*vect.z; 
    result.z = tens.xz*vect.x + tens.yz*vect.y + tens.zz*vect.z;
    return result;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates a 2D tensor 2D vector product
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  returns product
 *  @param[in] tens a 2D tensor structure
 *  @param[in] vect a 2D vector
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline SVECT2D tensor2d_vector2d_product(STENSOR2D tens, SVECT2D vect) {
    SVECT2D result; svect2d_init(&result);
    result.x = tens.xx*vect.x + tens.xy*vect.y;
    result.y = tens.xy*vect.x + tens.yy*vect.y;
    return result;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates a 2D anisotropic tensor 2D vector product
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  \returns product
 *  @param[in] tens a 2D anisotropic tensor structure
 *  @param[in] vect a 2D vector
 *  @param[in] vdir a 2D direction vector
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline SVECT2D tensor2d_vector2d_AI_product(STENSOR2D tens, SVECT2D vect, SVECT2D vdir) {
    SVECT2D result; svect2d_init(&result);
    result.x = tens.xx*vdir.x*(vdir.x*vect.x+vdir.y*vect.y) + tens.xy*vect.x; \
    result.y = tens.yy*vdir.y*(vdir.x*vect.x+vdir.y*vect.y) + tens.xy*vect.y;    
    return result;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     stores 3D local nodes in a vector given a global element "ie"
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] elem_nodes a 3D vector storing the local nodes
 *  @param[in] grid an AdH grid
 *  @param[in] ie the 3D global element ID
 *  @param[in] elem_dpl the local/elemental nodal displacements
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void get_elem3d_node_coordinates(SGRID *grid, int ie, double *elem_dpl, SVECT *elem_nodes) {
    int i;
    for (i=0; i<grid->elem3d[ie].nnodes; i++) {
        elem_nodes[i].x = grid->node[grid->elem3d[ie].nodes[i]].x;
        elem_nodes[i].y = grid->node[grid->elem3d[ie].nodes[i]].y;
        elem_nodes[i].z = grid->node[grid->elem3d[ie].nodes[i]].z + elem_dpl[i];  // factor in bed dpl later ...
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     stores 2D local nodes in a vector given a global element "ie"
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] elem_nodes a 3D vector storing the local nodes
 *  @param[in] grid an AdH grid
 *  @param[in] ie the 2D global element ID
 *  @param[in] elem_dpl the local/elemental nodal displacements
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void get_elem2d_node_coordinates(SGRID *grid, int ie, double *elem_dpl, SVECT *elem_nodes) {
    int i;
    for (i=0; i<grid->elem2d[ie].nnodes; i++) {
        elem_nodes[i].x = grid->node[grid->elem2d[ie].nodes[i]].x;
        elem_nodes[i].y = grid->node[grid->elem2d[ie].nodes[i]].y;
        elem_nodes[i].z = grid->node[grid->elem2d[ie].nodes[i]].z + elem_dpl[i];  // factor in bed dpl later ...
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     stores local double values from the global arrays
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] the local/elemental values
 *  @param[in] global the global values
 *  @param[in] nodes an array if integer global node IDs
 *  @param[in] nnodes the number of nodes on the element
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void global_to_local_dbl(double *global, double *local, int *nodes, int nnodes) {
    int i=0;
    for (i=0; i<nnodes; i++) {
        local[i] = global[nodes[i]];
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     stores local integer values from the global arrays
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] the local/elemental values
 *  @param[in] global the global values
 *  @param[in] nodes an array if integer global node IDs
 *  @param[in] nnodes the number of nodes on the element
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void global_to_local_int(int *global, int *local, int *nodes, int nnodes) {
    int i=0;
    for (i=0; i<nnodes; i++) {
        local[i] = global[nodes[i]];
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     stores local 2D vector values from the global arrays
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] the local/elemental values
 *  @param[in] global the global values
 *  @param[in] nodes an array if integer global node IDs
 *  @param[in] nnodes the number of nodes on the element
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void global_to_local_svect2d(SVECT2D *global, SVECT2D *local, int *nodes, int nnodes) {
    int i=0;
    for (i=0; i<nnodes; i++) {
        local[i] = global[nodes[i]];
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     stores local 3D vector values from the global arrays
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] the local/elemental values
 *  @param[in] global the global values
 *  @param[in] nodes an array if integer global node IDs
 *  @param[in] nnodes the number of nodes on the element
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void global_to_local_svect(SVECT *global, SVECT *local, int *nodes, int nnodes) {
    int i=0;
    for (i=0; i<nnodes; i++) {
        local[i] = global[nodes[i]];
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     returns the a 3D vector storing all elemental point locations, including midpoints
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] p the elemental node positions
 *  @param[in] nodes the vertex nodes
 *  @param[in] z the vertical position of each vertex node
 *  @param[in] nnodes the number of nodes on the element
 *  @param[in] nnodes_quad the number of vertex + midpoint nodes
 *  @param[in] nd_on_Edge a double point to the vertex nodes on each line segment of the element
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void elem_get_midpt_locations(SNODE *nodes, double *z, SVECT *p, int nnodes, int nnodes_quad, int **nd_on_Edge) {
    int i;
    for (i=0; i<nnodes; i++) {
        p[i].x = nodes[i].x;
        p[i].y = nodes[i].y;
        p[i].z = z[i];
    }
    
    int nd1, nd2;
    double xmid, ymid, zmid;
    for (i=nnodes; i<nnodes_quad; i++) {
        nd1 = nd_on_Edge[i-nnodes][0];
        nd2 = nd_on_Edge[i-nnodes][1];
        xmid = one_2 * ( nodes[nd1].x + nodes[nd2].x );
        ymid = one_2 * ( nodes[nd1].y + nodes[nd2].y );
        zmid = one_2 * ( z[nd1] + z[nd2] );
        p[i].x = xmid;
        p[i].y = ymid;
        p[i].z = zmid;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     returns the a 3D vector storing all elemental point locations, including midpoints
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] p the elemental node positions
 *  @param[in] nodes the vertex nodes
 *  @param[in] dpl the displacement of the node
 *  @param[in] nnodes the number of nodes on the element
 *  @param[in] nnodes_quad the number of vertex + midpoint nodes
 *  @param[in] nd_on_Edge a double point to the vertex nodes on each line segment of the element
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void elem_get_midpt_locations2(SNODE *nodes, double *dpl, SVECT *p, int nnodes, int nnodes_quad, int **nd_on_Edge) {
    int i;
    for (i=0; i<nnodes; i++) {
        p[i].x = nodes[i].x;
        p[i].y = nodes[i].y;
        p[i].z = nodes[i].z + dpl[i];
    }
    
    int nd1, nd2;
    double xmid, ymid, zmid;
    for (i=nnodes; i<nnodes_quad; i++) {
        nd1 = nd_on_Edge[i-nnodes][0];
        nd2 = nd_on_Edge[i-nnodes][1];
        xmid = one_2 * ( p[nd1].x + p[nd2].x );
        ymid = one_2 * ( p[nd1].y + p[nd2].y );
        zmid = one_2 * ( p[nd1].z + p[nd2].z );
        p[i].x = xmid;
        p[i].y = ymid;
        p[i].z = zmid;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     returns an array of relative velocities (Lagrange  - grid speed)
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] elem_rel_vel an elemental array of relative velocities
 *  @param[in] elem_vel the lagrangian velocities
 *  @param[in] elem_displacement the elemental vertical displacement at time = t
 *  @param[in] elem_old_displacement the elemental vertical displacement at time = t - dt
 *  @param[in] elem_older_displacement the elemental vertical displacement at time = t - 2*dt
 *  @param[in] dt the time-step
 *  @param[in] nnodes the number of nodes on the element
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void elem_get_relative_velocities(SVECT *elem_vel, double *elem_displacement, double *elem_old_displacement, double *elem_older_displacement, double dt, int nnodes, SVECT *elem_rel_vel) {
    double elem_dpl_32dt[nnodes]; sarray_init_dbl(elem_dpl_32dt, nnodes);  // displacement at t+(3/2)*dt
    get_second_order_array(elem_dpl_32dt, elem_displacement, elem_old_displacement, nnodes);
    
    double elem_dpl_12dt[nnodes]; sarray_init_dbl(elem_dpl_12dt, nnodes);  // displacement at t+(1/2)*dt
    get_second_order_array(elem_dpl_12dt, elem_old_displacement, elem_older_displacement, nnodes);
    
    int idof;
    double gs[nnodes]; sarray_init_dbl(gs,nnodes);
    for (idof=0; idof<nnodes; idof++) {
        gs[idof] = (elem_dpl_32dt[idof] - elem_dpl_12dt[idof])/dt;
    }
    
    for (idof=0; idof<nnodes; idof++) {
        elem_rel_vel[idof].x = elem_vel[idof].x;
        elem_rel_vel[idof].y = elem_vel[idof].y;
        elem_rel_vel[idof].z = (elem_vel[idof].z - gs[idof]);
    }
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

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     return two double arrays from a 3D AdH vector array
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] u a double array of x-velocities
 *  @param[out] v a double array of y-velocities
 *  @param[out] w a double array of y-velocities
 *  @param[in] vel the 3D vector array
 *  @param[in] nnodes the total number of elements in the 2D vector array
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void dumpVector(SVECT *vel, int nnodes, double *u, double *v, double *w) {
    int i;
    for (i=0; i<nnodes; i++) {
        u[i] = vel[i].x;
        v[i] = vel[i].y;
        w[i] = vel[i].z;
    }
}





