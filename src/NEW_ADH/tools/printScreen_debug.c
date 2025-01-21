#include "adh.h"
static int printFieldWidth = 30;
static int printPrecision  = 20;
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
void printScreen_rhs_3dof(char *string, int nnodes, int ie, int *nodes, double *elem_rhs) {
    int i;
    printf("\n***************************************\n%s\n",string);
    for (i=0; i<nnodes; i++) {
        printf("RHS contrib from %s[%d] ie: %d node: %-10d :: ",string,i+1,ie+1,nodes[i]+1);
        printf("c_eq: %*.*e \t x_eq: %*.*e \t y_eq: %*.*e\n",
            printFieldWidth,printPrecision,elem_rhs[i*3],
            printFieldWidth,printPrecision,elem_rhs[i*3+1],
            printFieldWidth,printPrecision,elem_rhs[i*3+2]);
    
    }
    //printf("**************************************\n");
}
