/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_assign_elem_mat_db.c This file collects functions that assign dirichlet bcs
 *         to elemental matrices before global assembly.                                    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Adds Dirichlet velocity BCs to the Newton matrix.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in,out] (DOF_3 *) elem_mat_u  the u-velocity coefficient tensor
 * @param[in,out] (DOF_3 *) elem_mat_v  the v-velocity coefficient tensor
 * @param[in,out] (DOF_3 *) elem_mat_h  the h-velocity coefficient tensor
 * @param[in]     (int) inode the local node number of the head boundary
 * @param[in]     (int) nnodes the total number of nodes on the element
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_assign_mom_db_dof3(int inode, int nnodes, DOF_3 *elem_mat_u, DOF_3 *elem_mat_v, DOF_3 *elem_mat_h) {
    
    int j;			/* counter */
    int index;	    /* location in tensor */
    int jend;       /* end for loop */
    
    /* zero the row for momentum equations */
    index = nnodes * inode;
    for(j = index, jend = index + nnodes; j < jend; j++) {
        elem_mat_u[j].x_eq = 0.;
        elem_mat_u[j].y_eq = 0.;
        elem_mat_v[j].x_eq = 0.;
        elem_mat_v[j].y_eq = 0.;
        elem_mat_h[j].x_eq = 0.;
        elem_mat_h[j].y_eq = 0.;
    }
    
    /* zero the column associated with each velocity component */
    for(j = inode; j < nnodes * nnodes; j += nnodes) {
        elem_mat_u[j].x_eq = 0.;
        elem_mat_u[j].y_eq = 0.;
        elem_mat_u[j].c_eq = 0.;
        elem_mat_v[j].x_eq = 0.;
        elem_mat_v[j].y_eq = 0.;
        elem_mat_v[j].c_eq = 0.;
    }
    
    /* put 1. on the diagonal */
    elem_mat_u[inode * (nnodes+1)].x_eq = 1.;
    elem_mat_v[inode * (nnodes+1)].y_eq = 1.;
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Adds Dirichlet velocity BCs to the Newton matrix.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in,out] (DOF_3 *) elem_mat_u  the u-velocity coefficient tensor
 * @param[in,out] (DOF_3 *) elem_mat_v  the v-velocity coefficient tensor
 * @param[in,out] (DOF_3 *) elem_mat_h  the h-velocity coefficient tensor
 * @param[in]     (int) inode the local node number of the head boundary
 * @param[in]     (int) nnodes the total number of nodes on the element
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_assign_mom_db_dof4(int inode, int nnodes, DOF_4 *elem_mat_u, DOF_4 *elem_mat_v, DOF_4 *elem_mat_w, DOF_4 *elem_mat_p) {
    
    int j;			/* counter */
    int index;	    /* location in tensor */
    int jend;       /* end for loop */
    
    /* zero the row for momentum equations */
    index = nnodes * inode;
    for(j = index, jend = index + nnodes; j < jend; j++) {
        elem_mat_u[j].x_eq = 0.;
        elem_mat_u[j].y_eq = 0.;
        elem_mat_u[j].z_eq = 0.;
        elem_mat_v[j].x_eq = 0.;
        elem_mat_v[j].y_eq = 0.;
        elem_mat_v[j].z_eq = 0.;
        elem_mat_w[j].x_eq = 0.;
        elem_mat_w[j].y_eq = 0.;
        elem_mat_w[j].z_eq = 0.;
        elem_mat_p[j].x_eq = 0.;
        elem_mat_p[j].y_eq = 0.;
        elem_mat_p[j].z_eq = 0.;
    }
    
    /* zero the column associated with each velocity component */
    for(j = inode; j < nnodes * nnodes; j += nnodes) {
        elem_mat_u[j].x_eq = 0.;
        elem_mat_u[j].y_eq = 0.;
        elem_mat_u[j].z_eq = 0.;
        elem_mat_u[j].c_eq = 0.;
        elem_mat_v[j].x_eq = 0.;
        elem_mat_v[j].y_eq = 0.;
        elem_mat_v[j].z_eq = 0.;
        elem_mat_v[j].c_eq = 0.;
        elem_mat_w[j].x_eq = 0.;
        elem_mat_w[j].y_eq = 0.;
        elem_mat_w[j].z_eq = 0.;
        elem_mat_w[j].c_eq = 0.;
    }
    
    /* put 1. on the diagonal */
    elem_mat_u[inode * (nnodes+1)].x_eq = 1.;
    elem_mat_v[inode * (nnodes+1)].y_eq = 1.;
    elem_mat_w[inode * (nnodes+1)].z_eq = 1.;
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Adds Dirichlet velocity BCs to the Newton matrix.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in,out] (DOF_3 *) elem_mat_u  the u-velocity coefficient tensor
 * @param[in,out] (DOF_3 *) elem_mat_v  the v-velocity coefficient tensor
 * @param[in,out] (DOF_3 *) elem_mat_h  the h-velocity coefficient tensor
 * @param[in]     (int) inode the local node number of the head boundary
 * @param[in]     (int) nnodes the total number of nodes on the element
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_assign_prs_db_dof4(int inode, int nnodes, DOF_4 *elem_mat_u, DOF_4 *elem_mat_v, DOF_4 *elem_mat_w, DOF_4 *elem_mat_p) {
    int j;            /* counter */
    int index;        /* location in tensor */
    int jend;         /* end for loop */
    
    /* zero the row for continuity equation */
    index = nnodes * inode;
    for(j = index, jend = index + nnodes; j < jend; j++) {
        elem_mat_u[j].c_eq = 0.;
        elem_mat_v[j].c_eq = 0.;
        elem_mat_w[j].c_eq = 0.;
        elem_mat_p[j].c_eq = 0.;
    }
    
    /* zero the column associated with the pressure */
    for(j = inode; j < nnodes * nnodes; j += nnodes) {
        elem_mat_p[j].x_eq = 0.;
        elem_mat_p[j].y_eq = 0.;
        elem_mat_p[j].z_eq = 0.;
        elem_mat_p[j].c_eq = 0.;
    }
    
    /* put 1. on the diagonal */
    elem_mat_p[inode * (nnodes+1)].c_eq = 1.;
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Adds Dirichlet velocity BCs to the Newton matrix.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in,out] (DOF_3 *) elem_mat_u  the u-velocity coefficient tensor
 * @param[in,out] (DOF_3 *) elem_mat_v  the v-velocity coefficient tensor
 * @param[in,out] (DOF_3 *) elem_mat_h  the h-velocity coefficient tensor
 * @param[in]     (int) inode the local node number of the head boundary
 * @param[in]     (int) nnodes the total number of nodes on the element
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_assign_dacont_db_dof3(int inode, int nnodes, DOF_3 *elem_mat_u, DOF_3 *elem_mat_v, DOF_3 *elem_mat_h) {
    int j;            /* counter */
    int index;        /* location in tensor */
    int jend;         /* end for loop */
    
    /* zero the row for continuity equation */
    index = nnodes * inode;
    for(j = index, jend = index + nnodes; j < jend; j++) {
        elem_mat_u[j].c_eq = 0.;
        elem_mat_v[j].c_eq = 0.;
        elem_mat_h[j].c_eq = 0.;
    }
    
    /* zero the column associated with the pressure */
    for(j = inode; j < nnodes * nnodes; j += nnodes) {
        elem_mat_h[j].x_eq = 0.;
        elem_mat_h[j].y_eq = 0.;
        elem_mat_h[j].c_eq = 0.;
    }
    
    /* put 1. on the diagonal */
    elem_mat_h[inode * (nnodes+1)].c_eq = 1.;
    
}
