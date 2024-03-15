/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_sw2_2d_elem_resid.c This file collections functions responsible for
 *          the 2D shallow water body contributions to an element in the Newton residual.              */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

// prototypes
void fe_sw2_temporal(int ie, SELEM_2D *elem2d, int nnodes, SVECT *elem_nds, double djac, double drying_lower_limit, double *elem_head, SVECT2D *elem_vel, double wd_factor, double dt_factor, DOF_3 *elem_rhs, char *string, int DEBUG, int DEBUG_LOCAL, int wd_flag);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 2D shallow water elemental residual.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] elem_rhs      the 2D elemental residual array
 * @param[in]  mod           a pointer to the model struct
 * @param[in]  ie            the elemental id
 * @param[in]  pertubation   the Newton pertubation
 * @param[in]  perturb_node  the node to be perturbed
 * @param[in]  perturb_var   the variable to be perturbed
 * @param[in]  perturb_sign  the direction of Newton perturbation
 * @param[in]  DEBUG         a debug flag
 * @param[in]  PRESSURE_FLAG  a flag turning pressure contributions on/off
 *
 *  \details Solves the body integals of the following weak, discrete body terms of the 2D shallow water equation: \n
 *  \f{eqnarray*}{
 *  \weakSWDAcont{e}{i}{h} \\
 *  \weakSWMxDD{e}{i}{h} \\
 *  \weakSWMyDD{e}{i}{h}
 *  \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
