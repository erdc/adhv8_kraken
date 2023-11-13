/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Calculates the 3D Navier Stokes residuals.
 * \author    Gaurav Savant, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \author    David Smith
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] mod (SMODEL *) a pointer to an AdH model struct
 *
 * \details   Solves the following weak, discrete equation: \n
 * \f{eqnarray*}{
 *   \weakSWDaContReducedKinematic{i} \\
 *   \weakSWMxDDD{e}{i}{h} \\
 *   \weakSWMxDDD{e}{i}{h}
 * \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

static int printFieldWidth = 12;
static int printPrecision  = 10;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_ns3_resid(SSUPER_MODEL *sm, int imod) {
#ifndef _PETSC

    SMODEL *mod = &(sm->submodel[imod]);
#ifdef _DEBUG
    assert(mod->flag.NS3_FLOW == ON);
    assert(mod->nsys == 4);
#endif
    
    int DEBUG = OFF;
    int DEBUG_PICKETS = OFF;
#ifdef _DEBUG
    if (mod->t_prev > DEBUG_TIME) debug.residual = ON;
    if (debug.residual == ON) {
        DEBUG = ON;
    }
    time_t time1;  time(&time1);
#endif
    
    int i, j, k, ie2d, ie3d, GnodeID, id;
    double value = 0.;
    
    // alias
    SNS_3D *ns3 = mod->ns->d3;
    SGRID *grid = mod->grid;
    
    DOF_4 elem_rhs[NDONPRISM];
    
    /* initialize the residual */
    if (mod->amICoupled == NO) {
        sarray_init_dbl(sm->residual, grid->nnodes * mod->max_nsys);
    }
    
    /* initialize the supg storage arrays */
    //for (i=0; i<MAX_NNODES_ON_ELEM3D; i++) {
    //    sarray_init_dbl(ns3->elem_rhs_supg_dacont[i], grid->nelems3d);
    //}
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS
    
    for (ie2d=0; ie2d<grid->nelems2d; ie2d++) {
        
        // intialize and copy element variables
        dof4_init_array(elem_rhs, NDONPRISM);
        fe_ns3_boundary_resid(mod, elem_rhs, ie2d, 0, UNSET_INT, PERTURB_NONE, 0, DEBUG);
        
        /* sets the Dirichlet boundary conditions */
        for (i=0; i<grid->elem2d[ie2d].nnodes; i++) {
            GnodeID = grid->elem2d[ie2d].nodes[i];
            
            if(mod->str_values[grid->node[GnodeID].string].pressure.bc_flag == BCT_PRS_DIR ||
               mod->str_values[grid->node[GnodeID].string].pressure.bc_flag == BCT_HSP_DIR) {
                elem_rhs[i].c_eq = 0.;
            }
            
            // force 0 pressure at surface for now
            //if (grid->node[GnodeID].bflag == 0 && mod->flag.MG == ON)
            //if(grid->node[GnodeID].bflag == 0)
            //        elem_rhs[i].c_eq = 0.;
            
            if(mod->str_values[grid->node[GnodeID].string].flow.bc_flag == BCT_VEL_DIR) {
                elem_rhs[i].x_eq = 0.;
                elem_rhs[i].y_eq = 0.;
                elem_rhs[i].z_eq = 0.;
            }
        }
        
#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
        if (DEBUG) {
            printf("2d element: %d \n",ie2d);
            dof4_printScreen_array("elem_rhs",elem_rhs, grid->elem2d[ie2d].nnodes,__LINE__,__FILE__,printFieldWidth,printPrecision);
        }
#endif
        // Avengers Assemble
        for (i = 0; i < grid->elem2d[ie2d].nnodes; i++) {
            j = mod->nsys * mod->fmap[grid->elem2d[ie2d].nodes[i]];
            sm->residual[j + 0] -= elem_rhs[i].x_eq;
            sm->residual[j + 1] -= elem_rhs[i].y_eq;
            sm->residual[j + 2] -= elem_rhs[i].z_eq;
            sm->residual[j + 3] -= elem_rhs[i].c_eq;
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 3D ELEMENT CONTRIBUTIONS
    
    for (ie3d=0; ie3d<grid->nelems3d; ie3d++) {
        
        dof4_init_array(elem_rhs, NDONPRISM);
        fe_ns3_body_resid(mod, elem_rhs, ie3d, 0., UNSET_INT, PERTURB_NONE, 0, DEBUG);
        
        /* sets the Dirichlet boundary conditions */
        for (i=0; i<grid->elem3d[ie3d].nnodes; i++) {
            GnodeID = grid->elem3d[ie3d].nodes[i];
            
            if(mod->str_values[grid->node[GnodeID].string].pressure.bc_flag == BCT_PRS_DIR ||
               mod->str_values[grid->node[GnodeID].string].pressure.bc_flag == BCT_HSP_DIR) {
                elem_rhs[i].c_eq = 0.;
            }
            if(mod->str_values[grid->node[GnodeID].string].flow.bc_flag == BCT_VEL_DIR) {
                elem_rhs[i].x_eq = 0.;
                elem_rhs[i].y_eq = 0.;
                elem_rhs[i].z_eq = 0.;
            }
            
            if (ROTATE == ON) {
                j = mod->nsys * grid->elem3d[ie3d].nodes[i];
                if (mod->bc_mask[j] == 2) {
                    elem_rhs[i].y_eq = elem_rhs[i].x_eq * ns3->tanvec[grid->elem3d[ie3d].nodes[i]].x + elem_rhs[i].y_eq * ns3->tanvec[grid->elem3d[ie3d].nodes[i]].y;
                    elem_rhs[i].x_eq = 0.;
                }
            }
            
            // force 0 pressure at surface for now
            //if (grid->node[GnodeID].bflag == 0 && mod->flag.MG == ON)
            //if(grid->node[GnodeID].bflag == 0)
            //    elem_rhs[i].c_eq = 0.;
        }
        
#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
        if (DEBUG) {
            printf("3d element: %d \n",ie3d);
            dof4_printScreen_array("elem_rhs",elem_rhs, grid->elem3d[ie3d].nnodes,__LINE__,__FILE__,printFieldWidth,printPrecision);
        }
#endif
        // Avengers Assemble
        for (i = 0; i < grid->elem3d[ie3d].nnodes; i++) {
            j = mod->nsys * mod->fmap[grid->elem3d[ie3d].nodes[i]];
            sm->residual[j + 0] -= elem_rhs[i].x_eq;
            sm->residual[j + 1] -= elem_rhs[i].y_eq;
            sm->residual[j + 2] -= elem_rhs[i].z_eq;
            sm->residual[j + 3] -= elem_rhs[i].c_eq;
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    DEBUG = OFF;
    if (mod->amICoupled == NO) {
#ifdef _DEBUG
        if (DEBUG == ON
#ifdef _MESSG
            && grid->smpi->myid == 0
#endif
            ) {
            printf("\n");
            printf("printing residual: fe_ns3_resid @ line %d:%s \n",__LINE__, __FILE__);
            for (i = 0; i < mod->grid->nnodes; i++) {
                j = i * mod->nsys;
                printf("i=%d \t {%20.10f,%20.10f,%20.10f} \t final residual = ", (i + 1),mod->grid->node[i].x,mod->grid->node[i].y,mod->grid->node[i].z);
                for (k = 0; k < mod->nsys; k++)
                    printf(" %30.20e ", sm->residual[j + k]);
                printf("\n");
            }
        }
#endif
    }
    
#ifdef _DEBUG
    if (DEBUG == ON) {
        time_t time2;  time(&time2);
        TIME_IN_NS3_RESID += difftime(time2,time1);
    }
#endif

#endif
}
