/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_ns3_load.c This file collections functions responsible for loading the 3D NS matrix */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

static int printFieldWidth = 12;
static int printPrecision  = 10;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// file prototypes
void perturb_ns3(SMODEL *mod, int nodes_on_elem, int ie, int dim, DOF_4 *elem_mat_u, DOF_4 *elem_mat_v, DOF_4 *elem_mat_w, DOF_4 *elem_mat_p, int DEBUG);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Forms the 3D Navier Stokes Newton Jacobi matrix.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \author    David Smith
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out]  mod a pointer to an AdH model where depth should be updated
 *
 *  \details Elemental matrix arrays - the "_u", for example, indicates the derivative with respect to u
 *   since these are DOF_3's there is one of these for each equation type.
 *   A 4-node element example the y-velocity derivative of the depth-averaged continuity matrix is: \n
 *  \f{eqnarray*}{ elem\_mat\_v[i].c\_eq =
 *    \begin{bmatrix}
 *    \deriv{R_{c}^0}{v^0} & \deriv{R_{c}^0}{v^1} & \deriv{R_{c}^0}{v^2} & \deriv{R_{c}^0}{v^3} \\
 *    \deriv{R_{c}^1}{v^0} & \deriv{R_{c}^1}{v^1} & \deriv{R_{c}^1}{v^2} & \deriv{R_{c}^1}{v^3} \\
 *    \deriv{R_{c}^2}{v^0} & \deriv{R_{c}^2}{v^1} & \deriv{R_{c}^2}{v^2} & \deriv{R_{c}^2}{v^3} \\
 *    \deriv{R_{c}^3}{v^0} & \deriv{R_{c}^3}{v^1} & \deriv{R_{c}^3}{v^2} & \deriv{R_{c}^3}{v^3}
 *    \end{bmatrix}
 *  \f}
 * \n where i=0 corresponds to \f$ \deriv{R_{c}^0}{v^0} \f$, i=1 to \f$ \deriv{R_{c}^0}{v^1} \f$, etc.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_ns3_load(SSUPER_MODEL *sm, int imod) {
 
    SMODEL *mod = &(sm->submodel[imod]);
    int i,j;
    
#ifdef _DEBUG
    assert(mod->flag.NS3_FLOW == ON);
    assert(mod->nsys == 4);
#endif
    
    int DEBUG = OFF;
    int DEBUG_MATRIX = OFF;
    int DEBUG_PICKETS = OFF;
#ifdef _DEBUG
    if (mod->t_prev > DEBUG_TIME) debug.load = ON;
    if (debug.load == ON) {
        DEBUG = ON;
    }
    if (debug.matrix == ON) {
        DEBUG_MATRIX = ON;
    }
    time_t time1;  time(&time1);
#endif
    
    // aliases
    SNS_3D *ns3 = mod->ns->d3;
    SGRID *grid = mod->grid;
    
    // local variables
    int ie2d, ie3d, nnodes, GnodeID;
    int nmatrix = NDONPRISM * NDONPRISM; // define using max 3D element size
    
    // elemental matrix stored as arrays
    DOF_4 elem_mat_p[nmatrix], elem_mat_u[nmatrix], elem_mat_v[nmatrix], elem_mat_w[nmatrix];
    
    // initialize matrix
    if (mod->amICoupled == NO) { // if this is a monolothic run, do this in the super load routine
        init_adh_matrix(grid->nnodes, mod->max_nsys_sq, sm->matrix, sm->diagonal);
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 3D ELEMENT CONTRIBUTIONS
    
    for (ie3d=0; ie3d<grid->nelems3d; ie3d++) {
        nnodes = grid->elem3d[ie3d].nnodes;
        
        dof4_init_array(elem_mat_p, nmatrix);
        dof4_init_array(elem_mat_u, nmatrix);
        dof4_init_array(elem_mat_v, nmatrix);
        dof4_init_array(elem_mat_w, nmatrix);
        
        perturb_ns3(mod, grid->elem3d[ie3d].nnodes, ie3d, 3, elem_mat_u, elem_mat_v, elem_mat_w, elem_mat_p, DEBUG);
        
        // matrix dirichlet bcs
        for(i=0; i<nnodes; i++) {
            GnodeID = grid->elem3d[ie3d].nodes[i];
            if(mod->str_values[grid->node[GnodeID].string].pressure.bc_flag == BCT_PRS_DIR ||
               mod->str_values[grid->node[GnodeID].string].pressure.bc_flag == BCT_HSP_DIR)
                fe_assign_prs_db_dof4(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_w, elem_mat_p);
            if(mod->str_values[grid->node[GnodeID].string].flow.bc_flag == BCT_VEL_DIR)
                fe_assign_mom_db_dof4(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_w, elem_mat_p);
            
            
            // force 0 presure on surface for now :: NOTE :: this is saying dpl = 0 too
            //if(grid->node[GnodeID].bflag == 0 && mod->flag.MG == ON)
            //if(grid->node[GnodeID].bflag == 0)
            //    fe_assign_prs_db_dof4(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_w, elem_mat_p);
        }
#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
        if (DEBUG) {
            printf("3d element: %d \n",ie3d);
            dof4_printScreen_array("elem_mat_p",elem_mat_p, nmatrix,__LINE__,__FILE__,printFieldWidth,printPrecision);
            dof4_printScreen_array("elem_mat_u",elem_mat_u, nmatrix,__LINE__,__FILE__,printFieldWidth,printPrecision);
            dof4_printScreen_array("elem_mat_v",elem_mat_v, nmatrix,__LINE__,__FILE__,printFieldWidth,printPrecision);
            dof4_printScreen_array("elem_mat_w",elem_mat_v, nmatrix,__LINE__,__FILE__,printFieldWidth,printPrecision);
        }
#endif
        
        // assemble global Jacobi Matrix
        fe_global_matrix_assemble_ns3(3, grid, mod->nsys, nnodes, ie3d, mod->fmap, elem_mat_u, elem_mat_v, elem_mat_w, elem_mat_p, sm->diagonal, sm->matrix);
    }
    
#ifdef _DEBUG
    if (DEBUG) {
        printScreen_matrix("FINAL ns3 matrix 3d element addtion", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
        //exit(-1);
    }
#endif
    //exit(-1);
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS
    
    for(ie2d = 0; ie2d < grid->nelems2d; ie2d++) {
        nnodes = grid->elem2d[ie2d].nnodes;
        
        dof4_init_array(elem_mat_p, nmatrix);
        dof4_init_array(elem_mat_u, nmatrix);
        dof4_init_array(elem_mat_v, nmatrix);
        dof4_init_array(elem_mat_w, nmatrix);
        
        perturb_ns3(mod, grid->elem2d[ie2d].nnodes, ie2d, 2, elem_mat_u, elem_mat_v, elem_mat_w, elem_mat_p, DEBUG);
        
        // matrix dirichlet bcs
        for(i = 0; i<grid->elem2d[ie2d].nnodes; i++) {
            GnodeID = grid->elem2d[ie2d].nodes[i];
            if(mod->str_values[grid->node[GnodeID].string].pressure.bc_flag == BCT_PRS_DIR ||
               mod->str_values[grid->node[GnodeID].string].pressure.bc_flag == BCT_HSP_DIR)
                fe_assign_prs_db_dof4(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_w, elem_mat_p);
            if(mod->str_values[grid->node[GnodeID].string].flow.bc_flag == BCT_VEL_DIR)
                fe_assign_mom_db_dof4(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_w, elem_mat_p);
            
            // force 0 presure on surface for now :: NOTE :: this is saying dpl = 0 too
            //if(grid->node[GnodeID].bflag == 0 && mod->flag.MG == ON)
            //if(grid->node[GnodeID].bflag == 0)
            //    fe_assign_prs_db_dof4(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_w, elem_mat_p);
        }
#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
        if (DEBUG) {
            printf("2d element: %d \n",ie2d);
            dof4_printScreen_array("elem_mat_p",elem_mat_p, nmatrix,__LINE__,__FILE__,printFieldWidth,printPrecision);
            dof4_printScreen_array("elem_mat_u",elem_mat_u, nmatrix,__LINE__,__FILE__,printFieldWidth,printPrecision);
            dof4_printScreen_array("elem_mat_v",elem_mat_v, nmatrix,__LINE__,__FILE__,printFieldWidth,printPrecision);
            dof4_printScreen_array("elem_mat_w",elem_mat_v, nmatrix,__LINE__,__FILE__,printFieldWidth,printPrecision);
        }
#endif
        
        fe_global_matrix_assemble_ns3(2, grid, mod->nsys, nnodes, ie2d, mod->fmap, elem_mat_u, elem_mat_v, elem_mat_w, elem_mat_p, sm->diagonal, sm->matrix);
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    //cjt :: not sure if this is needed
    //checks the diagonal for nonexistent entries
    if (mod->amICoupled == NO) { // if this is a monolothic run, do this in the super load routine
        double small_number = NOT_QUITE_SMALL;
        for (i = 0; i < grid->nnodes; i++) {
            j = i * mod->max_nsys_sq;
            if (fabs(sm->diagonal[j]) < small_number) {
                sm->diagonal[j] = small_number;
            }
            j = j + 4;
            if (fabs(sm->diagonal[j]) < small_number) {
                sm->diagonal[j] = small_number;
            }
            j = j + 4;
            if (fabs(sm->diagonal[j]) < small_number)
                sm->diagonal[j] = small_number;
        }
#ifdef _MESSG
        comm_update_double(sm->diagonal, mod->max_nsys_sq, mod->grid->smpi); // CJT :: questionable whether needed
#endif
        
#ifdef _DEBUG
        if (DEBUG_MATRIX == ON
#ifdef _MESSG
            && grid->smpi->myid == 0
#endif
            ) {
            printScreen_matrix("ns3 matrix", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
            exit(-1);
        }
#endif
    }
    
#ifdef _DEBUG
    if (DEBUG) {
        time_t time2;  time(&time2);
        TIME_IN_NS3_LOAD += difftime(time2,time1);
    }
#endif
    
    //exit(-1);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Adds perturbation to 3D velocity and pressure for each node in element for Newton Jacobian calculation.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] elem_mat_u stores the Jacobian, elemental matrix for the x-velocity perturbation
 *  @param[in,out] elem_mat_v stores the Jacobian, elemental matrix for the y-velocity perturbation
 *  @param[in,out] elem_mat_w stores the Jacobian, elemental matrix for the z-velocity perturbation
 *  @param[in,out] elem_mat_w stores the Jacobian, elemental matrix for the pressure perturbation
 *  @param[in] mod a pointer to an AdH model where depth should be updated
 *  @param[in] nodes_on_elem the number of nodes on the element
 *  @param[in] ie the element
 *  @param[in] dim the dimension of the resid
 *  @param[in] DEBUG a debug option
 *
 *  \note stuff
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void perturb_ns3(SMODEL *mod, int nodes_on_elem, int ie, int dim, DOF_4 *elem_mat_u, DOF_4 *elem_mat_v, DOF_4 *elem_mat_w, DOF_4 *elem_mat_p, int DEBUG) {
    
#ifdef _DEBUG
    assert(dim == 2 || dim == 3);
#endif
    
    int i, j, irow, GnodeID;
    double epsilon_u = 0., epsilon_v = 0., epsilon_2u = 0., epsilon_2v = 0., epsilon_d  = 0.;
    double epsilon_w = 0., epsilon_p = 0., epsilon_2w = 0., epsilon_2p = 0., epsilon_2d = 0.;
    double perturbation = sqrt(SMALL) * 0.01;
    
    DOF_4 elem_rhs_u_P[nodes_on_elem]; /* the residual resulting from a positive (P) perturbation of u */
    DOF_4 elem_rhs_u_M[nodes_on_elem]; /* the residual resulting from a negative (M) perturbation of u */
    DOF_4 elem_rhs_v_P[nodes_on_elem];
    DOF_4 elem_rhs_v_M[nodes_on_elem];
    DOF_4 elem_rhs_w_P[nodes_on_elem];
    DOF_4 elem_rhs_w_M[nodes_on_elem];
    DOF_4 elem_rhs_p_P[nodes_on_elem];
    DOF_4 elem_rhs_p_M[nodes_on_elem];
    
    // get maximum velocity over the element
    SVECT max_vel; svect_init_value(&max_vel,1.,1.,1.);
    double max_prs = 1., max_dpl = 1.;
    
    // cjt :: I'm not sure this should be done this way, since it changes the nodes perturbation value depending on what element it is solving
    for (i=0; i<nodes_on_elem; i++) {
        if (dim==3) {
            GnodeID = mod->grid->elem3d[ie].nodes[i];
        } else {
            GnodeID = mod->grid->elem2d[ie].nodes[i];
        }
        if (fabs(mod->ns->d3->vel[GnodeID].x) > max_vel.x) max_vel.x = fabs(mod->ns->d3->vel[GnodeID].x);
        if (fabs(mod->ns->d3->vel[GnodeID].y) > max_vel.y) max_vel.y = fabs(mod->ns->d3->vel[GnodeID].y);
        if (fabs(mod->ns->d3->vel[GnodeID].z) > max_vel.z) max_vel.z = fabs(mod->ns->d3->vel[GnodeID].z);
        if (fabs(mod->ns->d3->prs[GnodeID]) > max_prs) max_prs = fabs(mod->ns->d3->prs[GnodeID]);
        if (fabs(mod->ns->d3->displacement[GnodeID]) > max_dpl) max_dpl = fabs(mod->ns->d3->displacement[GnodeID]);
    }
    
    NUM_DIFF_EPSILON(epsilon_u, epsilon_2u, max_vel.x, perturbation);    // calculates epsilon and 2*epsilon
    NUM_DIFF_EPSILON(epsilon_v, epsilon_2v, max_vel.y, perturbation);    // calculates epsilon and 2*epsilon
    NUM_DIFF_EPSILON(epsilon_w, epsilon_2w, max_vel.z, perturbation);    // calculates epsilon and 2*epsilon
    NUM_DIFF_EPSILON(epsilon_p, epsilon_2p, max_prs, perturbation);      // calculates epsilon and 2*epsilon
    NUM_DIFF_EPSILON(epsilon_d, epsilon_2d, max_dpl, perturbation);      // calculates epsilon and 2*epsilon
    
    
    for (i=0; i<nodes_on_elem; i++) {
        if (dim == 2) { // BOUNDARY TERMS
            GnodeID = mod->grid->elem2d[ie].nodes[i];
            NUM_DIFF_EPSILON(epsilon_u, epsilon_2u, mod->ns->d3->vel[GnodeID].x, perturbation);    // calculates epsilon and 2*epsilon
            NUM_DIFF_EPSILON(epsilon_v, epsilon_2v, mod->ns->d3->vel[GnodeID].y, perturbation);    // calculates epsilon and 2*epsilon
            NUM_DIFF_EPSILON(epsilon_w, epsilon_2w, mod->ns->d3->vel[GnodeID].z, perturbation);    // calculates epsilon and 2*epsilon
            NUM_DIFF_EPSILON(epsilon_p, epsilon_2p, mod->ns->d3->prs[GnodeID], perturbation);      // calculates epsilon and 2*epsilon
            NUM_DIFF_EPSILON(epsilon_d, epsilon_2d, mod->ns->d3->displacement[GnodeID], perturbation);      // calculates epsilon and 2*epsilon
            
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) boundary perturbation of u-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +u for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_u);
#endif
            fe_ns3_boundary_resid(mod,elem_rhs_u_P,ie,epsilon_u,i,PERTURB_U,+1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) boundary perturbation of u-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -u for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_u);
#endif
            fe_ns3_boundary_resid(mod,elem_rhs_u_M,ie,epsilon_u,i,PERTURB_U,-1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) boundary perturbation of v-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +v for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_v);
#endif
            fe_ns3_boundary_resid(mod,elem_rhs_v_P,ie,epsilon_v,i,PERTURB_V,+1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) boundary perturbation of v-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -v for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_v);
#endif
            fe_ns3_boundary_resid(mod,elem_rhs_v_M,ie,epsilon_v,i,PERTURB_V,-1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) boundary perturbation of w-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +w for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_w);
#endif
            fe_ns3_boundary_resid(mod,elem_rhs_w_P,ie,epsilon_w,i,PERTURB_W,+1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) boundary perturbation of w-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -w for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_w);
#endif
            fe_ns3_boundary_resid(mod,elem_rhs_w_M,ie,epsilon_w,i,PERTURB_W,-1,DEBUG);
            
            if (mod->flag.MG == ON && mod->grid->node[GnodeID].bflag == 0) {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // (+) boundary perturbation of displacement +++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: pertubing +d for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_d);
#endif
                fe_ns3_boundary_resid(mod,elem_rhs_p_P,ie,epsilon_d,i,PERTURB_D,+1,DEBUG);
                
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // (-) boundary perturbation of pressure +++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: pertubing -d for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_d);
#endif
                fe_ns3_boundary_resid(mod,elem_rhs_p_M,ie,epsilon_d,i,PERTURB_D,-1,DEBUG);
            } else {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // (+) boundary perturbation of pressure +++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: pertubing +p for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_p);
#endif
                fe_ns3_boundary_resid(mod,elem_rhs_p_P,ie,epsilon_p,i,PERTURB_P,+1,DEBUG);
                
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // (-) boundary perturbation of pressure +++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: pertubing -p for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_p);
#endif
                fe_ns3_boundary_resid(mod,elem_rhs_p_M,ie,epsilon_p,i,PERTURB_P,-1,DEBUG);
            }
            
        } else if (dim == 3) { // BODY TERMS
            GnodeID = mod->grid->elem3d[ie].nodes[i];
            NUM_DIFF_EPSILON(epsilon_u, epsilon_2u, mod->ns->d3->vel[GnodeID].x, perturbation);    // calculates epsilon and 2*epsilon
            NUM_DIFF_EPSILON(epsilon_v, epsilon_2v, mod->ns->d3->vel[GnodeID].y, perturbation);    // calculates epsilon and 2*epsilon
            NUM_DIFF_EPSILON(epsilon_w, epsilon_2w, mod->ns->d3->vel[GnodeID].z, perturbation);    // calculates epsilon and 2*epsilon
            NUM_DIFF_EPSILON(epsilon_p, epsilon_2p, mod->ns->d3->prs[GnodeID], perturbation);      // calculates epsilon and 2*epsilon
            NUM_DIFF_EPSILON(epsilon_d, epsilon_2d, mod->ns->d3->displacement[GnodeID], perturbation);      // calculates epsilon and 2*epsilon
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) body perturbation of u-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +u for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_u);
#endif
            fe_ns3_body_resid(mod,elem_rhs_u_P,ie,epsilon_u,i,PERTURB_U,+1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) body perturbation of u-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -u for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_u);
#endif
            fe_ns3_body_resid(mod,elem_rhs_u_M,ie,epsilon_u,i,PERTURB_U,-1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) body perturbation of v-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +v for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_v);
#endif
            fe_ns3_body_resid(mod,elem_rhs_v_P,ie,epsilon_v,i,PERTURB_V,+1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) body perturbation of v-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -v for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_v);
#endif
            fe_ns3_body_resid(mod,elem_rhs_v_M,ie,epsilon_v,i,PERTURB_V,-1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) body perturbation of w-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +w for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_w);
#endif
            fe_ns3_body_resid(mod,elem_rhs_w_P,ie,epsilon_w,i,PERTURB_W,+1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) body perturbation of w-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -w for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_w);
#endif
            fe_ns3_body_resid(mod,elem_rhs_w_M,ie,epsilon_w,i,PERTURB_W,-1,DEBUG);
            
            if (mod->flag.MG == ON && mod->grid->node[GnodeID].bflag == 0) {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // (+) body perturbation of displacement +++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: pertubing +d for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_d);
#endif
                fe_ns3_body_resid(mod,elem_rhs_p_P,ie,epsilon_d,i,PERTURB_D,+1,DEBUG);
                
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // (-) body perturbation of displacement +++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: pertubing -d for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_d);
#endif
                fe_ns3_body_resid(mod,elem_rhs_p_M,ie,epsilon_d,i,PERTURB_D,-1,DEBUG);
            } else {
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // (+) body perturbation of pressure +++++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: pertubing +p for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_p);
#endif
                fe_ns3_body_resid(mod,elem_rhs_p_P,ie,epsilon_p,i,PERTURB_P,+1,DEBUG);
                
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // (-) body perturbation of pressure +++++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: pertubing -p for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_p);
#endif
                fe_ns3_body_resid(mod,elem_rhs_p_M,ie,epsilon_p,i,PERTURB_P,-1,DEBUG);
            }
            
        } else {
            
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> dimension (dim) must be either 2 or 3.");
            
        }
        
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // calculate Jacobian ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elem_matrix_deriv_4dof(i, nodes_on_elem, elem_rhs_u_P, elem_rhs_u_M, elem_mat_u, epsilon_2u);
        elem_matrix_deriv_4dof(i, nodes_on_elem, elem_rhs_v_P, elem_rhs_v_M, elem_mat_v, epsilon_2v);
        elem_matrix_deriv_4dof(i, nodes_on_elem, elem_rhs_w_P, elem_rhs_w_M, elem_mat_w, epsilon_2w);
        if (mod->flag.MG == ON && mod->grid->node[GnodeID].bflag == 0) {
            elem_matrix_deriv_4dof(i, nodes_on_elem, elem_rhs_p_P, elem_rhs_p_M, elem_mat_p, epsilon_2d);
        } else {
            elem_matrix_deriv_4dof(i, nodes_on_elem, elem_rhs_p_P, elem_rhs_p_M, elem_mat_p, epsilon_2p);
        }
        
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Vel BC residual replacement +++++++++++++++++++++++++++++++++++++++++++++++++
        // CJT :: I *THINK* THIS SHOULD JUST BE LATERAL
        ///*
        if (ROTATE == ON) {
            if (dim == 3) {
                for (j=0; j<nodes_on_elem; j++) {
                    GnodeID = mod->grid->elem3d[ie].nodes[j];
                    if (mod->bc_mask[GnodeID * mod->nsys] == 2) {
                        irow = j * nodes_on_elem + i;
                        elem_mat_u[irow].y_eq = mod->ns->d3->tanvec[GnodeID].x * elem_mat_u[irow].x_eq + mod->ns->d3->tanvec[GnodeID].y * elem_mat_u[irow].y_eq;
                        elem_mat_u[irow].x_eq = 0.;
                        elem_mat_v[irow].y_eq = mod->ns->d3->tanvec[GnodeID].x * elem_mat_v[irow].x_eq + mod->ns->d3->tanvec[GnodeID].y * elem_mat_v[irow].y_eq;
                        elem_mat_v[irow].x_eq = 0.;
                    }
                }
            }
            //exit(-1);
        }
        //*/
    }
}




