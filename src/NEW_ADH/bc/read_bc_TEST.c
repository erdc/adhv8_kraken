/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_TEST.c  This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file verifcation cases
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] *sm (SUPER_MODEL *)  a double pointer to an array of AdH supermodels
 * @param[in] **token (CHAR) a BC file line string token
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_bc_TEST(SMODEL_SUPER *mod, char **token) {
    int ng, n1, n2;
    
    if (strcmp(*token, "NB") == 0) {
        test_case_flag.on = ON;
        test_case_flag.nb = ON;
        
    } else if (strcmp(*token, "WD") == 0) {
        test_case_flag.on = ON;
        test_case_flag.wet_dry = ON;
        
    } else if (strcmp(*token, "DAM") == 0) {
        test_case_flag.on = ON;
        test_case_flag.dam2d = ON;
        test_case_flag.dam_location = get_next_token_dbl(token);
        test_case_flag.hl = get_next_token_dbl(token);
        test_case_flag.hr = get_next_token_dbl(token);
        
    } else if (strcmp(*token, "COR") == 0) {
        test_case_flag.node1 = get_next_token_int(token); - 1;
        test_case_flag.on = ON;
        test_case_flag.coriolis = ON;
        
    } else if (strcmp(*token, "XCOR") == 0) {
        test_case_flag.node1 = get_next_token_int(token); - 1;
        test_case_flag.node2 = get_next_token_int(token); - 1;
        printf("segment endpt 1: node[%d] :: {%20.10e %20.10e %20.10e}\n",
               test_case_flag.node1,
               mod->grid->node[test_case_flag.node1].x,
               mod->grid->node[test_case_flag.node1].y,
               mod->grid->node[test_case_flag.node1].z);
        printf("segment endpt 2: node[%d] :: {%20.10e %20.10e %20.10e}\n\n",
               test_case_flag.node2,
               mod->grid->node[test_case_flag.node2].x,
               mod->grid->node[test_case_flag.node2].y,
               mod->grid->node[test_case_flag.node2].z);
        test_case_flag.on = ON;
        test_case_flag.xcoriolis = ON;
        
    } else if (strcmp(*token, "YCOR") == 0) {
        test_case_flag.node1 = get_next_token_int(token); - 1;
        test_case_flag.node2 = get_next_token_int(token); - 1;
        printf("** analytic solution line segment found with endpoints: \n");
        printf("pt 1: node[%d] :: {%20.10e %20.10e %20.10e}\n",
               test_case_flag.node1,
               mod->grid->node[test_case_flag.node1].x,
               mod->grid->node[test_case_flag.node1].y,
               mod->grid->node[test_case_flag.node1].z);
        printf("pt 2: node[%d] :: {%20.10e %20.10e %20.10e}\n\n",
               test_case_flag.node2,
               mod->grid->node[test_case_flag.node2].x,
               mod->grid->node[test_case_flag.node2].y,
               mod->grid->node[test_case_flag.node2].z);
        
        
        test_case_flag.on = ON;
        test_case_flag.ycoriolis = ON;
        
    } else if (strcmp(*token, "PTEST") == 0) {
        test_case_flag.on = ON;
        test_case_flag.ptest = ON;
        
    } else if (strcmp(*token, "DIS") == 0) {
        test_case_flag.on = ON;
        test_case_flag.discharge = ON;
        
    } else if (strcmp(*token, "SLOSH") == 0) {
        test_case_flag.nodeID = get_next_token_int(token); - 1;
#ifndef _MESSG
        if (test_case_flag.nodeID < 0 || test_case_flag.nodeID > mod->grid->nnodes-1) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> The node ID for slosh model error calculation is not on grid.\n");
        }
#endif
        test_case_flag.a = get_next_token_dbl(token);;
        if (test_case_flag.a<0) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> the slosh disturbance amplitube cannot be negative.\n");
        }
        test_case_flag.L = get_next_token_dbl(token);;
        if (test_case_flag.L<0) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> the slosh grid length cannot be negative.\n");
        }
        test_case_flag.B = get_next_token_dbl(token);;
        if (test_case_flag.B<0) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> the slosh grid width cannot be negative.\n");
        }
        test_case_flag.H = get_next_token_dbl(token);;
        if (test_case_flag.H<0) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> the slosh water depth cannot be negative.\n");
        }
        if (test_case_flag.H/test_case_flag.L > 0.05) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> Your test case parameters violated hydrostatic approximation (H/L > 0.05).\n");
        }
        
        test_case_flag.on = ON;
        test_case_flag.slosh = ON;
        
    } else if (strcmp(*token, "SLSH23") == 0) {
        // printf("\n%s", subdata);
        test_case_flag.nodeID = get_next_token_int(token); - 1;
        // printf("\nnodeID = %i",test_case_flag.nodeID);
#ifndef _MESSG
        if (test_case_flag.nodeID < 0 || test_case_flag.nodeID > mod->grid->nnodes-1) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> The node ID for slosh model error calculation is not on grid.\n");
        }
#endif
        test_case_flag.a = get_next_token_dbl(token);;
        // printf("\na = %10.5e",test_case_flag.a);
        if (test_case_flag.a<0) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> the slosh disturbance amplitube cannot be negative.\n");
        }
        test_case_flag.L = get_next_token_dbl(token);;
        // printf("\nL = %10.5e",test_case_flag.L);
        if (test_case_flag.L<0) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> the slosh grid length cannot be negative.\n");
        }
        test_case_flag.B = get_next_token_dbl(token);;
        // printf("\nB = %10.5e",test_case_flag.B);
        if (test_case_flag.B<0) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> the slosh grid width cannot be negative.\n");
        }
        test_case_flag.H = get_next_token_dbl(token);;
        // printf("\nH = %10.5e",test_case_flag.H);
        if (test_case_flag.H<0) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> the slosh water depth cannot be negative.\n");
        }
        if (test_case_flag.H/test_case_flag.L > 0.05) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> Your test case parameters violated hydrostatic approximation (H/L > 0.05).\n");
        }
        
        test_case_flag.on = ON;
        test_case_flag.slosh2d3d = ON;
        
    } else if (strcmp(*token, "WNDH") == 0) {
        
        if (mod->grid->ndim==3) {
            // node ID
            test_case_flag.nodeID = get_next_token_int(token); - 1;
            test_case_flag.winds_Huang = ON;
#ifndef _MESSG
            if (test_case_flag.nodeID < 0 || test_case_flag.nodeID > mod->grid->nnodes-1) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> The node ID for wind model error calculation is not on grid.\n");
            }
            // make sure this is a surface node and get the surface node ID
            int check = NO;
            int nd = -1, inode = -1;
            ID_LIST_ITEM *ptr;
            for (inode=0; inode<mod->grid->nnodes_sur; inode++) {
                ptr = mod->grid->vertical_list[inode];
                nd = ptr->id;
                if (nd == test_case_flag.nodeID) {
                    test_case_flag.nodeID = inode;  // this is now a surface node in the surnode list
                    check = YES;
                    break;
                }
            }
            if (check == NO) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> The node ID for wind model error calculation is not a surface node.\n");
            }
#endif
        } else {
            test_case_flag.winds_and_waves = ON;
            test_case_flag.node1 = get_next_token_int(token); - 1;
            test_case_flag.node2 = get_next_token_int(token); - 1;
            printf("segment endpt 1: node[%d] :: {%20.10e %20.10e %20.10e}\n",
                   test_case_flag.node1,
                   mod->grid->node[test_case_flag.node1].x,
                   mod->grid->node[test_case_flag.node1].y,
                   mod->grid->node[test_case_flag.node1].z);
            printf("segment endpt 2: node[%d] :: {%20.10e %20.10e %20.10e}\n\n",
                   test_case_flag.node2,
                   mod->grid->node[test_case_flag.node2].x,
                   mod->grid->node[test_case_flag.node2].y,
                   mod->grid->node[test_case_flag.node2].z);
        }
        test_case_flag.on = ON;
        
    } else if (strcmp(*token, "SOURCE") == 0) {
        test_case_flag.on = ON;
        test_case_flag.water_source = ON;
        
    } else if (strcmp(*token, "LOCK") == 0) {
        test_case_flag.on = ON;
        test_case_flag.lock = ON;
        
        // barotropic transport
    } else if (strcmp(*token, "FRONT") == 0) {
        printf("Steep Front 3D Transport Test Active\n");
        if (debug.no_hydro == OFF){
            printf("NOTERM HYDRO card must be specified for this test case\n");
            exit(-1);}
        test_case_flag.on = ON;
        test_case_flag.steep_front = ON;
        test_case_flag.init_conc = get_next_token_dbl(token);;
        test_case_flag.nodeID = get_next_token_int(token); - 1;
        
        // baroclinic transport
    } else if (strcmp(*token, "SALT") == 0) {
        test_case_flag.on = ON;
        test_case_flag.salt = ON;
        
        
    } else if (strcmp(*token, "CSLOSH") == 0) {
        test_case_flag.nodeID = get_next_token_int(token); - 1;
        if (test_case_flag.nodeID < 0 || test_case_flag.nodeID > mod->grid->nnodes-1) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> The node ID for circular slosh model error calculation is not on grid.\n");
        }
        test_case_flag.a = get_next_token_dbl(token);;
        if (test_case_flag.a<0) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> the circular slosh disturbance amplitude cannot be negative.\n");
        }
        test_case_flag.R = get_next_token_dbl(token);;
        if (test_case_flag.R<0) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> the circular slosh grid radius cannot be negative.\n");
        }
        test_case_flag.H = get_next_token_dbl(token);;
        if (test_case_flag.H<0) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> the circular slosh water depth cannot be negative.\n");
        }
        if (test_case_flag.H/test_case_flag.R > 0.05) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> Your test case parameters violated hydrostatic approximation (H/L > 0.05).\n");
        }
        
        test_case_flag.on = ON;
        test_case_flag.cslosh = ON;
        
    } else if (strcmp(*token, "OF") == 0) {
        test_case_flag.on = ON;
        test_case_flag.outflow = ON;
        
    } else if (strcmp(*token, "TIDE") == 0) {
        ng = get_next_token_int(token); - 1;
        n1 = get_next_token_int(token); - 1;
        n2 = get_next_token_int(token); - 1;
        
        test_case_flag.on = ON;
        
        if (mod->grid->ndim == 2) {
            test_case_flag.tide2d = ON;
            test_case_flag.nodeID_2d = ng;
            test_case_flag.node1_2d = n1;
            test_case_flag.node2_2d = n2;
        }
        if (mod->grid->ndim == 3) {
            test_case_flag.tide3d = ON;
            test_case_flag.nodeID_3d = ng;
            test_case_flag.node1_3d = n1;
            test_case_flag.node2_3d = n2;
        }
        
        printf("** analytic solution line segment found with endpoints: \n");
        printf("pt 1: node[%d] :: {%20.10e %20.10e %20.10e}\n",n1,
               mod->grid->node[n1].x,
               mod->grid->node[n1].y,
               mod->grid->node[n1].z);
        printf("pt 2: node[%d] :: {%20.10e %20.10e %20.10e}\n\n",n2,
               mod->grid->node[n2].x,
               mod->grid->node[n2].y,
               mod->grid->node[n2].z);
        
        // diffusive wave - Hunter test case
    } else if (strcmp(*token, "DWEH") == 0) {
        test_case_flag.on = ON;
        test_case_flag.dwe_hunter = ON;
    }
}
