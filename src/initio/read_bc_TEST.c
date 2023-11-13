#include "global_header.h"

// these case will stop the code

void read_bc_TEST(SMODEL * mod, char *data) {

    char line[MAXLINE];         /* the input line */
    char *subdata = NULL;       /* the data after the second card is read */
    char *subsubdata = NULL;
    int ng, n1, n2;
    
    SIO info = *(mod->io);

    //printf("data: %s :: parsed: %d\n",data,parse_card(data,&subdata));

    switch (parse_card(data, &subdata)) {
        
        case CARD_NB:
            test_case_flag.on = ON;
            test_case_flag.nb = ON;
            break;

        case CARD_WD:
            test_case_flag.on = ON;
            test_case_flag.wet_dry = ON;
            break;

        case CARD_DAM:
            test_case_flag.on = ON;
            test_case_flag.dam2d = ON;
            test_case_flag.dam_location = read_dbl_field(info, &subdata);
            test_case_flag.hl = read_dbl_field(info, &subdata);
            test_case_flag.hr = read_dbl_field(info, &subdata);
            break;
            
        case CARD_COR:
            test_case_flag.node1 = read_int_field(info, &subdata) - 1;
            test_case_flag.on = ON;
            test_case_flag.coriolis = ON;
            break;

        case CARD_XCOR:
            test_case_flag.node1 = read_int_field(info, &subdata) - 1;
            test_case_flag.node2 = read_int_field(info, &subdata) - 1;
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
            break;

        case CARD_YCOR:
            test_case_flag.node1 = read_int_field(info, &subdata) - 1;
            test_case_flag.node2 = read_int_field(info, &subdata) - 1;
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
            break;

        case CARD_PTEST:
            test_case_flag.on = ON;
            test_case_flag.ptest = ON;
            break;
    
        case CARD_DIS:
            test_case_flag.on = ON;
            test_case_flag.discharge = ON;
            break;

        case CARD_SLOSH:
            test_case_flag.nodeID = read_int_field(info, &subdata) - 1;
#ifndef _MESSG
            if (test_case_flag.nodeID < 0 || test_case_flag.nodeID > mod->grid->nnodes-1) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> The node ID for slosh model error calculation is not on grid.\n");
            }
#endif
            test_case_flag.a = read_dbl_field(info, &subdata);
            if (test_case_flag.a<0) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> the slosh disturbance amplitube cannot be negative.\n");
            }
            test_case_flag.L = read_dbl_field(info, &subdata);
            if (test_case_flag.L<0) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> the slosh grid length cannot be negative.\n");
            }
            test_case_flag.B = read_dbl_field(info, &subdata);
            if (test_case_flag.B<0) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> the slosh grid width cannot be negative.\n");
            }
            test_case_flag.H = read_dbl_field(info, &subdata);
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
            break;

        case CARD_SLSH23:
            // printf("\n%s", subdata);
            test_case_flag.nodeID = read_int_field(info, &subdata) - 1;
            // printf("\nnodeID = %i",test_case_flag.nodeID);
#ifndef _MESSG
            if (test_case_flag.nodeID < 0 || test_case_flag.nodeID > mod->grid->nnodes-1) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> The node ID for slosh model error calculation is not on grid.\n");
            }
#endif
            test_case_flag.a = read_dbl_field(info, &subdata);
            // printf("\na = %10.5e",test_case_flag.a);
            if (test_case_flag.a<0) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> the slosh disturbance amplitube cannot be negative.\n");
            }
            test_case_flag.L = read_dbl_field(info, &subdata);
            // printf("\nL = %10.5e",test_case_flag.L);
            if (test_case_flag.L<0) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> the slosh grid length cannot be negative.\n");
            }
            test_case_flag.B = read_dbl_field(info, &subdata);
            // printf("\nB = %10.5e",test_case_flag.B);
            if (test_case_flag.B<0) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> the slosh grid width cannot be negative.\n");
            }
            test_case_flag.H = read_dbl_field(info, &subdata);
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
            break;

        case CARD_WNDH:
            
            if (mod->grid->ndim==3) {
                // node ID
                test_case_flag.nodeID = read_int_field(info, &subdata) - 1;
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
                test_case_flag.node1 = read_int_field(info, &subdata) - 1;
                test_case_flag.node2 = read_int_field(info, &subdata) - 1;
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
            break;

        case CARD_SOURCE:
            test_case_flag.on = ON;
            test_case_flag.water_source = ON;
            break;

        case CARD_LOCK:
            test_case_flag.on = ON;
            test_case_flag.lock = ON;
            break;

		// barotropic transport
		case CARD_FRONT:
			printf("Steep Front 3D Transport Test Active\n");
			if (debug.no_hydro == OFF){
				printf("NOTERM HYDRO card must be specified for this test case\n");
				exit(-1);}
			test_case_flag.on = ON;
			test_case_flag.steep_front = ON;
			test_case_flag.init_conc = read_dbl_field(info, &subdata);
			test_case_flag.nodeID = read_int_field(info, &subdata) - 1;
			break;

        // baroclinic transport
        case CARD_SALT:
            test_case_flag.on = ON;
            test_case_flag.salt = ON;
            break;


        case CARD_CSLOSH:
            test_case_flag.nodeID = read_int_field(info, &subdata) - 1;
            if (test_case_flag.nodeID < 0 || test_case_flag.nodeID > mod->grid->nnodes-1) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> The node ID for circular slosh model error calculation is not on grid.\n");
            }
            test_case_flag.a = read_dbl_field(info, &subdata);
            if (test_case_flag.a<0) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> the circular slosh disturbance amplitude cannot be negative.\n");
            }
            test_case_flag.R = read_dbl_field(info, &subdata);
            if (test_case_flag.R<0) {
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> the circular slosh grid radius cannot be negative.\n");
            }
            test_case_flag.H = read_dbl_field(info, &subdata);
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
            break;

        case CARD_OF:
            test_case_flag.on = ON;
            test_case_flag.outflow = ON;
            break;

        case CARD_TIDE:
            ng = read_int_field(info, &subdata) - 1;
            n1 = read_int_field(info, &subdata) - 1;
            n2 = read_int_field(info, &subdata) - 1;
            
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
            break;
            
        // diffusive wave - Hunter test case
        case CARD_DWEH:
            test_case_flag.on = ON;
            test_case_flag.dwe_hunter = ON;
            break;

        default:
            break;
    }
}
