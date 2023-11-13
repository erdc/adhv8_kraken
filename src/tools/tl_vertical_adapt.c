/* Looks like sigma stretching */
#include "global_header.h"

static int DEBUG = OFF;

//***********************************************************************//
//***********************************************************************//
//***********************************************************************//
// distributes the surface displacement perturbation down the column, and scales it
// -- (cjt) note :: distributing it proportionally down the column as in original code is less stable, so
// -- (cjt) note :: I have determined that setting the bottom node to 0 perturbation causes the convergence issues
void tl_get_dpl_perturbation(SGRID *grid, double perturbation, double *dpl, double *dpl_perturbation) {
    
    ID_LIST_ITEM *ptr = NULL;
    int nd=-1, iseg=-1, top_node=-1, bot_node=-1;
    double top_z = 0., bot_z = 0., elev_factor = 0., dum = 0., scalar = 1., surf_dpl = 0.;
    
    /* Loop over the surface nodes */
    for (iseg=0; iseg<grid->nnodes_sur; iseg++) {
        
        ptr = grid->vertical_list[iseg];
        top_node = ptr->id;
        bot_node = (ptr->prev)->id;
        top_z = grid->node[top_node].z;
        bot_z = grid->node[bot_node].z;
        surf_dpl = dpl[top_node]; // scale the perturbations by the surface displacement
        
        // distribute the perturbation fractionally down the column
        while(ptr->next != NULL) {
            nd = ptr->id;
            dpl_perturbation[nd] = 0.;
            elev_factor = 1;
            
            // scale the perturbation by surface displacement
            NUM_DIFF_EPSILON_GENERAL(dpl_perturbation[nd], dum, fabs(surf_dpl), perturbation); 
            
            // over-ride scaling for now, since for angle source is gives bunk results
            dpl_perturbation[nd] = perturbation;
            
            
            elev_factor = (grid->node[nd].z - bot_z)/(top_z - bot_z);
            dpl_perturbation[nd] *= elev_factor;
            //printf("dpl_perturbation[%d]: %20.10e\n",nd,dpl_perturbation[nd]);
            
            // fix bottom peturbation to 0 (already done if X elev_factor)
            //if (fabs(grid->node[nd].z - bot_z) < 1e-6) dpl_perturbation[nd] = 0.; // code fails to converge if this isn't here
            
            // fix all subsurfouce nodes
            //if (fabs(grid->node[nd].z - top_z) > 1e-6) dpl_perturbation[nd] = 0.;
            ptr=ptr->next;
        }
    }

    //exit(-1);
}

//***********************************************************************//
//***********************************************************************//
//***********************************************************************//
// adapts non-surface nodes in the column via their displacements
// accordian style (like old trunk)
void tl_vertical_adapt(SGRID *grid, double *dpl) {
    
    ID_LIST_ITEM *ptr = NULL;
    int nd=-1, iseg=-1, top_node=-1, bot_node=-1;
    double top_z = 0., bot_z = 0., elev_factor = 0.;
    double surf_dpl = 0.;
    
    /* Loop over the surface nodes */
    for (iseg=0; iseg<grid->nnodes_sur; iseg++) {
        ptr = grid->vertical_list[iseg];
        
        top_node = ptr->id;
        bot_node = (ptr->prev)->id;
        top_z = grid->node[top_node].z;
        bot_z = grid->node[bot_node].z;
        
        surf_dpl = dpl[top_node];
        //ptr=ptr->next; // leave the surface node alone (it's lagrangian)
        
        while(ptr->next != NULL) {
            nd = ptr->id;
            elev_factor = (grid->node[nd].z - bot_z)/(top_z - bot_z);
            dpl[nd] = surf_dpl * elev_factor; // to match trunk
            //printf("dpl[%d]: %30.20e top_z: %20.10f bot_z: %20.10f \n",nd,dpl[nd],top_z, bot_z);
            //dpl[nd] = 0.; // force internal points to not adapt
            ptr=ptr->next;
        }
    }
}

//*******************************************************************//
//*******************************************************************//
//*******************************************************************//

/* cjt - vertically adapt nodes by monitoring the gradient of a supplied funciton */
void tl_vertical_adapt_ALE(SGRID *grid, SSW_3D *sw3) { 
    
    ID_LIST_ITEM *ptr, *ptr_up, *ptr_down;
    int nd, nvert_nnodes = 0, count = 0, ivert = 0;
    int i, iseg;
    double ele_top, ele_bot, monitor_max = 0., new_dz_sum = 0.;
    double *monitor = NULL, *new_dz = NULL;
    double grid_dz = 0., dz1 = 0., dz2 = 0., dz3 = 0., z = 0., new_dz_norm = 0., z1 = 0., z2 = 0.;
    double elev_factor_old = 0.;
    
    double max_dz = 0.5; // this parameter means the maximum a dz (un-normalized) can change in one dt
 
    /* monitor function choice */
    for (i=0; i<grid->nnodes; i++) {
        sw3->darray[i] = sqrt(pow(sw3->vel[i].x,2) + pow(sw3->vel[i].y,2));
    }

    
#ifdef _DEBUG
    if (DEBUG) {
        printf("\ndebugging tl_adapt_elev_factors ************************** \n");
    }
#endif
    
    // averaging weights (more weight on adjacent nodes slow down adaption)
    double w1 = 0.2;
    double w2 = 0.6;
    double w3 = 0.2;
    
    // monitor coefficient
    double coef = 0.1;
    
    // // this is assumed initially constant, come back later for initial nonconstant
    ptr = grid->vertical_list[0];
    grid_dz = fabs(grid->node[ptr->id].z - grid->node[ptr->next->id].z);
    
    // minimum vertical spacing allowed
    double dz_min = grid_dz * .10; //0.5;// dz_min_percentage; //0.005;
    
#ifdef _DEBUG
    if (DEBUG) {
        printf("grid_dz: %20.10f \t dz_min: %20.10f \n",grid_dz, dz_min);
    }
#endif
    
    // ---------------------------------------------------
    
    for (iseg=0; iseg<grid->nnodes_sur; iseg++) {
        
        ptr = grid->vertical_list[iseg];
        
#ifdef _DEBUG
        if (DEBUG) {
            printf("top displacement: %20.10f\n",sw3->displacement[ptr->id]);
        }
#endif
        
        // count the vertical nodes in this column
        nvert_nnodes = 0;
        while (ptr->next != NULL) {
            ptr = ptr->next;
            nvert_nnodes++;
            
        }
#ifdef _DEBUG
        if (DEBUG) {
            printf("number of vertical nodes: %d \n",nvert_nnodes);
            ptr = grid->vertical_list[iseg];
             while(ptr->next != NULL) {
                printf("f[%d]: %20.10f\n",ivert,sw3->darray[ptr->id]);
                ptr=ptr->next;
            }
        }
#endif
        
        // allocate monitor array
        monitor = (double *) tl_alloc(sizeof(double), nvert_nnodes);
        new_dz = (double *) tl_alloc(sizeof(double), nvert_nnodes);
        
        // calculate monitor using finite differences
        count = 0.; monitor_max = 0.;
        ptr = grid->vertical_list[iseg];
        while (((ptr->next)->next)->next != NULL) {
            
            ptr = ptr->next; // node n
            
            ptr_up = ptr->prev;   // node n-1
            ptr_down = ptr->next; // node n+1
            
            count++;
            
            // the dz here should probably be on moving grid ...
            z1 = grid->node[ptr_up->id].z + sw3->displacement[ptr_up->id];
            z2 = grid->node[ptr_down->id].z + sw3->displacement[ptr_down->id];
            monitor[count] = coef * fabs((sw3->darray[ptr_up->id] - sw3->darray[ptr_down->id])/(z1-z2)); // estimate gradient
            if (fabs(monitor[count]) > monitor_max) monitor_max = fabs(monitor[count]);
#ifdef _DEBUG
            if (DEBUG) {
                printf("z1: %20.10f  z2: %20.10f f1: %20.10f  f2: %20.10f \n",z1, z2, sw3->darray[ptr_up->id], sw3->darray[ptr_down->id]);
            }
#endif

        }
        monitor[0] = 0.;
        monitor[nvert_nnodes-1] = 0.;
        
        if (count+2 != nvert_nnodes) {
            tl_error(">> BUG: count+2 != nvert_nnodes in file: tl_create_elev_factors.");
        }
        
        // normalize monitor (function gradients)
        if (monitor_max > 1e-6) {
            for (ivert=0; ivert<nvert_nnodes; ivert++) {
                monitor[ivert] /= monitor_max;
            }
        } else {
            // the function is constant down the column, we're done
            //return;
        }
        
#ifdef _DEBUG
        if (DEBUG) {
            for (ivert=0; ivert<nvert_nnodes; ivert++) {
                printf("normalized monitor[%d]: %20.10f\n",ivert,monitor[ivert]);
            }
        }
#endif
        
        // calculate new dz's
        new_dz[0] = grid_dz;
        new_dz[nvert_nnodes-1] = grid_dz;
        new_dz_sum = 0.;
        for (ivert=1; ivert<nvert_nnodes-1; ivert++) {
            dz1 = grid_dz * (1 - monitor[ivert-1] * max_dz);
            dz2 = grid_dz * (1 - monitor[ivert] * max_dz);
            dz3 = grid_dz * (1 - monitor[ivert+1] * max_dz);
            new_dz[ivert] = w1*dz1 + w2*dz2 + w3*dz3;
            if (new_dz[ivert] < dz_min) new_dz[ivert] = dz_min;
            new_dz_sum += new_dz[ivert];
        }
        
#ifdef _DEBUG
        if (DEBUG) {
            for (ivert=0; ivert<nvert_nnodes; ivert++) {
                printf("new un-normalized dz[%d]: %20.10f\n",ivert,new_dz[ivert]);
            }
        }
#endif
        
        // calculate new displacement factors
        ptr = grid->vertical_list[iseg];
        nd = ptr->id;
        ele_top = grid->node[nd].z;
        while(ptr->next != NULL) {
            nd = ptr->id;
            ele_bot = grid->node[nd].z;
            ptr = ptr->next;
        }
        
        new_dz_norm = (new_dz_sum / ((ele_top - ele_bot) - grid_dz));
#ifdef _DEBUG
        if (DEBUG) {
            printf("new dz norm: %0.10f\n",new_dz_norm);
        }
#endif
        
        ptr = grid->vertical_list[iseg];
        count = 0;
        z = grid->node[ptr->id].z;
        while(ptr->next != NULL) {
            nd = ptr->id;
            if (count != 0 && count != nvert_nnodes-1) { // surface node is displaced by solution
                z -= new_dz[count]/new_dz_norm; // normalize new dz's by column depth;
                sw3->displacement[nd] = .99*sw3->old_displacement[nd] + 0.01*(z - grid->node[nd].z); // weighted average to slow node movement
            } else {
                z = grid->node[ptr->id].z;
            }
            
#ifdef _DEBUG
            if (DEBUG) {
                printf("dpl: %20.10f old z + dpl: %20.10f \t z+dpl: %20.10f \n",sw3->displacement[nd],grid->node[nd].z + sw3->old_displacement[nd],grid->node[nd].z + sw3->displacement[nd]);
            }
#endif
            ptr=ptr->next;
            count ++;
        }
        // end of column
        monitor = tl_free(sizeof(double), nvert_nnodes, monitor);
        new_dz = tl_free(sizeof(double), nvert_nnodes, new_dz);
#ifdef _DEBUG
        if (DEBUG) {
            exit(-1);
        }
#endif
    }
    
#ifdef _DEBUG
    if (DEBUG) {
        printf("************************************************************ \n");
    }
#endif
    
}

