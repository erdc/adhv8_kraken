/* This routine loops over submodels to calculate the shallow water residuals */
#include "global_header.h"

/*****************************************************************************************/
/*****************************************************************************************/
static int DEBUG           = OFF;
static int DEBUG_INTERFACE = OFF; /* Only for coupled models */

/* Gajanan gkc WARNING :*/
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* KEEP THESE VALUES SAME AS IN FE_SW_HYBRID_WVEL_LOAD.C ! */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

static int is_wvel_dirichlet  = YES;

/*****************************************************************************************/
/*****************************************************************************************/

void fe_sw_hybrid_wvel_resid(SSUPER_MODEL *sm, int idum) {
#ifndef _PETSC
    int i, j, k, l, idof, ndof, nodeID_2D, isub,inode;
    
    SMODEL *mod, *mod2d, *mod3d;
    SINTERFACE_NODELIST *ndlist;
    assert(sm->nsys == 1);
    
    sarray_init_dbl(sm->residual, sm->total_nnodes * sm->max_nsys);
    
    // 3D SW Contribution
    for (i=0; i<sm->nsubmodels; i++){
        mod = &(sm->submodel[i]);
        if(mod->flag.SW3_FLOW){
            fe_wvel_resid(sm,i);
        }
    }
    
    
    // 2D SW Contribution distributed evenly down interfacial column
    int model_2d_id = UNSET_INT;
    int model_3d_id = UNSET_INT;
    if (is_wvel_dirichlet == NO) {
        
        for (j=0; j<sm->NumInterfaces; j++){
            if (sm->submodel[sm->interface[j].model_id[0]].flag.SW2_FLOW &&
                sm->submodel[sm->interface[j].model_id[1]].flag.SW3_FLOW){
                model_2d_id = 0;
                model_3d_id = 1;
            }
            else if (sm->submodel[sm->interface[j].model_id[0]].flag.SW3_FLOW &&
                     sm->submodel[sm->interface[j].model_id[1]].flag.SW2_FLOW){
                model_2d_id = 1;
                model_3d_id = 0;
            }
            else{
                continue; /* this interface is 2D-2D or 3D-3D. */
            }
            
            // What follows is only executed for 2D-3D coupling.
            //int *fmap = mod->fmap_wvel; // bugged, not used
            mod = &(sm->submodel[sm->interface[j].model_id[model_2d_id]]); /* 2d model alias */
#ifdef _DEBUG
            assert(mod);
#endif
            
            for (k=0; k<sm->interface[j].NumNodeColumns; k++){
                ndlist = &(sm->interface[j].nodelist[k]);
                assert(ndlist);
                for (l=0; l<ndlist->size[model_3d_id]; l++){
                    sm->residual[ndlist->couplednodes[model_3d_id][l]] += mod->sw->d2->dacontResid[ndlist->couplednodes[model_2d_id][0]]/((double) ndlist->size[model_3d_id]);
                }
            }
            
            mod = NULL;
        }
        
        
    } else {
        
        
        for (j=0; j<sm->NumInterfaces; j++){
            if (sm->submodel[sm->interface[j].model_id[0]].flag.SW2_FLOW &&
                sm->submodel[sm->interface[j].model_id[1]].flag.SW3_FLOW){
                model_2d_id = 0;
                model_3d_id = 1;
            } else if (sm->submodel[sm->interface[j].model_id[0]].flag.SW3_FLOW &&
                       sm->submodel[sm->interface[j].model_id[1]].flag.SW2_FLOW){
                model_2d_id = 1;
                model_3d_id = 0;
            }
            else{
                continue; /* this interface is 2D-2D or 3D-3D. */
            }
            
            int *fmap = mod->fmap_wvel;  // bugged, not used
            SINTERFACE_NODELIST *ndlist;
            mod3d = &(sm->submodel[sm->interface[j].model_id[model_3d_id]]); /* 2d model alias */
            assert(mod3d);
            
            
            for (k=0; k<sm->interface[j].NumNodeColumns; k++){
                ndlist = &(sm->interface[j].nodelist[k]);
                assert(ndlist);
                if (ndlist->size[model_3d_id] > 0){
                    int surfnode = ndlist->couplednodes[model_3d_id][0];
                    int botnode  = ndlist->couplednodes[model_3d_id][ndlist->size[model_3d_id]-1];
                    
                    assert(surfnode > -1 && surfnode < mod3d->grid->nnodes);
                    assert(botnode > -1 && botnode < mod3d->grid->nnodes);
                    
                    double zsurf = mod3d->grid->node[surfnode].z + mod3d->sw->d3->displacement[surfnode];
                    double zbot  = mod3d->grid->node[botnode].z  + mod3d->sw->d3->displacement[botnode];
                    double depth = zsurf - zbot;
                    
                    assert(depth > 0);
                    
                    // at sur :: dz/dt + u * dz/dx + v * dz/dy + w = 0 || Dz/dt = w
                    double dpl_32dt;  NUM_GET_SECONDORDER(dpl_32dt, mod3d->sw->d3->displacement[surfnode], mod3d->sw->d3->old_displacement[surfnode]);
                    double dpl_12dt;  NUM_GET_SECONDORDER(dpl_12dt, mod3d->sw->d3->old_displacement[surfnode], mod3d->sw->d3->older_displacement[surfnode]);
                    double w_sur = (dpl_32dt - dpl_12dt)/mod3d->dt;
                    
                    // at bed :: dz/dt + u * dz/dx + v * dz/dy - w = 0 || Dz/dt = w || u * dz/dx + v * dz/dy - w = 0 for no bed movement
                    double w_bed =  (mod3d->sw->d3->vel[botnode].x * mod3d->sw->d3->grad_bed[botnode].x + mod3d->sw->d3->vel[botnode].y * mod3d->sw->d3->grad_bed[botnode].y);
                    
                    for (i = 0; i<ndlist->size[model_3d_id]; i++) {
                        int inode = ndlist->couplednodes[model_3d_id][i];
                        double znode = mod3d->grid->node[inode].z + mod3d->sw->d3->displacement[inode];
                        double toplength = zsurf - znode;
                        double botlength = znode - zbot;
                        double topfactor = toplength / depth;
                        double botfactor = botlength / depth;
                        
                        assert(inode > -1 && inode < mod3d->grid->nnodes);
                        
                        sm->residual[inode] = mod3d->sw->d3->vel[inode].z - (w_sur * botfactor + w_bed * topfactor);
                        //printf("inode: %d \t w_sur: %20.10f \t w_bed: %20.10f \t botfactor: %20.10f \t topfactor %20.10f \t bonode vel: %20.10f %20.10f %20.10f \t grad_bed: %20.10f %20.10f\n",inode, w_sur, w_bed, topfactor, botfactor,mod3d->sw->d3->vel[botnode].x,mod3d->sw->d3->vel[botnode].y,mod3d->sw->d3->vel[botnode].z,mod3d->sw->d3->grad_bed[botnode].x,mod3d->sw->d3->grad_bed[botnode].y);
                    }
                }
            }
            mod3d = NULL;
        }
    }
    
    
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<sm->supersmpi_wvel->npes;proc_count++){
            if(proc_count==sm->supersmpi_wvel->myid){
                printf("***********myid %d",sm->supersmpi_wvel->myid);
#endif
                printScreen_resid("Final WVEL residual", sm->residual, sm->wvel_nnodes, sm->nsys, __LINE__, __FILE__);
#ifdef _MESSG
            }
            messg_barrier(MPI_COMM_WORLD);
        }
        exit(-1);
#endif
    }
#endif
#endif
}
