#include "extrusion.h"

void create_sw3_hot(SMODEL *mod, SGRID *grid3d, int *node_start, int *node_count) {
    
    int i,j,k,flag,icon;
	double pressure;
    
    open_output_file(&(mod->io->fout_hot), "3D hotstart file", FALSE);
    
    // Depths/displacements incorporated into the grid file
    //fprintf(mod->io->fout_hot.fp,"");
 
    // Pressures
        fprintf(mod->io->fout_hot.fp,"DATASET\nOBJTYPE \"mesh3d\"\nBEGSCL\nND %d\nNC %d\nNAME IP\nTS 0  0.0\n",grid3d->nnodes,grid3d->nelems3d);
        for (i=0; i<mod->grid->nnodes; i++) {

            for(j=node_start[i]; j<node_start[i]+node_count[i]; j++) {
                //printf("j: %d \t grid3d->node[j].z: %20.10f \t grid3d->node[node_start[i]].z: %20.10f \n",j,grid3d->node[j].z,grid3d->node[node_start[i]].z);
				pressure = - 1000. * 9.8 * (grid3d->node[j].z - grid3d->node[node_start[i]].z);
                fprintf(mod->io->fout_hot.fp,"%20.10f\n",pressure);
            }
        }
        fprintf(mod->io->fout_hot.fp,"ENDDS\n");

    // Velocities
    flag = OFF;
    for (i=0; i<mod->grid->nnodes; i++) {
        if (fabs(mod->sw->d2->vel[i].x) > SMALL8 || fabs(mod->sw->d2->vel[i].y) > SMALL8) flag = ON;
    }
    if (flag == ON) {
        fprintf(mod->io->fout_hot.fp,"DATASET\nOBJTYPE \"mesh3d\"\nBEGVEC\nND %d\nNC %d\nNAME IV\nTS 0  0.0\n",grid3d->nnodes,grid3d->nelems3d);
        for (i=0; i<mod->grid->nnodes; i++) {
            for(j=node_start[i]; j<node_start[i]+node_count[i]; j++) {
                fprintf(mod->io->fout_hot.fp,"%20.10f  %20.10f  %20.10f\n",mod->sw->d2->vel[i].x,mod->sw->d2->vel[i].y,DZERO);
            }
        }
        fprintf(mod->io->fout_hot.fp,"ENDDS\n");
    }
    
    // Old Velocities
    flag = OFF;
    for (i=0; i<mod->grid->nnodes; i++) {
        if (fabs(mod->sw->d2->old_vel[i].x) > SMALL8 || fabs(mod->sw->d2->old_vel[i].y) > SMALL8) flag = ON;
    }
    if (flag == ON) {
        fprintf(mod->io->fout_hot.fp,"DATASET\nOBJTYPE \"mesh3d\"\nBEGVEC\nND %d\nNC %d\nNAME IPV\nTS 0  0.0\n",grid3d->nnodes,grid3d->nelems3d);
        for (i=0; i<mod->grid->nnodes; i++) {
            for(j=node_start[i]; j<node_start[i]+node_count[i]; j++) {
                fprintf(mod->io->fout_hot.fp,"%20.10f  %20.10f  %20.10f\n",mod->sw->d2->old_vel[i].x,mod->sw->d2->old_vel[i].y,DZERO);
            }
        }
        fprintf(mod->io->fout_hot.fp,"ENDDS\n");
    }
    
    // Constituents
    for (icon=0; icon<mod->ntransport; icon++) {
        flag = OFF;
        for (i=0; i<mod->grid->nnodes; i++) {
            if (fabs(mod->con[icon].concentration[i]) > SMALL8) flag = ON;
        }
        if (flag == ON) {
            fprintf(mod->io->fout_hot.fp,"DATASET\nOBJTYPE \"mesh3d\"\nBEGSCL\nND %d\nNC %d\nNAME ICON %d \nTS 0  0.0\n",grid3d->nnodes,grid3d->nelems3d,icon);
            for (i=0; i<mod->grid->nnodes; i++) {
                for(j=node_start[i]; j<node_start[i]+node_count[i]; j++) {
                    fprintf(mod->io->fout_hot.fp,"%20.10f\n",mod->con[icon].concentration[i]);
                }
            }
            fprintf(mod->io->fout_hot.fp,"ENDDS\n");
        }
    }
}
