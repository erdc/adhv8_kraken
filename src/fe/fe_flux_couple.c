#include "global_header.h"

// cjt || FLUX_TAG

static int DEBUG = 0;
static int sw2_superModel_ID = UNSET_INT;
static int sw3_superModel_ID = UNSET_INT;

void calculate_supermodel_fluxes(int isuper_model, int isubmodel, SMODEL *mod) {
    
    DEBUG = 1;
    
    int ie, n1, n2, n3, macro;
    double h1, h2, h3, u1, u2, u3, v1, v2, v3, w1, w2, w3, nx, ny, nz, c1;
    double density = mod->density;
    
    total_mass_flux_from_2d = 0.;
    total_mass_flux_from_3d = 0.;
    
    if (DEBUG) printf("\n");
    if (mod->grid->ndim == 2) {
        //*****************************************************************//
        // SW 2D **********************************************************//
        
        sw2_superModel_ID = isuper_model;
        for (ie = 0; ie < mod->grid->nelems1d; ie++) {
            // only calculate at Tailwater/Elevation boundaries (assume this the interface for now)
            if (mod->str_values[mod->grid->elem1d[ie].string].ol_flow.bc_flag == BCT_FLUX_COUPLE) {printf("ADDING FLUX\n");
                macro = mod->grid->elem1d[ie].gid;
                
                n1 = mod->grid->elem1d[ie].nodes[0];
                n2 = mod->grid->elem1d[ie].nodes[1];
                
                h1 = mod->sw->d2->head[n1];
                h2 = mod->sw->d2->head[n2];
                
                u1 = mod->sw->d2->vel[n1].x;
                u2 = mod->sw->d2->vel[n2].x;
                
                v1 = mod->sw->d2->vel[n1].y;
                v2 = mod->sw->d2->vel[n2].y;
                
                nx = mod->grid->elem1d[ie].nrml.x;
                ny = mod->grid->elem1d[ie].nrml.y;
                
                // calculated SW 2D normal edge flux (edge normal depth-averaged mass flux  = int[ rho * h * (vbar dot n) ] )
                // note :: djac = ds = sqrt(diff(x,that)^2 + diff(y,that)^2) = dist
                c1 = (1./6.) * density * mod->grid->elem1d[ie].djac;
                coupled_normal_flux[isuper_model][isubmodel][ie] = c1 * ((2*h1 + h2)*nx*u1 + (h1 + 2*h2)*nx*u2 +
                                                                         (2*h1 + h2)*ny*v1 + (h1 + 2*h2)*ny*v2);
                
                total_mass_flux_from_2d += coupled_normal_flux[isuper_model][isubmodel][macro];
                
                if (DEBUG) printf("nelems1d: %d || coupled_normal_flux[iSuperModel=%d][iElem1d=%d]: %20.10f \n",mod->grid->nelems1d,isuper_model,ie,coupled_normal_flux[isuper_model][ie][macro]);
            }
        }
        if (DEBUG) printf("total_mass_flux: %20.10f \n",total_mass_flux_from_2d);
        
    } else if (mod->grid->ndim == 3) {
        //*****************************************************************//
        // SW 3D **********************************************************//
        
        sw3_superModel_ID = isuper_model;
        for (ie = 0; ie < mod->grid->nelems2d; ie++) {
            coupled_normal_flux[isuper_model][isubmodel][ie] = 0.;
            // only calculate at Tailwater/Elevation boundaries (assume this the interface for now)
            if (mod->str_values[mod->grid->elem2d[ie].string].ol_flow.bc_flag == BCT_FLUX_COUPLE) {
                macro = mod->grid->elem2d[ie].gid;
                n1 = mod->grid->elem2d[ie].nodes[0];
                n2 = mod->grid->elem2d[ie].nodes[1];
                n3 = mod->grid->elem2d[ie].nodes[2];
                
                u1 = mod->sw->d3->vel[n1].x;
                u2 = mod->sw->d3->vel[n2].x;
                u2 = mod->sw->d3->vel[n3].x;
                
                v1 = mod->sw->d3->vel[n1].y;
                v2 = mod->sw->d3->vel[n2].y;
                v3 = mod->sw->d3->vel[n3].y;
                
                w1 = mod->sw->d3->vel[n1].z;
                w2 = mod->sw->d3->vel[n2].z;
                w3 = mod->sw->d3->vel[n3].z;
                
                nx = mod->grid->elem2d[ie].nrml.x;
                ny = mod->grid->elem2d[ie].nrml.y;
                nz = mod->grid->elem2d[ie].nrml.z;
                
                // calculate SW 3D normal face fluxes (face normal mass flux  = int[ rho * (v dot n) ] )
                c1 = (1./3.) * density * mod->grid->elem2d[ie].djac3d;
                coupled_normal_flux[isuper_model][isubmodel][macro] = c1 * (nx*u1 + nx*u2 + nx*u3 + ny*v1 + ny*v2 + ny*v3 + nz*w1 + nz*w2 + nz*w3);
                
                //coupled_normal_flux[isuper_model][ie]
                
                
                total_mass_flux_from_3d += coupled_normal_flux[isuper_model][isubmodel][macro];
            }
        }
        
    }
    //if (DEBUG) exit(-1);
    
}


