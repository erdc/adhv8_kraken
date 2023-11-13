#include "global_header.h"

#define RK_TOP_ORDER 4
SVECT get_particle_value_on_elem_vector_2d(SGRID *grid, SPARTICLE *p, SVECT r, int ie, double tL, SVECT *vL, double tR, SVECT *vR, double t, int vel_flag);
SVECT get_particle_value_on_elem_vector_3d(SGRID *grid, SPARTICLE *p, SVECT r, int ie, double tL, SVECT *vL, double tR, SVECT *vR, double t, int vel_flag);
double get_particle_value_on_elem_scalar_2d(SGRID *grid, SVECT r, int ie, double tL, double *fL, double tR, double *fR, double t);
double get_particle_value_on_elem_scalar_3d(SGRID *grid, SVECT r, int ie, double tL, double *fL, double tR, double *fR, double t);
int project_vel_3D(SGRID *grid, int ie, int lnd, int ledge, int lface, SVECT r, SVECT *pv);


//SVECT get_point_velocity(SGRID *grid, int *ielem, SVECT *velocity, SVECT p, double t, double *, int);
void interpolate_velocity_field(SVECT *vL, SVECT *vR, double tL, double tR, double t, int nnodes, SVECT *v);

void error_mssg_point_not_in_elem(const char *file, const unsigned long line, int flag, SVECT r, int ie, SGRID *grid) {
    if (flag != 1) {
        printf("ERROR || FILE: %s || LINE: %lu || vector evaluation point is not on given element: %d\n",file,line,ie);
        printf("p: %20.15f %20.15f %20.15f\n",r.x,r.y,r.z);
        if (grid->ndim == 2) {
          int i,nd;
          for (i=0; i<3; i++) {
             nd = grid->elem2d[ie].nodes[i];
            printf("element node 1 || gid: %d ||  %20.15f %20.15f %20.15f\n",nd,grid->node[nd].x,grid->node[nd].y,grid->node[nd].z);
          }
        }
        if (grid->ndim == 3) {
          int i,nd;
          for (i=0; i<4; i++) {
             nd = grid->elem3d[ie].nodes[i];
            printf("element node 1 || gid: %d ||  %20.15f %20.15f %20.15f\n",nd,grid->node[nd].x,grid->node[nd].y,grid->node[nd].z);
          }
        }
        exit(-1);
    }
}

void printScreen_Search(SGRID *grid, int ip, char *STRING, double t, int ieo, int ie, SVECT ro, SVECT r, SVECT v, char *ivalue, int idum) {
    printf("-- %s  ",STRING);
    if (ivalue != NULL) {
        printf("%s: %d || ip: %d || t: %f || r(t0): {%f %f %f} || r(t): {%f %f %f} || ie(t0): %d || ie(t): %d || v(t): {%f %f %f}\n",
               ivalue,idum,ip,t,ro.x,ro.y,ro.z,r.x,r.y,r.z,ieo,ie,v.x,v.y,v.z);
    } else {
        printf("ip: %d || t: %f || r(t0): {%f %f %f} || r(t): {%f %f %f} || ie(t0): %d || ie(t): %d || v(t): {%f %f %f}\n",
               ip,t,ro.x,ro.y,ro.z,r.x,r.y,r.z,ieo,ie,v.x,v.y,v.z);
    }
    int i;
    if (grid->ndim == 3) {
        for (i=0; i<4; i++) {
            printf("ie %d node[%d] :: {%20.15e %20.15e %20.15f}\n",ie,i,
                   grid->node[grid->elem3d[ie].nodes[i]].x,
                   grid->node[grid->elem3d[ie].nodes[i]].y,
                   grid->node[grid->elem3d[ie].nodes[i]].z);
        }
    } else {
        for (i=0; i<3; i++) {
            printf("ie %d node[%d] :: {%20.15e %20.15e %20.15f}\n",ie,i,
                   grid->node[grid->elem2d[ie].nodes[i]].x,
                   grid->node[grid->elem2d[ie].nodes[i]].y,
                   grid->node[grid->elem2d[ie].nodes[i]].z);
        }
    }
}

void get_runge_kutta_coeffs(double rk_coeffs[RK_TOP_ORDER+1][RK_TOP_ORDER], double rk_dt_coeffs[RK_TOP_ORDER+1][RK_TOP_ORDER]) {
    // [5][4]
    
    // rk 1 - nothing
    rk_coeffs[0][0] = 0.0; rk_dt_coeffs[0][0] = 0.0;
    
    // rk 1 - Euler Method
    rk_coeffs[1][0] = 1; rk_dt_coeffs[1][0] = 0.0;
    
    // rk 2 - Explicit Midpoint Method
    rk_coeffs[2][0] = 0.0; rk_dt_coeffs[2][0] = 0.0;
    rk_coeffs[2][1] = 1.0; rk_dt_coeffs[2][1] = 0.5;
    
    // rk 3 - Heun's third order method
    rk_coeffs[3][0] = 1.0/4.0;  rk_dt_coeffs[3][0] = 0.0;
    rk_coeffs[3][1] = 0.0;      rk_dt_coeffs[3][1] = 1.0/3.0;
    rk_coeffs[3][2] = 3.0/4.0;  rk_dt_coeffs[3][2] = 2.0/3.0;
    
    // rk 4
    rk_coeffs[4][0] = 1.0/6.0; rk_dt_coeffs[4][0] = 0.0;
    rk_coeffs[4][1] = 1.0/3.0; rk_dt_coeffs[4][1] = 0.5;
    rk_coeffs[4][2] = 1.0/3.0; rk_dt_coeffs[4][2] = 0.5;
    rk_coeffs[4][3] = 1.0/6.0; rk_dt_coeffs[4][3] = 1.0;
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
SVECT update_simple(SVECT p, double dt, SVECT v, double c) {
    SVECT pnew;
    pnew.x = p.x + c * dt * v.x;
    pnew.y = p.y + c * dt * v.y;
    pnew.z = p.z + c * dt * v.z;
    return pnew;
}

void update_simple_replace(SVECT *p, double dt, SVECT v, double c) {
    SVECT pnew;
    p->x += c * dt * v.x;
    p->y += c * dt * v.y;
    p->z += c * dt * v.z;
}

int update_and_search(SGRID *grid, SVECT point_velocity, double dt, SVECT *p, int *p_elem) {
    int flag = -999;
    SVECT p_orig = *p,pi;
    
//   printf("point_velocity: {%20.10e, %20.10e, %20.10e} dt: %f\n",point_velocity.x,point_velocity.y,point_velocity.z,dt);
//    double weights[3] = {0., 0., 0.};
//    int lnd = UNSET_INT, ledge = UNSET_INT;
//     flag = compute_interpolate_weights_2D_triangle(grid, *p_elem, p->x, p->y, weights, &lnd, &ledge);
//     assert(flag==1);
    
    if (fabs(point_velocity.x) + fabs(point_velocity.y) + fabs(point_velocity.z) > 1e-10) {
        *p = update_simple(p_orig,dt,point_velocity,1.0);
        if (grid->ndim == 2) {
            flag = elementSearch2D(grid,p_orig,*p,p_elem,true,&pi);
        } else {
            flag = elementSearch3D(grid,p_orig,*p,p_elem,true,&pi);
        }
        if (flag == -2) {
            *p = pi;
            if (grid->ndim == 2) p->z = pi.z;
        }
    }
    
    return flag;
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
void time_update(SMODEL *mod, int rk_order) {
    
    
    int DEBUG = OFF;
    int i,rk,ip,ie,ie_orig,did_elem_change_flag,flag;
    char string[100];
    double t = 0.;
    SVECT point_velocity, r_orig, r, pi, k[rk_order];
    svect_init(&point_velocity);
    SGRID *grid = mod->grid;
    SPARTICLE *p = NULL;
    
    if (DEBUG == ON) printf("\n");
    
    // get Runge-Kutta coefficients
    assert(rk_order > 0);
    double rk_coeffs[RK_TOP_ORDER+1][RK_TOP_ORDER], rk_dt_coeffs[RK_TOP_ORDER+1][RK_TOP_ORDER];
    get_runge_kutta_coeffs(rk_coeffs,rk_dt_coeffs);
    
    // instead of updating dpl at each RK stage (cause maddening issues), use dpl at time + dt
    // to update the grid before stepping.  Not as accurate, but fine on a boundary.
    if (mod->grid->vertical_dpl_flag == ON) {
        t = mod->time + mod->dt;
        
        // update the whole grid befor the particle loop
        for (i=0; i<mod->grid->nnodes; i++) {
            mod->dpl[i] = interpolate1D(t,mod->tL_dpl,mod->tR_dpl,mod->dpl_tL[i],mod->dpl_tR[i]);
            mod->grid->node[i].z = mod->grid->node_t0[i].z + mod->dpl[i];
        }
        
        for (ip=0; ip<mod->np; ip++) {
            p = &mod->p[ip];
            did_elem_change_flag = 0;
            if (DEBUG == ON) printScreen_Search(grid,ip,"before checkElementAfterDisplacement",t,p->elemID,p->elemID,p->r,p->r,point_velocity,"change_flag",did_elem_change_flag);
            did_elem_change_flag = checkElementAfterDisplacement(mod->grid, &p->r, &p->elemID);
            if (DEBUG == ON) printScreen_Search(grid,ip,"after checkElementAfterDisplacement",t,p->elemID,p->elemID,p->r,p->r,point_velocity,"change_flag",did_elem_change_flag);
        }
    }
    
    // loop over particles
    for (ip=0; ip<mod->np; ip++) {
        p = &mod->p[ip];
        r_orig = p->r;
        ie_orig = p->elemID;
        
        // *****************************************************************************
        // check to see if the particle hits a boundary in this dt, if it does,
        // take a fraction step at current velocity to the boundary, then use projected
        // velocity for remaining dt
        r = r_orig;
        ie = ie_orig;
        if (grid->ndim == 2) {
            point_velocity = get_particle_value_on_elem_vector_2d(grid,p,r,ie,mod->tL,mod->vel_tL,mod->tR,mod->vel_tR,mod->time,1);
        } else {
            point_velocity = get_particle_value_on_elem_vector_3d(grid,p,r,ie,mod->tL,mod->vel_tL,mod->tR,mod->vel_tR,mod->time,1);
        }
        flag = UNSET_INT;
        flag = update_and_search(grid,point_velocity,mod->dt,&r,&ie); // returns r at boundary intersection
        
        if (flag == -2) {
            if (DEBUG == ON) {
                printf("\n*************************************************************\n");
                printf("Particle hit boundary during RK update ... using Euler method\n");
                printf("*************************************************************\n");
                if (DEBUG == ON) printScreen_Search(grid,ip,"EULER START",t,ie_orig,ie,r_orig,r,point_velocity,"flag",flag);
            }
            
            double dt_frac = 0.;
            if (fabs(point_velocity.x) > 1e-10) {
                dt_frac = (r.x - r_orig.x) / point_velocity.x;
            } else if (fabs(point_velocity.y) > 1e-10) {
                dt_frac = (r.y - r_orig.y) / point_velocity.y;
            } else if (fabs(point_velocity.z) > 1e-10) {
                dt_frac = (r.z - r_orig.z) / point_velocity.z;
            }
            if (fabs(dt_frac - mod->dt) > 1e-10 && fabs(dt_frac) > 1e-10) {
                if (DEBUG == ON) printf("moving dt_frac: %10.5f\n",dt_frac);
                if (DEBUG == ON) printScreen_Search(grid,ip,"EULER 1",t,ie_orig,ie,r_orig,r,point_velocity,"flag",flag);
                if (grid->ndim == 2) {
                    point_velocity = get_particle_value_on_elem_vector_2d(grid,p,r,ie,mod->tL,mod->vel_tL,mod->tR,mod->vel_tR,mod->time + dt_frac,1);
                } else {
                    point_velocity = get_particle_value_on_elem_vector_3d(grid,p,r,ie,mod->tL,mod->vel_tL,mod->tR,mod->vel_tR,mod->time + dt_frac,1);
                }
                flag = UNSET_INT;
                flag = update_and_search(grid,point_velocity,mod->dt - dt_frac,&r,&ie); // returns r at boundary intersection
                if (DEBUG == ON) printScreen_Search(grid,ip,"EULER 2",mod->dt - dt_frac,ie_orig,ie,r_orig,r,point_velocity,"flag",flag);
            }
            
        } else {
            
            if (DEBUG == ON) {
                printf("\n*************************************************************\n");
                printf("Particle in RK update\n");
                printf("*************************************************************\n");
            }
            
            r = r_orig;
            ie = ie_orig;
            
            svect_init(&point_velocity);
            for (rk=0; rk<rk_order; rk++) {
                
                t = mod->time + rk_dt_coeffs[rk_order][rk] * mod->dt;
                
                // Euler update using last k
                r = r_orig;
                ie = ie_orig;
                
                if (rk != 0) {
                    flag = UNSET_INT;
                    flag = update_and_search(grid,point_velocity,rk_dt_coeffs[rk_order][rk]*mod->dt,&r,&ie);
                    //assert(flag == -1); // cjt just added
                    assert(flag != -2); // we have already made sure it doesn't hit a boundary above
                }
                
                // now evaluate point velocity at this point
                if (grid->ndim == 2) {
                    
                    // call water quality first to update fields
                    if (mod->oxygen != NULL)
                        p->exp_oxygen   = get_particle_value_on_elem_scalar_2d(grid,r,ie,mod->tL_wq,mod->oxygen_tL,mod->tR_wq,mod->oxygen_tR,t);
                    if (mod->salinity != NULL)
                        p->exp_salinity = get_particle_value_on_elem_scalar_2d(grid,r,ie,mod->tL_wq,mod->salinity_tL,mod->tR_wq,mod->salinity_tR,t);
                    if (mod->sunlight != NULL)
                        p->exp_sunlight = get_particle_value_on_elem_scalar_2d(grid,r,ie,mod->tL_wq,mod->sunlight_tL,mod->tR_wq,mod->sunlight_tR,t);
                    
                    // now call velocity
                    point_velocity = get_particle_value_on_elem_vector_2d(grid,p,r,ie,mod->tL,mod->vel_tL,mod->tR,mod->vel_tR,t,1);
                    
                } else if (grid->ndim == 3) {
                    
                    // call water quality first to update fields
                    if (mod->oxygen != NULL)
                        p->exp_oxygen   = get_particle_value_on_elem_scalar_3d(grid,r,ie,mod->tL_wq,mod->oxygen_tL,mod->tR_wq,mod->oxygen_tR,t);
                    if (mod->salinity != NULL)
                        p->exp_salinity = get_particle_value_on_elem_scalar_3d(grid,r,ie,mod->tL_wq,mod->salinity_tL,mod->tR_wq,mod->salinity_tR,t);
                    if (mod->sunlight != NULL)
                        p->exp_sunlight = get_particle_value_on_elem_scalar_3d(grid,r,ie,mod->tL_wq,mod->sunlight_tL,mod->tR_wq,mod->sunlight_tR,t);
                    
                    // now call velocity
                    point_velocity = get_particle_value_on_elem_vector_3d(grid,p,r,ie,mod->tL,mod->vel_tL,mod->tR,mod->vel_tR,t,1);
                    
                }
                k[rk].x = point_velocity.x;
                k[rk].y = point_velocity.y;
                k[rk].z = point_velocity.z;
                
                if (DEBUG == ON) printScreen_Search(grid,ip,"RK",t,ie_orig,ie,r_orig,r,point_velocity,"STAGE",rk+1);
            }
            
            assert(flag != -2); // we have already made sure it doesn't hit a boundary above
            
            svect_init(&point_velocity);
            for (rk=0; rk<rk_order; rk++) {
                point_velocity.x += rk_coeffs[rk_order][rk] * k[rk].x;
                point_velocity.y += rk_coeffs[rk_order][rk] * k[rk].y;
                point_velocity.z += rk_coeffs[rk_order][rk] * k[rk].z;
            }
            
            r = r_orig;
            ie = ie_orig;
            flag = UNSET_INT;
            flag = update_and_search(grid,point_velocity,mod->dt,&r,&ie);
            if (DEBUG == ON) printScreen_Search(grid,ip,"RK FINAL",mod->time + mod->dt,ie_orig,ie,r_orig,r,point_velocity,NULL,0);
        }
        
        t = mod->time + mod->dt;
        p->r = r;
        p->elemID = ie;

        
        // update particle age
        p->age += mod->dt;
    }
    
    // update simulation time
    mod->time += mod->dt;
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
void interpolate_velocity_field(SVECT *vL, SVECT *vR, double tL, double tR, double t, int nnodes, SVECT *v) {
    int i;
    double factor = (t - tL) / (tR - tL);
    for (i=0; i<nnodes; i++) {
        v[i].x = vL[i].x + factor * (vR[i].x - vL[i].x);
        v[i].y = vL[i].y + factor * (vR[i].y - vL[i].y);
        v[i].z = vL[i].z + factor * (vR[i].z - vL[i].z);
    }
}


// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
int nodeExternalFaceVelProj(SGRID *grid, int ielem, int lnd, double *weights, SVECT *velocity, SVECT *v) {
    int i,j,k;
    svect_init(v);
    
    for (i=0; i<4; i++) {
        if (grid->elem3d[ielem].face_flag[i] == UNSET_INT) {
            
            // found an external boundary face on element
            int lnd1_ID = nd_on_fc[i][0];
            int lnd2_ID = nd_on_fc[i][1];
            int lnd3_ID = nd_on_fc[i][2];
            int nd1_ID = grid->elem3d[ielem].nodes[lnd1_ID];
            int nd2_ID = grid->elem3d[ielem].nodes[lnd2_ID];
            int nd3_ID = grid->elem3d[ielem].nodes[lnd3_ID];
            
            // is the node with partcle on this face?
            int gnd = grid->elem3d[ielem].nodes[lnd];
            //printf("gnd: %d || face nodes: %d %d %d\n",nd1_ID,nd2_ID,nd3_ID);
            if ((nd1_ID == gnd) || (nd2_ID == gnd) || (nd3_ID == gnd)) {
                
                SVECT nd1 = grid->node[nd1_ID];
                SVECT nd2 = grid->node[nd2_ID];
                SVECT nd3 = grid->node[nd3_ID];
                
                // External edge matched to external boundary face, replace node velocity on 3 external face nodes with projected velocities
                SVECT v_proj[4]; for (j=0; j<4; j++) v_proj[j] = velocity[ grid->elem3d[ielem].nodes[j] ];
                //printf("projecting nd %d using face nodes: %d %d %d\n",gnd,nd1_ID,nd2_ID,nd3_ID);
                project_velocity_3D(velocity[nd1_ID],velocity[nd2_ID],velocity[nd3_ID],nd1,nd2,nd3,&v_proj[lnd1_ID],&v_proj[lnd2_ID],&v_proj[lnd3_ID]);
                
                for (k=0; k<4; k++) {
                    v->x += v_proj[k].x * weights[k];
                    v->y += v_proj[k].y * weights[k];
                    v->z += v_proj[k].z * weights[k];
                }
                return 1;
            }
        }
    }
    return 0;
}



int edgeExternalFaceVelProj(SGRID *grid, int ielem, int ledge, double *weights, SVECT *velocity, SVECT *v) {
    int i,j,k;
    svect_init(v);
    
    for (i=0; i<4; i++) { // loop over faces
        if (grid->elem3d[ielem].face_flag[i] == UNSET_INT) {
            
            // found an external boundary face on element
            int lnd1_ID = nd_on_fc[i][0];
            int lnd2_ID = nd_on_fc[i][1];
            int lnd3_ID = nd_on_fc[i][2];
            int nd1_ID = grid->elem3d[ielem].nodes[lnd1_ID];
            int nd2_ID = grid->elem3d[ielem].nodes[lnd2_ID];
            int nd3_ID = grid->elem3d[ielem].nodes[lnd3_ID];
            
            // is the edge with particle on this face?
            int lnd1_ID_edge = grid->nd_on_TetEdge[ledge][0];
            int lnd2_ID_edge = grid->nd_on_TetEdge[ledge][1];
            int nd1_ID_edge = grid->elem3d[ielem].nodes[lnd1_ID_edge];
            int nd2_ID_edge = grid->elem3d[ielem].nodes[lnd2_ID_edge];
            
            if ((nd1_ID_edge == nd1_ID && nd2_ID_edge == nd2_ID) ||
                (nd1_ID_edge == nd2_ID && nd2_ID_edge == nd1_ID) ||
                (nd1_ID_edge == nd1_ID && nd2_ID_edge == nd3_ID) ||
                (nd1_ID_edge == nd3_ID && nd2_ID_edge == nd1_ID) ||
                (nd1_ID_edge == nd2_ID && nd2_ID_edge == nd3_ID) ||
                (nd1_ID_edge == nd3_ID && nd2_ID_edge == nd2_ID)) {
                SVECT nd1 = grid->node[nd1_ID];
                SVECT nd2 = grid->node[nd2_ID];
                SVECT nd3 = grid->node[nd3_ID];
                // External edge matched to external boundary face, replace node velocity on 3 external face nodes with projected velocities
                SVECT v_proj[4]; for (j=0; j<4; j++) v_proj[j] = velocity[ grid->elem3d[ielem].nodes[j] ];
                //printf("projecting edge with nodes {%d,%d} using face nodes: %d %d %d\n",nd1_ID_edge,nd2_ID_edge,nd1_ID,nd2_ID,nd3_ID);
                project_velocity_3D(velocity[nd1_ID],velocity[nd2_ID],velocity[nd3_ID],nd1,nd2,nd3,&v_proj[lnd1_ID],&v_proj[lnd2_ID],&v_proj[lnd3_ID]);
                for (k=0; k<4; k++) {
                    v->x += v_proj[k].x * weights[k];
                    v->y += v_proj[k].y * weights[k];
                    v->z += v_proj[k].z * weights[k];
                }
                return 1;
            }
        }
    }
    return 0;
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
// compute 2D velocity projection for two edge nodes at an external closed boundary
void project_velocity_2D(SVECT v1, SVECT v2, SVECT nd1, SVECT nd2, SVECT *v1_proj, SVECT *v2_proj) {
    double dx = nd2.x - nd1.x;
    double dy = nd2.y - nd1.y;
    
    // first makes sure velocity is only projecting if it exits the grid (same direction as normal)
    double nx = dy, ny = -dx; //printf("nx: %f \t ny: %f\n",nx,ny); // may need to flip these signs
    double l12,dum;
    
    if (v1.x * nx + v1.y * ny > 0) {
        l12 = pow(dx*dx + dy*dy,0.5);
        dum = (v1.x * dx + v1.y * dy) / (l12  * l12);
        v1_proj->x = dum * dx;
        v1_proj->y = dum * dy;
    } else {
        v1_proj->x = v1.x;
        v1_proj->y = v1.y;
    }
    
    if (v2.x * nx + v2.y * ny > 0) {
        l12 = pow(dx*dx + dy*dy,0.5);
        dum = (v2.x * dx + v2.y * dy) / (l12  * l12);
        v2_proj->x = dum * dx;
        v2_proj->y = dum * dy;
    } else {
        v2_proj->x = v2.x;
        v2_proj->y = v2.y;
    }
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
// compute 2D velocity projection for two edge nodes at an external closed boundary
void project_point_velocity_2D(SVECT *v, SVECT nd1, SVECT nd2) {
    double dx = nd2.x - nd1.x;
    double dy = nd2.y - nd1.y;
    
    // first makes sure velocity is only projecting if it exits the grid (same direction as normal)
    double nx = dy, ny = -dx; //printf("nx: %f \t ny: %f\n",nx,ny); // may need to flip these signs
    double l12,dum;
    
    if (v->x * nx + v->y * ny > 0) {
        l12 = pow(dx*dx + dy*dy,0.5);
        dum = (v->x * dx + v->y * dy) / (l12  * l12);
        v->x = dum * dx;
        v->y = dum * dy;
    }
}
// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
// compute 3D velocity projection for three face nodes at an external closed boundary
void project_velocity_3D(SVECT v1, SVECT v2, SVECT v3, SVECT nd1, SVECT nd2, SVECT nd3,
                         SVECT *v1_proj, SVECT *v2_proj, SVECT *v3_proj) {
    
    SVECT v21 = svect_subtract(nd2,nd1);
    SVECT v31 = svect_subtract(nd3,nd1);
    SVECT nrml = svect_cross(v21,v31);
    double nrml_2 = pow(svect_mag(nrml),2);
    double dum = 0.0, dotp = 0.0;
    
    dotp = svect_dotp(v1,nrml);
    if (dotp > 0) { // only project for flows out of grid
        dum = dotp / nrml_2;
        v1_proj->x = v1.x - dum * nrml.x;
        v1_proj->y = v1.y - dum * nrml.y;
        v1_proj->z = v1.z - dum * nrml.z;
        //printf("node 1 projected from v: {%f,%f,%f} --> {%f,%f,%f}\n",nd1,v1.x,v1.y,v1.z,v1_proj->x,v1_proj->y,v1_proj->z);
    } else {
        //printf("node 1 NOT projected\n",nd1);
        v1_proj->x = v1.x;
        v1_proj->y = v1.y;
        v1_proj->z = v1.z;
    }
    
    dotp = svect_dotp(v2,nrml);
    if (dotp > 0) { // only project for flows out of grid
        dum = dotp  / nrml_2;
        v2_proj->x = v2.x - dum * nrml.x;
        v2_proj->y = v2.y - dum * nrml.y;
        v2_proj->z = v2.z - dum * nrml.z;
        //printf("node 2 projected from v: {%f,%f,%f} --> {%f,%f,%f}\n",nd2,v2.x,v2.y,v2.z,v2_proj->x,v2_proj->y,v2_proj->z);
    } else {
        //printf("node 2 NOT projected\n",nd2);
        v2_proj->x = v2.x;
        v2_proj->y = v2.y;
        v2_proj->z = v2.z;
    }
    
    dotp = svect_dotp(v3,nrml);
    if (dotp > 0) { // only project for flows out of grid
        dum = dotp / nrml_2;
        v3_proj->x = v3.x - dum * nrml.x;
        v3_proj->y = v3.y - dum * nrml.y;
        v3_proj->z = v3.z - dum * nrml.z;
        //printf("node 3 projected from v: {%f,%f,%f} --> {%f,%f,%f}\n",nd3,v3.x,v3.y,v3.z,v3_proj->x,v3_proj->y,v3_proj->z);
    } else {
        //printf("node 3 NOT projected\n",nd3);
        v3_proj->x = v3.x;
        v3_proj->y = v3.y;
        v3_proj->z = v3.z;
    }
}

void project_point_velocity_3D(SVECT *v, SVECT nd1, SVECT nd2, SVECT nd3) {
    
    SVECT v21 = svect_subtract(nd2,nd1);
    SVECT v31 = svect_subtract(nd3,nd1);
    SVECT nrml = svect_cross(v21,v31);
    double nrml_2 = pow(svect_mag(nrml),2);
    double dum = 0.0, dotp = 0.0;
    //printf("normal: %f %f %f\n",nrml.x,nrml.y,nrml.z);
    
    dotp = svect_dotp(*v,nrml);
    if (dotp > 0) { // only project for flows out of grid
        dum = dotp / nrml_2;
        v->x = v->x - dum * nrml.x;
        v->y = v->y - dum * nrml.y;
        v->z = v->z - dum * nrml.z;
        //printf("node 1 projected from v: {%f,%f,%f} --> {%f,%f,%f}\n",nd1,v1.x,v1.y,v1.z,v1_proj->x,v1_proj->y,v1_proj->z);
    }
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
SVECT get_particle_value_on_elem_vector_2d(SGRID *grid, SPARTICLE *p, SVECT r, int ie, double tL, SVECT *vL, double tR, SVECT *vR, double t, int vel_flag) {
    
    int i,j,ie2;
    int proj_flag = 0;
    int break_flag = 0;
    SVECT pv; // the return point vector
    
    // interpolate velocity to RK time at triangle vertices
    SVECT v[3]; svect_init(v);
    for (i=0; i<3; i++) {
        v[i].x = interpolate1D(t,tL,tR,vL[grid->elem2d[ie].nodes[i]].x,vR[grid->elem2d[ie].nodes[i]].x);
        v[i].y = interpolate1D(t,tL,tR,vL[grid->elem2d[ie].nodes[i]].y,vR[grid->elem2d[ie].nodes[i]].y);
    }
    
    // interpolate in space from triangle vertices to point
    double weights[3];
    int lnd = UNSET_INT, ledge = UNSET_INT, flag = UNSET_INT;
    flag = compute_interpolate_weights_2D_triangle(grid, ie, r.x, r.y, weights, &lnd, &ledge);
    error_mssg_point_not_in_elem(__FILE__,__LINE__,flag,r,ie,grid);
    svect_init(&pv);
    for (i=0; i<3; i++) {
        pv.x += v[i].x * weights[i];
        pv.y += v[i].y * weights[i];
    }
    
    // behavioral velocity adjustment
    if (p->ptype == OYSTER && vel_flag == 1) {
        if (p->isActive == true) {
            soyster_get_behavorial_velocity(p);
            pv.x += p->v_behavoir.x;
            pv.y += p->v_behavoir.y;
        }
    }
    
    // if the p lies on external edge or node, project velocity (relies on proper grid node numbering!)
    int nd1_ID = UNSET_INT, nd2_ID = UNSET_INT;
    SVECT nd1, nd2, vedge;
    if (lnd != UNSET_INT && grid->node_flag[grid->elem2d[ie].nodes[lnd]] == 1) {
        // particle found on external node, find an edge in the direction of the velocity and use it to project
        for (i=0; i<grid->nc_nelems[ grid->elem2d[ie].nodes[lnd] ]; i++) {
            ie2 = grid->nc_elems[ grid->elem2d[ie].nodes[lnd] ][i];
            for (j=0;j<3;j++) { // loop over edges
                if (grid->elem2d[ie2].edge_flag[j] == UNSET_INT) { // found an external edge
                    int nd1_ID = grid->elem2d[ie2].nodes[grid->nd_on_TriEdge[j][0]];
                    int nd2_ID = grid->elem2d[ie2].nodes[grid->nd_on_TriEdge[j][1]];
                    if (grid->elem2d[ie].nodes[lnd] == grid->elem2d[ie2].nodes[grid->nd_on_TriEdge[j][0]]) {
                        nd1 = grid->node[nd1_ID];
                        nd2 = grid->node[nd2_ID];
                    } else {
                        nd1 = grid->node[nd2_ID];
                        nd2 = grid->node[nd1_ID];
                    }
                    vedge.x = nd2.x - nd1.x;
                    vedge.y = nd2.y - nd1.y;
                    if (svect_dotp(vedge,pv)>-1e-12) {
                        // found an external edge in the direction of point vel
                        project_point_velocity_2D(&pv,grid->node[nd1_ID],grid->node[nd2_ID]);
                        proj_flag = 1;
                        break_flag = 1;
                        break;
                    }
                }
            }
            if (break_flag == 1) break;
        }
    } else if (ledge != UNSET_INT && grid->elem2d[ie].edge_flag[ledge] == UNSET_INT) {
        // particle found on external edge
        nd1 = grid->node[grid->elem2d[ie].nodes[grid->nd_on_TriEdge[ledge][0]]];
        nd2 = grid->node[grid->elem2d[ie].nodes[grid->nd_on_TriEdge[ledge][1]]];
        project_point_velocity_2D(&pv,nd1,nd2);
        proj_flag = 1;
    }
//    if (proj_flag == 1) {
//        printf("point velocity after projection: %20.15e %20.15e %20.15e\n",pv.x,pv.y,pv.z);
//        //exit(-1);
//    }
    
    return pv;
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
SVECT get_particle_value_on_elem_vector_3d(SGRID *grid, SPARTICLE *p, SVECT r, int ie, double tL, SVECT *vL, double tR, SVECT *vR, double t, int vel_flag) {
    
    int DEBUG = OFF;
    int i,j,ie2;
    int proj_flag = 0;
    int break_flag = 0;
    SVECT pv; // the return point vector

    if (DEBUG == ON) printf("***** in get_particle_value_on_elem_vector_3d\n");
    
    // interpolate velocity to RK time at triangle vertices
    SVECT v[4]; svect_init(v);
    for (i=0; i<4; i++) {
        v[i].x = interpolate1D(t,tL,tR,vL[ grid->elem3d[ie].nodes[i] ].x,vR[ grid->elem3d[ie].nodes[i] ].x);
        v[i].y = interpolate1D(t,tL,tR,vL[ grid->elem3d[ie].nodes[i] ].y,vR[ grid->elem3d[ie].nodes[i] ].y);
        v[i].z = interpolate1D(t,tL,tR,vL[ grid->elem3d[ie].nodes[i] ].z,vR[ grid->elem3d[ie].nodes[i] ].z);
    }
    
    // interpolate in space from triangle vertices to point
    double weights[4];
    int lnd = UNSET_INT, ledge = UNSET_INT, lface = UNSET_INT, flag = UNSET_INT;
    flag = compute_interpolate_weights_3D_tetrahedron(grid, ie, r.x, r.y, r.z, weights, &lnd, &ledge, &lface);
    error_mssg_point_not_in_elem(__FILE__,__LINE__,flag,r,ie,grid);
    svect_init(&pv);
    for (i=0; i<4; i++) {
        pv.x += v[i].x * weights[i];
        pv.y += v[i].y * weights[i];
        pv.z += v[i].z * weights[i];
    }
    
    if (p->ptype == OYSTER && vel_flag == 1) {
        // check water quality effect on oyster larvae (must happen before velocity adjustment in case oyster dies)
        if (p->isActive == true) {
            soyster_check_wq_exposure(p);
        }
        if (p->isActive == true) {
            // behavioral velocity adjustment
            soyster_get_behavorial_velocity(p);
            pv.x += p->v_behavoir.x;
            pv.y += p->v_behavoir.y;
            pv.z += p->v_behavoir.z;
        } else {
            // just add fall velocity
            soyster_get_fall_velocity(p);
            pv.z += p->w_fall;
        }
    }
    
    proj_flag = project_vel_3D(grid,ie,lnd,ledge,lface,r,&pv);
    
    if (proj_flag == 1) {
        if (DEBUG == ON) printf("point velocity after projection: %20.15e %20.15e %20.15e\n",pv.x,pv.y,pv.z);
        //exit(-1);
    }
    
    if (DEBUG == ON) printf("***** laeving get_particle_value_on_elem_vector_3d with pv: %f %f %f\n",pv.x,pv.y,pv.z);
    
    return pv;
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------

int project_vel_3D(SGRID *grid, int ie, int lnd, int ledge, int lface, SVECT r, SVECT *pv) {
    
    int DEBUG = OFF;
    int i,j,ie2;
    int proj_flag  = 0;
    int break_flag = 0;
    
    // if the r lies on external node, edge or face, project velocity
    int nd1_ID = UNSET_INT, nd2_ID = UNSET_INT, nd3_ID = UNSET_INT;
    SVECT vel_store, ptest;
    if (lnd != UNSET_INT && grid->node_flag[grid->elem3d[ie].nodes[lnd]] == 1) {
        if (DEBUG == ON) printf("particle found on external node, find a face in the direction of the velocity and use it to project\n");
        vel_store = *pv;
        for (i=0; i<grid->nc_nelems[ grid->elem3d[ie].nodes[lnd] ]; i++) {
            ie2 = grid->nc_elems[ grid->elem3d[ie].nodes[lnd] ][i];
            for (j=0;j<4;j++) { // loop over faces
                if (grid->elem3d[ie2].face_flag[j] == UNSET_INT) { // found an external face
                    nd1_ID = grid->elem3d[ie2].nodes[nd_on_fc[j][0]];
                    nd2_ID = grid->elem3d[ie2].nodes[nd_on_fc[j][1]];
                    nd3_ID = grid->elem3d[ie2].nodes[nd_on_fc[j][2]];
                    // project velocity onto face
                    *pv = vel_store;
                    project_point_velocity_3D(pv,grid->node[nd1_ID],grid->node[nd2_ID],grid->node[nd3_ID]);
                    // move forward in time a tiny bit
                    ptest = update_simple(r,1e-4,*pv,1);
                    // now test point for this element
                    int flag = UNSET_INT, lnd = UNSET_INT, ledge = UNSET_INT, lface = UNSET_INT;
                    double weights[4] = {0., 0., 0., 0.};
                    flag = compute_interpolate_weights_3D_tetrahedron(grid, ie2, ptest.x, ptest.y, ptest.z, weights, &lnd, &ledge, &lface);
                    if (flag == 1) { // we have found a face the projected velocity vector crosses
                        proj_flag = 1;
                        break_flag = 1;
                        break;
                    }
                }
            }
            if (break_flag == 1) break;
        }
    } else if (ledge != UNSET_INT && grid->elem3d[ie].edge_flag[ledge] == 1) {
        if (DEBUG == ON) printf("particle found on external edge, find a face in the direction of the velocity and use it to project\n");
        vel_store = *pv;
        for (i=0; i<grid->elem3d[ie].nedge2elem[ledge]; i++) {
            ie2 = grid->elem3d[ie].edge2elem[ledge][i];
            for (j=0;j<4;j++) { // loop over faces
                if (grid->elem3d[ie2].face_flag[j] == UNSET_INT) { // found an external face
                    nd1_ID = grid->elem3d[ie2].nodes[nd_on_fc[j][0]]; //printf("face node 1: %f %f %f\n",grid->node[nd1_ID].x,grid->node[nd1_ID].y,grid->node[nd1_ID].z);
                    nd2_ID = grid->elem3d[ie2].nodes[nd_on_fc[j][1]]; //printf("face node 2: %f %f %f\n",grid->node[nd2_ID].x,grid->node[nd2_ID].y,grid->node[nd2_ID].z);
                    nd3_ID = grid->elem3d[ie2].nodes[nd_on_fc[j][2]]; //printf("face node 3: %f %f %f\n",grid->node[nd3_ID].x,grid->node[nd3_ID].y,grid->node[nd3_ID].z);
                    // project velocity onto face
                    *pv = vel_store;
                    project_point_velocity_3D(pv,grid->node[nd1_ID],grid->node[nd2_ID],grid->node[nd3_ID]);
                    // move forward in time a tiny bit (dt * 1e-4)
                    ptest = update_simple(r,1e-4,*pv,1);
                    // now test point for this element
                    int flag = UNSET_INT, lnd = UNSET_INT, ledge = UNSET_INT, lface = UNSET_INT;
                    double weights[4] = {0., 0., 0., 0.};
                    flag = compute_interpolate_weights_3D_tetrahedron(grid, ie2, ptest.x, ptest.y, ptest.z, weights, &lnd, &ledge, &lface);
                    //printf("flag: %d || ie: %d face nodes: %d %d %d || pvel: %f %f %f ||ptest: %f %f %f || pvel: %f %f %f\n",flag,ie,nd1_ID,nd2_ID,nd3_ID,p_vel.x,p_vel.y,p_vel.z,ptest.x,ptest.y,ptest.z,point_velocity.x,point_velocity.y,point_velocity.z);
                    if (flag == 1) { // our test point with projected velocity lays on this face, use it
                        proj_flag = 1;
                        break_flag = 1;
                        break;
                    }
                }
            }
            if (break_flag == 1) break;
        }
        //printf("original velocity: %f %f %f || point velocity: %f %f %f\n",vel_store.x,vel_store.y,vel_store.z,point_velocity.x,point_velocity.y,point_velocity.z);
        //exit(-1);
    } else if (lface != UNSET_INT && grid->elem3d[ie].face_flag[lface] == UNSET_INT) {
        if (DEBUG == ON) printf("Euler_Update :: particle found on external face\n");
        SVECT nd1, nd2, nd3;
        nd1 = grid->node[grid->elem3d[ie].nodes[nd_on_fc[lface][0]]];
        nd2 = grid->node[grid->elem3d[ie].nodes[nd_on_fc[lface][1]]];
        nd3 = grid->node[grid->elem3d[ie].nodes[nd_on_fc[lface][2]]];
        project_point_velocity_3D(pv,nd1,nd2,nd3);
        proj_flag = 1;
    }
    
    return proj_flag;
}


// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
double get_particle_value_on_elem_scalar_2d(SGRID *grid, SVECT r, int ie, double tL, double *fL, double tR, double *fR, double t) {
    
    int i;
    
    // interpolate velocity to RK time at triangle vertices
    double f[3];
    for (i=0; i<3; i++) {
        f[i] = interpolate1D(t,tL,tR,fL[grid->elem2d[ie].nodes[i]],fR[grid->elem2d[ie].nodes[i]]);
    }
    
    // interpolate in space from triangle vertices to point
    double weights[3];
    int lnd = UNSET_INT, ledge = UNSET_INT, flag = UNSET_INT;
    flag = compute_interpolate_weights_2D_triangle(grid, ie, r.x, r.y, weights, &lnd, &ledge);
    error_mssg_point_not_in_elem(__FILE__,__LINE__,flag,r,ie,grid);
    
    double value = 0.0;
    for (i=0; i<3; i++) {value += f[i] * weights[i];}
    return value;
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
double get_particle_value_on_elem_scalar_3d(SGRID *grid, SVECT r, int ie, double tL, double *fL, double tR, double *fR, double t) {
    
    int i;
    
    // interpolate velocity to RK time at triangle vertices
    double f[4];
    for (i=0; i<4; i++) {
        f[i] = interpolate1D(t,tL,tR,fL[grid->elem3d[ie].nodes[i]],fR[grid->elem3d[ie].nodes[i]]);
    }
    
    // interpolate in space from triangle vertices to point
    double weights[4];
    int lnd = UNSET_INT, ledge = UNSET_INT, lface = UNSET_INT, flag = UNSET_INT;
    flag = compute_interpolate_weights_3D_tetrahedron(grid, ie, r.x, r.y, r.z, weights, &lnd, &ledge, &lface);
    error_mssg_point_not_in_elem(__FILE__,__LINE__,flag,r,ie,grid);
    
    double value = 0.0;
    for (i=0; i<4; i++) {value += f[i] * weights[i];}
    return value;
}
