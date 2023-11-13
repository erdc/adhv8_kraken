#include "global_header.h"

// file prototypes
double get_function_quad(double *lshape, double *f, int nnodes);
//SVECT get_velocity_quad(SQUAD_PT qp, SVECT *v);
SVECT get_velocity2d_quad(SQUAD_PT qp, SVECT *v);
//double get_strong_continuity_quad(SQUAD_PT qp, SVECT *v, int nnodes);
//SVECT2D get_strong_momentum_quad(SQUAD_PT qp, SSW_3D_ELEM *elem, SMAT_SW mat, SVECT elem_vel_avg, SVECT elem_vel_rel_avg, double gravity, double tau_temporal, double dt, int DEBUG);

//*****************************************************************//
//*****************************************************************//
// Galerkin DA-CONT SW3D integration
void fe_sw3_int_continuity_quad(SQUAD *quad, double *elem_rhs, SVECT *v, double dt, int nnodes) {
    int iqp=0, idof=0;
    double constant = 0.;
    SVECT v_quad;
    
    //    // this is the correct way to do it, but results are slightly different than trunk
    //    for (iqp=0; iqp<nqp; iqp++) {
    //        v_quad = get_velocity_quad(qp, iqp, v); // evaluate velocity at quadrature point (can store and send them from ssw_3d_elem.c)
    //        constant = -dt * qp[iqp].djac * qp[iqp].w;
    //        for (idof=0; idof<nnodes; idof++) {
    //            elem_rhs[idof] += constant * svect_dotp(qp[iqp].grad_shp[idof], v_quad);
    //        }
    //    }
    
    
    // this is a linear integration, so use a 1 point quadrature rule
    double lshape[nnodes], weight = 1./6.;
    elem3d_get_tet_local_shape(0.25, 0.25, 0.25, lshape);
    v_quad.x = v[0].x*lshape[0] + v[1].x*lshape[1] + v[2].x*lshape[2]  + v[3].x*lshape[3];
    v_quad.y = v[0].y*lshape[0] + v[1].y*lshape[1] + v[2].y*lshape[2]  + v[3].y*lshape[3];
    v_quad.z = v[0].z*lshape[0] + v[1].z*lshape[1] + v[2].z*lshape[2]  + v[3].z*lshape[3];
    constant = weight * (-dt) * quad[1].pt[0].djac;
    for (idof=0; idof<nnodes; idof++) {
        elem_rhs[idof] += constant * svect_dotp(quad[1].pt[0].grad_shp[idof], v_quad);
    }
    
    
}


////*****************************************************************//
////*****************************************************************//
//// Petrov-Galerkin DA-CONT SW3D integration
//void fe_sw3_int_dacontPG_quad(SQUADRATURE *qp, int nqp, DOF_3 *elem_rhs, SVECT *v, SVECT2D da_vel_avg, double supg_coeff, double dt, int nnodes) {
//    int iqp=0, idof=0;
//    SVECT v_quad;
//
//    for (iqp=0; iqp<nqp; iqp++) {
//
//        strong_cont_quad = get_strong_continuity_quad(qp, iqp, v);
//        for (idof=0; idof<nnodes; idof++) {
//            elem_rhs[idof].c_eq += dt * qp[iqp].djac * qp[iqp].w * strong_cont_quad * supg_coeff *
//                (da_vel_avg.x * qp[iqp].grad_shp[idof].x + da_vel_avg.y * qp[iqp].grad_shp[idof].y);
//        }
//
//    }
//}

//*****************************************************************//
//*****************************************************************//
// Petrov-Galerkin FULL CONT SW3D integration
void fe_sw3_int_contPG_quad(SQUAD *quad, double *elem_rhs, SVECT *v, double supg_coeff, double dt, int nnodes) {
    int iqp=0, idof=0;
    double strong_cont_quad = 0.;
    //    for (iqp=0; iqp<nqp; iqp++) {
    //        strong_cont_quad = get_strong_continuity_quad(qp, iqp, v, nnodes);
    //        for (idof=0; idof<nnodes; idof++) {
    //            elem_rhs[idof] += dt * qp[iqp].djac * qp[iqp].w * strong_cont_quad * supg_coeff * qp[iqp].grad_shp[idof].z;
    //        }
    //    }
    
    // this is linear
    iqp = 0;
    strong_cont_quad = get_strong_continuity_quad(quad[1].pt[iqp], v, nnodes);
    for (idof=0; idof<nnodes; idof++) {
        elem_rhs[idof] += -dt * quad[1].pt[0].djac * (1./6.) * strong_cont_quad * supg_coeff * quad[1].pt[0].grad_shp[idof].z;
    }
}

//*****************************************************************//
//*****************************************************************//
// Convective SW3D Integration :: quadratic :: integral{ u * (grad(phi_i) dot vrel) }dOmega
void fe_sw3_int_convection_quad(SQUAD *quad, DOF_3 * elem_rhs, SVECT *v, SVECT *vrel, double dt, int nnodes) {
    int iqp=0, idof=0;
    SVECT v_quad;
    SVECT vrel_quad;
    
    int order = 2; // quadrature order
    double constant = -dt; // off by a factor of 5
    for (iqp=0; iqp<quad[order].n; iqp++) {
        
        printf("SHOULD NOT BE HERE  !!!!!!!! UNCOMMENT BELOW\n\n\n\n\n");
        //v_quad = get_velocity_quad(quad[order].pt[iqp], v); // evaluate velocity at quadrature point
        //vrel_quad = get_velocity_quad(quad[order].pt[iqp], vrel); // evaluate relative velocity at quadrature point
        
        for (idof=0; idof<nnodes; idof++) {
            elem_rhs[idof].x_eq += constant * quad[order].pt[iqp].djac * quad[order].pt[iqp].w * svect_dotp(quad[order].pt[iqp].grad_shp[idof], svect_scale(vrel_quad,v_quad.x));
            elem_rhs[idof].y_eq += constant * quad[order].pt[iqp].djac * quad[order].pt[iqp].w * svect_dotp(quad[order].pt[iqp].grad_shp[idof], svect_scale(vrel_quad,v_quad.y));
        }
        
        //                // grouping - gives more accurate results for angle_NB, but more nonlinear iterations
        //                // evaluate u * vrel at quadrature point
        //                vrel_quad.x = v[0].x*vrel[0].x*qp[iqp].lshape[0] + v[1].x*vrel[1].x*qp[iqp].lshape[1] + v[2].x*vrel[2].x*qp[iqp].lshape[2] + v[3].x*vrel[3].x*qp[iqp].lshape[3];
        //                vrel_quad.y = v[0].x*vrel[0].y*qp[iqp].lshape[0] + v[1].x*vrel[1].y*qp[iqp].lshape[1] + v[2].x*vrel[2].y*qp[iqp].lshape[2] + v[3].x*vrel[3].y*qp[iqp].lshape[3];
        //                vrel_quad.z = v[0].x*vrel[0].z*qp[iqp].lshape[0] + v[1].x*vrel[1].z*qp[iqp].lshape[1] + v[2].x*vrel[2].z*qp[iqp].lshape[2] + v[3].x*vrel[3].z*qp[iqp].lshape[3];
        //                for (idof=0; idof<nnodes; idof++) {
        //                    elem_rhs[idof].x_eq += constant * qp[iqp].w * svect_dotp(qp[iqp].grad_shp[idof], vrel_quad);
        //                }
        //                // evaluate v * vrel at quadrature point
        //                vrel_quad.x = v[0].y*vrel[0].x*qp[iqp].lshape[0] + v[1].y*vrel[1].x*qp[iqp].lshape[1] + v[2].y*vrel[2].x*qp[iqp].lshape[2] + v[3].y*vrel[3].x*qp[iqp].lshape[3];
        //                vrel_quad.y = v[0].y*vrel[0].y*qp[iqp].lshape[0] + v[1].y*vrel[1].y*qp[iqp].lshape[1] + v[2].y*vrel[2].y*qp[iqp].lshape[2] + v[3].y*vrel[3].y*qp[iqp].lshape[3];
        //                vrel_quad.z = v[0].y*vrel[0].z*qp[iqp].lshape[0] + v[1].y*vrel[1].z*qp[iqp].lshape[1] + v[2].y*vrel[2].z*qp[iqp].lshape[2] + v[3].y*vrel[3].z*qp[iqp].lshape[3];
        //                for (idof=0; idof<nnodes; idof++) {
        //                    elem_rhs[idof].y_eq += constant * qp[iqp].w * svect_dotp(qp[iqp].grad_shp[idof], vrel_quad);
        //                }
        
        
        
    }
}

//*****************************************************************//
//*****************************************************************//
// Temporal SW3D integration :: quadratic :: integral{ u * phi_i }dOmega
void fe_sw3_int_time_quad(SQUAD *quad, DOF_3 *elem_rhs, SVECT *v, double constant, int nnodes) {
    int iqp=0, idof=0;
    SVECT v_quad;
    
    int order = 2; // quadrature order
    for (iqp=0; iqp<quad[order].n; iqp++) {
       
                printf("SHOULD NOT BE HERE  !!!!!!!! UNCOMMENT BELOW\n\n\n\n\n");
     //   v_quad = get_velocity_quad(quad[order].pt[iqp], v); // evaluate velocity at quadrature point
        for (idof=0; idof<nnodes; idof++) {
            elem_rhs[idof].x_eq += constant * quad[order].pt[iqp].djac * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] * v_quad.x;
            elem_rhs[idof].y_eq += constant * quad[order].pt[iqp].djac * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] * v_quad.y;
        }
    }
    
}

//*****************************************************************//
//*****************************************************************//
// Pressure SW3D integration :: quadratic
void fe_sw3_int_prs_quad(SQUAD *quad, DOF_3 *elem_rhs, double *prs, double constant, int nnodes) {
    int iqp=0, idof=0;
    int order = 1; // quadrature order
    
    // constant = -gravity * dt
    
    // calculate avg pressure over element (note: multiplied then divided by jac, so left out)
    // original code uses this as constant over element
    double sum1 = prs[0] + prs[1] + prs[2] + prs[3];
    double sum2 = prs[4] + prs[5] + prs[6] + prs[7] + prs[8] + prs[9];
    double avg_pressure = 0.05 * (-sum1 + 4. * sum2);
    
    //    double prs_quad = 0.;
    //    for (iqp=0; iqp<nqp; iqp++) {
    //
    //        // evaluate pressure at quad point
    //        //prs_quad = 0.;
    //        //for (idof=0; idof<10; idof++) {
    //        //    prs_quad += prs[idof]*qp[iqp].lshape_quad[idof];
    //        //}
    //        prs_quad = avg_pressure;
    //
    //        for (idof=0; idof<nnodes; idof++) {
    //            elem_rhs[idof].x_eq += constant * qp[iqp].djac * qp[iqp].w * qp[iqp].grad_shp[idof].x * prs_quad;
    //            elem_rhs[idof].y_eq += constant * qp[iqp].djac * qp[iqp].w * qp[iqp].grad_shp[idof].y * prs_quad;
    //        }
    //    }
    
    
    // if we're assuming pressure is constant and equal to average on element, then this greatly simplifies too ...
    for (idof=0; idof<nnodes; idof++) {
        elem_rhs[idof].x_eq += constant * avg_pressure * quad[order].pt[iqp].djac * (1./6.) * quad[order].pt[iqp].grad_shp[idof].x; // here 1/6 is quadrature constant weight
        elem_rhs[idof].y_eq += constant * avg_pressure * quad[order].pt[iqp].djac * (1./6.) * quad[order].pt[iqp].grad_shp[idof].y; // here 1/6 is quadrature constant weight
    }
}

//*****************************************************************//
//*****************************************************************//
// Diffusion SW3D integration :: quadratic :: integral{ u * (grad(phi_i) dot vrel) }dOmega
void fe_sw3_int_diff_quad(SQUAD *quad, DOF_3 *elem_rhs, SVECT *v, SMAT_SW *mat, SVECT elem_vel_avg, double elem_dep_avg, double *elem_density, double viscosity, double g, double dt, int nnodes) {
    
    int iqp=0, idof=0;
    double ddensity_dz_quad;
    SVECT grad_u_quad, grad_v_quad;
    SVECT dy_vis_x;   svect_init(&dy_vis_x);
    SVECT dy_vis_y;   svect_init(&dy_vis_y);
    
    /* Compute requirements for the turbulence models */
    /* Mode 0 = Smagorinkski */
    /* Mode 1 = MY 2 */
    /* Mode 2 = MY 2.5 */
    /* Mode 3 = MY 3 */
    /* Mode 4 = K-e  GSavant */
    /* Modes are open ended, we should be able to include other models using the transport capabilities and Operator-Splitting */
    double evxz = 0.; double evyz = 0.; double tur_viscosity = 0.;  double dist_above_bed = 0.;
    double v_mag_actual = 0.; double avg_den = 0.; double avg_factor = 0.;
    
    STENSOR ev;
    stensor_copy(&ev, mat->ev);  /* copy the eddy viscosity to the local element */
    
    int order = 2; // quadrature order
    
    for (iqp=0; iqp<quad[order].n; iqp++) {
        
        // evaluate velocity gradients at quad point
        svect_init(&grad_u_quad); svect_init(&grad_v_quad);
        ddensity_dz_quad = 0.;
        for (idof=0; idof<nnodes; idof++) {
            grad_u_quad.x += v[idof].x * quad[order].pt[iqp].grad_shp[idof].x;
            grad_u_quad.y += v[idof].x * quad[order].pt[iqp].grad_shp[idof].y;
            grad_u_quad.z += v[idof].x * quad[order].pt[iqp].grad_shp[idof].z;
            grad_v_quad.x += v[idof].y * quad[order].pt[iqp].grad_shp[idof].x;
            grad_v_quad.y += v[idof].y * quad[order].pt[iqp].grad_shp[idof].y;
            grad_v_quad.z += v[idof].y * quad[order].pt[iqp].grad_shp[idof].z;
            ddensity_dz_quad += elem_density[idof] * quad[order].pt[iqp].grad_shp[idof].z;
        }
        
        if (mat->turbulence_model_xy == 0) {
            tur_viscosity = tur_smag(grad_u_quad.x, grad_u_quad.y, grad_v_quad.x, grad_v_quad.y, quad[order].pt[iqp].djac, mat->smag_coeff);
        }
        if (mat->turbulence_model_xy == 1) {
            v_mag_actual = svect_mag_safe(elem_vel_avg);
            tur_viscosity = tur_ws(v_mag_actual, elem_dep_avg, mat->smag_coeff);
        }
        if (mat->turbulence_model_z == 1) {
            double avg_den=0., avg_factor = 0.;
            avg_den = sarray_avg_dbl(elem_density, nnodes);
            dist_above_bed = elem_dep_avg;
            evxz = tur_MY_2(g, elem_dep_avg, dist_above_bed, grad_u_quad.z, grad_v_quad.z, ddensity_dz_quad, avg_den, mat->supression_func, 0);
            evyz = tur_MY_2(g, elem_dep_avg, dist_above_bed, grad_u_quad.z, grad_v_quad.z, ddensity_dz_quad, avg_den, mat->supression_func, 1);
            evxz /= avg_den;
            evyz /= avg_den;
        }
        if (mat->turbulence_model_z <= 0) {
            evxz = 0.;
            evyz = 0.;
        }
        
        /* Store elemental eddy viscosity -- HOW DO I DO THIS WITH QUADRATURE??? */
        //grid->hyd_eddy[elem3d->id] = evxz + (ev.xz + ev.zz + ev.yz)/3.;
        //printf("Eddy Visc %lf for node %d\n", grid->hyd_eddy[elem3d->id], elem3d->id+1);
        
        double multiply = 1.;
        if (ddensity_dz_quad > 0) multiply = 10.;
        
        dy_vis_x.x = multiply * (ev.xx + viscosity + tur_viscosity) * grad_u_quad.x;
        dy_vis_x.y = multiply * (ev.xy + viscosity + tur_viscosity) * grad_u_quad.y;
        dy_vis_x.z = multiply * (ev.xz + viscosity + evxz) * grad_u_quad.z;
        dy_vis_y.x = multiply * (ev.xy + viscosity + tur_viscosity) * grad_v_quad.x;
        dy_vis_y.y = multiply * (ev.yy + viscosity + tur_viscosity) * grad_v_quad.y;
        dy_vis_y.z = multiply * (ev.yz + viscosity + evyz) * grad_v_quad.z;
        
        for (idof=0; idof<nnodes; idof++) {
            elem_rhs[idof].x_eq += dt * quad[order].pt[iqp].djac * quad[order].pt[iqp].w * svect_dotp(quad[order].pt[iqp].grad_shp[idof], dy_vis_x);
            elem_rhs[idof].y_eq += dt * quad[order].pt[iqp].djac * quad[order].pt[iqp].w * svect_dotp(quad[order].pt[iqp].grad_shp[idof], dy_vis_y);
            
        }
        
    }
}

//*****************************************************************//
//*****************************************************************//
// Coriolis SW3D integration :: linear
void fe_sw3_int_coriolis_quad(SQUAD *quad, DOF_3 *elem_rhs, SVECT *v, double coriolis, double dt, int nnodes) {
    int iqp=0, idof=0;
    int order = 1; // quadrature order
    
    SVECT v_quad;
    double angular_speed = get_coriolis_angular_speed(coriolis);
    //printf("angular speed: %20.10f coriolis: %20.10f \n",angular_speed, coriolis);
    for (iqp=0; iqp<quad[order].n; iqp++) {
        for (idof=0; idof<nnodes; idof++) {
    
                    printf("SHOULD NOT BE HERE  !!!!!!!! UNCOMMENT BELOW\n\n\n\n\n");
          //  v_quad = get_velocity_quad(quad[order].pt[iqp], v); // evaluate velocity at quadrature point (can store and send them from ssw_3d_elem.c)
            elem_rhs[idof].x_eq += -dt * quad[order].pt[iqp].djac * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] * angular_speed * v_quad.y;
            elem_rhs[idof].y_eq += +dt * quad[order].pt[iqp].djac * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] * angular_speed * v_quad.x;
        }
        
    }
}

//*****************************************************************//
//*****************************************************************//
// Petrov-Galerkin SW3D momentum and cont. contribution to momentum integration :: linear
void fe_sw3_int_supg_quad(SQUAD *quad, DOF_3 *elem_rhs, SSW_3D_ELEM *elem, SMAT_SW mat, double gravity, double tau_temporal,
                          double supg_momC_coef, double supg_momXY_coef, double supg_cont_coef, double dt, int DEBUG) {
    
    int i, iqp=0, idof=0;
    double strong_con_quad = 0., constant = 0.;
    SVECT2D strong_mom_quad; svect2d_init(&strong_mom_quad);
    int nnodes3d = elem->elem3d->nnodes;
    double supg_cont_store[nnodes3d]; sarray_init_dbl(supg_cont_store, nnodes3d);
    
    // get grid velocities
    double dpl[nnodes3d];
    double dpl_old[nnodes3d];
    elem_get_tposition(dpl_old, elem->old_displacement, elem->older_displacement, tau_temporal, dt, dt, nnodes3d); // get gs at t(i+1/2)
    elem_get_tposition(dpl, elem->displacement, elem->old_displacement, tau_temporal, dt, dt, nnodes3d); // get gs at t(i+3/2)
    SVECT grid_vel[nnodes3d]; svect_init_array(grid_vel, nnodes3d);
    for (i = 0; i < nnodes3d; i++) {
        grid_vel[i].z = (dpl[i] - dpl_old[i])/dt;
        //grid_vel[i].z = dpl[i]/dt - dpl_old[i]/dt;
    }
    
    // get relative velocities
    SVECT rel_vel[nnodes3d];  svect_init_array(rel_vel, nnodes3d);
    svect_subtract_array2(rel_vel, elem->vel, grid_vel, nnodes3d); // relative velocity
    
    // elementally averaged hydro and relative velocity
    SVECT elem_vel_avg;  svect_init(&elem_vel_avg);
    SVECT elem_vel_rel_avg; svect_init(&elem_vel_rel_avg);
    elem_vel_avg = svect_average_array(elem->vel, nnodes3d);
    elem_vel_rel_avg = svect_average_array(rel_vel, nnodes3d);
    
    // elementally averaged depth-average velocity
    SVECT2D da_vel_avg; svect2d_init(&da_vel_avg);
    da_vel_avg = svect2d_average_array(elem->depth_avg_vel, elem->elem2d_sur->nnodes);
    
    DOF_3 rhs_mom_supg[nnodes3d]; dof3_init_array(rhs_mom_supg, nnodes3d);
    DOF_3 rhs_dacont_supg[nnodes3d]; dof3_init_array(rhs_dacont_supg, nnodes3d);
    
    //    for (iqp=0; iqp<nqp; iqp++) {
    //
    //        strong_mom_quad = get_strong_momentum_quad(qp, iqp, elem, elem_vel_avg, elem_vel_rel_avg, gravity, tau_temporal, dt, nnodes3d, DEBUG);
    //        strong_con_quad = get_strong_continuity_quad(qp, iqp, elem->vel, nnodes3d);
    //        constant = dt * qp[iqp].djac * qp[iqp].w;
    //
    //#ifdef _DEBUG
    //        if (DEBUG) {
    //            printf("strong_mom_quad: (%20.10f,%20.10f) strong_con_quad: %20.10f dt*djac: %20.10f\n",strong_mom_quad.x,strong_mom_quad.y, strong_con_quad, dt * qp[iqp].djac);
    //            //printf("elem_vel_rel_avg: %20.10f %20.10f\n",elem_vel_rel_avg.x,elem_vel_rel_avg.y);
    //        }
    //#endif
    //
    //        for (idof=0; idof<nnodes3d; idof++) {
    //
    //            // depth-averaged continuity supg
    //            rhs_dacont_supg[idof].c_eq += constant* supg_cont_coef * (da_vel_avg.x * qp[iqp].grad_shp[idof].x + da_vel_avg.y * qp[iqp].grad_shp[idof].y)* strong_con_quad;
    //
    //            // depth-averaged continuity supg contribution to momentum   (THIS GIVES NANS FOR ANGLE_SOURCE - I think because momC is big)
    //            rhs_dacont_supg[idof].x_eq += constant * supg_momC_coef * qp[iqp].grad_shp[idof].x * strong_con_quad;
    //            rhs_dacont_supg[idof].y_eq += constant * supg_momC_coef * qp[iqp].grad_shp[idof].y * strong_con_quad;
    //
    //            // momentum supg
    //            rhs_mom_supg[idof].x_eq += constant * supg_momXY_coef * svect_dotp(elem_vel_rel_avg, qp[iqp].grad_shp[idof]) * strong_mom_quad.x;
    //            rhs_mom_supg[idof].y_eq += constant * supg_momXY_coef * svect_dotp(elem_vel_rel_avg, qp[iqp].grad_shp[idof]) * strong_mom_quad.y;
    //
    //            // momentum supg contribution to depth-averaged continuity
    //            rhs_mom_supg[idof].c_eq += constant * cont_coef * (qp[iqp].grad_shp[idof].x * strong_mom_quad.x + qp[iqp].grad_shp[idof].y * strong_mom_quad.y);
    //
    //
    //        }
    //
    //    }
    
    // These are all constants, by design
    int order = 1; // quadrature order
    
    strong_mom_quad = get_strong_momentum_quad(quad[order].pt[iqp], elem, mat, elem_vel_avg, elem_vel_rel_avg, gravity, tau_temporal, dt, DEBUG);
    strong_con_quad = get_strong_continuity_quad(quad[order].pt[iqp], elem->vel, nnodes3d);
    constant = dt * elem->elem3d->djac; //qp[iqp].djac * (1./6.); // the 1/6 is order=1 quad weight
    for (idof=0; idof<nnodes3d; idof++) {
        
        // depth-averaged continuity supg
        rhs_dacont_supg[idof].c_eq += constant* supg_cont_coef * (da_vel_avg.x * quad[order].pt[iqp].grad_shp[idof].x + da_vel_avg.y * quad[order].pt[iqp].grad_shp[idof].y)* strong_con_quad;
        
        // depth-averaged continuity supg contribution to momentum   (THIS GIVES NANS FOR ANGLE_SOURCE - I think because momC is big)
        rhs_dacont_supg[idof].x_eq += constant * supg_momC_coef * quad[order].pt[iqp].grad_shp[idof].x * strong_con_quad;
        rhs_dacont_supg[idof].y_eq += constant * supg_momC_coef * quad[order].pt[iqp].grad_shp[idof].y * strong_con_quad;
        
        // momentum supg
        rhs_mom_supg[idof].x_eq += constant * supg_momXY_coef * svect_dotp(elem_vel_rel_avg, quad[order].pt[iqp].grad_shp[idof]) * strong_mom_quad.x;
        rhs_mom_supg[idof].y_eq += constant * supg_momXY_coef * svect_dotp(elem_vel_rel_avg, quad[order].pt[iqp].grad_shp[idof]) * strong_mom_quad.y;
        
        // momentum supg contribution to depth-averaged continuity
        rhs_mom_supg[idof].c_eq += constant * supg_cont_coef * (quad[order].pt[iqp].grad_shp[idof].x * strong_mom_quad.x + quad[order].pt[iqp].grad_shp[idof].y * strong_mom_quad.y);
        
    }
    
    for (idof=0; idof<nnodes3d; idof++) {
        supg_cont_store[idof] = rhs_dacont_supg[idof].c_eq + rhs_mom_supg[idof].c_eq;
        elem_rhs[idof].c_eq += supg_cont_store[idof];
        elem_rhs[idof].x_eq += rhs_dacont_supg[idof].x_eq + rhs_mom_supg[idof].x_eq;
        elem_rhs[idof].y_eq += rhs_dacont_supg[idof].y_eq + rhs_mom_supg[idof].y_eq;
    }
    
#ifdef _DEBUG
    if (DEBUG) {
        printf("strong_mom_quad: (%20.10f,%20.10f) strong_con_quad: %20.10f dt*djac: %20.10f\n",strong_mom_quad.x,strong_mom_quad.y, strong_con_quad, dt * quad[1].pt[iqp].djac);
        //printf("elem_vel_rel_avg: %20.10f %20.10f\n",elem_vel_rel_avg.x,elem_vel_rel_avg.y);
        printf("supg_momC_coef: %20.10f supg_momXY_coef: %20.10f supg_cont_coef: %20.10f \n",supg_momC_coef,supg_momXY_coef,supg_cont_coef);
        rhs_3dof("depth-averaged continuity supg", nnodes3d, elem->elem3d->id, elem->elem3d->nodes, rhs_dacont_supg);
        rhs_3dof("momentum supg", nnodes3d, elem->elem3d->id, elem->elem3d->nodes, rhs_mom_supg);
    }
#endif
    
}

//*****************************************************************//
//*****************************************************************//
// Return strong momentum equation evaluated at Gauss Point
SVECT2D get_strong_momentum_quad(SQUAD_PT qp, SSW_3D_ELEM *elem, SMAT_SW mat, SVECT elem_vel_avg, SVECT elem_vel_rel_avg, double gravity, double tau_temporal, double dt, int DEBUG) {
    
    // note :: elemental averages are good, since they are constant in element for any element
    // note :: gradients are not, since for prisms, they will not be constant
    
    int i, idof;
    int nnodes3d = elem->elem3d->nnodes;
    
    // gradients at quadrature points (constant for tet, linear for bilinear prism)
    SVECT grad_u_quad; svect_init(&grad_u_quad);
    SVECT grad_v_quad; svect_init(&grad_v_quad);
    for (idof=0; idof<nnodes3d; idof++) {
        grad_u_quad.x += elem->vel[idof].x*qp.grad_shp[idof].x;
        grad_u_quad.y += elem->vel[idof].x*qp.grad_shp[idof].y;
        grad_u_quad.z += elem->vel[idof].x*qp.grad_shp[idof].z;
        grad_v_quad.x += elem->vel[idof].y*qp.grad_shp[idof].x;
        grad_v_quad.y += elem->vel[idof].y*qp.grad_shp[idof].y;
        grad_v_quad.z += elem->vel[idof].y*qp.grad_shp[idof].z;
    }
    
    // evaluate pressure at quad point (note, this really should loop over all 10 basis functions here ...)
    SVECT grad_p_quad; svect_init(&grad_p_quad);
    //if (elem->PRS_FLAG == ON) {
    for (idof=0; idof<nnodes3d; idof++) {
        grad_p_quad.x += elem->pressure[idof]*qp.grad_shp[idof].x;
        grad_p_quad.y += elem->pressure[idof]*qp.grad_shp[idof].y;
    }
    //}
    
    // calculate velocities at time t+3/2*dt and t+1/2*dt
    SVECT elem_vel_new[nnodes3d]; svect_init_array(elem_vel_new, nnodes3d);
    SVECT elem_vel_old[nnodes3d]; svect_init_array(elem_vel_old, nnodes3d);
    elem_get_tposition_vect(elem_vel_new, elem->vel, elem->old_vel, tau_temporal, dt, dt, nnodes3d);
    elem_get_tposition_vect(elem_vel_old, elem->old_vel, elem->older_vel, tau_temporal, dt, dt, nnodes3d);
    SVECT avg_new_vel; svect_init(&avg_new_vel);
    SVECT avg_old_vel; svect_init(&avg_old_vel);
    avg_new_vel = svect_average_array(elem_vel_new, nnodes3d);
    avg_old_vel = svect_average_array(elem_vel_old, nnodes3d);
    
    double angular_speed = 0.;
    //if (elem->CORIOLIS_FLAG == ON) {
    angular_speed = get_coriolis_angular_speed(mat.coriolis);
    //}
    
#ifdef _DEBUG
    if (DEBUG) {
        //printf("elem_vel_rel_avg: %20.10f %20.10f %20.10f dt: %20.10f avg_new_vel.x: %20.10f avg_old_vel.x: %20.10f elem_vel_avg.y: %20.10f angular_speed: %20.10f grad_p_quad.x: %20.10f dot: %20.10f grad_u %20.10f %20.10f %20.10f shape.x: %20.10f\n",elem_vel_rel_avg.x, elem_vel_rel_avg.y, elem_vel_rel_avg.z, dt,avg_new_vel.x,avg_old_vel.x,elem_vel_avg.y,angular_speed,grad_p_quad.x,svect_dotp(elem_vel_rel_avg, grad_u_quad),grad_u_quad.x, grad_u_quad.y, grad_u_quad.z,qp[iqp].grad_shp[0].x);
        
        //svect_printScreen_array("elem_vel", elem->vel, "elem_vel", nnodes, __LINE__, __FILE__);
        //svect_printScreen(elem_vel_rel_avg, "elem_vel_rel_avg");
        //svect_printScreen(grad_u_quad, "grad_u_quad");
        //svect_printScreen_array("grad_shape", qp[iqp].grad_shp, "grad_shp", nnodes, __LINE__, __FILE__);
    }
#endif
    
    SVECT2D momentum; svect2d_init(&momentum);
    momentum.x = (avg_new_vel.x - avg_old_vel.x)/dt + (svect_dotp(elem_vel_rel_avg, grad_u_quad) - elem_vel_avg.y * angular_speed + gravity * grad_p_quad.x);
    momentum.y = (avg_new_vel.y - avg_old_vel.y)/dt + (svect_dotp(elem_vel_rel_avg, grad_v_quad) + elem_vel_avg.x * angular_speed + gravity * grad_p_quad.y);
    
    //    // cjt :: this way is MUCH more stable ....
    //    momentum.x = avg_new_vel.x/dt - avg_old_vel.x/dt + svect_dotp(elem_vel_rel_avg, grad_u_quad) - elem_vel_avg.y * angular_speed;
    //    momentum.y = avg_new_vel.y/dt - avg_old_vel.y/dt + svect_dotp(elem_vel_rel_avg, grad_v_quad) + elem_vel_avg.x * angular_speed;
    //
    //    //if (elem->PRS_FLAG) {
    //        momentum.x += gravity * grad_p_quad.x;
    //        momentum.y += gravity * grad_p_quad.y;
    //    //}
    
    
    return (momentum);
}


//*****************************************************************//
//*****************************************************************//
// Return strong continuity equation evaluated at Gauss Point
double get_strong_continuity_quad(SQUAD_PT qp, SVECT *v, int nnodes) {
    
    double du_dx_quad = v[0].x*qp.grad_shp[0].x + v[1].x*qp.grad_shp[1].x + v[2].x*qp.grad_shp[2].x + v[3].x*qp.grad_shp[3].x;
    double dv_dy_quad = v[0].y*qp.grad_shp[0].y + v[1].y*qp.grad_shp[1].y + v[2].y*qp.grad_shp[2].y + v[3].y*qp.grad_shp[3].y;
    double dw_dz_quad = v[0].z*qp.grad_shp[0].z + v[1].z*qp.grad_shp[1].z + v[2].z*qp.grad_shp[2].z + v[3].z*qp.grad_shp[3].z;
    
    // evaluate continuity eq. at gauss point
    double cont_quad = du_dx_quad + dv_dy_quad + dw_dz_quad;
    
    return cont_quad;
}

////*****************************************************************//
////*****************************************************************//
//// Return the 3d velocity vector at Gauss Point
//SVECT get_velocity_quad(SQUAD_PT qp, SVECT *v) {
//    SVECT v_quad; svect_init(&v_quad);
//    v_quad.x = v[0].x*qp.lshape[0] + v[1].x*qp.lshape[1] + v[2].x*qp.lshape[2] + v[3].x*qp.lshape[3];
//    v_quad.y = v[0].y*qp.lshape[0] + v[1].y*qp.lshape[1] + v[2].y*qp.lshape[2] + v[3].y*qp.lshape[3];
//    v_quad.z = v[0].z*qp.lshape[0] + v[1].z*qp.lshape[1] + v[2].z*qp.lshape[2] + v[3].z*qp.lshape[3];
//    return v_quad;
//}

// Return 3d velocity on a 2d element
SVECT get_velocity2d_quad(SQUAD_PT qp, SVECT *v) {
    SVECT v_quad; svect_init(&v_quad);
    v_quad.x = v[0].x*qp.lshape[0] + v[1].x*qp.lshape[1] + v[2].x*qp.lshape[2];
    v_quad.y = v[0].y*qp.lshape[0] + v[1].y*qp.lshape[1] + v[2].y*qp.lshape[2];
    v_quad.z = v[0].z*qp.lshape[0] + v[1].z*qp.lshape[1] + v[2].z*qp.lshape[2];
    return v_quad;
}

double get_function_quad(double *lshape, double *f, int nnodes) {
    int inode = 0.;
    double fout = 0.;
    for (inode=0; inode<nnodes; inode++) {
        fout += f[inode]*lshape[inode];
        

        if(solv_isnan(fout)) {
            printf("f: %20.10f lshape: %20.10f\n",f[inode], lshape[inode]);
        }
    }
    return fout;
}

//*****************************************************************//
//*****************************************************************//
// Return the djac at quadrature points
void get_djac_quad_3d(SQUAD *quad, double elem_djac) {
    int iorder, iqp;
    for (iorder=1; iorder<=MAX_QUAD_ORDER; iorder++) {
        for (iqp=0; iqp<quad[iorder].n; iqp++) {
            quad[iorder].pt[iqp].djac = elem_djac * 6.; // djac in elem is really area
        }
    }
}

void get_djac_quad_2d(SQUAD *quad, double elem_djac2d,  double elem_djac2d_3d, double elem_djac3d_fixed) {
    int iorder, iqp;
    for (iorder=1; iorder<=MAX_QUAD_ORDER; iorder++) {
        for (iqp=0; iqp<quad[iorder].n; iqp++) {
            quad[iorder].pt[iqp].djac2d = elem_djac2d * 2.;
            quad[iorder].pt[iqp].djac2d_3d = elem_djac2d_3d * 2.;
            quad[iorder].pt[iqp].djac3d_fixed = elem_djac3d_fixed * 2.;
        }
    }
}

////*****************************************************************//
////*****************************************************************//
//// Return the djac and grad(phi_i) at quadrature points :: later this will need to be adapted for prisms
//void get_elem_djac_grad_quad_3d(SQUAD *quad, SELEM_3D *elem3d, double elem_djac[][MAX_QUAD_POINTS], SVECT elem_grad_shp[][MAX_QUAD_POINTS][MAX_NNODES_ON_ELEM3D]) {
//    int iorder, iqp, idof;
//    for (iorder=1; iorder<=MAX_QUAD_ORDER; iorder++) {
//        for (iqp=0; iqp<quad[iorder].n; iqp++) {
//            elem_djac[iorder][iqp] = elem3d->djac * 6.; // djac in elem is really area
//            for (idof=0; idof<elem3d->nnodes; idof++) {
//                elem_grad_shp[iorder][iqp][idof] = elem3d->grad_shp[idof];
//            }
//        }
//    }
//}

//*****************************************************************//
//*****************************************************************//
// Return the 3d velocity vector at Gauss Point
void get_djac_grad_quad_3d(SQUAD *quad, double elem_djac, SVECT *elem_grad_shp, int nnodes3d) {
    int iorder, iqp;
    for (iorder=1; iorder<=MAX_QUAD_ORDER; iorder++) {
        for (iqp=0; iqp<quad[iorder].n; iqp++) {
            quad[iorder].pt[iqp].djac = elem_djac * 6.; // djac in elem is really area
            svect_copy_array(quad[iorder].pt[iqp].grad_shp, elem_grad_shp, nnodes3d);
            
            //        // below is not working for some reason ...
            //        elem3d_get_tet_local_shape_gradients(qp[i].xhat, qp[i].yhat, qp[i].zhat, qp[i].grad_shp);
            //
            //        for (j=0; j<4; j++) {
            //            qp[i].grad_shp[j].x /= qp[i].djac;
            //            qp[i].grad_shp[j].y /= qp[i].djac;
            //            qp[i].grad_shp[j].z /= qp[i].djac;
            //
            //            printf("new grad shp: %20.10f \n",qp[0].grad_shp[j].x);
            //            printf("old grad shp: %20.10f \n",elem_grad_shp[j].x);
            //        }
            //        exit(-1);
        }
    }
}

//*****************************************************************//
//*****************************************************************//
// Convective SW3D Integration :: 2D NB integration :: normal_discharge (Natural BC), normal_vel (Dirichlet BC)
void fe_sw3_2d_int_NB(SQUAD *quad, DOF_3 *elem_rhs, SSW_3D_ELEM *elem, double normal_discharge, double *normal_vel, double dt, int nnodes_quad, int *mask, double gravity) {
    int iqp=0, idof=0, order=1;
    double temp = 0., drag = 0., elem_prs_quad=0., elem_v_quad_mag = 0., elem_nrml_vel_quad = 0., term1 = 0., term2 = 0.;
    SVECT elem_v_quad; svect_init(&elem_v_quad);
    
    int nnodes = elem->elem2d->nnodes;
    
    //    int iorder;
    //    for (iorder=1; iorder<4; iorder++) {
    //        for (iqp=0; iqp<quad[iorder].n; iqp++) {
    //            printf("iorder: %d  iqp: %d  djac2d: %20.10f djac2d_3d: %20.10f djac3d_fixed: %20.10f\n",
    //                   iorder, iqp, quad[iorder].pt[iqp].djac2d,quad[iorder].pt[iqp].djac2d_3d,quad[iorder].pt[iqp].djac3d_fixed);
    //        }
    //    }
    //svect_printScreen_array("velocity",elem->vel, "velocity",3,__LINE__,__FILE__);
    //printf("discharge: %20.10f\n",normal_discharge);
    //sarray_printScreen_dbl(normal_vel, nnodes, "normal_vel");
    
    for (idof=0; idof<nnodes; idof++) {
        
        // depth-averaged continuity additions
        if (mask[idof * nnodes] == YES) {
            // NON-projected BC imposed, quadratic integrand
            order = 2;
            for (iqp=0; iqp<quad[order].n; iqp++) {
                elem_nrml_vel_quad = get_function_quad(quad[order].pt[iqp].lshape, normal_vel, nnodes); // evaluate velocity at quadrature point
                elem_rhs[idof].c_eq += dt * quad[order].pt[iqp].djac3d_fixed * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] * elem_nrml_vel_quad; // quadratic
            }
            //printf("1: elem_rhs[idof].c_eq: %20.10f\n",elem_rhs[idof].c_eq);
        } else {
            // Natural normal flow BC imposed, linear integrand
            order=1; iqp = 0;
            elem_rhs[idof].c_eq += dt * quad[order].pt[iqp].djac3d_fixed * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] * normal_discharge; // linear
            //printf("2: elem_rhs[idof].c_eq: %20.10f\n",elem_rhs[idof].c_eq);
        }
        
        // momentum additions (in projected orthogonal coordinates)
        if (mask[idof * nnodes] != 2) {
            order=3;
            for (iqp=0; iqp<quad[order].n; iqp++) {
                elem_prs_quad = get_function_quad(quad[order].pt[iqp].lshape_quad, elem->pressure, nnodes_quad); // this needs to be quadratric basis!!
                elem_v_quad = get_velocity2d_quad(quad[order].pt[iqp], elem->vel);
                elem_v_quad_mag = sqrt(elem_v_quad.x*elem_v_quad.x + elem_v_quad.y*elem_v_quad.y);
                temp = dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof]; // linear
                
                elem_rhs[idof].x_eq += temp * elem->elem2d->nrml.x * gravity * elem_prs_quad; // cubic
                elem_rhs[idof].y_eq += temp * elem->elem2d->nrml.y * gravity * elem_prs_quad; // cubic
                elem_rhs[idof].x_eq += temp * (1./2.) * drag * elem_v_quad_mag * elem_v_quad.x; // cubic
                elem_rhs[idof].y_eq += temp * (1./2.) * drag * elem_v_quad_mag * elem_v_quad.y; // cubic
            }
            //printf("1: elem_rhs[idof].x_eq: %20.10f\n",elem_rhs[idof].x_eq);
        } else {
            order=3;
            for (iqp=0; iqp<quad[order].n; iqp++) {

                elem_prs_quad = get_function_quad(quad[order].pt[iqp].lshape_quad, elem->pressure, nnodes_quad); // quadratic pressure at quadrature pt
                //elem_prs_quad = get_function_quad(quad[order].pt[iqp].lshape, elem->pressure, nnodes); // linear pressure at quadrature pt
                
                elem_nrml_vel_quad = get_function_quad(quad[order].pt[iqp].lshape, normal_vel, nnodes); // evaluate velocity at quadrature pt
                elem_v_quad = get_velocity2d_quad(quad[order].pt[iqp], elem->vel);
                elem_v_quad_mag = sqrt(elem_v_quad.x*elem_v_quad.x + elem_v_quad.y*elem_v_quad.y);
                term1 = elem->elem2d->nrml.x * elem->tanvec[idof].x + elem->elem2d->nrml.y * elem->tanvec[idof].y; // really this should be linearly integrated
                term2 = elem_v_quad.x * elem->tanvec[idof].x + elem_v_quad.y * elem->tanvec[idof].y; // really this should be linearly integrated
                
                temp = dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] ; // linear
                elem_rhs[idof].x_eq += temp * elem_nrml_vel_quad;  // quadratic
                elem_rhs[idof].x_eq -= temp * normal_discharge;  // linear
                elem_rhs[idof].y_eq += temp * gravity * elem_prs_quad * term1;  // cubic // issue!!!!
                elem_rhs[idof].y_eq += temp * (1./2.) * drag * elem_v_quad_mag * term2; // cubic
            }
        }
        
        
        //printf("nnodes_quad: %d :: elem_prs_quad: %20.10f \n",nnodes_quad,elem_prs_quad);
        //printf("nrml: %20.10f %20.10f :: elem_vquad: %20.10f %20.10f \n",elem->elem2d->nrml.x,elem->elem2d->nrml.y,elem_v_quad.x,elem_v_quad.y);
        //printf("2: elem_rhs[idof].y_eq: %20.10f term1: %20.10f term2: %20.10f \n",elem_rhs[idof].y_eq, term1, term2);
        //printf("2: elem_rhs[idof].x_eq: %20.10f\n",elem_rhs[idof].x_eq);
        //Is_DoubleArray_Inf_or_NaN(&elem_rhs[idof].y_eq, 1, __FILE__, __LINE__);
    }
    
    
//    // testing different order of operations for pressure term
//    order = 3;
//    for (iqp=0; iqp<quad[order].n; iqp++) {
//        for (idof=0; idof<nnodes; idof++) {
//            elem_prs_quad = get_function_quad(quad[order].pt[iqp].lshape_quad, elem->pressure, nnodes_quad); // this needs to be quadratric basis!!
//            //elem_prs_quad = get_function_quad(quad[order].pt[iqp].lshape, elem->pressure, nnodes);
//            temp = dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] ; // linear
//            term1 = elem->elem2d->nrml.x * elem->tanvec[idof].x + elem->elem2d->nrml.y * elem->tanvec[idof].y; // really this should be linearly integrated
//            elem_rhs[idof].y_eq += temp * gravity * elem_prs_quad * term1;  // cubic // issue!!!!
//        }
//        
//    }
    
    
    
}

//*****************************************************************//
//*****************************************************************//
// Convective WVEL SW3D Integration :: 2D NB integration :: normal_discharge (Natural BC), normal_vel (Dirichlet BC)
void fe_sw3_2d_int_WVEL_NB(SQUAD *quad, double *elem_rhs, SSW_3D_ELEM *elem, double normal_discharge, double *normal_vel, double dt, int *mask) {
    int iqp=0, idof=0, order=1;
    double temp = 0., drag = 0., elem_v_quad_mag = 0., elem_nrml_vel_quad = 0.;
    SVECT elem_v_quad; svect_init(&elem_v_quad);
    
    int nnodes = elem->elem2d->nnodes;
    
    for (idof=0; idof<nnodes; idof++) {
        
        // depth-averaged continuity additions
        if (mask[idof * nnodes] == YES) {
            // NON-projected BC imposed, quadratic integrand
            order = 2;
            for (iqp=0; iqp<quad[order].n; iqp++) {
                elem_nrml_vel_quad = get_function_quad(quad[order].pt[iqp].lshape, normal_vel, nnodes); // evaluate velocity at quadrature point
                elem_rhs[idof] += dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] * elem_nrml_vel_quad; // quadratic
            }
        } else {
            // Natural normal flow BC imposed, linear integrand
            order=1; iqp = 0;
            elem_rhs[idof] += dt * quad[order].pt[iqp].djac3d_fixed * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] * normal_discharge; // linear
        }
    }
}

//*****************************************************************//
//*****************************************************************//
// Convective SW3D Integration :: 2D NB integration :: Outflow
void fe_sw3_2d_int_OUTFLOW(SQUAD *quad, DOF_3 *elem_rhs, SSW_3D_ELEM *elem, double *normal_vel, double dt, int nnodes_quad, int *mask, double gravity) {
    int iqp=0, idof=0, order=0;
    double elem_prs_quad=0., elem_nrml_vel_quad=0., temp = 0.;
    SVECT elem_v_quad; svect_init(&elem_v_quad);
    
    int nnodes = elem->elem2d->nnodes;
    
    for (idof=0; idof<nnodes; idof++) {
        
        // Depth-Averaged Continuity Equation
        order = 2;
        for (iqp=0; iqp<quad[order].n; iqp++) {
            elem_nrml_vel_quad = get_function_quad(quad[order].pt[iqp].lshape, normal_vel, nnodes); // evaluate velocity at quadrature point
            elem_rhs[idof].c_eq += dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] * elem_nrml_vel_quad; // quadratic
        }
        
        // Momentum Equation (in projected orthogonal coordinates)
        if (mask[idof * nnodes] != 2) {
            order=3;
            for (iqp=0; iqp<quad[order].n; iqp++) {
                elem_prs_quad = get_function_quad(quad[order].pt[iqp].lshape_quad, elem->pressure, nnodes_quad); // this needs to be quadratric basis!!
                elem_v_quad = get_velocity2d_quad(quad[order].pt[iqp], elem->vel);
                elem_nrml_vel_quad = get_function_quad(quad[order].pt[iqp].lshape, normal_vel, nnodes); // evaluate velocity at quadrature point
                
                temp = dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof]; // linear
                
                elem_rhs[idof].x_eq += temp * elem->elem2d->nrml.x * gravity * elem_prs_quad; // cubic
                elem_rhs[idof].y_eq += temp * elem->elem2d->nrml.y * gravity * elem_prs_quad; // cubic
                elem_rhs[idof].x_eq += temp * elem_nrml_vel_quad * elem_v_quad.x; // cubic
                elem_rhs[idof].y_eq += temp * elem_nrml_vel_quad * elem_v_quad.y; // cubic
            }
        } else {
            order=3;
            for (iqp=0; iqp<quad[order].n; iqp++) {
                elem_prs_quad = get_function_quad(quad[order].pt[iqp].lshape_quad, elem->pressure, nnodes_quad); // this needs to be quadratric basis!!
                elem_nrml_vel_quad = get_function_quad(quad[order].pt[iqp].lshape, normal_vel, nnodes); // evaluate velocity at quadrature point
                elem_v_quad = get_velocity2d_quad(quad[order].pt[iqp], elem->vel);
                
                temp = dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof]; // linear
                
                elem_rhs[idof].y_eq += temp * gravity * elem_prs_quad * (elem->tanvec[idof].x * elem->elem2d->nrml.x + elem->tanvec[idof].y * elem->elem2d->nrml.y); // cubic
                elem_rhs[idof].y_eq += temp * elem_nrml_vel_quad * (elem->tanvec[idof].x * elem_v_quad.x + elem->tanvec[idof].y * elem_v_quad.y); // cubic
            }
        }
        
    }
}

//*****************************************************************//
//*****************************************************************//
// Convective WVEL SW3D Integration :: 2D NB integration :: Outflow
void fe_sw3_2d_int_WVEL_OUTFLOW(SQUAD *quad, double *elem_rhs, SSW_3D_ELEM *elem, double *normal_vel, double dt) {
    int iqp=0, idof=0, order=2;
    double elem_nrml_vel_quad=0., temp = 0.;
    SVECT elem_v_quad; svect_init(&elem_v_quad);
    
    int nnodes = elem->elem2d->nnodes;
    
    for (idof=0; idof<nnodes; idof++) {
        for (iqp=0; iqp<quad[order].n; iqp++) {
            elem_nrml_vel_quad = get_function_quad(quad[order].pt[iqp].lshape, normal_vel, nnodes); // evaluate velocity at quadrature point
            elem_rhs[idof] += dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] * elem_nrml_vel_quad; // quadratic
        }
    }
}

//*****************************************************************//
//*****************************************************************//
// SW3D Integration :: 2D NB integration :: BED
void fe_sw3_2d_int_BED(SQUAD *quad, DOF_3 *elem_rhs, SSW_3D_ELEM *elem, double dt, int nnodes_quad, int *mask, double gravity, double drag) {
    int iqp=0, idof=0, order=0;
    double elem_prs_quad=0., temp = 0., elem_v_quad_mag = 0.;
    SVECT elem_v_quad; svect_init(&elem_v_quad);
    
    int nnodes = elem->elem2d->nnodes;
    
    for (idof=0; idof<nnodes; idof++) {
        
        // momentum additions
        if (mask[idof * nnodes] != 2) {
            order=3;
            for (iqp=0; iqp<quad[order].n; iqp++) {
                elem_prs_quad = get_function_quad(quad[order].pt[iqp].lshape_quad, elem->pressure, nnodes_quad); // this needs to be quadratric basis!!
                elem_v_quad = get_velocity2d_quad(quad[order].pt[iqp], elem->vel);
                elem_v_quad_mag = sqrt(elem_v_quad.x*elem_v_quad.x + elem_v_quad.y*elem_v_quad.y);
                temp = dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof]; // linear
                
                elem_rhs[idof].x_eq += temp * elem->elem2d->nrml.x * gravity * elem_prs_quad; // cubic
                elem_rhs[idof].y_eq += temp * elem->elem2d->nrml.y * gravity * elem_prs_quad; // cubic
                elem_rhs[idof].x_eq += temp * (1./2.) * drag * elem_v_quad_mag * elem_v_quad.x; // cubic
                elem_rhs[idof].y_eq += temp * (1./2.) * drag * elem_v_quad_mag * elem_v_quad.y; // cubic
            }
        } else { //(in projected orthogonal coordinates)
            order=3;
            for (iqp=0; iqp<quad[order].n; iqp++) {
                elem_prs_quad = get_function_quad(quad[order].pt[iqp].lshape_quad, elem->pressure, nnodes_quad); // this needs to be quadratric basis!!
                elem_v_quad = get_velocity2d_quad(quad[order].pt[iqp], elem->vel);
                elem_v_quad_mag = sqrt(elem_v_quad.x*elem_v_quad.x + elem_v_quad.y*elem_v_quad.y);
                
                temp = dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] ; // linear
                elem_rhs[idof].y_eq += temp * gravity * elem_prs_quad * (elem->elem2d->nrml.x * elem->tanvec[idof].x + elem->elem2d->nrml.y * elem->tanvec[idof].y);  // cubic
                elem_rhs[idof].y_eq += temp * (1./2.) * drag * elem_v_quad_mag * (elem->tanvec[idof].x * elem_v_quad.x + elem->tanvec[idof].y * elem_v_quad.y); // cubic
            }
        }
    }
}

// depth-averaged continuity additions (only when bed displacement happens from sediment)
// NEEDS TESTING!
void fe_sw3_2d_int_BED_DPL(SQUAD *quad, DOF_3 *elem_rhs, double *elem_eta, double *elem_eta_old, int nnodes) {
    int iqp=0, idof=0, order=2; // quadratic
    double elem_eta_quad = 0., elem_eta_old_quad = 0., temp=0.;
    for (idof=0; idof<nnodes; idof++) {
        for (iqp=0; iqp<quad[order].n; iqp++) {
            elem_eta_quad = get_function_quad(quad[order].pt[iqp].lshape, elem_eta, nnodes);
            elem_eta_old_quad = get_function_quad(quad[order].pt[iqp].lshape, elem_eta_old, nnodes);
            temp = quad[order].pt[iqp].djac2d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof];
            
            elem_rhs[idof].c_eq += temp * (elem_eta_quad - elem_eta_old_quad);
        }
    }
}


//*****************************************************************//
//*****************************************************************//
// SW3D Integration :: 2D FRS DPL & WATER SOURCE integration
void fe_sw3_2d_int_FRS(SQUAD *quad, DOF_3 *elem_rhs, double *elem_eta, double *elem_eta_old, double source, double dt, int nnodes) {
    int iqp=0, idof=0, order=2; // quadratic
    double elem_eta_quad = 0., elem_eta_old_quad = 0., temp=0.;
    
    for (idof=0; idof<nnodes; idof++) {
        for (iqp=0; iqp<quad[order].n; iqp++) {
            elem_eta_quad = get_function_quad(quad[order].pt[iqp].lshape, elem_eta, nnodes);
            elem_eta_old_quad = get_function_quad(quad[order].pt[iqp].lshape, elem_eta_old, nnodes);
            temp = quad[order].pt[iqp].djac2d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof];
            
            elem_rhs[idof].c_eq += temp * (elem_eta_quad - elem_eta_old_quad + source * dt);
        }
    }
}

// SW3D Integration :: 2D FRS STRESS integration
void fe_sw3_2d_int_FRS_STRESS(SQUAD *quad, DOF_3 *elem_rhs, SVECT2D *elem_tanvec, double dt, SVECT2D stress, int *mask, int nnodes) {
    int iqp=0, idof=0, order=1; // Surfaces stresses are treated as constant over an element, so this is linear integration
    double temp=0;
    
    for (idof=0; idof<nnodes; idof++) {
        // momentum additions
        if (mask[idof * nnodes] != 2) {
            for (iqp=0; iqp<quad[order].n; iqp++) {
                temp = dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof]; // linear
                
                elem_rhs[idof].x_eq -= temp * stress.x; // linear
                elem_rhs[idof].y_eq -= temp * stress.y; // linear
            }
        } else { //(in projected orthogonal coordinates)
            for (iqp=0; iqp<quad[order].n; iqp++) {
                temp = dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] ; // linear
                elem_rhs[idof].y_eq -= temp * (stress.x * elem_tanvec[idof].x + stress.y * elem_tanvec[idof].y);  // linear
            }
        }
    }
}

//*****************************************************************//
//*****************************************************************//
// WVEL SW3D Integration :: 2D FRS DPL & WATER SOURCE integration
void fe_sw3_2d_int_WVEL_FRS(SQUAD *quad, double *elem_rhs, double *elem_eta, double *elem_eta_old, double source, double dt, double nrml_z, int nnodes) {
    int iqp=0, idof=0, order=2; // quadratic
    double elem_eta_quad = 0., elem_eta_old_quad = 0., temp=0.;
    
    for (idof=0; idof<nnodes; idof++) {
        for (iqp=0; iqp<quad[order].n; iqp++) {
            elem_eta_quad = get_function_quad(quad[order].pt[iqp].lshape, elem_eta, nnodes);
            elem_eta_old_quad = get_function_quad(quad[order].pt[iqp].lshape, elem_eta_old, nnodes);
            temp = quad[order].pt[iqp].djac2d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof];
            
            elem_rhs[idof] += temp * nrml_z * (elem_eta_quad - elem_eta_old_quad + source * dt);
        }
    }
}


////*****************************************************************//
////*****************************************************************//
// Convective SW3D Integration :: 2D NB integration :: Outflow
void fe_sw3_2d_int_TAILWATER(SQUAD *quad, DOF_3 *elem_rhs, SSW_3D_ELEM *elem, double *normal_vel, double dt, int nnodes_quad, int *mask, double gravity, double *dprs) {
    int iqp=0, idof=0, order=0;
    double elem_prs_quad=0., elem_dprs_quad=0., elem_nrml_vel_quad=0., temp = 0.;
    SVECT elem_v_quad; svect_init(&elem_v_quad);
    
    int nnodes = elem->elem2d->nnodes;
    
    //**********************************
    // Water Elevation/Pressure ********
    
    for (idof=0; idof<nnodes; idof++) {
        
        if (mask[idof * nnodes] != 2) {
            order=3;
            for (iqp=0; iqp<quad[order].n; iqp++) {
                elem_prs_quad = get_function_quad(quad[order].pt[iqp].lshape_quad, elem->pressure, nnodes_quad); // original pressure
                elem_dprs_quad = get_function_quad(quad[order].pt[iqp].lshape_quad, dprs, nnodes_quad); // pressure difference applying new water elevation
                temp = dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof]; // linear
                
                elem_rhs[idof].x_eq += temp * elem->elem2d->nrml.x * gravity * (elem_prs_quad + elem_dprs_quad); // cubic
                elem_rhs[idof].y_eq += temp * elem->elem2d->nrml.y * gravity * (elem_prs_quad + elem_dprs_quad); // cubic
            }
        } else { // rotate equations
            order=3;
            for (iqp=0; iqp<quad[order].n; iqp++) {
                elem_prs_quad = get_function_quad(quad[order].pt[iqp].lshape_quad, elem->pressure, nnodes_quad); // original pressure
                elem_dprs_quad = get_function_quad(quad[order].pt[iqp].lshape_quad, dprs, nnodes_quad); // pressure difference applying new water elevation
                temp = dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof]; // linear
                
                elem_rhs[idof].y_eq += temp * gravity * (elem_prs_quad + elem_dprs_quad) * (elem->elem2d->nrml.x * elem->tanvec[idof].x + elem->elem2d->nrml.y * elem->tanvec[idof].y);  // cubic
            }
        }
    }
    
    
    //**********************************
    // Outflow *************************
    
    // Depth-Averaged Continuity Equation
    order = 2;
    for (idof=0; idof<nnodes; idof++) {
        for (iqp=0; iqp<quad[order].n; iqp++) {
            elem_nrml_vel_quad = get_function_quad(quad[order].pt[iqp].lshape, normal_vel, nnodes); // evaluate velocity at quadrature point
            elem_rhs[idof].c_eq += dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof] * elem_nrml_vel_quad; // quadratic
        }
    }
    
    // Momentum
    if (sarray_sum_dbl(normal_vel, nnodes) > 0) {
        for (idof=0; idof<nnodes; idof++) {
            
            if (mask[idof * nnodes] != 2) {
                order=3;
                for (iqp=0; iqp<quad[order].n; iqp++) {
                    elem_v_quad = get_velocity2d_quad(quad[order].pt[iqp], elem->vel);
                    elem_nrml_vel_quad = get_function_quad(quad[order].pt[iqp].lshape, normal_vel, nnodes); // evaluate velocity at quadrature point
                    
                    temp = dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof]; // linear
                    elem_rhs[idof].x_eq += temp * elem_nrml_vel_quad * elem_v_quad.x; // cubic
                    elem_rhs[idof].y_eq += temp * elem_nrml_vel_quad * elem_v_quad.y; // cubic
                }
            } else { // rotate equations
                order=3;
                for (iqp=0; iqp<quad[order].n; iqp++) {
                    elem_nrml_vel_quad = get_function_quad(quad[order].pt[iqp].lshape, normal_vel, nnodes); // evaluate velocity at quadrature point
                    elem_v_quad = get_velocity2d_quad(quad[order].pt[iqp], elem->vel);
                    
                    temp = dt * quad[order].pt[iqp].djac2d_3d * quad[order].pt[iqp].w * quad[order].pt[iqp].lshape[idof]; // linear
                    elem_rhs[idof].y_eq += temp * elem_nrml_vel_quad * (elem->tanvec[idof].x * elem_v_quad.x + elem->tanvec[idof].y * elem_v_quad.y); // cubic
                }
            }
        }
    }
}



