#include "global_header.h" 


// *************************************************************************
// *************************************************************************
// Analytic Integration on A Line Segment

// ------------------------------------------------------------------------
// integral(integral(constant phi * (vh dot n), yhat, -1,1), xhat, -1, 1)
void inline integrate_line_phi_VHdotN(SVECT2D *v, double *h, SVECT2D nrml, double s, double djac, double *integral) {
    
    double const_contrib0;  /* the contributions */
    double sum, sum_u, sum_v, sum_h;
    sum_h = (h[0] + h[1]);
    sum_u = (v[0].x + v[1].x);
    sum_v = (v[0].y + v[1].y);
    const_contrib0 = djac * s / 12.;
    
    /* contributions to the continuity equation */
    integral[0] = const_contrib0 * (sum_h * sum_u + 2. * h[0] * v[0].x) * nrml.x;
    integral[1] = const_contrib0 * (sum_h * sum_u + 2. * h[1] * v[1].x) * nrml.x;
    integral[0] += const_contrib0 * (sum_h * sum_v + 2. * h[0] * v[0].y) * nrml.y;
    integral[1] += const_contrib0 * (sum_h * sum_v + 2. * h[1] * v[1].y) * nrml.y;
    
}

// ------------------------------------------------------------------------
// integral(integral(constant phi, yhat, -1,1), xhat, -1, 1)
void inline integrate_line_phi(double s, double djac, double *integral) {
    
    double const_contrib = one_2 * djac * s;
    integral[0] += const_contrib;
    integral[1] += const_contrib;
    
}



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Analytic Integration on Quadrilaterals
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    

// ------------------------------------------------------------------------
// integral(integral(dphi_x * fx * gx * hx + dphi_y * fy * gy * hy, yhat, -1,1), xhat, -1, 1) || 2D convection type integrals
// note :: the following does NOT assume an elementally averaged depth, triangle AdH does
// note :: quadrature for a 5th-order integrand on quads = 4 evals X 5 vars X 8 expansion evals X 4 evals/basis = 640 operations;
// note :: quadrature for a 5th-order integrand on quads w/ grouping of 1 or 2 vars = 9 evals X 5 vars X 8 expansion evals X 4 evals/basis = 1440 operations;
// note :: total operation count here :: 929 operations
// note :: moral of the story, even if one combo of variables can be grouped, use quadrature.  Otherwise, analytic saves operations.
//void inline integrate_quad_gradPhi_dot_Vfh(SVECT *nd, double constant, double *fx, double *gx, double *hx, double *fy, double *gy, double *hy, double *integral) {
//    
//    double con;
//    double dy1, dy2, dx1, dx2, py1, py2, py3, py4, py5, py6, py7, py8, px1, px2, px3, px4, px5, px6, px7, px8;
//    double fyt1, fyt2, fyt3, fyt4, fyt5, fyt6, fyt7, fyt8, fyt9, fyt10, fyt11, fyt12, fyt13, fyt14;
//    double fxt1, fxt2, fxt3, fxt4, fxt5, fxt6, fxt7, fxt8, fxt9, fxt10, fxt11, fxt12, fxt13, fxt14;
//    double fyt15, fyt16, fyt17, fyt18, fyt19, fyt20, fyt21, fyt22, fyt23, fyt24, fyt25, fyt26, fyt27, fyt28;
//    double fxt15, fxt16, fxt17, fxt18, fxt19, fxt20, fxt21, fxt22, fxt23, fxt24, fxt25, fxt26, fxt27, fxt28;
//    double ty1, ty2, ty3, ty4, ty5, ty6, ty7, ty8, ty9, ty10, ty11, ty12, ty13, ty14, ty15, ty16, ty17, ty18, ty19, ty20;
//    double tx1, tx2, tx3, tx4, tx5, tx6, tx7, tx8, tx9, tx10, tx11, tx12, tx13, tx14, tx15, tx16, tx17, tx18, tx19, tx20;
//    double hty1, hty2, hty3, hty4, hty5, hty6;
//    double htx1, htx2, htx3, htx4, htx5, htx6;
//    
//    // operation count = 8 X 4 = 32;
//    double x1 = nd[0].x, x2 = nd[1].x, x3 = nd[2].x, x4 = nd[3].x;
//    double y1 = nd[0].y, y2 = nd[1].y, y3 = nd[2].y, y4 = nd[3].y;
//    double fx1 = fx[0], fx2 = fx[1], fx3 = fx[2], fx4 = fx[3];
//    double gx1 = gx[0], gx2 = gx[1], gx3 = gx[2], gx4 = gx[3];
//    double hx1 = hx[0], hx2 = hx[1], hx3 = hx[2], hx4 = hx[3];
//    double fy1 = fy[0], fy2 = fy[1], fy3 = fy[2], fy4 = fy[3];
//    double gy1 = gy[0], gy2 = gy[1], gy3 = gy[2], gy4 = gy[3];
//    double hy1 = hy[0], hy2 = hy[1], hy3 = hy[2], hy4 = hy[3];
//
//    // operation count = 67
//    dy1 = (fy2 - fy4); dx1 = (fx2 - fx4);
//    dy2 = (fy1 - fy3); dx2 = (fx1 - fx3);
//    py1 = (3*fy1 + 6*fy2 + fy3); px1 = (3*fx1 + 6*fx2 + fx3);
//    py2 = (3*fy1 + fy3 + 6*fy4); px2 = (3*fx1 + fx3 + 6*fx4);
//    py3 = (fy1 + 6*fy2 + 3*fy3); px3 = (fx1 + 6*fx2 + 3*fx3);
//    py4 = (fy1 + 3*fy3 + 6*fy4); px4 = (fx1 + 3*fx3 + 6*fx4);
//    py5 = (6*fy1 + 3*fy2 + fy4); px5 = (6*fx1 + 3*fx2 + fx4);
//    py6 = (6*fy1 + fy2 + 3*fy4); px6 = (6*fx1 + fx2 + 3*fx4);
//    py7 = (3*fy2 + 6*fy3 + fy4); px7 = (3*fx2 + 6*fx3 + fx4);
//    py8 = (fy2 + 6*fy3 + 3*fy4); px8 = (fx2 + 6*fx3 + 3*fx4);
//    
//    // operation count = 28 * 2 * 6 = 336
//    fyt1 =  (12*fy1 + 4*fy2 + fy3 + 3*fy4);     fxt1 =  (12*fx1 + 4*fx2 + fx3 + 3*fx4);
//    fyt2 =  (4*fy1 + 4*fy2 + fy3 + fy4);        fxt2 =  (4*fx1 + 4*fx2 + fx3 + fx4);
//    fyt3 =  (3*fy1 + 3*fy2 + 2*fy3 + 2*fy4);    fxt3 =  (3*fx1 + 3*fx2 + 2*fx3 + 2*fx4);
//    fyt4 =  (9*fy1 + 3*fy2 + 2*fy3 + 6*fy4);    fxt4 =  (9*fx1 + 3*fx2 + 2*fx3 + 6*fx4);
//    fyt5 =  (4*fy1 + 12*fy2 + 3*fy3 + fy4);     fxt5 =  (4*fx1 + 12*fx2 + 3*fx3 + fx4);
//    fyt6 =  (3*fy1 + 9*fy2 + 6*fy3 + 2*fy4);    fxt6 =  (3*fx1 + 9*fx2 + 6*fx3 + 2*fx4);
//    fyt7 =  (2*fy1 + 6*fy2 + 9*fy3 + 3*fy4);    fxt7 =  (2*fx1 + 6*fx2 + 9*fx3 + 3*fx4);
//    fyt8 =  (2*fy1 + 2*fy2 + 3*fy3 + 3*fy4);    fxt8 =  (2*fx1 + 2*fx2 + 3*fx3 + 3*fx4);
//    fyt9 =  (6*fy1 + 2*fy2 + 3*fy3 + 9*fy4);    fxt9 =  (6*fx1 + 2*fx2 + 3*fx3 + 9*fx4);
//    fyt10 = (6*fy1 + 27*fy2 + 6*fy3 + fy4);     fxt10 = (6*fx1 + 27*fx2 + 6*fx3 + fx4);
//    fyt11 = (6*fy1 + fy2 + 6*fy3 + 27*fy4);     fxt11 = (6*fx1 + fx2 + 6*fx3 + 27*fx4);
//    fyt12 = (12*fy1 + 3*fy2 + fy3 + 4*fy4);     fxt12 = (12*fx1 + 3*fx2 + fx3 + 4*fx4);
//    fyt13 = (9*fy1 + 6*fy2 + 2*fy3 + 3*fy4);    fxt13 = (9*fx1 + 6*fx2 + 2*fx3 + 3*fx4);
//    fyt14 = (3*fy1 + 2*fy2 + 2*fy3 + 3*fy4);    fxt14 = (3*fx1 + 2*fx2 + 2*fx3 + 3*fx4);
//    fyt15 = (4*fy1 + fy2 + fy3 + 4*fy4);        fxt15 = (4*fx1 + fx2 + fx3 + 4*fx4);
//    fyt16 = (6*fy1 + 9*fy2 + 3*fy3 + 2*fy4);    fxt16 = (6*fx1 + 9*fx2 + 3*fx3 + 2*fx4);
//    fyt17 = (2*fy1 + 3*fy2 + 3*fy3 + 2*fy4);    fxt17 = (2*fx1 + 3*fx2 + 3*fx3 + 2*fx4);
//    fyt18 = (2*fy1 + 3*fy2 + 9*fy3 + 6*fy4);    fxt18 = (2*fx1 + 3*fx2 + 9*fx3 + 6*fx4);
//    fyt19 = (3*fy1 + 2*fy2 + 6*fy3 + 9*fy4);    fxt19 = (3*fx1 + 2*fx2 + 6*fx3 + 9*fx4);
//    fyt20 = (4*fy1 + fy2 + 3*fy3 + 12*fy4);     fxt20 = (4*fx1 + fx2 + 3*fx3 + 12*fx4);
//    fyt21 = (fy1 + fy2 + 4*fy3 + 4*fy4);        fxt21 = (fx1 + fx2 + 4*fx3 + 4*fx4);
//    fyt22 = (3*fy1 + fy2 + 4*fy3 + 12*fy4);     fxt22 = (3*fx1 + fx2 + 4*fx3 + 12*fx4);
//    fyt23 = (fy1 + 3*fy2 + 12*fy3 + 4*fy4);     fxt23 = (fx1 + 3*fx2 + 12*fx3 + 4*fx4);
//    fyt24 = (fy1 + 4*fy2 + 4*fy3 + fy4);        fxt24 = (fx1 + 4*fx2 + 4*fx3 + fx4);
//    fyt25 = (fy1 + 6*fy2 + 27*fy3 + 6*fy4);     fxt25 = (fx1 + 6*fx2 + 27*fx3 + 6*fx4);
//    fyt26 = (3*fy1 + 12*fy2 + 4*fy3 + fy4);     fxt26 = (3*fx1 + 12*fx2 + 4*fx3 + fx4);
//    fyt27 = (fy1 + 4*fy2 + 12*fy3 + 3*fy4);     fxt27 = (fx1 + 4*fx2 + 12*fx3 + 3*fx4);
//    fyt28 = (27*fy1 + 6*fy2 + fy3 + 6*fy4);     fxt28 = (27*fx1 + 6*fx2 + fx3 + 6*fx4);
//
//    // operation count = 20 * 2 * 9 = 360
//    ty1 =  (3*fyt2*gy1 + 3*fyt5*gy2 + fyt6*gy3 + fyt3*gy4);      tx1 =  (3*fxt2*gx1 + 3*fxt5*gx2 + fxt6*gx3 + fxt3*gx4);
//    ty2 =  (fyt3*gy1 + fyt6*gy2 + fyt7*gy3 + fyt8*gy4);          tx2 =  (fxt3*gx1 + fxt6*gx2 + fxt7*gx3 + fxt8*gx4);
//    ty3 =  (fyt4*gy1 + fyt3*gy2 + fyt8*gy3 + fyt9*gy4);          tx3 =  (fxt4*gx1 + fxt3*gx2 + fxt8*gx3 + fxt9*gx4);
//    ty4 =  (3*fyt1*gy1 + 3*fyt2*gy2 + fyt3*gy3 + fyt4*gy4);      tx4 =  (3*fxt1*gx1 + 3*fxt2*gx2 + fxt3*gx3 + fxt4*gx4);
//    ty5 =  (3*fyt12*gy1 + fyt13*gy2 + fyt14*gy3 + 3*fyt15*gy4);  tx5 =  (3*fxt12*gx1 + fxt13*gx2 + fxt14*gx3 + 3*fxt15*gx4);
//    ty6 =  (fyt14*gy1 + fyt17*gy2 + fyt18*gy3 + fyt19*gy4);      tx6 =  (fxt14*gx1 + fxt17*gx2 + fxt18*gx3 + fxt19*gx4);
//    ty7 =  (fyt13*gy1 + fyt16*gy2 + fyt17*gy3 + fyt14*gy4);      tx7 =  (fxt13*gx1 + fxt16*gx2 + fxt17*gx3 + fxt14*gx4);
//    ty8 =  (fyt16*gy1 + 3*fyt26*gy2 + 3*fyt24*gy3 + fyt17*gy4);  tx8 =  (fxt16*gx1 + 3*fxt26*gx2 + 3*fxt24*gx3 + fxt17*gx4);
//    ty9 =  (fyt17*gy1 + 3*fyt24*gy2 + 3*fyt27*gy3 + fyt18*gy4);  tx9 =  (fxt17*gx1 + 3*fxt24*gx2 + 3*fxt27*gx3 + fxt18*gx4);
//    ty10 = (3*fyt15*gy1 + fyt14*gy2 + fyt19*gy3 + 3*fyt20*gy4);  tx10 = (3*fxt15*gx1 + fxt14*gx2 + fxt19*gx3 + 3*fxt20*gx4);
//    ty11 = (fyt8*gy1 + fyt7*gy2 + 3*fyt23*gy3 + 3*fyt21*gy4);    tx11 = (fxt8*gx1 + fxt7*gx2 + 3*fxt23*gx3 + 3*fxt21*gx4);
//    ty12 = (fyt9*gy1 + fyt8*gy2 + 3*fyt21*gy3 + 3*fyt22*gy4);    tx12 = (fxt9*gx1 + fxt8*gx2 + 3*fxt21*gx3 + 3*fxt22*gx4);
//    ty13 = (fyt28*gy1 + py5*gy2 + dy2*gy3 + py6*gy4);            tx13 = (fxt28*gx1 + px5*gx2 + dx2*gx3 + px6*gx4);
//    ty14 = (3*dy1*gy1 + py1*gy2 + dy1*gy3 - py2*gy4);            tx14 = (3*dx1*gx1 + px1*gx2 + dx1*gx3 - px2*gx4);
//    ty15 = (dy1*gy1 + py3*gy2 + 3*dy1*gy3 - py4*gy4);            tx15 = (dx1*gx1 + px3*gx2 + 3*dx1*gx3 - px4*gx4);
//    ty16 = (py1*gy1 + fyt10*gy2 + py3*gy3 + dy1*gy4);            tx16 = (px1*gx1 + fxt10*gx2 + px3*gx3 + dx1*gx4);
//    ty17 = (py2*gy1 - dy1*gy2 + py4*gy3 + fyt11*gy4);            tx17 = (px2*gx1 - dx1*gx2 + px4*gx3 + fxt11*gx4);
//    ty18 = (py5*gy1 + 3*dy2*gy2 - py7*gy3 + dy2*gy4);            tx18 = (px5*gx1 + 3*dx2*gx2 - px7*gx3 + dx2*gx4);
//    ty19 = (dy2*gy1 - py7*gy2 - fyt25*gy3 - py8*gy4);            tx19 = (dx2*gx1 - px7*gx2 - fxt25*gx3 - px8*gx4);
//    ty20 = (py6*gy1 + dy2*gy2 - py8*gy3 + 3*dy2*gy4);            tx20 = (px6*gx1 + dx2*gx2 - px8*gx3 + 3*dx2*gx4);
//   
//    // operation count = 6 * 2 * 7 = 84
//    hty1 = (ty4*hy1 + ty1*hy2 + ty2*hy3 + ty3*hy4);      htx1 = (tx4*hx1 + tx1*hx2 + tx2*hx3 + tx3*hx4);
//    hty2 = (ty14*hy1 + ty16*hy2 + ty15*hy3 - ty17*hy4);  htx2 = (tx14*hx1 + tx16*hx2 + tx15*hx3 - tx17*hx4);
//    hty3 = (ty5*hy1 + ty7*hy2 + ty6*hy3 + ty10*hy4);     htx3 = (tx5*hx1 + tx7*hx2 + tx6*hx3 + tx10*hx4);
//    hty4 = (ty7*hy1 + ty8*hy2 + ty9*hy3 + ty6*hy4);      htx4 = (tx7*hx1 + tx8*hx2 + tx9*hx3 + tx6*hx4);
//    hty5 = (ty13*hy1 + ty18*hy2 + ty19*hy3 + ty20*hy4);  htx5 = (tx13*hx1 + tx18*hx2 + tx19*hx3 + tx20*hx4)
//    hty6 = (ty3*hy1 + ty2*hy2 + ty11*hy3 + ty12*hy4);    htx6 = (tx3*hx1 + tx2*hx2 + tx11*hx3 + tx12*hx4);
//    
//    // operation count = 4 * 12 + 2 = 50
//    con = (1./720.) * constant;
//    integral[0] = con * (-hty1*x2 + hty2*x3 + hty3*x4 + htx1*y2 - htx2*y3 - htx3*y4);
//    integral[1] = con * (+hty1*x1 - hty4*x3 - hty5*x4 - htx1*y1 + htx4*y3 + htx5*y4);
//    integral[2] = con * (-hty2*x1 + hty4*x2 - hty6*x4 + htx2*y1 - htx4*y2 + htx6*y4);
//    integral[3] = con * (-hty3*x1 + hty5*x2 + hty6*x3 + htx3*y1 - htx5*y2 - htx6*y3);
//
//}
//
//
//
//
//
//
//
//
//
//
//// ------------------------------------------------------------------------
//// integral(integral( constant * (dphi_i_dx * h * sigma_x + dphi_i_dy * h * sigma_y), yhat, -1,1), xhat, -1, 1)
//void inline integrate_quad_diffusion_quadrature(SVECT *nd, SQUAD *quad, double constant, double ev_st, double vst_i, double vst_j, double ev_tr, double *h, double *u, double *v, double *integral_x, double *integral_y) {
//    
//    int i, iqp, quad_order = 4; // quadrature order
//    double qp_h = 0., qp_sigma_xx = 0., qp_sigma_xy = 0., qp_sigma_yx = 0., qp_sigma_yy = 0.;
//    SVECT qp_grad_u, qp_grad_v;
//    SQUAD_PT *qp = NULL;
//    
//    sarray_init_dbl(integral_x,NDONQUAD);
//    sarray_init_dbl(integral_y,NDONQUAD);
//    for (iqp=0; iqp<quad[quad_order].n; iqp++) {
//        qp = &(quad[quad_order].pt[iqp]);
//        
//        qp->djac = elem2d_get_quad_djac(qp->xhat, qp->yhat, qp->zhat, nd); // evaluate prism djac at quadrature point
//        elem2d_get_quad_shape_gradients(qp->xhat, qp->yhat, nd, qp->grad_shp); // evaluate shape gradients at quadrature points
//        
//        // evaluate depth at quadrature point
//        qp_h = SQUAD_get_function(qp, h, NDONQUAD);
//
//        // evaluate grad(u) and grad(v) at quadrature points
//        svect_init(&qp_grad_u); svect_init(&qp_grad_v);
//        for (i=0; i<NDONQUAD; i++) {
//            qp_grad_u.x += u[i] * qp->grad_shp[i].x;  qp_grad_v.x += v[i] * qp->grad_shp[i].x;
//            qp_grad_u.y += u[i] * qp->grad_shp[i].y;  qp_grad_v.z += v[i] * qp->grad_shp[i].z;
//        }
//        
//        // sigma_xx = 2 ev.xx du_dx
//        // sigma_xy = ev.xy (du_dy + dv_dx)
//        // sigma_yx = sigma_xy
//        // sigma_yy = 2 ev.yy dv_dy
//        qp_sigma_xx = 2. * ev_st * vst_i * (vst_i * qp_grad_u.x + vst_j * qp_grad_u.y) + 2. * ev_tr * qp_grad_u.x;
//        qp_sigma_xy = ev_st * vst_j * (vst_j * qp_grad_u.y + vst_i * qp_grad_u.x + vst_i * qp_grad_v.x + vst_j * qp_grad_v.y) + ev_tr * (qp_grad_u.y + qp_grad_v.x);
//        qp_sigma_yx = ev_st * vst_i * (vst_j * qp_grad_u.y + vst_i * qp_grad_u.x + vst_i * qp_grad_v.x + vst_j * qp_grad_v.y) + ev_tr * (qp_grad_u.y + qp_grad_v.x);
//        qp_sigma_yy = 2. * ev_st * vst_j * (vst_j * qp_grad_v.y  + vst_i * qp_grad_v.x) + 2. * ev_tr * qp_grad_v.y ;
//        
//        double t1 = constant * qp->djac * qp->w * qp_h;
//        for (i=0; i<NDONQUAD; i++) {
//            integral_x[i] += t1 * (qp->grad_shp[i].x * qp_sigma_xx + qp->grad_shp[i].y * qp_sigma_xy);
//            integral_y[i] += t1 * (qp->grad_shp[i].x * qp_sigma_xy + qp->grad_shp[i].y * qp_sigma_yy);
//        }
//    }
//}



    
// *************************************************************************
// *************************************************************************
// Analytic Integration on Triangular Prisms




// ------------------------------------------------------------------------
// integral(integral(cosntant * dphi_i/d{x,y} * df/d{x,y}, yhat, 0,1-xhat), xhat, 0, 1), zhat, -1, 1)
// uses quadrature, cannot directly integrate in sage
// total operation count = alot
void inline integrate_triPrism_dphiXdfg_quadrature(SVECT *nd, SQUAD *quad, double constant, SVECT constant_f, SVECT constant_g, double *f, double *g, double *integral_x, double *integral_y) {

    int i, iqp, quad_order = 4; // quadrature order
    double qp_f = 0., qp_g = 0.;
    SQUAD_PT *qp = NULL;
    
    sarray_init_dbl(integral_x,NDONPRISM);
    sarray_init_dbl(integral_y,NDONPRISM);
    for (iqp=0; iqp<quad[quad_order].n; iqp++) {
        qp = &(quad[quad_order].pt[iqp]);
        
        qp->djac = elem3d_get_prism_djac(qp->xhat, qp->yhat, qp->zhat, nd); // evaluate prism djac at quadrature point
        elem3d_get_triprism_shape_gradients(qp->xhat, qp->yhat, qp->zhat, nd, qp->grad_shp); // evaluate shape gradients at quadrature points
        qp_f = SQUAD_get_function(qp, f, NDONPRISM); // evaluate function at quadrature point
        qp_g = SQUAD_get_function(qp, g, NDONPRISM); // evaluate function at quadrature point
        
        // evalute grad(f) at quadrature points
        SVECT qp_grad_f; svect_init(&qp_grad_f);
        SVECT qp_grad_g; svect_init(&qp_grad_g);
        for (i=0; i<NDONPRISM; i++) {
            qp_grad_f.x += f[i] * qp->grad_shp[i].x;
            qp_grad_f.y += f[i] * qp->grad_shp[i].y;
            qp_grad_f.z += f[i] * qp->grad_shp[i].z;
            qp_grad_g.x += g[i] * qp->grad_shp[i].x;
            qp_grad_g.y += g[i] * qp->grad_shp[i].y;
            qp_grad_g.z += g[i] * qp->grad_shp[i].z;
        }
        qp_grad_f.x *= constant_f.x;
        qp_grad_f.y *= constant_f.y;
        qp_grad_f.z *= constant_f.z;
        qp_grad_g.x *= constant_g.x;
        qp_grad_g.y *= constant_g.y;
        qp_grad_g.z *= constant_g.z;
        
        double t1 = constant * qp->djac * qp->w;
        for (i=0; i<NDONPRISM; i++) {
            integral_x[i] += t1 * svect_dotp(qp->grad_shp[i], qp_grad_f);
            integral_y[i] += t1 * svect_dotp(qp->grad_shp[i], qp_grad_g);
        }
    }
}

// ------------------------------------------------------------------------
// integral(integral(constant * dphi_i/d{x,y} * f, yhat, 0,1-xhat), xhat, 0, 1), zhat, -1, 1)
// uses quadrature, cannot directly integrate in sage
// total operation count = alot
void inline integrate_triPrism_dphiXf_quadrature(SVECT *nd, SQUAD *quad, double constant, double *f, double *integral_x, double *integral_y, int quad_order) {
    
    int i, iqp;
    double qp_f = 0.;
    SQUAD_PT *qp = NULL;
    
    sarray_init_dbl(integral_x,NDONPRISM);
    sarray_init_dbl(integral_y,NDONPRISM);
    for (iqp=0; iqp<quad[quad_order].n; iqp++) {
        qp = &(quad[quad_order].pt[iqp]);
        
         // evaluate prism djac at quadrature point :: 50 operations
        qp->djac = elem3d_get_prism_djac(qp->xhat, qp->yhat, qp->zhat, nd);
        
        // evaluate shape gradients at quadrature points :: 415 operations
        elem3d_get_triprism_shape_gradients(qp->xhat, qp->yhat, qp->zhat, nd, qp->grad_shp);
        
        // evaluate function at quadrature point :: 14 operations
        qp_f = SQUAD_get_function(qp, f, NDONPRISM);
        
        // 16 operations
        double t1 = constant * qp_f * qp->djac * qp->w;
        for (i=0; i<NDONPRISM; i++) {
            integral_x[i] += t1 * qp->grad_shp[i].x;
            integral_y[i] += t1 * qp->grad_shp[i].y;
        }
        
        // total operation count = about 500 * 8 quadrature evaluations = 4000
    }
}









