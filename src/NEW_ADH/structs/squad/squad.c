#include "adh.h"

static int DEBUG = OFF;

/***********************************************************/
/***********************************************************/
/***********************************************************/
// on parent element ::
//      - gets quadrature points
//      - gets element gradients at quadrature points
//      - gets element determinate at quadrature points
// the variables are the same for all elements and do not change in time

/***********************************************************/
/***********************************************************/

int squad_segment_alloc_init(SQUAD **quad_seg) {
    // [-1:1]
    
    int iorder, iqp, ndof = 0;
    double a=0., b=0., c=0., d=0., e=0.;
    
    (*quad_seg) = (SQUAD *) tl_alloc(sizeof(SQUAD), MAX_QUAD_ORDER+1);
    SQUAD *qp = (*quad_seg); // alias
    
    qp[0].n = 0;
    qp[0].pt = NULL; // not used
    
    // ***************************************************************
    // triangle ******************************************************
    ndof = 2; // # of basis functions on element
    
    // linear
    qp[1].n = 1;
    squad_pt_alloc_array(&(qp[1].pt),qp[1].n,ndof);
    squad_pt_init_array(qp[1].pt,qp[1].n,ndof);
    qp[1].pt[0].xhat = 0.0;   qp[1].pt[0].w = 2.;
    
    // quadratic
    qp[2].n = 2;
    squad_pt_alloc_array(&(qp[2].pt),qp[2].n,ndof);
    squad_pt_init_array(qp[2].pt,qp[2].n,ndof);
    qp[2].pt[0].xhat = +sqrt(1./3.);  qp[2].pt[0].w = 1.;
    qp[2].pt[1].xhat = -sqrt(1./3.);  qp[2].pt[1].w = 1.;
    
    // cubic
    qp[3].n = 3;
    squad_pt_alloc_array(&(qp[3].pt),qp[3].n,ndof);
    squad_pt_init_array(qp[3].pt,qp[3].n,ndof);
    qp[3].pt[0].xhat = 0.0;           qp[3].pt[0].w =  8./9.;
    qp[3].pt[1].xhat = +sqrt(3./5.);  qp[3].pt[1].w =  5./9.;
    qp[3].pt[2].xhat = -sqrt(3./5.);  qp[3].pt[2].w =  5./9.;
    
    // quartic
    qp[4].n = 4;
    squad_pt_alloc_array(&(qp[4].pt),qp[4].n,ndof);
    squad_pt_init_array(qp[4].pt,qp[4].n,ndof);
    double tp = sqrt( (3./7.) + (2./7.) * sqrt((6./5.)) );
    double tm = sqrt( (3./7.) - (2./7.) * sqrt((6./5.)) );
    qp[4].pt[0].xhat = +tm;  qp[4].pt[0].w = (18. + sqrt(30))/36.;
    qp[4].pt[1].xhat = -tm;  qp[4].pt[1].w = (18. + sqrt(30))/36.;
    qp[4].pt[2].xhat = +tp;  qp[4].pt[2].w = (18. - sqrt(30))/36.;
    qp[4].pt[3].xhat = -tp;  qp[4].pt[3].w = (18. - sqrt(30))/36.;
    
    double xh = 0.0;
    for (iorder=0; iorder<=MAX_QUAD_ORDER; iorder++) {
        for (iqp=0; iqp<qp[iorder].n; iqp++) {
            xh = qp[iorder].pt[iqp].xhat ;
            qp[iorder].pt[iqp].lshape[0] = 0.5 * ( 1 - xh );
            qp[iorder].pt[iqp].lshape[1] = 0.5 * ( 1 + xh );
            
            // CJT NOT SURE ABOUT THIS
            qp[iorder].pt[iqp].lshape_quad[0] = 0.5 * xh * (xh - 1);
            qp[iorder].pt[iqp].lshape_quad[1] = 1 - xh*xh;
            qp[iorder].pt[iqp].lshape_quad[2] = 0.5 * xh * (xh + 1);
        }
    }
    
#ifdef _DEBUG
    tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    return 0;
}

/***********************************************************/
/***********************************************************/

int squad_triangle_alloc_init(SQUAD **quad_tri) {
    
    int iorder, iqp, ndof = 0;
    double a=0., b=0., c=0., d=0., e=0.;
    
    (*quad_tri) = (SQUAD *) tl_alloc(sizeof(SQUAD), MAX_QUAD_ORDER+1);
    SQUAD *qp_tri = (*quad_tri); // alias
    
    qp_tri[0].n = 0;
    qp_tri[0].pt = NULL; // not used
    
    // ***************************************************************
    // triangle ******************************************************
    ndof = 3;
    
    // linear
    qp_tri[1].n = 1;
    squad_pt_alloc_array(&(qp_tri[1].pt),qp_tri[1].n,ndof);
    squad_pt_init_array(qp_tri[1].pt,qp_tri[1].n,ndof);
    qp_tri[1].pt[0].xhat = 1./3.;  qp_tri[1].pt[0].yhat = 1./3.;  qp_tri[1].pt[0].w = 1./2.;
    
    // quadratic
    qp_tri[2].n = 3;
    squad_pt_alloc_array(&(qp_tri[2].pt),qp_tri[2].n,ndof);
    squad_pt_init_array(qp_tri[2].pt,qp_tri[2].n,ndof);
    qp_tri[2].pt[0].xhat = 0.;     qp_tri[2].pt[0].yhat = 1./2.;  qp_tri[2].pt[0].w = 1./6.;
    qp_tri[2].pt[1].xhat = 1./2.;  qp_tri[2].pt[1].yhat = 0.;     qp_tri[2].pt[1].w = 1./6.;
    qp_tri[2].pt[2].xhat = 1./2.;  qp_tri[2].pt[2].yhat = 1./2.;  qp_tri[2].pt[2].w = 1./6.;
    
    // cubic
    qp_tri[3].n = 4;
    squad_pt_alloc_array(&(qp_tri[3].pt),qp_tri[3].n,ndof);
    squad_pt_init_array(qp_tri[3].pt,qp_tri[3].n,ndof);
    qp_tri[3].pt[0].xhat = 1./3.;  qp_tri[3].pt[0].yhat = 1./3.;   qp_tri[3].pt[0].w = -27./96.;
    qp_tri[3].pt[1].xhat = 1./5.;  qp_tri[3].pt[1].yhat = 3./5.;   qp_tri[3].pt[1].w =  25./96.;
    qp_tri[3].pt[2].xhat = 1./5.;  qp_tri[3].pt[2].yhat = 1./5.;   qp_tri[3].pt[2].w =  25./96.;
    qp_tri[3].pt[3].xhat = 3./5.;  qp_tri[3].pt[3].yhat = 1./5.;   qp_tri[3].pt[3].w =  25./96.;
    
    // quartic
    qp_tri[4].n = 7;
    squad_pt_alloc_array(&(qp_tri[4].pt),qp_tri[4].n,ndof);
    squad_pt_init_array(qp_tri[4].pt,qp_tri[4].n,ndof);
    qp_tri[4].pt[0].xhat = 0.;     qp_tri[4].pt[0].yhat = 0.;     qp_tri[4].pt[0].w = 1./40.;
    qp_tri[4].pt[1].xhat = 1./2.;  qp_tri[4].pt[1].yhat = 0.;     qp_tri[4].pt[1].w = 1./15.;
    qp_tri[4].pt[2].xhat = 1.;     qp_tri[4].pt[2].yhat = 0.;     qp_tri[4].pt[2].w = 1./40.;
    qp_tri[4].pt[3].xhat = 1./2.;  qp_tri[4].pt[3].yhat = 1./2.;  qp_tri[4].pt[3].w = 1./15.;
    qp_tri[4].pt[4].xhat = 0.;     qp_tri[4].pt[4].yhat = 1.;     qp_tri[4].pt[4].w = 1./40.;
    qp_tri[4].pt[5].xhat = 0.;     qp_tri[4].pt[5].yhat = 1./2.;  qp_tri[4].pt[5].w = 1./15.;
    qp_tri[4].pt[6].xhat = 1./3.;  qp_tri[4].pt[6].yhat = 1./3.;  qp_tri[4].pt[6].w = 9./40.;
    
    for (iorder=0; iorder<=MAX_QUAD_ORDER; iorder++) {
        for (iqp=0; iqp<qp_tri[iorder].n; iqp++) {
            // calculate local element shape values at quad points
            selem2d_get_triangle_local_shape(qp_tri[iorder].pt[iqp].xhat, qp_tri[iorder].pt[iqp].yhat, qp_tri[iorder].pt[iqp].zhat, qp_tri[iorder].pt[iqp].lshape);
            
            // calculate local element quadratic shape values at quad points
            selem2d_get_triangle_local_shape_quad(qp_tri[iorder].pt[iqp].xhat, qp_tri[iorder].pt[iqp].yhat, qp_tri[iorder].pt[iqp].zhat, qp_tri[iorder].pt[iqp].lshape_quad);
        }
    }
    
#ifdef _DEBUG
    //tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    return 0;
}

/***********************************************************/
/***********************************************************/

int squad_rectangle_alloc_init(SQUAD **quad_rect) {
    
    int iorder, iqp, ndof = 0;
    double a=0., b=0., c=0., d=0., e=0., wa, wb;
    
    (*quad_rect) = (SQUAD *) tl_alloc(sizeof(SQUAD), MAX_QUAD_ORDER+1);
    SQUAD *qp_rect = (*quad_rect); // alias
    
    qp_rect[0].n = 0;
    qp_rect[0].pt = NULL; // not used
    
    // ***************************************************************
    // quadrilateral *************************************************
    ndof = 4;
    
    // linear
    qp_rect[1].n = 1;
    squad_pt_alloc_array(&(qp_rect[1].pt),qp_rect[1].n,ndof);
    squad_pt_init_array(qp_rect[1].pt,qp_rect[1].n,ndof);
    qp_rect[1].pt[0].xhat = 0.;  qp_rect[1].pt[0].yhat = 0.;  qp_rect[1].pt[0].w = 4.;
    
    // quadratic
    qp_rect[2].n = 4;
    squad_pt_alloc_array(&(qp_rect[2].pt),qp_rect[2].n,ndof);
    squad_pt_init_array(qp_rect[2].pt,qp_rect[2].n,ndof);
    a = 1./sqrt(3.);
    qp_rect[2].pt[0].xhat = -a;  qp_rect[2].pt[0].yhat = -a;  qp_rect[2].pt[0].w = 1.;
    qp_rect[2].pt[1].xhat =  a;  qp_rect[2].pt[1].yhat = -a;  qp_rect[2].pt[1].w = 1.;
    qp_rect[2].pt[2].xhat =  a;  qp_rect[2].pt[2].yhat =  a;  qp_rect[2].pt[2].w = 1.;
    qp_rect[2].pt[3].xhat = -a;  qp_rect[2].pt[3].yhat =  a;  qp_rect[2].pt[3].w = 1.;
    
    // cubic
    qp_rect[3].n = 4;
    squad_pt_alloc_array(&(qp_rect[3].pt),qp_rect[3].n,ndof);
    squad_pt_init_array(qp_rect[3].pt,qp_rect[3].n,ndof);
    a = 1./sqrt(3.);
    qp_rect[3].pt[0].xhat = -a;  qp_rect[3].pt[0].yhat = -a;  qp_rect[3].pt[0].w = 1.;
    qp_rect[3].pt[1].xhat =  a;  qp_rect[3].pt[1].yhat = -a;  qp_rect[3].pt[1].w = 1.;
    qp_rect[3].pt[2].xhat =  a;  qp_rect[3].pt[2].yhat =  a;  qp_rect[3].pt[2].w = 1.;
    qp_rect[3].pt[3].xhat = -a;  qp_rect[3].pt[3].yhat =  a;  qp_rect[3].pt[3].w = 1.;
    
    // 4th/5th order
    qp_rect[4].n = 9;
    squad_pt_alloc_array(&(qp_rect[4].pt),qp_rect[4].n,ndof);
    squad_pt_init_array(qp_rect[4].pt,qp_rect[4].n,ndof);
    a = sqrt(3./5.); b = 5./9.; c = 8./9.;
    qp_rect[4].pt[0].xhat =  a;  qp_rect[4].pt[0].yhat =  a;   qp_rect[4].pt[0].w = b*b;
    qp_rect[4].pt[1].xhat =  0;  qp_rect[4].pt[1].yhat =  a;   qp_rect[4].pt[1].w = b*c;
    qp_rect[4].pt[2].xhat = -a;  qp_rect[4].pt[2].yhat =  a;   qp_rect[4].pt[2].w = b*b;
    qp_rect[4].pt[3].xhat =  a;  qp_rect[4].pt[3].yhat =  0;   qp_rect[4].pt[3].w = b*c;
    qp_rect[4].pt[4].xhat =  0;  qp_rect[4].pt[4].yhat =  0;   qp_rect[4].pt[4].w = c*c;
    qp_rect[4].pt[5].xhat = -a;  qp_rect[4].pt[5].yhat =  0;   qp_rect[4].pt[5].w = b*c;
    qp_rect[4].pt[6].xhat =  a;  qp_rect[4].pt[6].yhat = -a;   qp_rect[4].pt[6].w = b*b;
    qp_rect[4].pt[7].xhat =  0;  qp_rect[4].pt[7].yhat = -a;   qp_rect[4].pt[7].w = b*c;
    qp_rect[4].pt[8].xhat = -a;  qp_rect[4].pt[8].yhat = -a;   qp_rect[4].pt[8].w = b*b;

    
    for (iorder=0; iorder<=MAX_QUAD_ORDER; iorder++) {
        for (iqp=0; iqp<qp_rect[iorder].n; iqp++) {
            // calculate local element shape values at quad points
            selem2d_get_quadrilateral_local_shape(qp_rect[iorder].pt[iqp].xhat, qp_rect[iorder].pt[iqp].yhat, qp_rect[iorder].pt[iqp].zhat, qp_rect[iorder].pt[iqp].lshape);
            
            // calculate local element quadratic shape values at quad points
            selem2d_get_quadrilateral_local_shape_quad(qp_rect[iorder].pt[iqp].xhat, qp_rect[iorder].pt[iqp].yhat, qp_rect[iorder].pt[iqp].zhat, qp_rect[iorder].pt[iqp].lshape_quad);
        }
    }
    
#ifdef _DEBUG
    tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    return 0;
}

/***********************************************************/
/***********************************************************/

int squad_tetrahedron_alloc_init(SQUAD **quad_tet) {
    
    int iorder, iqp, ndof = 0;
    double a=0., b=0., c=0., d=0., e=0.;
    
    (*quad_tet) = (SQUAD *) tl_alloc(sizeof(SQUAD), MAX_QUAD_ORDER+1);
    SQUAD *qp_tet = (*quad_tet); // alias
    
    qp_tet[0].n = 0;
    qp_tet[0].pt = NULL; // not used
    
    // ***************************************************************
    // tetrahedron ***************************************************
    ndof = 4;
    
    // linears
    qp_tet[1].n = 1;
    squad_pt_alloc_array(&(qp_tet[1].pt),qp_tet[1].n,ndof);
    squad_pt_init_array(qp_tet[1].pt,qp_tet[1].n,ndof);
    qp_tet[1].pt[0].xhat = 1./4.;  qp_tet[1].pt[0].yhat = 1./4.;  qp_tet[1].pt[0].zhat = 1./4.;  qp_tet[1].pt[0].w = 1./6.;
    
    // quadratic
    qp_tet[2].n = 4;
    squad_pt_alloc_array(&(qp_tet[2].pt),qp_tet[2].n,ndof);
    squad_pt_init_array(qp_tet[2].pt,qp_tet[2].n,ndof);
    b = (5. + 3*sqrt(5.))/20.;
    a = (5. - sqrt(5.))/20.;
    qp_tet[2].pt[0].xhat = a;  qp_tet[2].pt[0].yhat = a;  qp_tet[2].pt[0].zhat = a;  qp_tet[2].pt[0].w = 1./24.;
    qp_tet[2].pt[1].xhat = a;  qp_tet[2].pt[1].yhat = b;  qp_tet[2].pt[1].zhat = a;  qp_tet[2].pt[1].w = 1./24.;
    qp_tet[2].pt[2].xhat = a;  qp_tet[2].pt[2].yhat = a;  qp_tet[2].pt[2].zhat = b;  qp_tet[2].pt[2].w = 1./24.;
    qp_tet[2].pt[3].xhat = b;  qp_tet[2].pt[3].yhat = a;  qp_tet[2].pt[3].zhat = a;  qp_tet[2].pt[3].w = 1./24.;
    
    // cubic
    qp_tet[3].n = 5;
    squad_pt_alloc_array(&(qp_tet[3].pt),qp_tet[3].n,ndof);
    squad_pt_init_array(qp_tet[3].pt,qp_tet[3].n,ndof);
    a = 1./4.;  b = 1./2.; c = 1./6.;
    qp_tet[3].pt[0].xhat = a;  qp_tet[3].pt[0].yhat = a;  qp_tet[3].pt[0].zhat = a;  qp_tet[3].pt[0].w = -4./30.;
    qp_tet[3].pt[1].xhat = b;  qp_tet[3].pt[1].yhat = c;  qp_tet[3].pt[1].zhat = c;  qp_tet[3].pt[1].w = 9./120.;
    qp_tet[3].pt[2].xhat = c;  qp_tet[3].pt[2].yhat = b;  qp_tet[3].pt[2].zhat = c;  qp_tet[3].pt[2].w = 9./120.;
    qp_tet[3].pt[3].xhat = c;  qp_tet[3].pt[3].yhat = c;  qp_tet[3].pt[3].zhat = b;  qp_tet[3].pt[3].w = 9./120.;
    qp_tet[3].pt[4].xhat = c;  qp_tet[3].pt[4].yhat = c;  qp_tet[3].pt[4].zhat = c;  qp_tet[3].pt[4].w = 9./120.;
    
    // quartic
    qp_tet[4].n = 11;
    squad_pt_alloc_array(&(qp_tet[4].pt),qp_tet[4].n,ndof);
    squad_pt_init_array(qp_tet[4].pt,qp_tet[4].n,ndof);
    a = (1. + sqrt((5./14.)))/4.;
    b = (1. - sqrt((5./14.)))/4.;
    c = 1./4.;
    d = 11./14.;
    e = 1./14.;
    qp_tet[4].pt[0].xhat  = c;  qp_tet[4].pt[0].yhat  = c;  qp_tet[4].pt[0].zhat  = c;  qp_tet[4].pt[0].w  = -74./5625.;
    qp_tet[4].pt[1].xhat  = d;  qp_tet[4].pt[1].yhat  = e;  qp_tet[4].pt[1].zhat  = e;  qp_tet[4].pt[1].w  = 343./45000.;
    qp_tet[4].pt[2].xhat  = e;  qp_tet[4].pt[2].yhat  = d;  qp_tet[4].pt[2].zhat  = e;  qp_tet[4].pt[2].w  = 343./45000.;
    qp_tet[4].pt[3].xhat  = e;  qp_tet[4].pt[3].yhat  = e;  qp_tet[4].pt[3].zhat  = d;  qp_tet[4].pt[3].w  = 343./45000.;
    qp_tet[4].pt[4].xhat  = e;  qp_tet[4].pt[4].yhat  = e;  qp_tet[4].pt[4].zhat  = e;  qp_tet[4].pt[4].w  = 343./45000.;
    qp_tet[4].pt[5].xhat  = a;  qp_tet[4].pt[5].yhat  = a;  qp_tet[4].pt[5].zhat  = b;  qp_tet[4].pt[5].w  = 56./2250.;
    qp_tet[4].pt[6].xhat  = a;  qp_tet[4].pt[6].yhat  = b;  qp_tet[4].pt[6].zhat  = a;  qp_tet[4].pt[6].w  = 56./2250.;
    qp_tet[4].pt[7].xhat  = a;  qp_tet[4].pt[7].yhat  = b;  qp_tet[4].pt[7].zhat  = b;  qp_tet[4].pt[7].w  = 56./2250.;
    qp_tet[4].pt[8].xhat  = b;  qp_tet[4].pt[8].yhat  = a;  qp_tet[4].pt[8].zhat  = a;  qp_tet[4].pt[8].w  = 56./2250.;
    qp_tet[4].pt[9].xhat  = b;  qp_tet[4].pt[9].yhat  = a;  qp_tet[4].pt[9].zhat  = b;  qp_tet[4].pt[9].w  = 56./2250.;
    qp_tet[4].pt[10].xhat = b;  qp_tet[4].pt[10].yhat = b;  qp_tet[4].pt[10].zhat = a;  qp_tet[4].pt[10].w = 56./2250.;
    
    for (iorder=0; iorder<=MAX_QUAD_ORDER; iorder++) {
        for (iqp=0; iqp<qp_tet[iorder].n; iqp++) {
            // calculate local element shape values at quad points
            selem3d_get_tet_local_shape(qp_tet[iorder].pt[iqp].xhat, qp_tet[iorder].pt[iqp].yhat, qp_tet[iorder].pt[iqp].zhat, qp_tet[iorder].pt[iqp].lshape);
            
            // calculate local element quadratic shape values at quad points
            selem3d_get_tet_local_shape_quad(qp_tet[iorder].pt[iqp].xhat, qp_tet[iorder].pt[iqp].yhat, qp_tet[iorder].pt[iqp].zhat, qp_tet[iorder].pt[iqp].lshape_quad);
        }
    }
    
#ifdef _DEBUG
    if (DEBUG) tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    return 0;
}

/***********************************************************/
/***********************************************************/

int squad_triprism_alloc_init(SQUAD **quad_prism) {
    
    int iorder, iqp, ndof = 0;
    double three_5 = 3./5.;
    
    (*quad_prism) = (SQUAD *) tl_alloc(sizeof(SQUAD), MAX_QUAD_ORDER+1);
    SQUAD *qp_prism = (*quad_prism); // alias
    
    // # --------------------------------------------------------
    // # line (order = 2*n - 1) ---------------------------------
    // # order = 1
    double xquad_line_1_1 = 0, wquad_line_1_1 = 2.;
    // # order = 2 - 3
    double xquad_line_23_1 = +sqrt(one_3), wquad_line_23_1 = 1.;
    double xquad_line_23_2 = -sqrt(one_3), wquad_line_23_2 = 1.;
    // # order = 4 - 5
    double xquad_line_45_1 = 0,              wquad_line_45_1 = 8./9.;
    double xquad_line_45_2 = +sqrt(three_5), wquad_line_45_2 = 5./9.;
    double xquad_line_45_3 = -sqrt(three_5), wquad_line_45_3 = 5./9.;
    
    // # ---------------------------------------------------------
    // # triangle ------------------------------------------------
    // # order = 1
    double xquad_tri_1_1 = one_3, yquad_tri_1_1 = one_3, wquad_tri_1_1 = one_2;
    // # order = 2
    double xquad_tri_2_1 = one_2, yquad_tri_2_1 = one_2, wquad_tri_2_1 = one_6;
    double xquad_tri_2_2 = 0,     yquad_tri_2_2 = one_2, wquad_tri_2_2 = one_6;
    double xquad_tri_2_3 = one_2, yquad_tri_2_3 = 0,     wquad_tri_2_3 = one_6;
    // # order = 3
    double xquad_tri_3_1 = one_3,   yquad_tri_3_1 = one_3,   wquad_tri_3_1 = -27./96.;
    double xquad_tri_3_2 = one_5,   yquad_tri_3_2 = one_5,   wquad_tri_3_2 =  25./96.;
    double xquad_tri_3_3 = three_5, yquad_tri_3_3 = one_5,   wquad_tri_3_3 =  25./96.;
    double xquad_tri_3_4 = one_5,   yquad_tri_3_4 = three_5, wquad_tri_3_4 =  25./96.;
    // # order = 4
    double a4 = 0.445948490915965;
    double b4 = 0.091576213509771;
    double c4 = 0.111690794839005;
    double d4 = 0.054975871827661;
    double xquad_tri_4_1 = a4,     yquad_tri_4_1 = a4,     wquad_tri_4_1 = c4;
    double xquad_tri_4_2 = 1-2*a4, yquad_tri_4_2 = a4,     wquad_tri_4_2 = c4;
    double xquad_tri_4_3 = a4,     yquad_tri_4_3 = 1-2*a4, wquad_tri_4_3 = c4;
    double xquad_tri_4_4 = b4,     yquad_tri_4_4 = b4,     wquad_tri_4_4 = d4;
    double xquad_tri_4_5 = 1-2*b4, yquad_tri_4_5 = b4,     wquad_tri_4_5 = d4;
    double xquad_tri_4_6 = b4,     yquad_tri_4_6 = 1-2*b4, wquad_tri_4_6 = d4;
    
//    double xquad_tri_4_1 = 0.,     yquad_tri_4_1 = 0.,     wquad_tri_4_1 = 1./40.;
//    double xquad_tri_4_2 = 1./2.,  yquad_tri_4_2 = 0.,     wquad_tri_4_2 = 1./15.;
//    double xquad_tri_4_3 = 1.,     yquad_tri_4_3 = 0.,     wquad_tri_4_3 = 1./40.;
//    double xquad_tri_4_4 = 1./2.,  yquad_tri_4_4 = 1./2.,  wquad_tri_4_4 = 1./15.;
//    double xquad_tri_4_5 = 0.,     yquad_tri_4_5 = 1.,     wquad_tri_4_5 = 1./40.;
//    double xquad_tri_4_6 = 0.,     yquad_tri_4_6 = 1./2.,  wquad_tri_4_6 = 1./15.;
//    double xquad_tri_4_7 = 1./3.,  yquad_tri_4_7 = 1./3.,  wquad_tri_4_7 = 9./40.;
    
    qp_prism[0].n = 0;
    qp_prism[0].pt = NULL; // not used
    
    // ***************************************************************
    // triangular prism **********************************************
    ndof = 6;
    
    // linear
    qp_prism[1].n = 1;
    squad_pt_alloc_array(&(qp_prism[1].pt),qp_prism[1].n,ndof);
    squad_pt_init_array(qp_prism[1].pt,qp_prism[1].n,ndof);
    qp_prism[1].pt[0].xhat = xquad_tri_1_1;  qp_prism[1].pt[0].yhat = yquad_tri_1_1;  qp_prism[1].pt[0].zhat = xquad_line_1_1;  qp_prism[1].pt[0].w = wquad_line_1_1*wquad_tri_1_1;
    
    // quadratic
    qp_prism[2].n = 6;
    squad_pt_alloc_array(&(qp_prism[2].pt),qp_prism[2].n,ndof);
    squad_pt_init_array(qp_prism[2].pt,qp_prism[2].n,ndof);
    qp_prism[2].pt[0].xhat = xquad_tri_2_1;  qp_prism[2].pt[0].yhat = yquad_tri_2_1;  qp_prism[2].pt[0].zhat = xquad_line_23_1;  qp_prism[2].pt[0].w = wquad_line_23_1*wquad_tri_2_1;
    qp_prism[2].pt[1].xhat = xquad_tri_2_2;  qp_prism[2].pt[1].yhat = yquad_tri_2_2;  qp_prism[2].pt[1].zhat = xquad_line_23_1;  qp_prism[2].pt[1].w = wquad_line_23_1*wquad_tri_2_2;
    qp_prism[2].pt[2].xhat = xquad_tri_2_3;  qp_prism[2].pt[2].yhat = yquad_tri_2_3;  qp_prism[2].pt[2].zhat = xquad_line_23_1;  qp_prism[2].pt[2].w = wquad_line_23_1*wquad_tri_2_3;
    
    qp_prism[2].pt[3].xhat = xquad_tri_2_1;  qp_prism[2].pt[3].yhat = yquad_tri_2_1;  qp_prism[2].pt[3].zhat = xquad_line_23_2;  qp_prism[2].pt[3].w = wquad_line_23_1*wquad_tri_2_1;
    qp_prism[2].pt[4].xhat = xquad_tri_2_2;  qp_prism[2].pt[4].yhat = yquad_tri_2_2;  qp_prism[2].pt[4].zhat = xquad_line_23_2;  qp_prism[2].pt[4].w = wquad_line_23_1*wquad_tri_2_2;
    qp_prism[2].pt[5].xhat = xquad_tri_2_3;  qp_prism[2].pt[5].yhat = yquad_tri_2_3;  qp_prism[2].pt[5].zhat = xquad_line_23_2;  qp_prism[2].pt[5].w = wquad_line_23_1*wquad_tri_2_3;
    
    // cubic
    qp_prism[3].n = 8;
    squad_pt_alloc_array(&(qp_prism[3].pt),qp_prism[3].n,ndof);
    squad_pt_init_array(qp_prism[3].pt,qp_prism[3].n,ndof);
    qp_prism[3].pt[0].xhat = xquad_tri_3_1;  qp_prism[3].pt[0].yhat = yquad_tri_3_1;  qp_prism[3].pt[0].zhat = xquad_line_23_1;  qp_prism[3].pt[0].w = wquad_line_23_1*wquad_tri_3_1;
    qp_prism[3].pt[1].xhat = xquad_tri_3_2;  qp_prism[3].pt[1].yhat = yquad_tri_3_2;  qp_prism[3].pt[1].zhat = xquad_line_23_1;  qp_prism[3].pt[1].w = wquad_line_23_1*wquad_tri_3_2;
    qp_prism[3].pt[2].xhat = xquad_tri_3_3;  qp_prism[3].pt[2].yhat = yquad_tri_3_3;  qp_prism[3].pt[2].zhat = xquad_line_23_1;  qp_prism[3].pt[2].w = wquad_line_23_1*wquad_tri_3_3;
    qp_prism[3].pt[3].xhat = xquad_tri_3_4;  qp_prism[3].pt[3].yhat = yquad_tri_3_4;  qp_prism[3].pt[3].zhat = xquad_line_23_1;  qp_prism[3].pt[3].w = wquad_line_23_1*wquad_tri_3_4;
    
    qp_prism[3].pt[4].xhat = xquad_tri_3_1;  qp_prism[3].pt[4].yhat = yquad_tri_3_1;  qp_prism[3].pt[4].zhat = xquad_line_23_2;  qp_prism[3].pt[4].w = wquad_line_23_2*wquad_tri_3_1;
    qp_prism[3].pt[5].xhat = xquad_tri_3_2;  qp_prism[3].pt[5].yhat = yquad_tri_3_2;  qp_prism[3].pt[5].zhat = xquad_line_23_2;  qp_prism[3].pt[5].w = wquad_line_23_2*wquad_tri_3_2;
    qp_prism[3].pt[6].xhat = xquad_tri_3_3;  qp_prism[3].pt[6].yhat = yquad_tri_3_3;  qp_prism[3].pt[6].zhat = xquad_line_23_2;  qp_prism[3].pt[6].w = wquad_line_23_2*wquad_tri_3_3;
    qp_prism[3].pt[7].xhat = xquad_tri_3_4;  qp_prism[3].pt[7].yhat = yquad_tri_3_4;  qp_prism[3].pt[7].zhat = xquad_line_23_2;  qp_prism[3].pt[7].w = wquad_line_23_2*wquad_tri_3_4;
    
    // quartic
    qp_prism[4].n = 18;
    squad_pt_alloc_array(&(qp_prism[4].pt),qp_prism[4].n,ndof);
    squad_pt_init_array(qp_prism[4].pt,qp_prism[4].n,ndof);
    qp_prism[4].pt[0].xhat  = xquad_tri_4_1; qp_prism[4].pt[0].yhat  = yquad_tri_4_1; qp_prism[4].pt[0].zhat  = xquad_line_45_1; qp_prism[4].pt[0].w  = wquad_line_45_1*wquad_tri_4_1;
    qp_prism[4].pt[1].xhat  = xquad_tri_4_2; qp_prism[4].pt[1].yhat  = yquad_tri_4_2; qp_prism[4].pt[1].zhat  = xquad_line_45_1; qp_prism[4].pt[1].w  = wquad_line_45_1*wquad_tri_4_2;
    qp_prism[4].pt[2].xhat  = xquad_tri_4_3; qp_prism[4].pt[2].yhat  = yquad_tri_4_3; qp_prism[4].pt[2].zhat  = xquad_line_45_1; qp_prism[4].pt[2].w  = wquad_line_45_1*wquad_tri_4_3;
    qp_prism[4].pt[3].xhat  = xquad_tri_4_4; qp_prism[4].pt[3].yhat  = yquad_tri_4_4; qp_prism[4].pt[3].zhat  = xquad_line_45_1; qp_prism[4].pt[3].w  = wquad_line_45_1*wquad_tri_4_4;
    qp_prism[4].pt[4].xhat  = xquad_tri_4_5; qp_prism[4].pt[4].yhat  = yquad_tri_4_5; qp_prism[4].pt[4].zhat  = xquad_line_45_1; qp_prism[4].pt[4].w  = wquad_line_45_1*wquad_tri_4_5;
    qp_prism[4].pt[5].xhat  = xquad_tri_4_6; qp_prism[4].pt[5].yhat  = yquad_tri_4_6; qp_prism[4].pt[5].zhat  = xquad_line_45_1; qp_prism[4].pt[5].w  = wquad_line_45_1*wquad_tri_4_6;

    qp_prism[4].pt[6].xhat  = xquad_tri_4_1; qp_prism[4].pt[6].yhat  = yquad_tri_4_1; qp_prism[4].pt[6].zhat  = xquad_line_45_2; qp_prism[4].pt[6].w  = wquad_line_45_2*wquad_tri_4_1;
    qp_prism[4].pt[7].xhat  = xquad_tri_4_2; qp_prism[4].pt[7].yhat  = yquad_tri_4_2; qp_prism[4].pt[7].zhat  = xquad_line_45_2; qp_prism[4].pt[7].w  = wquad_line_45_2*wquad_tri_4_2;
    qp_prism[4].pt[8].xhat  = xquad_tri_4_3; qp_prism[4].pt[8].yhat  = yquad_tri_4_3; qp_prism[4].pt[8].zhat  = xquad_line_45_2; qp_prism[4].pt[8].w  = wquad_line_45_2*wquad_tri_4_3;
    qp_prism[4].pt[9].xhat  = xquad_tri_4_4; qp_prism[4].pt[9].yhat  = yquad_tri_4_4; qp_prism[4].pt[9].zhat  = xquad_line_45_2; qp_prism[4].pt[9].w  = wquad_line_45_2*wquad_tri_4_4;
    qp_prism[4].pt[10].xhat = xquad_tri_4_5; qp_prism[4].pt[10].yhat = yquad_tri_4_5; qp_prism[4].pt[10].zhat = xquad_line_45_2; qp_prism[4].pt[10].w = wquad_line_45_2*wquad_tri_4_5;
    qp_prism[4].pt[11].xhat = xquad_tri_4_6; qp_prism[4].pt[11].yhat = yquad_tri_4_6; qp_prism[4].pt[11].zhat = xquad_line_45_2; qp_prism[4].pt[11].w = wquad_line_45_2*wquad_tri_4_6;
    
    qp_prism[4].pt[12].xhat = xquad_tri_4_1; qp_prism[4].pt[12].yhat = yquad_tri_4_1; qp_prism[4].pt[12].zhat = xquad_line_45_3; qp_prism[4].pt[12].w = wquad_line_45_3*wquad_tri_4_1;
    qp_prism[4].pt[13].xhat = xquad_tri_4_2; qp_prism[4].pt[13].yhat = yquad_tri_4_2; qp_prism[4].pt[13].zhat = xquad_line_45_3; qp_prism[4].pt[13].w = wquad_line_45_3*wquad_tri_4_2;
    qp_prism[4].pt[14].xhat = xquad_tri_4_3; qp_prism[4].pt[14].yhat = yquad_tri_4_3; qp_prism[4].pt[14].zhat = xquad_line_45_3; qp_prism[4].pt[14].w = wquad_line_45_3*wquad_tri_4_3;
    qp_prism[4].pt[15].xhat = xquad_tri_4_4; qp_prism[4].pt[15].yhat = yquad_tri_4_4; qp_prism[4].pt[15].zhat = xquad_line_45_3; qp_prism[4].pt[15].w = wquad_line_45_3*wquad_tri_4_4;
    qp_prism[4].pt[16].xhat = xquad_tri_4_5; qp_prism[4].pt[16].yhat = yquad_tri_4_5; qp_prism[4].pt[16].zhat = xquad_line_45_3; qp_prism[4].pt[16].w = wquad_line_45_3*wquad_tri_4_5;
    qp_prism[4].pt[17].xhat = xquad_tri_4_6; qp_prism[4].pt[17].yhat = yquad_tri_4_6; qp_prism[4].pt[17].zhat = xquad_line_45_3; qp_prism[4].pt[17].w = wquad_line_45_3*wquad_tri_4_6;
    
    
//    qp_prism[4].n = 21;
//    squad_pt_alloc_array(&(qp_prism[4].pt),qp_prism[4].n,ndof);
//    squad_pt_init_array(qp_prism[4].pt,qp_prism[4].n,ndof);
//    qp_prism[4].pt[0].xhat  = xquad_tri_4_1; qp_prism[4].pt[0].yhat  = yquad_tri_4_1; qp_prism[4].pt[0].zhat  = xquad_line_45_1; qp_prism[4].pt[0].w  = wquad_line_45_1*wquad_tri_4_1;
//    qp_prism[4].pt[1].xhat  = xquad_tri_4_2; qp_prism[4].pt[1].yhat  = yquad_tri_4_2; qp_prism[4].pt[1].zhat  = xquad_line_45_1; qp_prism[4].pt[1].w  = wquad_line_45_1*wquad_tri_4_2;
//    qp_prism[4].pt[2].xhat  = xquad_tri_4_3; qp_prism[4].pt[2].yhat  = yquad_tri_4_3; qp_prism[4].pt[2].zhat  = xquad_line_45_1; qp_prism[4].pt[2].w  = wquad_line_45_1*wquad_tri_4_3;
//    qp_prism[4].pt[3].xhat  = xquad_tri_4_4; qp_prism[4].pt[3].yhat  = yquad_tri_4_4; qp_prism[4].pt[3].zhat  = xquad_line_45_1; qp_prism[4].pt[3].w  = wquad_line_45_1*wquad_tri_4_4;
//    qp_prism[4].pt[4].xhat  = xquad_tri_4_5; qp_prism[4].pt[4].yhat  = yquad_tri_4_5; qp_prism[4].pt[4].zhat  = xquad_line_45_1; qp_prism[4].pt[4].w  = wquad_line_45_1*wquad_tri_4_5;
//    qp_prism[4].pt[5].xhat  = xquad_tri_4_6; qp_prism[4].pt[5].yhat  = yquad_tri_4_6; qp_prism[4].pt[5].zhat  = xquad_line_45_1; qp_prism[4].pt[5].w  = wquad_line_45_1*wquad_tri_4_6;
//    qp_prism[4].pt[6].xhat  = xquad_tri_4_7; qp_prism[4].pt[6].yhat  = yquad_tri_4_7; qp_prism[4].pt[6].zhat  = xquad_line_45_1; qp_prism[4].pt[6].w  = wquad_line_45_1*wquad_tri_4_7;
//    
//    qp_prism[4].pt[7].xhat  = xquad_tri_4_1; qp_prism[4].pt[7].yhat  = yquad_tri_4_1; qp_prism[4].pt[7].zhat  = xquad_line_45_2; qp_prism[4].pt[7].w  = wquad_line_45_2*wquad_tri_4_1;
//    qp_prism[4].pt[8].xhat  = xquad_tri_4_2; qp_prism[4].pt[8].yhat  = yquad_tri_4_2; qp_prism[4].pt[8].zhat  = xquad_line_45_2; qp_prism[4].pt[8].w  = wquad_line_45_2*wquad_tri_4_2;
//    qp_prism[4].pt[9].xhat  = xquad_tri_4_3; qp_prism[4].pt[9].yhat  = yquad_tri_4_3; qp_prism[4].pt[9].zhat  = xquad_line_45_2; qp_prism[4].pt[9].w  = wquad_line_45_2*wquad_tri_4_3;
//    qp_prism[4].pt[10].xhat = xquad_tri_4_4; qp_prism[4].pt[10].yhat = yquad_tri_4_4; qp_prism[4].pt[10].zhat = xquad_line_45_2; qp_prism[4].pt[10].w = wquad_line_45_2*wquad_tri_4_4;
//    qp_prism[4].pt[11].xhat = xquad_tri_4_5; qp_prism[4].pt[11].yhat = yquad_tri_4_5; qp_prism[4].pt[11].zhat = xquad_line_45_2; qp_prism[4].pt[11].w = wquad_line_45_2*wquad_tri_4_5;
//    qp_prism[4].pt[12].xhat = xquad_tri_4_6; qp_prism[4].pt[12].yhat = yquad_tri_4_6; qp_prism[4].pt[12].zhat = xquad_line_45_2; qp_prism[4].pt[12].w = wquad_line_45_2*wquad_tri_4_6;
//    qp_prism[4].pt[13].xhat = xquad_tri_4_7; qp_prism[4].pt[13].yhat = yquad_tri_4_7; qp_prism[4].pt[13].zhat = xquad_line_45_2; qp_prism[4].pt[13].w = wquad_line_45_2*wquad_tri_4_7;
//    
//    qp_prism[4].pt[14].xhat = xquad_tri_4_1; qp_prism[4].pt[14].yhat = yquad_tri_4_1; qp_prism[4].pt[14].zhat = xquad_line_45_3; qp_prism[4].pt[14].w = wquad_line_45_3*wquad_tri_4_1;
//    qp_prism[4].pt[15].xhat = xquad_tri_4_2; qp_prism[4].pt[15].yhat = yquad_tri_4_2; qp_prism[4].pt[15].zhat = xquad_line_45_3; qp_prism[4].pt[15].w = wquad_line_45_3*wquad_tri_4_2;
//    qp_prism[4].pt[16].xhat = xquad_tri_4_3; qp_prism[4].pt[16].yhat = yquad_tri_4_3; qp_prism[4].pt[16].zhat = xquad_line_45_3; qp_prism[4].pt[16].w = wquad_line_45_3*wquad_tri_4_3;
//    qp_prism[4].pt[17].xhat = xquad_tri_4_4; qp_prism[4].pt[17].yhat = yquad_tri_4_4; qp_prism[4].pt[17].zhat = xquad_line_45_3; qp_prism[4].pt[17].w = wquad_line_45_3*wquad_tri_4_4;
//    qp_prism[4].pt[18].xhat = xquad_tri_4_5; qp_prism[4].pt[18].yhat = yquad_tri_4_5; qp_prism[4].pt[18].zhat = xquad_line_45_3; qp_prism[4].pt[18].w = wquad_line_45_3*wquad_tri_4_5;
//    qp_prism[4].pt[19].xhat = xquad_tri_4_6; qp_prism[4].pt[19].yhat = yquad_tri_4_6; qp_prism[4].pt[19].zhat = xquad_line_45_3; qp_prism[4].pt[19].w = wquad_line_45_3*wquad_tri_4_6;
//    qp_prism[4].pt[20].xhat = xquad_tri_4_7; qp_prism[4].pt[20].yhat = yquad_tri_4_7; qp_prism[4].pt[20].zhat = xquad_line_45_3; qp_prism[4].pt[20].w = wquad_line_45_3*wquad_tri_4_7;
    
    for (iorder=0; iorder<=MAX_QUAD_ORDER; iorder++) {
        for (iqp=0; iqp<qp_prism[iorder].n; iqp++) {
            // calculate local element shape values at quad points
            selem3d_get_triprism_local_shape(qp_prism[iorder].pt[iqp].xhat, qp_prism[iorder].pt[iqp].yhat, qp_prism[iorder].pt[iqp].zhat, qp_prism[iorder].pt[iqp].lshape);
            
            // calculate local element quadratic shape values at quad points
            selem3d_get_triprism_local_shape_quad(qp_prism[iorder].pt[iqp].xhat, qp_prism[iorder].pt[iqp].yhat, qp_prism[iorder].pt[iqp].zhat, qp_prism[iorder].pt[iqp].lshape_quad);
        }
    }
    
#ifdef _DEBUG
    tl_check_all_pickets(__FILE__,__LINE__);
#endif
    
    return 0;
}

/***********************************************************/
/***********************************************************/

void squad_free(SQUAD *quad, int ndof) {
    
    int i;
    for (i=1; i<=MAX_QUAD_ORDER; i++) {
        squad_pt_free_array(quad[i].pt,quad[i].n,ndof);
    }
    quad = (SQUAD *) tl_free(sizeof(SQUAD), MAX_QUAD_ORDER+1, quad);
}

/***********************************************************/
/***********************************************************/
void squad_pt_alloc_array(SQUAD_PT **quad_pts, int nqp, int ndof) {
    (*quad_pts) = (SQUAD_PT *) tl_alloc(sizeof(SQUAD_PT), nqp);
    int iqp;
    for (iqp=0; iqp<nqp; iqp++) {
        (*quad_pts)[iqp].lshape = (double *) tl_alloc(sizeof(double), ndof);
        (*quad_pts)[iqp].lshape_quad = (double *) tl_alloc(sizeof(double), 15); // for now, max it
        (*quad_pts)[iqp].lgrad_shp = (SVECT *) tl_alloc(sizeof(SVECT), ndof);
        (*quad_pts)[iqp].grad_shp  = (SVECT *) tl_alloc(sizeof(SVECT), ndof);
    }
}

/***********************************************************/
/***********************************************************/
void squad_pt_free_array(SQUAD_PT *quad_pts, int nqp, int ndof) {
    int iqp;
    for (iqp=0; iqp<nqp; iqp++) {
        quad_pts[iqp].lshape = (double *) tl_free(sizeof(double), ndof, quad_pts[iqp].lshape);
        quad_pts[iqp].lshape_quad = (double *) tl_free(sizeof(double), 15, quad_pts[iqp].lshape_quad); // for now max it
        quad_pts[iqp].lgrad_shp = (SVECT *) tl_free(sizeof(SVECT), ndof, quad_pts[iqp].lgrad_shp);
        quad_pts[iqp].grad_shp  = (SVECT *) tl_free(sizeof(SVECT), ndof, quad_pts[iqp].grad_shp);
    }
    quad_pts = (SQUAD_PT *) tl_free(sizeof(SQUAD_PT), nqp, quad_pts);
}

/***********************************************************/
/***********************************************************/
void squad_pt_init(SQUAD_PT *qp, int ndof) {
    qp->xhat = 0.;
    qp->yhat = 0.;
    qp->zhat = 0.;
    qp->djac = 0.;
    qp->djac2d = 0.;
    qp->djac2d_3d = 0.;
    qp->djac3d_fixed = 0.;
    svect_init_array(qp->lgrad_shp, ndof);
    svect_init_array(qp->grad_shp, ndof);
}

/***********************************************************/
/***********************************************************/
void squad_pt_init_array(SQUAD_PT *qp, int nqp, int ndof) {
    int iquad;
    for (iquad=0; iquad<nqp; iquad++) {
        squad_pt_init(&(qp[iquad]), ndof);
    }
}

/***********************************************************/
/***********************************************************/
void squad_printScreen(SQUAD *qp, int flag, int ndof, char *geometry) {
    int iorder, iquad, idof;
    printf("%s :: QUADRATURE POINTS *******************\n",geometry);
    for (iorder=0; iorder<=MAX_QUAD_ORDER; iorder++) {
        printf("** order: %d\n",iorder);
        for (iquad=0; iquad<qp[iorder].n; iquad++) {
            printf("---- point %d: xhat: %10.5f yhat: %10.5f zhat: %10.5f w: %10.5f djac: %20.10e\n",
                   iquad,
                   qp[iorder].pt[iquad].xhat,
                   qp[iorder].pt[iquad].yhat,
                   qp[iorder].pt[iquad].zhat,
                   qp[iorder].pt[iquad].w,
                   qp[iorder].pt[iquad].djac);
            if (flag) {
                for (idof=0; idof<ndof; idof++) {
                    printf("      dof %d: local shape: %10.5f || local shape gradients {dx,dy,dz}: %10.5f %10.5f %10.5f || shape gradients {dx,dy,dz}: %10.5f %10.5f %10.5f \n",
                           idof,
                           qp[iorder].pt[iquad].lshape[idof],
                           qp[iorder].pt[iquad].lgrad_shp[idof].x,
                           qp[iorder].pt[iquad].lgrad_shp[idof].y,
                           qp[iorder].pt[iquad].lgrad_shp[idof].z,
                           qp[iorder].pt[iquad].grad_shp[idof].x,
                           qp[iorder].pt[iquad].grad_shp[idof].y,
                           qp[iorder].pt[iquad].grad_shp[idof].z);
                }
            }
        }
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/
// quadrature evaluations of fields

SVECT squad_get_svect(SQUAD_PT *qp, SVECT *v, int nnodes) {
    SVECT v_quad; svect_init(&v_quad);
    int inode = 0.;
    for (inode=0; inode<nnodes; inode++) {
        v_quad.x += v[inode].x * qp->lshape[inode];
        v_quad.y += v[inode].y * qp->lshape[inode];
        v_quad.z += v[inode].z * qp->lshape[inode];
    }
    return v_quad;
}

SVECT2D squad_get_svect2d(SQUAD_PT *qp, SVECT2D *v, int nnodes) {
    SVECT2D v_quad; svect2d_init(&v_quad);
    int inode = 0.;
    for (inode=0; inode<nnodes; inode++) {
        v_quad.x += v[inode].x * qp->lshape[inode];
        v_quad.y += v[inode].y * qp->lshape[inode];
    }
    return v_quad;
}

double squad_get_function(SQUAD_PT *qp, double *f, int nnodes) {
    int inode = 0.;
    double fout = 0.;
    for (inode=0; inode<nnodes; inode++) {
        fout += f[inode] * qp->lshape[inode];
    }
    return fout;
}

double squad_get_quadFunction(SQUAD_PT *qp, double *f, int nnodes) {
    int inode = 0.;
    double fout = 0.;
    for (inode=0; inode<nnodes; inode++) {
        fout += f[inode] * qp->lshape_quad[inode];
    }
    return fout;
}
