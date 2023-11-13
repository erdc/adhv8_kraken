#include "global_header.h"

void test_geometry() {

    int iresult;
    bool result = false;
    SVECT line1_p1, line1_p2, line2_p1, line2_p2, pi;
    SVECT p1, p2, p3;
    
    // ---------------------------------------------------------------------------------------
    // testing  point on line ----------------------------------------------------------------
    printf("-------- TEST :: 2D point on line || ");
    
    // test endpoint
    p1.x = 50.0; p1.y = 50.0;
    p2.x = 0.0; p2.y = 0.0;
    pi.x = 0; pi.y = 0;
    assert(pointOnLine2D(p1,p2,pi) == true);
    
    // test inside
    p1.x = 50.0; p1.y = 50.0;
    p2.x = 0.0; p2.y = 0.0;
    pi.x = 25; pi.y = 25;
    assert(pointOnLine2D(p1,p2,pi) == true);
    
    // test outside
    p1.x = 0.5; p1.y = 0.5;
    p2.x = 0.0; p2.y = 0.0;
    pi.x = 1; pi.y = 1;
    assert(pointOnLine2D(p1,p2,pi) == false);
    
    printf("PASSED\n");
    
    printf("-------- TEST :: 3D point on line || ");
    
    // test endpoint
    p1.x = 50.0; p1.y = 50.0; p1.z = 50.0;
    p2.x = 0.0;  p2.y = 0.0;  p2.z = 0.0;
    pi.x = 0; pi.y = 0; pi.z = 0;
    assert(pointOnLine(p1,p2,pi) == true);
    
    // test inside
    p1.x = 50.0; p1.y = 50.0; p1.z = 50.0;
    p2.x = 0.0;  p2.y = 0.0;  p2.z = 0.0;
    pi.x = 25; pi.y = 25; pi.z = 25;
    assert(pointOnLine(p1,p2,pi) == true);
    
    // test outside
    p1.x = 0.5; p1.y = 0.5; p1.z = 0.5;
    p2.x = 0.0; p2.y = 0.0; p2.z = 0.5;
    pi.x = 1; pi.y = 1; pi.z = 1;
    assert(pointOnLine(p1,p2,pi) == false);
    
    printf("PASSED\n");
    
    // ---------------------------------------------------------------------------------------
    // testing 2D line segment/segment intersections -----------------------------------------
    printf("-------- TEST :: 2D line segment/segment intersection || ");
    
    // lines are colinear and do not overlap
    line1_p1.x = 0.0;  line1_p1.y = 0.0; line1_p1.z = 0.0;
    line1_p2.x = 1.0;  line1_p2.y = 1.0; line1_p2.z = 0.0;
    line2_p1.x =-2.0;  line2_p1.y = -2.0; line2_p1.z = 0.0;
    line2_p2.x =-3.0;  line2_p2.y = -3.0; line2_p2.z = 0.0;
    iresult = intersect2D_SegSeg(line1_p1,line1_p2,line2_p1,line2_p2,&pi);
    //printf("line intersect result: %d || intersection point :: {%f %f %f} \n",result,pi.x,pi.y,pi.z);
    assert(iresult == 0);
    
    // lines intersect but segments do not
    line1_p1.x = 1.0;  line1_p1.y = 1.0; line1_p1.z = 0.0;
    line1_p2.x = 4.0;  line1_p2.y = 4.0; line1_p2.z = 0.0;
    line2_p1.x = 1.0;  line2_p1.y = 8.0; line2_p1.z = 0.0;
    line2_p2.x = 2.0;  line2_p2.y = 3.0; line2_p2.z = 0.0;
    iresult = intersect2D_SegSeg(line1_p1,line1_p2,line2_p1,line2_p2,&pi);
    //printf("line intersect result: %d || intersection point :: {%f %f %f} \n",result,pi.x,pi.y,pi.z);
    assert(iresult == 0);
    
    // lines intersect
    line1_p1.x = 1.0;  line1_p1.y = 1.0; line1_p1.z = 0.0;
    line1_p2.x = 4.0;  line1_p2.y = 4.0; line1_p2.z = 0.0;
    line2_p1.x = 1.0;  line2_p1.y = 3.0; line2_p1.z = 0.0;
    line2_p2.x = 3.0;  line2_p2.y = 1.0; line2_p2.z = 0.0;
    iresult = intersect2D_SegSeg(line1_p1,line1_p2,line2_p1,line2_p2,&pi);
    //printf("line intersect result: %d || intersection point :: {%f %f %f} \n",result,pi.x,pi.y,pi.z);
    assert(iresult == 1 && fabs(pi.x - 2.0) < 1e-12 && fabs(pi.y - 2.0) < 1e-12);
    
    printf("PASSED\n");
    
    // ---------------------------------------------------------------------------------------
    // testing 3D line segment/segment intersections -----------------------------------------
    printf("-------- TEST :: 3D line segment/segment intersection || ");
    
    // lines are colinear and do not overlap
    line1_p1.x = 0.0;  line1_p1.y = 0.0; line1_p1.z = 0.0;
    line1_p2.x = 1.0;  line1_p2.y = 1.0; line1_p2.z = 1.0;
    line2_p1.x =-2.0;  line2_p1.y = -2.0; line2_p1.z = -2.0;
    line2_p2.x =-3.0;  line2_p2.y = -3.0; line2_p2.z = -3.0;
    iresult = intersect3D_SegSeg(line1_p1,line1_p2,line2_p1,line2_p2,&pi);
    //printf("line intersect result: %d || intersection point :: {%f %f %f} \n",result,pi.x,pi.y,pi.z);
    assert(iresult == 0);
    
    // lines are colinear and overlap
    line1_p1.x = 0.0;  line1_p1.y = 0.0; line1_p1.z = 0.0;
    line1_p2.x = 1.0;  line1_p2.y = 1.0; line1_p2.z = 1.0;
    line2_p1.x =-2.0;  line2_p1.y = -2.0; line2_p1.z = -2.0;
    line2_p2.x = 0.5;  line2_p2.y = 0.5; line2_p2.z = 0.5;
    iresult = intersect3D_SegSeg(line1_p1,line1_p2,line2_p1,line2_p2,&pi);
    //printf("line intersect result: %d || intersection point :: {%f %f %f} \n",result,pi.x,pi.y,pi.z);
    assert(iresult == -1);
    
    // lines are colinear and intersect at edge point
    line1_p1.x = 0.0;  line1_p1.y = 0.0; line1_p1.z = 0.0;
    line1_p2.x = 1.0;  line1_p2.y = 1.0; line1_p2.z = 1.0;
    line2_p1.x =-1.0;  line2_p1.y = -1.0; line2_p1.z = -1.0;
    line2_p2.x = 0;  line2_p2.y = 0; line2_p2.z = 0;
    iresult = intersect3D_SegSeg(line1_p1,line1_p2,line2_p1,line2_p2,&pi);
    //printf("line intersect result: %d || intersection point :: {%f %f %f} \n",result,pi.x,pi.y,pi.z);
    assert(iresult == -1);
    
    // lines are parallel
    line1_p1.x = 0.0;  line1_p1.y = 0.0; line1_p1.z = 0.0;
    line1_p2.x = 1.0;  line1_p2.y = 1.0; line1_p2.z = 1.0;
    line2_p1.x = 1.0;  line2_p1.y = 0.0; line2_p1.z = 0.0;
    line2_p2.x = 2.0;  line2_p2.y = 1.0; line2_p2.z = 1.0;
    iresult = intersect3D_SegSeg(line1_p1,line1_p2,line2_p1,line2_p2,&pi);
    //printf("line intersect result: %d || intersection point :: {%f %f %f} \n",iresult,pi.x,pi.y,pi.z);
    assert(iresult == 0);
    
    // lines intersect by segments do not!
    line1_p1.x = 0.0;  line1_p1.y = 0.0; line1_p1.z = 0.0;
    line1_p2.x = 1.0;  line1_p2.y = 0.0; line1_p2.z = 0.0;
    line2_p1.x = 1.5;  line2_p1.y = -0.5; line2_p1.z = 0.0;
    line2_p2.x = 1.5;  line2_p2.y =  0.5; line2_p2.z = 0.0;
    iresult = intersect3D_SegSeg(line1_p1,line1_p2,line2_p1,line2_p2,&pi);
    //printf("line intersect result: %d || intersection point :: {%f %f %f} \n",result,pi.x,pi.y,pi.z);
    assert(iresult == 0);
    
    // line intersects
    line1_p1.x = 0.0;  line1_p1.y = 0.0; line1_p1.z = 0.0;
    line1_p2.x = 1.0;  line1_p2.y = 0.0; line1_p2.z = 0.0;
    line2_p1.x = 0.5;  line2_p1.y = -0.5; line2_p1.z = 0.0;
    line2_p2.x = 0.5;  line2_p2.y =  0.5; line2_p2.z = 0.0;
    iresult = intersect3D_SegSeg(line1_p1,line1_p2,line2_p1,line2_p2,&pi);
    //printf("line intersect result: %d || intersection point :: {%f %f %f} \n",result,pi.x,pi.y,pi.z);
    assert(iresult == 1 && fabs(pi.x - 0.5) < 1e-12 && fabs(pi.y - 0.0) < 1e-12 && fabs(pi.z - 0.0) < 1e-12);
    
    // line intersects
    line1_p1.x = 5.0;  line1_p1.y = 5.0; line1_p1.z = 4.0;
    line1_p2.x = 10.0;  line1_p2.y = 10.0; line1_p2.z = 6.0;
    line2_p1.x = 5.0;  line2_p1.y =  5.0; line2_p1.z = 5.0;
    line2_p2.x = 10.0;  line2_p2.y =  10.0; line2_p2.z = 3.0;
    iresult = intersect3D_SegSeg(line1_p1,line1_p2,line2_p1,line2_p2,&pi);
    //printf("line intersect result: %d || intersection point :: {%f %f %f} \n",result,pi.x,pi.y,pi.z);
    assert(iresult == 1 && fabs(pi.x - 6.25) < 1e-12 && fabs(pi.y - 6.25) < 1e-12 && fabs(pi.z - 4.5) < 1e-12);
    
    line1_p1.x = 6.0;  line1_p1.y = 8.0; line1_p1.z = 4.0;
    line1_p2.x = 12.0;  line1_p2.y = 15.0; line1_p2.z = 4.0;
    line2_p1.x = 6.0;  line2_p1.y =  8.0; line2_p1.z = 2.0;
    line2_p2.x = 12.0;  line2_p2.y =  15.0; line2_p2.z = 6.0;
    iresult = intersect3D_SegSeg(line1_p1,line1_p2,line2_p1,line2_p2,&pi);
    //printf("line intersect result: %d || intersection point :: {%f %f %f} \n",result,pi.x,pi.y,pi.z);
    assert(iresult == 1 && fabs(pi.x - 9) < 1e-12 && fabs(pi.y - (23./2.)) < 1e-12 && fabs(pi.z - 4) < 1e-12);
    
    printf("PASSED\n");
    
    // ---------------------------------------------------------------------------------------
    // testing line segment through face intersection ----------------------------------------
    printf("-------- TEST :: 3D line segment/triangle intersection || ");
    
    // intersection
    line1_p1.x = 0.5;  line1_p1.y = 0.25;  line1_p1.z = 0.0;
    line1_p2.x = 0.5;  line1_p2.y = 0.25;  line1_p2.z = 1.0;
    p1.x = 0.0; p1.y = 0.0; p1.z = 0.0;
    p2.x = 1.0; p2.y = 1.0; p2.z = 1.0;
    p3.x = 1.0; p3.y = 0.0; p3.z = 0.0;
    iresult = intersect3D_SegTriangle(line1_p1,line1_p2,p1,p2,p3,&pi);
    //printf("line/face intersect result: %d || intersection point :: {%f %f %f} \n",iresult,pi.x,pi.y,pi.z);
    assert(iresult == 1 && fabs(pi.x - 0.5) < 1e-12 && fabs(pi.y - 0.25) < 1e-12 && fabs(pi.z - 0.25) < 1e-12);
    
    
    // intersection on an edge
    line1_p1.x = 0.5;  line1_p1.y = 0.5;  line1_p1.z = 0.0;
    line1_p2.x = 0.5;  line1_p2.y = 0.5;  line1_p2.z = 1.0;
    p1.x = 0.0; p1.y = 0.0; p1.z = 0.0;
    p2.x = 1.0; p2.y = 1.0; p2.z = 1.0;
    p3.x = 1.0; p3.y = 0.0; p3.z = 0.0;
    iresult = intersect3D_SegTriangle(line1_p1,line1_p2,p1,p2,p3,&pi);
    //printf("line/face intersect result: %d || intersection point :: {%f %f %f} \n",iresult,pi.x,pi.y,pi.z);
    assert(iresult == 1 && fabs(pi.x - 0.5) < 1e-12 && fabs(pi.y - 0.5) < 1e-12 && fabs(pi.z - 0.5) < 1e-12);
    
    // intersection at a node
    line1_p1.x = 0.0;  line1_p1.y = 0.0;  line1_p1.z = 0.0;
    line1_p2.x = 0.5;  line1_p2.y = 0.5;  line1_p2.z = 1.0;
    p1.x = 0.0; p1.y = 0.0; p1.z = 0.0;
    p2.x = 1.0; p2.y = 1.0; p2.z = 1.0;
    p3.x = 1.0; p3.y = 0.0; p3.z = 0.0;
    iresult = intersect3D_SegTriangle(line1_p1,line1_p2,p1,p2,p3,&pi);
    //printf("line/face intersect result: %d || intersection point :: {%f %f %f} \n",iresult,pi.x,pi.y,pi.z);
    assert(iresult == 1 && fabs(pi.x - 0) < 1e-12 && fabs(pi.y - 0) < 1e-12 && fabs(pi.z - 0) < 1e-12);
    
    // colinear
    // intersection at a node
    line1_p1.x = 0.0;  line1_p1.y = 0.0;  line1_p1.z = 0.0;
    line1_p2.x = 0.5;  line1_p2.y = 0.5;  line1_p2.z = 0.5;
    p1.x = 0.0; p1.y = 0.0; p1.z = 0.0;
    p2.x = 1.0; p2.y = 1.0; p2.z = 1.0;
    p3.x = 1.0; p3.y = 0.0; p3.z = 0.0;
    iresult = intersect3D_SegTriangle(line1_p1,line1_p2,p1,p2,p3,&pi);
    //printf("line/face intersect result: %d || intersection point :: {%f %f %f} \n",iresult,pi.x,pi.y,pi.z);
    assert(iresult == 2);
    
    // intersection on an edge
    line1_p1.x = 1.0;  line1_p1.y = 0.5;  line1_p1.z = 0.0;
    line1_p2.x = 1.0;  line1_p2.y = 0.5;  line1_p2.z = 0.5;
    p1.x = 0.0; p1.y = 0.0; p1.z = 0.0;
    p2.x = 1.0; p2.y = 1.0; p2.z = 1.0;
    p3.x = 1.0; p3.y = 0.0; p3.z = 0.0;
    iresult = intersect3D_SegTriangle(line1_p1,line1_p2,p1,p2,p3,&pi);
    //printf("line/face intersect result: %d || intersection point :: {%f %f %f} \n",iresult,pi.x,pi.y,pi.z);
    assert(iresult == 1 && fabs(pi.x - 1) < 1e-12 && fabs(pi.y - 0.5) < 1e-12 && fabs(pi.z - 0.5) < 1e-12);
    
    // line intersects but segment doesnt
    line1_p1.x = 0.5;  line1_p1.y = 0.5;  line1_p1.z = 0.0;
    line1_p2.x = 0.5;  line1_p2.y = 0.5;  line1_p2.z = 0.15;
    p1.x = 0.0; p1.y = 0.0; p1.z = 0.0;
    p2.x = 1.0; p2.y = 1.0; p2.z = 1.0;
    p3.x = 1.0; p3.y = 0.0; p3.z = 0.0;
    iresult = intersect3D_SegTriangle(line1_p1,line1_p2,p1,p2,p3,&pi);
    //printf("line/face intersect result: %d || intersection point :: {%f %f %f} \n",iresult,pi.x,pi.y,pi.z);
    assert(iresult == 0);
    
    printf("PASSED\n");
    

    // TEST POINT IN POLYGON
    spolygon_test();
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
bool pointOnLine2D(SVECT p1, SVECT p2, SVECT p) {
    p1.z = 0; p2.z = 0; p.z = 0;
    double AB = sqrt((p2.x-p1.x)*(p2.x-p1.x)+(p2.y-p1.y)*(p2.y-p1.y)+(p2.z-p1.z)*(p2.z-p1.z));
    double AP = sqrt((p.x-p1.x)*(p.x-p1.x)+(p.y-p1.y)*(p.y-p1.y)+(p.z-p1.z)*(p.z-p1.z));
    double PB = sqrt((p2.x-p.x)*(p2.x-p.x)+(p2.y-p.y)*(p2.y-p.y)+(p2.z-p.z)*(p2.z-p.z));
    if(fabs(AB - (AP + PB)) < 1e-12) {
        return true;
    } else return false;
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
bool pointOnLine(SVECT p1, SVECT p2, SVECT p) {
    double AB = sqrt((p2.x-p1.x)*(p2.x-p1.x)+(p2.y-p1.y)*(p2.y-p1.y)+(p2.z-p1.z)*(p2.z-p1.z));
    double AP = sqrt((p.x-p1.x)*(p.x-p1.x)+(p.y-p1.y)*(p.y-p1.y)+(p.z-p1.z)*(p.z-p1.z));
    double PB = sqrt((p2.x-p.x)*(p2.x-p.x)+(p2.y-p.y)*(p2.y-p.y)+(p2.z-p.z)*(p2.z-p.z));
    if(fabs(AB - (AP + PB)) < 1e-12) {
        return true;
    } else return false;
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
// intersect2D_SegSeg(): find the 2D intersection of two lines *segments*
//    Input:  line 1 points (line1_p1, line1_p2) and line 2 points (line2_p1, line2_p2)
//    Output: *x,*y = intersection point (when it exists)
//    Return:  0 = no intersection
//             1 = intersection
//            -1 = line segments are colinear
int intersect2D_SegSeg(SVECT p1,SVECT p2,SVECT p3,SVECT p4, SVECT *pi) {
    
    pi->x = -99999999.0; pi->y = -99999999.0; pi->z = -99999999.0;
    
    // Line a1x + b1y = c1
    double a1 = p2.y - p1.y;
    double b1 = p1.x - p2.x;
    double c1 = a1 * p1.x + b1 * p1.y;
    
    // Line a2x + b2y = c2
    double a2 = p4.y - p3.y;
    double b2 = p3.x - p4.x;
    double c2 = a2 * p3.x + b2 * p3.y;
    
    double determinant = a1*b2 - a2*b1;
    
    if (fabs(determinant)<1e-12) {
        // krs :: switch pointOnLine --> pointOnLine2D
        if (pointOnLine2D(p3,p4,p1) || pointOnLine2D(p3,p4,p2)) {
            return -1; // segments are colinear
        }
        else {
            //printf("here\n");
            return 0; // segments are parallel or colinear and disjoint
        }
    } else {
        pi->x = (b2*c1 - b1*c2)/determinant;
        pi->y = (a1*c2 - a2*c1)/determinant;
        pi->z = p1.z; // doesn't matter
        
        if (pointOnLine2D(p1,p2,*pi) && pointOnLine2D(p3,p4,*pi)) {
            return 1;
        } else {
            //printf("here || p1: %f %f %f || p2: %f %f %f || nd1: %f %f %f || nd2: %f %f %f || pi: %f %f %f\n",p1.x,p1.y,p1.z,p2.x,p2.y,p2.z,p3.x,p3.y,p3.z,p4.x,p4.y,p4.z,pi->x,pi->y,pi->z);
            return 0;
        }
    }
    
    return 0;
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
//   Calculate the line segment PaPb that is the shortest route between
//   two lines P1P2 and P3P4. Calculate also the values of mua and mub where
//      Pa = P1 + mua (P2 - P1)
//      Pb = P3 + mub (P4 - P3)
//    Return:  0 = no intersection
//             1 = intersection
//            -1 = line segments are colinear and overlap
//            -2 = one of the segments is just a point!
int intersect3D_SegSeg(SVECT p1,SVECT p2,SVECT p3,SVECT p4, SVECT *pi) {
    SVECT p13,p43,p21;
    double d1343,d4321,d1321,d4343,d2121;
    double numer,denom;
    double mua, mub;
    
    pi->x = -99999999.0; pi->y = -99999999.0; pi->z = -99999999.0;
    
    p13.x = p1.x - p3.x;
    p13.y = p1.y - p3.y;
    p13.z = p1.z - p3.z;
    
    p43.x = p4.x - p3.x;
    p43.y = p4.y - p3.y;
    p43.z = p4.z - p3.z;
    if (fabs(p43.x) < 1e-12 && fabs(p43.y) < 1e-12 && fabs(p43.z) < 1e-12)
        return(-2); // 2nd segment is just a point
    p21.x = p2.x - p1.x;
    p21.y = p2.y - p1.y;
    p21.z = p2.z - p1.z;
    if (fabs(p21.x) < 1e-12 && fabs(p21.y) < 1e-12 && fabs(p21.z) < 1e-12)
        return(-2); // 1rst segment is just a point
    
    d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
    d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
    d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
    d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
    d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;
    
    denom = d2121 * d4343 - d4321 * d4321;
    if (fabs(denom) < 1e-12) {
        if (pointOnLine(p1,p2,p3) || pointOnLine(p1,p2,p4)) return -1; // segments are colinear
        else return 0; // segments are parallel
    }
    numer = d1343 * d4321 - d1321 * d4343;
    
    mua = numer / denom;
    mub = (d1343 + d4321 * (mua)) / d4343;
    
    SVECT pa = svect_add(p1,svect_scale(p21,mua)); // pa = p1 * mua * p21
    SVECT pb = svect_add(p3,svect_scale(p43,mub)); // pb = p3 * mua * p43
    
    if (svect_mag(svect_subtract(pa,pb)) < 1e-12) {
        // 2 lines intersect, now check to see if point is on both segments
        //printf("%f %f %f %d %d\n",pa.x,pa.y,pa.z,pointOnLine(p1,p2,pa),pointOnLine(p2,p3,pa));
        if (pointOnLine(p1,p2,pa) && pointOnLine(p3,p4,pa)) {
            *pi = pa;
            return 1;
        }
    }
    
    return 0;
}


// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
// intersect3D_PointTriangle(): find the 3D intersection of a point with a triangle
//    Input:  a triangles 3 nodes (tri_node0,1,2) and a point
//    Return:  0 = intersection
//             1 = no intersection
int isPointTriangle_3D(SVECT point, SVECT triV1, SVECT triV2, SVECT triV3) {
    //To find out if the point is inside the triangle, we will first find the area
    //of the entire triangle. After we find the area of the entire triangle, we will
    //use the point as a vertex, and create 3 more triangles using the three vertices
    //which make up the first triangle. We will find the area of the three triangles we
    //make using the point, then add the three triangle areas together. If the area
    //of the three triangles added up is the same as the first triangles area, the point
    //is inside the triangle. If the area of the three triangles added up is greater than
    //the area of the first triangle, the point is outside the triangle.

    //Find area of first triangle
    double distX = triV1.x - triV2.x;
    double distY = triV1.y - triV2.y;
    double distZ = triV1.z - triV2.z;

    double edgeLength1 = sqrt(distX*distX + distY*distY + distZ*distZ);

    distX = triV1.x - triV3.x;
    distY = triV1.y - triV3.y;
    distZ = triV1.z - triV3.z;

    double edgeLength2 = sqrt(distX*distX + distY*distY + distZ*distZ);

    distX = triV2.x - triV3.x;
    distY = triV2.y - triV3.y;
    distZ = triV2.z - triV3.z;

    double edgeLength3 = sqrt(distX*distX + distY*distY + distZ*distZ);

    double s = (edgeLength1 + edgeLength2 + edgeLength3)/2.0;

    double mainTriArea = sqrt(s*(s-edgeLength1)*(s-edgeLength2)*(s-edgeLength3));

    //Find areas of the three triangles created with the point

    double smallTriArea[3] = {0.0, 0.0, 0.0};
    SVECT triVert[4];
    triVert[0] = triV1;
    triVert[1] = triV2;
    triVert[2] = triV3;
    triVert[3] = triV1; //When i=2, i+1 will be triV1

    //Find 3 triangle areas using the plane intersecting point
    for(int i = 0; i < 3; i++)
    {
        distX = point.x - triVert[i].x;
        distY = point.y - triVert[i].y;
        distZ = point.z - triVert[i].z;

        edgeLength1 = sqrt(distX*distX + distY*distY + distZ*distZ);

        distX = point.x - triVert[i+1].x;
        distY = point.y - triVert[i+1].y;
        distZ = point.z - triVert[i+1].z;

        edgeLength2 = sqrt(distX*distX + distY*distY + distZ*distZ);

        distX = triVert[i].x - triVert[i+1].x;
        distY = triVert[i].y - triVert[i+1].y;
        distZ = triVert[i].z - triVert[i+1].z;

        edgeLength3 = sqrt(distX*distX + distY*distY + distZ*distZ);

        s = (edgeLength1 + edgeLength2 + edgeLength3)/2.0;

        smallTriArea[i] = sqrt(s*(s-edgeLength1)*(s-edgeLength2)*(s-edgeLength3));
    }

    double totalSmallTriArea = smallTriArea[0] + smallTriArea[1] + smallTriArea[2];

    //Compare the three smaller triangle areas with the main triangles area
    //Subtract a small value from totalSmallTriArea to make up for inacuracy
    if(mainTriArea >= (totalSmallTriArea - 1e-6))
    {
        return 1;
    }
    
    return 0; // point not found in element, on element nodes or on element edges
}

// ----------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------
// intersect3D_RayTriangle(): find the 3D intersection of a ray with a triangle
//    Input:  a triangles 3 nodes (tri_node0,1,2), two points of a line segment/ray (p0,p1)
//    Output: *I = intersection point (when it exists)
//    Return: -1 = triangle is degenerate (a segment or point)
//             0 =  disjoint (no intersect)
//             1 =  intersect in unique point I1
//             2 =  are in the same plane
int intersect3D_SegTriangle(SVECT p0, SVECT p1, SVECT tri_node0, SVECT tri_node1, SVECT tri_node2, SVECT *I) {
    SVECT    u, v, n;              // triangle vectors
    SVECT    dir, w0, w;           // ray vectors
    double   r, a, b;              // params to calc ray-plane intersect
    
    // get triangle edge vectors and plane normal
    u = svect_subtract(tri_node1, tri_node0);
    v = svect_subtract(tri_node2, tri_node0);
    n = svect_cross(u, v);
    
    // ray direction vector
    dir = svect_subtract(p1, p0);
    w0 = svect_subtract(p0, tri_node0);
    
    a = -svect_dotp(n,w0);
    b = svect_dotp(n,dir);
    if (fabs(b) < 1e-10) {     // ray is  parallel to triangle plane
        //printf("\n in plane\n");
        //if (a == 0) {                 // ray lies in triangle plane
        if (fabs(a) < 1e-10) { // ray lies in triangle plane}
            // does p1 fall on plane?
            if (isPointTriangle_3D(p1,tri_node0,tri_node1,tri_node2) == 1) {
//                printf("\n The ray lays in the the plane and p1 does intesect\\n");
//                printf("p0: %10.5e %10.5e %10.5e\n",p0.x,p0.y,p0.z);
//                printf("p1: %10.5e %10.5e %10.5e\n",p1.x,p1.y,p1.z);
//                printf("tri_node0: %10.5e %10.5e %10.5e\n",tri_node0.x,tri_node0.y,tri_node0.z);
//                printf("tri_node1: %10.5e %10.5e %10.5e\n",tri_node1.x,tri_node1.y,tri_node1.z);
//                printf("tri_node2: %10.5e %10.5e %10.5e\n",tri_node2.x,tri_node2.y,tri_node2.z);
                return 2; // yes!
            } else {
//                printf("\n The ray lays in the the plane but p1 does not intesect :: fabs(b): %20.15f\n",fabs(b));
//                printf("p0: %10.5e %10.5e %10.5e\n",p0.x,p0.y,p0.z);
//                printf("p1: %10.5e %10.5e %10.5e\n",p1.x,p1.y,p1.z);
//                printf("tri_node0: %10.5e %10.5e %10.5e\n",tri_node0.x,tri_node0.y,tri_node0.z);
//                printf("tri_node1: %10.5e %10.5e %10.5e\n",tri_node1.x,tri_node1.y,tri_node1.z);
//                printf("tri_node2: %10.5e %10.5e %10.5e\n",tri_node2.x,tri_node2.y,tri_node2.z);
                return 0; // no!
            }
        } else {
//            printf("\n The ray outside of triangle plane :: fabs(b): %20.15f\n",fabs(b));
//            printf("p0: %10.5e %10.5e %10.5e\n",p0.x,p0.y,p0.z);
//            printf("p1: %10.5e %10.5e %10.5e\n",p1.x,p1.y,p1.z);
//            printf("tri_node0: %10.5e %10.5e %10.5e\n",tri_node0.x,tri_node0.y,tri_node0.z);
//            printf("tri_node1: %10.5e %10.5e %10.5e\n",tri_node1.x,tri_node1.y,tri_node1.z);
//            printf("tri_node2: %10.5e %10.5e %10.5e\n",tri_node2.x,tri_node2.y,tri_node2.z);
            return 0;              // ray disjoint from plane
        }
    }
    
    // get intersect point of ray with triangle plane
    r = a / b;
    if (r < -1e-12)  {
        // ray goes away from triangle
//        printf("\n The ray goes away from triangle :: r: %20.15f\n",r);
//        printf("p0: %10.5e %10.5e %10.5e\n",p0.x,p0.y,p0.z);
//        printf("p1: %10.5e %10.5e %10.5e\n",p1.x,p1.y,p1.z);
//        printf("tri_node0: %10.5e %10.5e %10.5e\n",tri_node0.x,tri_node0.y,tri_node0.z);
//        printf("tri_node1: %10.5e %10.5e %10.5e\n",tri_node1.x,tri_node1.y,tri_node1.z);
//        printf("tri_node2: %10.5e %10.5e %10.5e\n",tri_node2.x,tri_node2.y,tri_node2.z);
        return 0; // => no intersect
    }
    
    // for a segment, also test if (r > 1.0) => no intersect
    I->x = p0.x + r * dir.x;
    I->y = p0.y + r * dir.y;
    I->z = p0.z + r * dir.z;
    
    // is I inside T?
    double uu, uv, vv, wu, wv, D;
    uu = svect_dotp(u,u);
    uv = svect_dotp(u,v);
    vv = svect_dotp(v,v);
    
    w.x = I->x - tri_node0.x;
    w.y = I->y - tri_node0.y;
    w.z = I->z - tri_node0.z;
    
    wu = svect_dotp(w,u);
    wv = svect_dotp(w,v);
    
    D = uv * uv - uu * vv;
    
    // get and test parametric coords
    double s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0 || s > 1.0) {         // I is outside T
//        printf("\n The ray  intersection is outside of T :: r: %20.15f\n",r);
//        printf("p0: %10.5e %10.5e %10.5e\n",p0.x,p0.y,p0.z);
//        printf("p1: %10.5e %10.5e %10.5e\n",p1.x,p1.y,p1.z);
//        printf("tri_node0: %10.5e %10.5e %10.5e\n",tri_node0.x,tri_node0.y,tri_node0.z);
//        printf("tri_node1: %10.5e %10.5e %10.5e\n",tri_node1.x,tri_node1.y,tri_node1.z);
//        printf("tri_node2: %10.5e %10.5e %10.5e\n",tri_node2.x,tri_node2.y,tri_node2.z);
        return 0;
    }
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  { // I is outside T
//        printf("\n The ray intersection is outside of T :: r: %20.15f\n",r);
//        printf("p0: %10.5e %10.5e %10.5e\n",p0.x,p0.y,p0.z);
//        printf("p1: %10.5e %10.5e %10.5e\n",p1.x,p1.y,p1.z);
//        printf("tri_node0: %10.5e %10.5e %10.5e\n",tri_node0.x,tri_node0.y,tri_node0.z);
//        printf("tri_node1: %10.5e %10.5e %10.5e\n",tri_node1.x,tri_node1.y,tri_node1.z);
//        printf("tri_node2: %10.5e %10.5e %10.5e\n",tri_node2.x,tri_node2.y,tri_node2.z);
        return 0;
    }

    
    // for a segment, make sure the intersection point is on the segment
    if (pointOnLine(p0,p1,*I) == false) {
//        printf("\n The ray segment intersection is outside of T :: r: %20.15f\n",r);
//        printf("p0: %10.5e %10.5e %10.5e\n",p0.x,p0.y,p0.z);
//        printf("p1: %10.5e %10.5e %10.5e\n",p1.x,p1.y,p1.z);
//        printf("tri_node0: %10.5e %10.5e %10.5e\n",tri_node0.x,tri_node0.y,tri_node0.z);
//        printf("tri_node1: %10.5e %10.5e %10.5e\n",tri_node1.x,tri_node1.y,tri_node1.z);
//        printf("tri_node2: %10.5e %10.5e %10.5e\n",tri_node2.x,tri_node2.y,tri_node2.z);
        return 0;
    }
    
//    printf("\n The ray segment intersection is inside of T :: r: %20.15f\n",r);
//    printf("p0: %10.5e %10.5e %10.5e\n",p0.x,p0.y,p0.z);
//    printf("p1: %10.5e %10.5e %10.5e\n",p1.x,p1.y,p1.z);
//    printf("tri_node0: %10.5e %10.5e %10.5e\n",tri_node0.x,tri_node0.y,tri_node0.z);
//    printf("tri_node1: %10.5e %10.5e %10.5e\n",tri_node1.x,tri_node1.y,tri_node1.z);
//    printf("tri_node2: %10.5e %10.5e %10.5e\n",tri_node2.x,tri_node2.y,tri_node2.z);
    
    return 1;                       // I is in T
}



