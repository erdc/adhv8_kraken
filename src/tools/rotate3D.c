#include "global_header.h"


//*****************************************************************************//
//*****************************************************************************//
// CJT :: this routine calculates a 3D rotation matrix about an arbitrary axis
// inputs  ::   angle  - the angle of rotation
//            {u,v,w}  - an axis of rotation
// outputs :: rotationMatrix - a 4X4 3D rotation matrix
//*****************************************************************************//
void getRotationMatrix(double angle, double u, double v, double w, double **rotationMatrix) {
    
    double L = (u*u + v * v + w * w);
    angle = angle * PI / 180.0; //converting to radian value
    double u2 = u * u;
    double v2 = v * v;
    double w2 = w * w;
    
    rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
    rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
    rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
    rotationMatrix[0][3] = 0.0;
    
    rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
    rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
    rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
    rotationMatrix[1][3] = 0.0;
    
    rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
    rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
    rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
    rotationMatrix[2][3] = 0.0;
    
    rotationMatrix[3][0] = 0.0;
    rotationMatrix[3][1] = 0.0;
    rotationMatrix[3][2] = 0.0;
    rotationMatrix[3][3] = 1.0;
}

//*****************************************************************************//
//*****************************************************************************//
// CJT :: this routine calculates a 3D rotation of an array of cartesion points
// inputs  :: {x,y,z}  - an array of input points
//              angle  - the angle of rotation
//            {u,v,w}  - an axis of rotation
// outputs :: {xnew, ynew, xnew} -  an array of the original points rotated
//*****************************************************************************//
void rotatePoints3D(double angle, double u, double v, double w, SVECT *nodes_in, SVECT *nodes_out, int nnodes) {
    
    int i;
    double position_in[4], position_out[4];
    
    double **rotationMatrix;
    rotationMatrix = (double **) tl_alloc(sizeof(double *), nnodes);
    for (i=0; i<nnodes; i++) {
        rotationMatrix[i] = (double *) tl_alloc(sizeof(double), nnodes);
    }
    getRotationMatrix(angle, u, v, w, rotationMatrix);
    
    for (i=0; i<nnodes; i++) {
        position_in[0] = nodes_in[i].x;
        position_in[1] = nodes_in[i].y;
        position_in[2] = nodes_in[i].z;
        position_in[3] = 1.0;
        sarray_vector_matrix_multiply(position_in, rotationMatrix, position_out, 4);
        nodes_out[i].x = position_out[0];
        nodes_out[i].y = position_out[1];
        nodes_out[i].z = position_out[2];
        //node_out[inode].x = position_out[3];
    }
    
    for (i=0; i<nnodes; i++) {
        //printf("x: %20.10f \t y: %20.10f \t z: %20.10f\n",nodes_out[i].x,nodes_out[i].y,nodes_out[i].z);
        rotationMatrix[i] = (double *) tl_free(sizeof(double), nnodes, rotationMatrix[i]);
    }
    rotationMatrix = (double **) tl_free(sizeof(double *), nnodes, rotationMatrix);
}
