#include "adh.h"

void stensor_init(STENSOR *tensor) {
    tensor->xx = 0.;
    tensor->xy = 0.;
    tensor->yy = 0.;
    tensor->yz = 0.;
    tensor->zz = 0.;
}

void stensor2d_init(STENSOR2D *tensor) {
    tensor->xx = 0.;
    tensor->xy = 0.;
    tensor->yy = 0.;
}

void stensor3d_init(STENSOR3D *tensor) {
    tensor->xx = 0.; tensor->xy = 0.; tensor->xz = 0.;
    tensor->yx = 0.; tensor->yy = 0.; tensor->yz = 0.;
    tensor->zx = 0.; tensor->zy = 0.; tensor->zz = 0.;
}

double stensor3d_max(STENSOR3D tensor) {
    double max = tensor.xx;
    if (tensor.xy > max) max = tensor.xy;
    if (tensor.xz > max) max = tensor.xz;
    if (tensor.yx > max) max = tensor.yx;
    if (tensor.yy > max) max = tensor.yy;
    if (tensor.yz > max) max = tensor.yz;
    if (tensor.zx > max) max = tensor.zx;
    if (tensor.zy > max) max = tensor.zy;
    if (tensor.zz > max) max = tensor.zz;
    return max;
}

void stensor2d_add(STENSOR2D *new_tensor, STENSOR2D d1, STENSOR2D d2) {
    new_tensor->xx = d1.xx + d2.xx;
    new_tensor->xy = d1.xy + d2.xy;
    new_tensor->yy = d1.yy + d2.yy;
}

void stensor2d_add_replace(STENSOR2D *d1, STENSOR2D d2) {
    d1->xx += d2.xx;
    d1->xy += d2.xy;
    d1->yy += d2.yy;
}

void stensor2dai_init(STENSOR2DAI *tensor) {
    tensor->xx = 0.;
    tensor->xy = 0.;
    tensor->yx = 0.;
    tensor->yy = 0.;
}

void stensor_copy(STENSOR *to, STENSOR from) {
    to->xx = from.xx;
    to->yy = from.yy;
    to->zz = from.zz;
    to->xy = from.xy;
    to->xz = from.xz;
    to->yz = from.yz;
}

void stensor2dai_printScreen(STENSOR2DAI t) {
    printf("xx: %20.10e \t yy: %20.10e \t xy: %20.10e \t yx: %20.10e \n",t.xx,t.yy,t.xy,t.yx);
}


