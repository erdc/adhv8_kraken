
#ifndef H_STENSOR_
#define H_STENSOR_

/**************************************/
/**************************************/
typedef struct {
      double xx, xy, xz, yy, yz, zz;    /* the components of the tensor */
} STENSOR;

/**************************************/
/**************************************/
typedef struct{
      double xx, xy, yy;                /* the components of the tensor */
} STENSOR2D;

/**************************************/
/**************************************/
typedef struct{
      double xx, xy, yx, yy;            /* the components of the tensor */
} STENSOR2DAI;

/**************************************/
/**************************************/
typedef struct{
    double xx, xy, xz, yx, yy, yz, zx, zy, zz;
} STENSOR3D;

/*********************************************************/
/* struct methods -------------------------------------- */

void stensor_init(STENSOR *);
void stensor2d_init(STENSOR2D *);
void stensor3d_init(STENSOR3D *);
double stensor3d_max(STENSOR3D tensor);
void stensor2d_add(STENSOR2D *new_tensor, STENSOR2D d1, STENSOR2D d2);
void stensor2d_add_replace(STENSOR2D *d1, STENSOR2D d2);
void stensor2dai_init(STENSOR2DAI *);
void stensor_copy(STENSOR *to, STENSOR from);
void stensor2dai_printScreen(STENSOR2DAI t);

#endif
