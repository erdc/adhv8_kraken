#ifndef H_SGW_H_
#define H_SGW_H_

/***********************************************************/
/***********************************************************/
/***********************************************************/

#ifdef _ADH_GROUNDWATER
typedef struct
{
  double saturation[MAX_NNODES_ON_ELEM3D];   /* Saturation at each node */
  double old_saturation[MAX_NNODES_ON_ELEM3D];   /* Old Saturation at each node */
} ELEMENT_3D_DATA;              /* Data in the element */

/* groundwater flow parameters */
typedef struct {

 //this is just head now, migrated up to sm 
 //double *gw_phead;    /* pressure from the current time step at time t_{n} */
 //double *old_gw_phead;    /* pressure from the previous time step at time t_{n-1} */
 //double *older_gw_phead;  /* pressure from the time step before last t_{n-2} */
 double *predict_gw_phead;    /* pressure predicted by forward Euler */
 double *gw_density;  /* water density */
 double *error; /*nodal error */
 ELEMENT_3D_DATA *elem_3d_data;   /* the nodal data by element */
 SVECT *elem_gw_flux; /* elemental darcy fluxes */
 /* general use array to keep from reallocating after adaption */
 double *darray;
 int *iarray;
 void *vwork;
 int vwork_size;
  
} SGW;

#endif


#endif 
