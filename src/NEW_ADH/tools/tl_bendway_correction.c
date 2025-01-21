/*  tl_bendway_correction.c */

/*  tl_bendway_correction is called by fe_sw2_elem_resid.c to calculate the streamwise acceleration */
/*  vector for the purpuse of determining the bendway correction factor.   */
/*  The source of the equation is from Bob Benard's "Depth-Average Numerical Modeling   */
/*  for Curved Channels" report, dated September 1991.                                  */

/*  READ ME IS REQUIRED !!!!!!!!!!!!!!  */
/*  READ ME IS REQUIRED !!!!!!!!!!!!!!  */
/*  READ ME IS REQUIRED !!!!!!!!!!!!!!  */
/*  READ ME IS REQUIRED !!!!!!!!!!!!!!  */
/*  This function is based on tau_s being constant, which in reality is not true.  
   If bendway correction functionality is not what is expected, it may be due to 
   this assumption in coding.  Review the code and make adjustments according to 
   equation 21 in Bob Bernard's "Depth-Average Numerical Modeling for Curved Channels" 
   report, dated September 1991. This, in turn, affects equation 10.  In this case, the 
   variable h is averaged, but really varies-as well as vorticity.  */

#include "adh.h"

SVECT2D tl_bendway_correction(
  SVECT2D * grad_shp,
  SVECT2D * elem_vel,
  double *elem_head,
  double *elem_c,
  double prop1,
  double prop2,
  double drying_lower_limit,
  double roughness,
  double fluid_density,
  int routine_flag,		/* flag indicator indicating what routine is call  bwc  */
  int ie
)
{
  double radius;		/* radius of curvature term  */
  double abs_u_mag;		/* absolute magnitude of the "u" velocity vector */
  double u;			/* x component of u_mag */
  double v;			/* y component of u_mag */
  double u_x, u_y, v_x, v_y;	/* spatial derivitives of each component */
  double c_f;			/* friction factor  */
  double h;			/* average depth    */
  double tau_s;			/* average tau sub s variable */
  SVECT2D stream_accel_vect;	/* Stream wise acceleration vector */
  double normal_derive;		/* the normal derivative of the depth-integrated shear
				   stress created by secondary flow  */
  double denominator;
 // double elem_c[NDPRFC];	/* the local old vorticity in the form of concentrations */
  double avg_c;			/* average vorticity  */
  double c_limit;		/* the limiting vorticity */
  double tau[NDPRFC];
  double cf[NDPRFC];
  double h_vor[NDPRFC];

  /* ******************************************************************   */
  /*  The average x and y components of the velocity vector             */

  u = (elem_vel[0].x + elem_vel[1].x + elem_vel[2].x) / 3;
  v = (elem_vel[0].y + elem_vel[1].y + elem_vel[2].y) / 3;

  /* ******************************************************************** */
  /* The derivatives with respect to x and y components of the velocity vector  */

  u_x =
    elem_vel[0].x * grad_shp[0].x + elem_vel[1].x * grad_shp[1].x +
    elem_vel[2].x * grad_shp[2].x;
  u_y =
    elem_vel[0].x * grad_shp[0].y + elem_vel[1].x * grad_shp[1].y +
    elem_vel[2].x * grad_shp[2].y;
  v_x =
    elem_vel[0].y * grad_shp[0].x + elem_vel[1].y * grad_shp[1].x +
    elem_vel[2].y * grad_shp[2].x;
  v_y =
    elem_vel[0].y * grad_shp[0].y + elem_vel[1].y * grad_shp[1].y +
    elem_vel[2].y * grad_shp[2].y;

  /* ***************************************************************** */
  /* Calculation of the radius of curvature                            */

  abs_u_mag = sqrt(u * u + v * v);

  /* if the velocity magnitude is near zero then  */
  /*   (0) The streamwise acceleration of the depth-average flow is zero  */
  /*   (1) The decay and source for vorticity are zero too    */

  if (abs_u_mag < pow(NOT_QUITE_SMALL,0.333))
    {
      stream_accel_vect.x = 0;
      stream_accel_vect.y = 0;

      return (stream_accel_vect);
    }
  denominator = ((u * v * (v_y - u_x)) + (u * u * v_x) - (v * v * u_y));

  if (denominator > -NOT_QUITE_SMALL && denominator < 0.0)
    denominator = -NOT_QUITE_SMALL;
  if (denominator < NOT_QUITE_SMALL && denominator >= 0.0)
    denominator = NOT_QUITE_SMALL;
  radius = (abs_u_mag * abs_u_mag * abs_u_mag) / denominator;

  /* ***************************************************************** */
  /* Calculation of the average depth or h                           */

  /* if any of the depths on a particular element is zero  */
  /*   (0) The acceleration due to secondary effects are zero  */
  /*   (1) The source & decay for vorticity are zero   */

 /* if((elem_head[0] <= 0) || (elem_head[1] <= 0) || (elem_head[2] <= 0))
    {
      stream_accel_vect.x = 0;
      stream_accel_vect.y = 0;
      return (stream_accel_vect);
    } */

   h_vor[0] = elem_head[0];
   if (elem_head[0] < drying_lower_limit) h_vor[0] =  drying_lower_limit;
   h_vor[1] = elem_head[1];
   if (elem_head[1] < drying_lower_limit) h_vor[1] =  drying_lower_limit;
   h_vor[2] = elem_head[2];
   if (elem_head[2] < drying_lower_limit) h_vor[2] =  drying_lower_limit;

  h = (h_vor[0] + h_vor[1] + h_vor[2]) / 3;

 /* now limit the radius to be no smalle than the average depth */

  if (fabs(radius) < h){
   if (radius >= 0) radius = h;
   if (radius < 0)  radius = -h;
  }

  /* ***************************************************************** */
  /* Calculation of C_f                                                */

  c_f = 0.5 * roughness;
  cf[0] = 0.5 * roughness;
  cf[1] = 0.5 * roughness;
  cf[2] = 0.5 * roughness;

  /* ***************************************************************** */
  /* Routine flag = 0 is being called by hydraulic routines while it = 1 is */
  /* called by transport routines.                                     */

  if(routine_flag == 0)
    {

      /* ***************************************************************** */
      /* Calculation of Average Concentration                              */

     //  ELEM2D_GET_LOCAL(concentration[transport.vor], elem_c, elem2d[ie].nodes);

 /* calculate a limiting vorticity, by assimng that the magnitude of the */
 /* near-bed trnsverse velocity should not exceed 2 times the streamwise velocity */

      c_limit = 0.25 * abs_u_mag * sqrt (6. * prop1) /
                (6. * (h + 1.e-8));

      if (fabs(elem_c[0]) > c_limit){
       if (elem_c[0] >= 0) elem_c[0] =  c_limit;
       if (elem_c[0] < 0)  elem_c[0] = -c_limit;
      }
      if (fabs(elem_c[1]) > c_limit){
       if (elem_c[1] >= 0) elem_c[1] =  c_limit;
       if (elem_c[1] < 0)  elem_c[1] = -c_limit;
      }
      if (fabs(elem_c[2]) > c_limit){
       if (elem_c[2] >= 0) elem_c[2] =  c_limit;
       if (elem_c[2] < 0)  elem_c[2] = -c_limit;
      }

      avg_c = (elem_c[0] + elem_c[1] + elem_c[2]) / 3;

      /* ***************************************************************** */
      /* Calculation of Tau_s                                              */
      /* avg_c is the old vorticity   */

      tau_s = fluid_density * h * avg_c * abs_u_mag * sqrt(c_f);

      tau[0] =
	fluid_density * h_vor[0] * elem_c[0] *
	(sqrt(elem_vel[0].x * elem_vel[0].x + elem_vel[0].y * elem_vel[0].y)) *
	(sqrt(cf[0]));
      tau[1] =
	fluid_density * h_vor[1] * elem_c[1] *
	(sqrt(elem_vel[1].x * elem_vel[1].x + elem_vel[1].y * elem_vel[1].y)) *
	(sqrt(cf[1]));
      tau[2] =
	fluid_density * h_vor[2] * elem_c[2] *
	(sqrt(elem_vel[2].x * elem_vel[2].x + elem_vel[2].y * elem_vel[2].y)) *
	(sqrt(cf[2]));

      /* ****************************************************************     */
      /* Normal unit vector to the velocity vector dot product with               
         gradient (h*tau_s) is the normal derivative of the depth-integrated
         shear stress created by secondary flow.  Known here as normal_derive     */

      normal_derive =
	(((h_vor[0] * grad_shp[0].x + h_vor[1] * grad_shp[1].x +
	   h_vor[2] * grad_shp[2].x) * v * tau_s + (tau[0] * grad_shp[0].x +
							tau[1] * grad_shp[1].x +
							tau[2] * grad_shp[2].x) * h * v) -
	 ((h_vor[0] * grad_shp[0].y + h_vor[1] * grad_shp[1].y +
	   h_vor[2] * grad_shp[2].y) * u * tau_s + (tau[0] * grad_shp[0].y +
							tau[1] * grad_shp[1].y +
							tau[2] * grad_shp[2].y) * h * u)) /
	abs_u_mag;

      /* *****************************************************************        */
      /* Streamwise Acceleration (S) unit vector                          */
      /*  This is returning the streamwise acceleration vector terms for the      */
      /*  momentum equations.  It is being called by fe_sw2_elem_resid.c          */

      stream_accel_vect.x =
	h * (1 / fluid_density) * (u / abs_u_mag) * ((normal_derive / h) + 2 * tau_s / radius);
      stream_accel_vect.y =
	h * (1 / fluid_density) * (v / abs_u_mag) * ((normal_derive / h) + 2 * tau_s / radius);
	//  printf("X %lf Y %lf DEN %lf U %lf H %lf tau %lf rad %lf ABS_U %lf C_f %lf C0 %lf c1 %lf c2 %lf\n", stream_accel_vect.x,  stream_accel_vect.y, fluid_density,u,  h, tau_s, radius, abs_u_mag, c_f, elem_c[0], elem_c[1],  elem_c[2]);
      return (stream_accel_vect);
    }
  else
    {
      /*  This is returning the production and destruction terms for the vorticity        */
      /*  equation, x and y respectively. It is being called by fe_sw2trns_elem_resid.c */
		
		stream_accel_vect.x =	(prop1 * sqrt(c_f) * abs_u_mag * abs_u_mag) / (radius *
										  h * (1 +
										       (9 *
											h *
											h) /
										       (radius
											*
											radius)));
      stream_accel_vect.y = -prop2 * sqrt(c_f) * abs_u_mag / h;
	 // printf("X %lf Y %lf \n", stream_accel_vect.x, stream_accel_vect.y);

      return (stream_accel_vect);
    }

}
