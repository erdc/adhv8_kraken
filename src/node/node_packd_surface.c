/* this routine packs the double information for the given node in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void node_packd_surface(void *buffer,   /* the void buffer to be packed */
                int ind,         /* the node to be packed */
                SSW *sw,
                SGRID *grid 

    )
{
  int ibuff = 0;                /* current counter in the buffer */
  double *dbuffer;              /* the double version of the buffer */
  int i,j,k;

  /* cast the buffer */
  dbuffer = (double *) buffer;

    dbuffer[ibuff++] = sw->d3->depth[ind];
    dbuffer[ibuff++] = sw->d3->old_depth[ind];
    //dbuffer[ibuff++] = sw->d3->bed_elevation[ind];
    dbuffer[ibuff++] = sw->d3->depth_avg_vel[ind].x;
    dbuffer[ibuff++] = sw->d3->depth_avg_vel[ind].y;
    dbuffer[ibuff++] = sw->d3->old_depth_avg_vel[ind].x;
    dbuffer[ibuff++] = sw->d3->old_depth_avg_vel[ind].y;
    dbuffer[ibuff++] = sw->d3->surface_vel[ind].x;
    dbuffer[ibuff++] = sw->d3->surface_vel[ind].y;
    //dbuffer[ibuff++] = sw->d3->bottom_vel[ind].x;
    //dbuffer[ibuff++] = sw->d3->bottom_vel[ind].y;
      
    if (sw->d3->waves != NULL) {
      dbuffer[ibuff++] = sw->d3->waves[ind].rads.xx;
      dbuffer[ibuff++] = sw->d3->waves[ind].rads.yy;
      dbuffer[ibuff++] = sw->d3->waves[ind].rads.xy;
      dbuffer[ibuff++] = sw->d3->waves[ind].stress.x;
      dbuffer[ibuff++] = sw->d3->waves[ind].stress.y;
    }

    if (sw->d3->winds != NULL ) {
      dbuffer[ibuff++] = sw->d3->winds[ind].stress.x;
      dbuffer[ibuff++] = sw->d3->winds[ind].stress.y;    
    }
}
#else
void node_packd_surface(void)
{
}
#endif
