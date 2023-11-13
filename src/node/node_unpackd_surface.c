/* this routine unpacks the double information for the given node in the given buffer */

#include "global_header.h"

#ifdef _MESSG
void node_unpackd_surface(void *buffer, /* the void buffer to be unpacked */
                  int ind,       /* the node to be unpacked */
                  SSW *sw
    )
{
  int ibuff = 0;                /* current counter in the buffer */
  double *dbuffer;              /* the double version of the buffer */
  int i,j,k;
  /* cast the buffer */
  dbuffer = (double *) buffer;


  /* sets the buffer */
    sw->d3->depth[ind] = dbuffer[ibuff++];
    sw->d3->old_depth[ind] = dbuffer[ibuff++];
    //sw->d3->bed_elevation[ind] = dbuffer[ibuff++];
    sw->d3->depth_avg_vel[ind].x = dbuffer[ibuff++];
    sw->d3->depth_avg_vel[ind].y = dbuffer[ibuff++];
    sw->d3->old_depth_avg_vel[ind].x = dbuffer[ibuff++];
    sw->d3->old_depth_avg_vel[ind].y = dbuffer[ibuff++];
    sw->d3->surface_vel[ind].x = dbuffer[ibuff++];
    sw->d3->surface_vel[ind].y = dbuffer[ibuff++];
    //sw->d3->bottom_vel[ind].x = dbuffer[ibuff++];
    //sw->d3->bottom_vel[ind].y = dbuffer[ibuff++];

    if (sw->d3->waves != NULL) {
      sw->d3->waves[ind].rads.xx = dbuffer[ibuff++];
      sw->d3->waves[ind].rads.yy = dbuffer[ibuff++];
      sw->d3->waves[ind].rads.xy = dbuffer[ibuff++];
      sw->d3->waves[ind].stress.x = dbuffer[ibuff++];
      sw->d3->waves[ind].stress.y = dbuffer[ibuff++];
    }


    if (sw->d3->winds) {
      sw->d3->winds[ind].stress.x = dbuffer[ibuff++];
     sw->d3->winds[ind].stress.y = dbuffer[ibuff++];
    }
    
}
#else
void node_unpackd_surface(void)
{
}
#endif
