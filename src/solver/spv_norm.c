/*!
   \file spv_norm.c
   \brief Computes the norm of the given sparse (SPV) vector 
 */

#include "global_header.h"

void spv_norm(
  SPARSE_VECT sv,		/* the sparse vector */
  double *dnorm,		/* the norm */
  int p				/* the number of equations being solved */
)
{
  double temp0, temp1, temp2, temp3;	/* temporary values */
  int i = 0;			/* loop counter */
  int id0, id11, id12, id13, id14, id21, id22, id23, id24, id31, id32, id33, id34, id41, id42, id43, id44;	/* indices for multiple equation multiplication */

  /* split over the number of equations being solved */
  if(p == 1)
    {
      /* initialize the return */
      dnorm[0] = DZERO;

      /* loops over the values and computes the norm */
      for(i = 0; i < sv.size; i++)
	{
	  temp0 = fabs(sv.value[i]);
	  if(temp0 > dnorm[0])
	    dnorm[0] = temp0;
	}
    }
  else if(p == 4)
    {
      /* initialize the return */
      dnorm[0] = 0.0;
      dnorm[1] = 0.0;
      dnorm[2] = 0.0;
      dnorm[3] = 0.0;

      /* loops over the values and computes the norm */
      for(i = 0, id0 = 0; i < sv.size; i++, id0 += 16)
	{
	  id11 = id0;
	  id12 = id0 + 1;
	  id13 = id0 + 2;
	  id14 = id0 + 3;
	  id21 = id0 + 4;
	  id22 = id0 + 5;
	  id23 = id0 + 6;
	  id24 = id0 + 7;
	  id31 = id0 + 8;
	  id32 = id0 + 9;
	  id33 = id0 + 10;
	  id34 = id0 + 11;
	  id41 = id0 + 12;
	  id42 = id0 + 13;
	  id43 = id0 + 14;
	  id44 = id0 + 15;
	  temp0 = fabs(sv.value[id11]);
	  temp1 = fabs(sv.value[id21]);
	  temp2 = fabs(sv.value[id31]);
	  temp3 = fabs(sv.value[id41]);
	  if(temp0 > dnorm[0])
	    dnorm[0] = temp0;
	  if(temp1 > dnorm[1])
	    dnorm[1] = temp1;
	  if(temp2 > dnorm[2])
	    dnorm[2] = temp2;
	  if(temp3 > dnorm[3])
	    dnorm[3] = temp3;
	  temp0 = fabs(sv.value[id12]);
	  temp1 = fabs(sv.value[id22]);
	  temp2 = fabs(sv.value[id32]);
	  temp3 = fabs(sv.value[id42]);
	  if(temp0 > dnorm[0])
	    dnorm[0] = temp0;
	  if(temp1 > dnorm[1])
	    dnorm[1] = temp1;
	  if(temp2 > dnorm[2])
	    dnorm[2] = temp2;
	  if(temp3 > dnorm[3])
	    dnorm[3] = temp3;
	  temp0 = fabs(sv.value[id13]);
	  temp1 = fabs(sv.value[id23]);
	  temp2 = fabs(sv.value[id33]);
	  temp3 = fabs(sv.value[id43]);
	  if(temp0 > dnorm[0])
	    dnorm[0] = temp0;
	  if(temp1 > dnorm[1])
	    dnorm[1] = temp1;
	  if(temp2 > dnorm[2])
	    dnorm[2] = temp2;
	  if(temp3 > dnorm[3])
	    dnorm[3] = temp3;
	  temp0 = fabs(sv.value[id14]);
	  temp1 = fabs(sv.value[id24]);
	  temp2 = fabs(sv.value[id34]);
	  temp3 = fabs(sv.value[id44]);
	  if(temp0 > dnorm[0])
	    dnorm[0] = temp0;
	  if(temp1 > dnorm[1])
	    dnorm[1] = temp1;
	  if(temp2 > dnorm[2])
	    dnorm[2] = temp2;
	  if(temp3 > dnorm[3])
	    dnorm[3] = temp3;
	}
    }
  else if(p == 3)
    {
      /* initialize the return */
      dnorm[0] = 0.0;
      dnorm[1] = 0.0;
      dnorm[2] = 0.0;

      /* loops over the values and computes the norm */
      for(i = 0, id0 = 0; i < sv.size; i++, id0 += p * p)
	{
	  id11 = id0;
	  id12 = id0 + 1;
	  id13 = id0 + 2;
	  id21 = id0 + 3;
	  id22 = id0 + 4;
	  id23 = id0 + 5;
	  id31 = id0 + 6;
	  id32 = id0 + 7;
	  id33 = id0 + 8;
	  temp0 = fabs(sv.value[id11]);
	  temp1 = fabs(sv.value[id21]);
	  temp2 = fabs(sv.value[id31]);
	  if(temp0 > dnorm[0])
	    dnorm[0] = temp0;
	  if(temp1 > dnorm[1])
	    dnorm[1] = temp1;
	  if(temp2 > dnorm[2])
	    dnorm[2] = temp2;
	  temp0 = fabs(sv.value[id12]);
	  temp1 = fabs(sv.value[id22]);
	  temp2 = fabs(sv.value[id32]);
	  if(temp0 > dnorm[0])
	    dnorm[0] = temp0;
	  if(temp1 > dnorm[1])
	    dnorm[1] = temp1;
	  if(temp2 > dnorm[2])
	    dnorm[2] = temp2;
	  temp0 = fabs(sv.value[id13]);
	  temp1 = fabs(sv.value[id23]);
	  temp2 = fabs(sv.value[id33]);
	  if(temp0 > dnorm[0])
	    dnorm[0] = temp0;
	  if(temp1 > dnorm[1])
	    dnorm[1] = temp1;
	  if(temp2 > dnorm[2])
	    dnorm[2] = temp2;
	}
    }
  else if(p == 2)
    {
      /* initialize the return */
      dnorm[0] = 0.0;
      dnorm[1] = 0.0;

      /* loops over the values and computes the norm */
      for(i = 0, id0 = 0; i < sv.size; i++, id0 += p * p)
	{
	  id11 = id0;
	  id12 = id0 + 1;
	  id21 = id0 + 2;
	  id22 = id0 + 3;
	  temp0 = fabs(sv.value[id11]);
	  temp1 = fabs(sv.value[id21]);
	  if(temp0 > dnorm[0])
	    dnorm[0] = temp0;
	  if(temp1 > dnorm[1])
	    dnorm[1] = temp1;
	  temp0 = fabs(sv.value[id12]);
	  temp1 = fabs(sv.value[id22]);
	  if(temp0 > dnorm[0])
	    dnorm[0] = temp0;
	  if(temp1 > dnorm[1])
	    dnorm[1] = temp1;
	}
    }

  else
    tl_error("Spv_norm has not been written for the desired system of equations.");
}
