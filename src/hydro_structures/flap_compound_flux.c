/* This routine calculates the
   flap gate flow  and passes the flows
   to the code as a boundary*/
/* GSavant 10/25/2008*/
/* Based on Swammee 1992, journal of Irrigation and Drainage Engr*/

#include "global_header.h"

void flap_compound_flux(int ifp, SNODE *node, int nnodes, double *ol_head, SVECT2D *ol_vel,
	                    int nflap, STR_VALUE *str_values, SFLAP_C *flap, double gravity, SGRID *grid)
{
  int i, j,ireceive0, isend0, ierr_code, WEIR_TAGU=9997, WEIR_TAGD=9998;
  double Upstream[3],Downstream[3],Zd,Vd,Zu,Vu,max, ireceive0_d, isend0_d;
  double isend_a[3], ireceive_a[3];
  double widthd,widthu,widthd_t,widthu_t,multf,direction;
  double Hvu,Hvd,Htu,Htd,He1,He2,L1,L2,Q1,Q2,s1,S2,Qr1,Qr2,dx,dy,Hcd;  /* Qr1 is submergence correction*/
  int up1 = 0;
  int down1 = 0;
  double upflux, getflux;
  double flp_mass, buoyancy, ang, rise;
  double del_l, A, B, C, D, E, F, Gg, H; /* Flap gate parameters */

  Upstream[0] = Upstream[1]= Upstream[2] = Downstream[0] = Downstream[1] = Downstream[2] = 0.0;

  for(i = 0; i < nnodes; i++)
   {


    if((node[i].edge_string == flap[ifp].egsu) || (node[i].string == flap[ifp].egsu))
    {
	  /*if (SW2_FLOW)*/
	  {
	  if (ol_head[i] > 0)
	   {
      Upstream[0]+= ol_head[i]+node[i].z;
      Upstream[1]+= sqrt(ol_vel[i].x*ol_vel[i].x+ol_vel[i].y*ol_vel[i].y);
	  Upstream[2]+= ol_head[i];
      up1 += 1 ;
	  }
	  }
	 /* TEST OUT  else if (SW3_FLOW)
	  {
      Upstream[0]+= displacement[i]+node[i].z;
      Upstream[1]+= sqrt(vel[i].x*vel[i].x+vel[i].y*vel[i].y+vel[i].z*vel[i].z);
	  Upstream[2]+= displacement[i];
      up1 += 1 ;
	  } TEST OUT */
    }
    if(((node[i].edge_string == flap[ifp].egsd) || (node[i].string == flap[ifp].egsd)) && (node[i].edge_string != flap[ifp].egsu) && (node[i].string != flap[ifp].egsu))
    {
	 /* if (SW2_FLOW)*/
	  {
	  if (ol_head[i] > 0)
	   {
      Downstream[0]+=  ol_head[i]+node[i].z;
      Downstream[1]+=  sqrt(ol_vel[i].x*ol_vel[i].x+ol_vel[i].y*ol_vel[i].y);
	  Downstream[2]+= ol_head[i];
      down1 += 1;
	  }
	  }
	/* TEST OUT  else if (SW3_FLOW)
	  {
      Downstream[0]+=  displacement[i]+node[i].z;
      Downstream[1]+=  sqrt(vel[i].x*vel[i].x+vel[i].y*vel[i].y+vel[i].z*vel[i].z);
	  Downstream[2]+= displacement[i];
      down1 += 1;
	  printf("\nDS REF STRING NODE %d WSEL %G\n",i+1, Downstream[0]); 
	  } TEST OUT */
	}
  }
 #ifdef _MESSG
  MPI_Allreduce(&up1, &ireceive0, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM);
  MPI_Allreduce(&down1, &isend0, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM);
  up1 =  ireceive0;
  down1 = isend0;
  MPI_Allreduce (Upstream, ireceive_a, 3, MPI_DOUBLE, MPI_SUM, grid->smpi->ADH_COMM);
  MPI_Allreduce (Downstream, isend_a, 3, MPI_DOUBLE, MPI_SUM, grid->smpi->ADH_COMM);
  Upstream[0] = ireceive_a[0];
  Upstream[1] = ireceive_a[1];
  Upstream[2] = ireceive_a[2];
  Downstream[0] = isend_a[0];
  Downstream[1] = isend_a[1];
  Downstream[2] = isend_a[2];
 #endif

  if (up1 > 0)
  {
  Upstream[0] = Upstream[0]/up1;
  Upstream[1] = Upstream[1]/up1;
  Upstream[2] = Upstream[2]/up1;
  }
  if (down1 > 0)
  {
  Downstream[0] = Downstream[0]/down1;
  Downstream[1] = Downstream[1]/down1;
  Downstream[2] = Downstream[2]/down1;
  }

 /* if (flap[ifp].calc == 1)
  {*/
	/*if (SW2_FLOW)*/
	{
    
    if ((Upstream[0] > Downstream[0]) && (Upstream[2] > 0) && (Downstream[2] > 0))
      upflux = flap[ifp].A * pow((Upstream[0] - Downstream[0]), flap[ifp].B) + flap[ifp].C * (Upstream[0] - Downstream[0]) + flap[ifp].D;
    else
      upflux = 0.;
	}
	/* TEST OUT if (SW3_FLOW)
	{
		if (Upstream[0] > Downstream[0])
          upflux = flap[ifp].A * pow((Upstream[0] - Downstream[0]), flap[ifp].B) + flap[ifp].C * (Upstream[0] - Downstream[0]) + flap[ifp].D;
	} TEST OUT */

  /*}*/
  /*else if (flap[ifp].calc == 2)
  {
	  flp_mass = flap[ifp].density * (flap[ifp].length * flap[ifp].width);

	  buoyancy = density * gravity * Downstream[0];
	
	  ang = atan(density * gravity * (Upstream[0] - Downstream[0])/(flp_mass * gravity - buoyancy));
	
	  del_l = (flap[ifp].length/cos(ang)) - flap[ifp].length;
	
	  rise = del_l * cos(ang);
	 
	  A = pow(gravity * Upstream[0],0.5);
	  B = pow(((Upstream[0] - rise)/(Upstream[0] + 15 * rise)),0.072);
	  C = pow((Upstream[0] - Downstream[0]),0.7);
	  D = pow(Downstream[0]/rise,0.72);
	  E = 0.81 * D * Downstream[0];
	  F = 0.32 * pow((E - Upstream[0]),0.7);
	  Gg = pow((Upstream[0] - Downstream[0]),0.7);
	  H = pow((F+Gg),(-1));
	
	  upflux = 0.864 * rise * A * B * C * H;
  } */

#ifdef _MESSG
  MPI_Allreduce (&upflux,&getflux,1,MPI_DOUBLE,MPI_SUM,grid->smpi->ADH_COMM);
#endif


  /* if (SW2_FLOW)*/
  {
  flap[ifp].fluxu= upflux/flap[ifp].GG ;
  flap[ifp].fluxd= -1 * flap[ifp].fluxu ;
  }
 /* TEST OUT else if (SW3_FLOW)
  {
	  for (i =0; i < nstring; i ++)
	  {
		  if (str_values[i].ol_flow.bc_flag = BCT_FLAPU)
			  for (j =0; j < nhydraulicstructure; j ++)
				  if (hydraulic_structure[j].up == i)
				  {
					  hydraulic_structure[j].flow = upflux;
					
  
				  }
	  }
   } TEST OUT */

  //if (myid <= 0)
	 // printf("\nThe Flap flow to the downstream is %lf \n", flap[ifp].fluxu);

}
