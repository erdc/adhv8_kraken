/* This routine calculates the
   flow over weirs (ogee type) and passes the flows
   to the code as a boundary*/
/* GSavant 10/25/2008*/

#include "global_header.h"

void weir_compound_flux(int iwn, SNODE *node, int nnodes, double *ol_head, SVECT2D *ol_vel, 
	                    int nweir, STR_VALUE *str_values, SWEIR_C *weir, double gravity, SGRID *grid)
 {
  int i,j, ireceive0, isend0, ierr_code, WEIR_TAGU=9997, WEIR_TAGD=9998;
  double Upstream[3],Downstream[3],Zd,Vd,Zu,Vu,max, ireceive0_d, isend0_d;
  double isend_a[3], ireceive_a[3];
  double widthd,widthu,widthd_t,widthu_t,multf,direction;
  double Hvu,Hvd,Htu,Htd,He1,He2,L1,L2,Q1,Q2,s1,S2,Qr1,Qr2,dx,dy,Hcd;  /* Qr1 is submergence correction*/
  int up1 = 0;
  int down1 = 0;
  double upflux, getflux;
   Upstream[0] = Upstream[1]= Upstream[2] = Downstream[0] = Downstream[1] = Downstream[2] = 0.0;
  for(i = 0; i < nnodes; i++)
   {


    if((node[i].edge_string == weir[iwn].egsu) || (node[i].string == weir[iwn].egsu))
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
	/* TEST OUT if (SW3_FLOW)
	{
      Upstream[0]+= displacement[i]+node[i].z;
      Upstream[1]+= sqrt(vel[i].x*vel[i].x+vel[i].y*vel[i].y+vel[i].z*vel[i].z);
	  Upstream[2]+= displacement[i];
      up1 += 1 ;
	} TEST OUT */
    }
    if(((node[i].edge_string == weir[iwn].egsd) || (node[i].string == weir[iwn].egsd)) && (node[i].edge_string != weir[iwn].egsu) && (node[i].string != weir[iwn].egsu))
    {
	  /*if (SW2_FLOW)*/
	  {
	  if (ol_head[i] > 0)
	   {
      Downstream[0]+=  ol_head[i]+node[i].z;
      Downstream[1]+=  sqrt(ol_vel[i].x*ol_vel[i].x+ol_vel[i].y*ol_vel[i].y);
	  Downstream[2]+= ol_head[i];
      down1 += 1;
	  }
	}
	 /* TEST OUT  	if (SW3_FLOW)
	{
      Downstream[0]+=  displacement[i]+node[i].z;
      Downstream[1]+=  sqrt(vel[i].x*vel[i].x+vel[i].y*vel[i].y+vel[i].z*vel[i].z);
	  Downstream[2]+= displacement[i];
      down1 += 1;
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

  if (up1 != 0){
  Upstream[0] = Upstream[0]/up1;
  Upstream[1] /= up1;
  Upstream[2] /=up1;}
  if (down1 != 0){
  Downstream[0] /= down1;
  Downstream[1] /= down1;
  Downstream[2] /= down1;}





  widthu_t=widthu;
  widthd_t=widthd;


  Zd=Downstream[0];
  Vd=Downstream[1];
  Zu=Upstream[0];
  Vu=Upstream[1];

  weir[iwn].cd=3.0;



	  L1 = weir[iwn].w1;

  Hcd = Zu - weir[iwn].z1;
  /*L2=weir[iwn].w2+He2*weir[iwn].sl2;*/
  /*weir[iwn].cd = 0.5 + 0.1 * pow((pow(Hcd/L1,5) + 1500 * pow(Hcd/L1,13))/(1 + 1000 * pow(Hcd/L1,3)),0.1);*/ /* Swamee , 1988*/


  weir[iwn].cd = (2./3.)*(0.611 + 0.075 * (Zu - weir[iwn].z1)/weir[iwn].hw)*pow(2. * gravity,0.5);
  



  if (Zu > weir[iwn].z1 && Upstream[0] > Downstream[0])
    Q1 = weir[iwn].cd * pow(Zu - weir[iwn].z1,1.5) * L1;
  else
	Q1 = 0.;


  if ( Upstream[0] > weir[iwn].z1 && Downstream[0] > weir[iwn].z1 && Upstream[0] > Downstream[0])
  Qr1 = pow( 1 - pow((Downstream[0] - weir[iwn].z1)/(Upstream[0] - weir[iwn].z1),1.5),0.385);
  else Qr1 = 1;



  upflux = Q1;
#ifdef _MESSG
  MPI_Allreduce (&upflux,&getflux,1,MPI_DOUBLE,MPI_SUM,grid->smpi->ADH_COMM);
#endif


  weir[iwn].fluxu= upflux/L1;
  weir[iwn].fluxd= -upflux/L1;




  /* if Upstream WSEL < Weir Crest elevation */
  /* if (SW2_FLOW)*/
  {
  if ((Upstream[0] < weir[iwn].z1) || (Upstream[2] <= 0))
  {
  weir[iwn].fluxu= 0;
  weir[iwn].fluxd= 0;
  }
  }
 /* TEST OUT if (SW3_FLOW)
  {
	  	  for (i =0; i < nstring; i ++)
	  {
		  if (str_values[i].ol_flow.bc_flag = BCT_WEIRU)
			  for (j =0; j < nhydraulicstructure; j ++)
				  if (hydraulic_structure[j].up == i)
				  {
					  hydraulic_structure[j].flow = upflux;
					
				  }
	  }
  } TEST OUT */
 // if (myid <= 0)
 // printf("\nThe Flux to the downstream (+ve going out) is %lf Weir Coeff is %lf The Upstream WSE is %lf Downstream is %lf Submergence Modifier is %lf\n",weir[iwn].fluxu, weir[iwn].cd, Upstream[0], Downstream[0],Qr1);

}
