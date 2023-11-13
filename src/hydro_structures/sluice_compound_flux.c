/* THis routine calculates flow through a sluice gate
based on a sluice equation presented in Swamee " Sluice-Gate Discharge Equations 
JHE, Vol 118, No. 1, 1992 */

#include "global_header.h"

void sluice_compound_flux(int iwn, SNODE *node, int nnodes, double *ol_head, SVECT2D *ol_vel, 
	                    int nsluice, STR_VALUE *str_values, SSERIES *series_head, SSLUICE_C *sluice, double gravity, SGRID *grid)
{
  int i,j, ireceive0, isend0, ierr_code, WEIR_TAGU=9997, WEIR_TAGD=9998;
  double Upstream[3],Downstream[3],Zd,Vd,Zu,Vu,max, ireceive0_d, isend0_d;
  double isend_a[3], ireceive_a[3];
  double widthd,widthu,widthd_t,widthu_t,multf,direction;
  double Cd, p1, p2, p3, p4, p5, value;
  double Hvu,Hvd,Htu,Htd,He1,He2,L1,L2,Q1,Q2,s1,S2,Qr1,Qr2,dx,dy,Hcd;  /* Qr1 is submergence correction*/
  int up1 = 0;
  int down1 = 0;
  int submerged, isers;
  double upflux, getflux;
   Upstream[0] = Upstream[1]= Upstream[2] = Downstream[0] = Downstream[1] = Downstream[2] = 0.0;
  for(i = 0; i < nnodes; i++)
   {


    if((node[i].edge_string == sluice[iwn].egsu) || (node[i].string == sluice[iwn].egsu))
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
    if(((node[i].edge_string == sluice[iwn].egsd) || (node[i].string == sluice[iwn].egsd)) && (node[i].edge_string != sluice[iwn].egsu) && (node[i].string != sluice[iwn].egsu))
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



  L1 = sluice[iwn].a;  /* Sluice Length or width cross channel*/
  isers = sluice[iwn].opening;
  value = sseries_get_value(isers, series_head, 0);

  /* DECIDE IF SUBMERGED FLOW */
  if ((Upstream[2] > Downstream[2]) && ( Upstream[2] < (0.81 * Downstream[2] * pow(Downstream[2]/value,0.72))))
	  submerged = 1;
  if ((Upstream[2] > (0.81 * Downstream[2] * pow(Downstream[2]/value,0.72))) || (Upstream[2] == (0.81 * Downstream[2] * pow(Downstream[2]/value,0.72))))
	  submerged = 0;
  if (submerged == 0)
	  Cd = 0.611 * pow((Upstream[2] - value)/(Upstream[2] + 15. * value),0.072);
 else if (submerged == 1)
  {
	  p1 = 0.611 * pow((Upstream[2] - value)/(Upstream[2] + 15. * value),0.072);
	  p2 = pow(Upstream[2] - Downstream[2], 0.7);
	  p3 = 0.32 * pow(0.81 * Downstream[2] * pow(Downstream[2]/value,0.72) - Upstream[2], 0.7);
	  p4 = pow(Upstream[2] - Downstream[2], 0.7);
	  p5 = pow(p3 + p4,-1.);
	  Cd = p1 * p2 * p5;
  }
else
   Cd = 0.;

Q1 = Cd * value *pow(2. * gravity * Upstream[2], 0.5);  /* Flow per unit length */



  upflux = Q1;
#ifdef _MESSG
  MPI_Allreduce (&upflux,&getflux,1,MPI_DOUBLE,MPI_SUM,grid->smpi->ADH_COMM);
#endif


  sluice[iwn].fluxu= upflux;
  sluice[iwn].fluxd= -upflux;



  /* if Upstream WSEL < Weir Crest elevation */
  /* if (SW2_FLOW)*/
 
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
//  if (myid <= 0)
//  printf("\nThe Flux to the downstream (+ve going out) is %lf Coeff is %lf The Upstream WSE is %lf Downstream is %lf Sumberged is %d (0: Free-Flow 1: Submerged) \n",sluice[iwn].fluxu, Cd, Upstream[0], Downstream[0], submerged);

}
