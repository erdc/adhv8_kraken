/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the GW  elemental residual.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] elem_rhs      the 3D elemental residual array
 * @param[in]  mod           a pointer to the model struct
 * @param[in]  ie            the elemental id
 * @param[in]  perturbation   the Newton perturbation
 * @param[in]  perturb_node  the node to be perturbed
 * @param[in]  perturb_var   the variable to be perturbed
 * @param[in]  perturb_sign  the direction of Newton perturbation
 * @param[in]  DEBUG         a debug flag
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_gw_body_resid(SMODEL *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
  int DEBUG_LOCAL = OFF;
    
#ifdef _DEBUG
    time_t time1;  time(&time1);
#endif
    
  int i;
  double weight;
    
  // aliases
  SGW *gw = mod->sgw;
  SGRID *grid = mod->grid;
  SELEM_3D *elem3d = &(grid->elem3d[ie]); 
  SVECT *grad_shp = elem3d->grad_shp;
  SVECT darcy_flux;
  int nnodes = elem3d->nnodes;
  double alpha = mod->tau_temporal;
  double dt = mod->dt;
  double djac = elem3d->djac;
  int imat = elem3d->mat;
  SMAT_GW *gwm = mod->mat[imat].gw;
  STENSOR k_total;
  SSERIES * psk_sat; SSERIES * psk_kr; 
  double sat,kr,sat_slope,kr_slope;
  double ref_density = mod->density;
  
  int isTet = NO; if (nnodes == NDONTET) isTet = YES;
  assert(isTet);
  assert(djac > SMALL);
  assert(alpha > -1E-6 && alpha <= 1.0);
  assert(imat >= 0); assert(imat < mod->nmat);
  assert(perturb_var == UNSET_INT || perturb_var == PERTURB_H);
  assert(gwm);
  assert(gwm->isat >= 0); assert(gwm->isat < 2*mod->nmat);
  assert(gwm->ikr >= 0); assert(gwm->ikr < 2*mod->nmat);
  psk_sat = sseries_search(gwm->isat,NULL,mod->series_gw_psk_head);
  assert(psk_sat);
  psk_kr = sseries_search(gwm->ikr,NULL,mod->series_gw_psk_head);
  assert(psk_kr);
  
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  // GRID VARIABLES
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
  double node_z[nnodes];
  SNODE nodes[nnodes];
  for (i=0; i<nnodes; i++) {
    snode_copy(&(nodes[i]), grid->node[elem3d->nodes[i]]);
    node_z[i] = nodes[i].z;
  }
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  // INDEPENDENT VARIABLES
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  double elem_phead[nnodes];
  global_to_local_dbl(gw->gw_phead, elem_phead, elem3d->nodes, nnodes);
  if (perturb_var == PERTURB_H) {
    //mwf debug 
#ifdef _DEBUG
    if (DEBUG_LOCAL==ON) printf("gw body resid perturb_node= %d sign= %d total perturb= %g val was \n %30.20f",perturb_node,perturb_sign,
	   perturb_sign * perturbation,
	   elem_phead[perturb_node]);
#endif
    
    elem_phead[perturb_node] += perturb_sign * perturbation;

  }
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  // DEPENDENT VARIABLES
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  double elem_old_phead[nnodes], elem_older_phead[nnodes];
  global_to_local_dbl(gw->old_gw_phead, elem_old_phead, elem3d->nodes, nnodes);
  global_to_local_dbl(gw->older_gw_phead, elem_older_phead, elem3d->nodes, nnodes);

  double elem_density[nnodes];
  global_to_local_dbl(gw->gw_density,elem_density, elem3d->nodes,nnodes);
  
  /* added a check for no-flow elements (rocks, mines, etc) */
  /* Treat as solid.  Return no contribution to the elemental residual */
  double mag_k = gwm->k.xx + gwm->k.yy + gwm->k.zz;
  if (gwm->porosity < NOT_QUITE_SMALL || mag_k < NOT_QUITE_SMALL) {
    return;
  }
  /**********************************************************************/
  /* get the material properties */
  double s_s = gwm->s_s;
  //VT_3D_TCOPY(gwm->k, k_total);

  /**********************************************************************/
  /* constitutive evaluations */
  double elem_sat[nnodes], elem_sat_old[nnodes]; 
  for (i=0; i < nnodes; i++) {
    elem_sat[i] = sgw_eval_sat(psk_sat,elem_phead[i],&elem_sat[i],&sat_slope);
    elem_sat_old[i] = gw->elem_3d_data[ie].old_saturation[i];
  }
  //double elem_kr = sgw_eval_kr_elem(psk_kr,nnodes,elem_phead);
  ///*mwf hack looking at nonlinearities 
  //elem_kr = 1.0;
  //for (i=0; i < nnodes; i++) {elem_sat[i]=1.;elem_sat_old[i]=1.;}
  //*/
  ///* scales by the relative permeability */
  //VT_3D_TSCALE(k_total, elem_kr);
  //double avg_dens = 0.25*(elem_density[0]+elem_density[1]+elem_density[2]+elem_density[3]);
  //
  //SVECT gradient_z;
  //gradient_z.x = 0.; gradient_z.y = 0.; gradient_z.z = 1.;
  //
  ///**********************************************************************/
  ///* solution gradient */
  //SVECT grad_phead;
  //grad_phead.x=0.; grad_phead.y=0.; grad_phead.z=0.;
  //for (i=0; i < nnodes; i++) {
  //  grad_phead.x += elem_phead[i]*grad_shp[i].x;
  //  grad_phead.y += elem_phead[i]*grad_shp[i].y;
  //  grad_phead.z += elem_phead[i]*grad_shp[i].z;
  //}
  //
  ///**********************************************************************/
  ///* elemental flux */
  //SVECT k_grad_phead,k_grad_z,darcy_flux;
  //VT_3D_TENS_VECT_PROD(k_grad_phead,k_total,grad_phead);
  //VT_3D_TENS_VECT_PROD(k_grad_z,k_total,gradient_z);

  ///*mwf not allowing density dependence for now */
  //assert(fabs(avg_dens/ref_density-1.0) <= SMALL);
  //darcy_flux.x = -(k_grad_phead.x + avg_dens/ref_density*k_grad_z.x);
  //darcy_flux.y = -(k_grad_phead.y + avg_dens/ref_density*k_grad_z.y);
  //darcy_flux.z = -(k_grad_phead.z + avg_dens/ref_density*k_grad_z.z);


  sgw_eval_flux(nnodes,psk_kr,gwm->k,ref_density,elem_density,elem_phead,grad_shp,&darcy_flux);

  /*mwf debug 
  printf("###fe_gw_body_resid elem %d djac= %g perturbation = %30.10e elem_kr= %g###\n",ie,djac,perturbation,elem_kr);
  for (i=0; i < nnodes; i++) {
    printf("---node %d : %d phead= %g \n",i,elem3d->nodes[i],elem_phead[i]);
    printf("---grad_shp = %g %g %g \n",grad_shp[i].x,grad_shp[i].y,grad_shp[i].z);
  }
  printf("---grad_phead = %g %g %g \n",grad_phead.x,grad_phead.y,grad_phead.z);
  */
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   *                                FINITE ELEMENT INTEGRATIONS
   *==========================================================================================*/
  double rhs[nnodes];
  sarray_init_dbl(elem_rhs,nnodes);

//  for (i=0; i < nnodes; i++) {
//    if (grid->smpi->myid ==0 && grid->node[elem3d->nodes[i]].gid == 49) printf("body resid elem %d node %d ph_new=%g ph_old=%g sat_new=%g sat_old=%g darcy flow: %f %f %f\n",ie,i,elem_phead[i],elem_old_phead[i],elem_sat[i],elem_sat_old[i],darcy_flux.x,darcy_flux.y,darcy_flux.z);
//  }

  // velocity terms midpoint quadrature
  weight = djac * dt;
  for (i=0; i < nnodes; i++) {
    elem_rhs[i] = -VT_3D_DOT(darcy_flux,grad_shp[i]) * weight;
    //mwf debug 
#ifdef _DEBUG
    if (DEBUG_LOCAL==ON) printf("fe_gw_body_resid ie=%d i=%d v.dw= %g\n",ie,i,elem_rhs[i]); 
#endif
    
  }
  //mwf debug
#ifdef _DEBUG
  if (DEBUG_LOCAL==ON) printf("fe_gw_body_resid ie=%d dt=%g\n",ie,dt);
#endif

  //source terms midpoint quadrature
  double test_val = 0.25;
  weight = dt * djac ;
  for (i=0; i < nnodes; i++) {
    elem_rhs[i] -= gwm->water_vol * test_val * weight;
  }
  //time term vertex quadrature
  weight = djac * 0.25;
  for (i=0; i < nnodes; i++) {
    /*mwf debug 
    printf("\telem_time[%d] weight = %g s_s=%g elem_sat=%g elem_sat_old=%g poro=%g elem_phead=%g elem_old_phead=%g\n",i,weight,
	   s_s,elem_sat[i],elem_sat_old[i],gwm->porosity,elem_phead[i],elem_old_phead[i]);
    
    printf("\telem_time[%d] = %30.10e\n",i,s_s * elem_sat[i] * (elem_phead[i]-elem_old_phead[i])  +
      gwm->porosity * (elem_sat[i] - elem_sat_old[i]) );
    */
    elem_rhs[i]+= s_s * elem_sat[i] * (elem_phead[i]-elem_old_phead[i]) * weight +
      gwm->porosity * (elem_sat[i] - elem_sat_old[i]) * weight;
    /*mwf hack 
    elem_rhs[i] = (elem_phead[i] - node_z[i]) * weight;
    if (perturb_var == PERTURB_H) {
      printf("body resid ie=%d i=%d weight=%g elem_phead=%30.20f perturb= %g node_z=%g res=%30.20f\n",ie,i,weight,elem_phead[i],perturbation,node_z[i],elem_rhs[i]);
    }
    */
  }

    
#ifdef _DEBUG
    time_t time2;  time(&time2);
    TIME_IN_GW_BODY_RESID += difftime(time2,time1);
#endif
    
  return;
}

