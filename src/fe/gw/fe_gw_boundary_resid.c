/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes the 2D element residual for the groundwater model.
 * @param[out] elem_rhs      the 2D elemental residual array
 * @param[in]  mod           a pointer to the model struct
 * @param[in]  ie            the elemental id
 * @param[in]  perturbation   the Newton perturbation
 * @param[in]  perturb_node  the node to be perturbed
 * @param[in]  perturb_var   the variable to be perturbed
 * @param[in]  perturb_sign  the direction of Newton perturbation
 * @param[in]  DEBUG         a debug flag
 *
 
TODO:
  - Incorporate metfile, surface exchange terms
    - options for direct calculation from code
    - area or area3d?
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_gw_boundary_resid(SMODEL *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    
//    printf("myid: %d  ie: %d  elem3did: %d string: %d\n",mod->grid->smpi->myid,ie,mod->grid->elem2d[ie].id_3d,mod->grid->elem2d[ie].string);
//    tl_check_all_pickets(__FILE__,__LINE__);
//    MPI_Barrier(MPI_COMM_WORLD); // tl_error("test");
    
  int i, DEBUG_LOCAL = OFF;
    
#ifdef _DEBUG
    time_t time1;  time(&time1);
#endif
    
  double weight;
    
  // aliases
  SGW *gw = mod->sgw;
  SGRID *grid = mod->grid;
  SELEM_2D *elem2d = &(grid->elem2d[ie]); 
  int nnodes = elem2d->nnodes;
  double dt = mod->dt;
  double djac = elem2d->djac3d;
  int string = elem2d->string;
  STR_VALUE *str_values = mod->str_values;

  assert(nnodes == NDONTRI); // assume only triangular faces for now
  assert(djac > 0.); //
    
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  // INDEPENDENT VARIABLES
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  double elem_phead[nnodes];
  global_to_local_dbl(gw->gw_phead, elem_phead, elem2d->nodes, nnodes);
  if (perturb_var == PERTURB_H) elem_phead[perturb_node] += perturb_sign * perturbation;

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  // DEPENDENT VARIABLES
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   *                                FINITE ELEMENT INTEGRATIONS
   *==========================================================================================*/
  sarray_init_dbl(elem_rhs, MAX_NNODES_ON_ELEM3D);
    
  if (string > NORMAL) {
    if (str_values[string].flow.bc_flag == BCT_VEL_NEU) {
      int isers = str_values[string].flow.isigma;
      double flux = -1.0 * sseries_get_value(isers, mod->series_head, 0);
      /* midpoint integration */
      double test_val = 1./3.;
      weight = djac;
      for (i=0; i < nnodes; i++) {
        elem_rhs[i] += dt * flux * test_val * weight;
      }
    } /*neumann flux */
  } /*string > NORMAL */
    

  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   *                                DIRICHLET BOUNDARY CONDITIONS
   *--------------------------------------------------------------------------------------------
   * Zeros matrix rows/columns when Dirichlet boundaries are used. \n
   * \note
   *********************************************************************************************/
  for (i=0; i<nnodes; i++) {
    string = mod->grid->node[ elem2d->nodes[i] ].string;
    if (string > NORMAL) {
      if (str_values[string].flow.bc_flag == BCT_PRS_DIR) {
        int isers = str_values[string].flow.iu_0;
        //mwf debug 
#ifdef _DEBUG
        if (DEBUG_LOCAL==ON) printf("in boundary residi ie=%d i=%d string= %d BCT_PRS_DIR sol=%g val=%g \n",ie,i,string,elem_phead[i],
	                                 sseries_get_value(isers, mod->series_head, 0));
#endif

        elem_rhs[i] = 0.0;
      }
      else if (str_values[string].flow.bc_flag == BCT_DIR) {
        elem_rhs[i] = 0.0;
      }
    }
  }
  
    //tl_check_all_pickets(__FILE__,__LINE__);
    //MPI_Barrier(MPI_COMM_WORLD); // tl_error("test");
    
    
#ifdef _DEBUG
    time_t time2;  time(&time2);
    TIME_IN_GW_BOUNDARY_RESID += difftime(time2,time1);
#endif
    
  return;
}
