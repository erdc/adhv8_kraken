/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     The FE engine for computing the next time-step of the 3D SW model. Returns
 *            NO if the nonlinear step fails, YES if successful.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] mod (SMODEL *) a pointer to a 3D SW wave model struct
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int fe_gw(SMODEL *mod) {
    
  assert(mod->sgw);
  assert(mod->grid);
  
  // alias
  SGW *gw = mod->sgw;
  int nnodes = mod->grid->nnodes;
    
  // prep solutions
  mod->old_dt = mod->dt;
  int i, j, inode = 0;
  for (inode = 0; inode < nnodes; inode++) {
    gw->older_gw_phead[inode]  = gw->old_gw_phead[inode];
    gw->old_gw_phead[inode]  = gw->gw_phead[inode];
  }
  int ie = 0;
  for (ie = 0; ie < mod->grid->nelems3d; ie++) {
    for (i=0; i < mod->grid->elem3d[ie].nnodes; i++) {
      gw->elem_3d_data[ie].old_saturation[i] = gw->elem_3d_data[ie].saturation[i];
    }
  }
    
  // calculate density 
  tl_density_calculator_metric(mod->density, NULL, 1.0, NULL, 1.0, mod->grid->nnodes, mod->sgw->gw_density, 3);
  //
  
  // Gajanan gkc :: moved to fe_main.c
//  mod->nsys=1;
//  mod->nsys_sq = 1;
//  mod->solver_info.refresh=YES;
//  mod->solver_info.PRN_NEWTON = GW;
//  mod->solver_info.LINEAR_PROBLEM = NO;
//#ifdef _DEBUG
//  /*mwf debug */
//  //mod->solver_info.max_nonlin_linesearch_cuts = 5;
//  //printf("in fe_gw setting max linesearches = %d\n",mod->solver_info.max_nonlin_linesearch_cuts);
//#endif
//
//  /** */
//  if (fe_newton(mod, mod->grid, fe_gw_init, fe_gw_update, fe_gw_resid, fe_gw_load, fe_gw_inc) == NO) {
//    return (NO);
//  }
//
//  sgw_evaluate_element_fluxes(mod);
//    
//  // reset to defaults
//  mod->solver_info.PRN_NEWTON = OFF;
//  
//  // check total grid mass
//  if (screen_output.grid_mass_error == ON) {
//        //mod->grid_mass_error = tl_find_grid_mass_error_tet_gw(mod->initial_grid_mass, mod->grid, gw->gw_phead, gw->elem_3d_data, mod->series_head, mod->str_values, mod->dt);
//        //printf("  grid_mass_error: %10.5e",mod->grid_mass_error);
//  }
//    
//  // calculate approximate continuity error over elements for adaption
//  //gw_3d_calculate_elem_error(gw, mod->grid, mod->mat, mod->dt);
    
  return YES;
}
