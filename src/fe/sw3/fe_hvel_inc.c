/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Increments the SW 3D HVEL solutions after a Newton iterate.
 * \author    Charlie Berger, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Gary Brown, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout]  mod (SMODEL *) pointer to the model struct
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_hvel_inc(SSUPER_MODEL *sm, int imod) {

    int DEBUG = OFF;
    
    // aliases
    SMODEL *mod = &(sm->submodel[imod]);
    SGRID *grid = mod->grid;
    SSW_3D *sw3 = mod->sw->d3;
    
#ifdef _PETSC
    PetscScalar const *values;
    PetscInt ierr;

    ierr = VecGetArrayRead(sm->sol,&(values));

    // adds the increment to the solution by looping over the surface nodes, and traverse the vertical line of nodes below
    int i = 0, j = 0, nd = 0, top_node = 0, i3 = 0;
    double surf_dpl = 0.;
    ID_LIST_ITEM *ptr;
    for(i = 0; i < grid->nnodes_sur; i++) {
        // Get the first node in the vertical list (which is the surface node) 
        ptr = grid->vertical_list[i];
        nd = ptr->id;
        top_node=nd;
        // Find the increment for the surface displacement
        i3 = mod->fmap[nd] * 3;
        surf_dpl = values[i3+2];
        sw3->displacement[nd] += surf_dpl; 
        // Now, increment the velocities and the displacement for each node
        while(ptr->next != NULL) {
            nd = ptr->id;
            i3 = mod->fmap[nd] * 3;
            sw3->vel[nd].x += values[i3];
            sw3->vel[nd].y += values[i3+1];
            ptr = ptr->next;
        }
    }
    ierr = VecRestoreArrayRead(sm->sol,&(values));
#else
    // adds the increment to the solution by looping over the surface nodes, and traverse the vertical line of nodes below
    int i = 0, j = 0, nd = 0, top_node = 0, i3 = 0;
    double surf_dpl = 0.;
    ID_LIST_ITEM *ptr;
    for(i = 0; i < grid->nnodes_sur; i++) {
        // Get the first node in the vertical list (which is the surface node) 
        ptr = grid->vertical_list[i];
        nd = ptr->id;
        top_node=nd;
        // Find the increment for the surface displacement
        i3 = mod->fmap[nd] * 3;
        surf_dpl = sm->sol[i3+2];
        sw3->displacement[nd] += surf_dpl; 
        // Now, increment the velocities and the displacement for each node
        while(ptr->next != NULL) {
            nd = ptr->id;
            i3 = mod->fmap[nd] * 3;
            sw3->vel[nd].x += sm->sol[i3];
            sw3->vel[nd].y += sm->sol[i3+1];
            ptr = ptr->next;
        }
    }

#ifdef _DEBUG
    if (DEBUG) {
      printScreen_dble_array("hvel_inc :: displacement", sw3->displacement, grid->nnodes, __LINE__, __FILE__);   printf("\n");
      svect_printScreen_array("hvel_inc",sw3->vel,"velocity",grid->nnodes,__LINE__,__FILE__); printf("\n");
      printScreen_dble_array("hvel_inc :: dpl_perturbation", sw3->dpl_perturbation, grid->nnodes, __LINE__, __FILE__); printf("\n");
      printScreen_dble_array("hvel_inc :: node pressure", sw3->prs, grid->nnodes, __LINE__, __FILE__); printf("\n");
      printScreen_dble_array("hvel_inc :: node pressure +", sw3->prs_plus, grid->nnodes, __LINE__, __FILE__); printf("\n");
      printScreen_dble_array("hvel_inc :: node pressure -", sw3->prs_minus, grid->nnodes, __LINE__, __FILE__); printf("\n");

      ///*
      printf("+sol :: mdpt pressures\n");
      for (i=0; i<grid->num_midpts; i++) {
          for (j=0; j<5; j++) {
              printf("%20.10e",grid->midpt_list[i]->value[j]);
          }
          printf("\n");
      }
      //*/
    }        
#endif

#endif
}
