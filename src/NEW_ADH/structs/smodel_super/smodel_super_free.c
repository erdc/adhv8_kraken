#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees and AdH grid
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] pgrid          (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_super_free_array(SMODEL_SUPER *sm, int nSuper) {
    int i;
    
    for (i=0; i<nSuper; i++) {smodel_super_free(&(sm[i]));}
    sm = (SMODEL_SUPER *) tl_free(sizeof(SMODEL_SUPER), nSuper, sm);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees and AdH grid
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] pgrid          (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_super_free(SMODEL_SUPER *sm) {
    /* free file struct */
    //needs some work!
    //if(sm->io != NULL) sio_free(sm->io, mod->ntransport, mod->nlayers, mod->nsed);

    //free solution vectors
    sm->sol = (double *) tl_free(sizeof(double), *(sm->ndofs), sm->sol);
    sm->sol_old = (double *) tl_free(sizeof(double), *(sm->ndofs), sm->sol_old);
    sm->sol_older = (double *) tl_free(sizeof(double), *(sm->ndofs), sm->sol_older);


    //free dofmap
    if(sm->dof_map_local!=NULL){
        sm->dof_map_local = (int *) tl_free(sizeof(int), sm->grid->nnodes, sm->dof_map_local);
    }

    //free smat physics
    if (sm->elem1d_physics_mat!=NULL){
        smat_physics_free_array(sm->elem1d_physics_mat,sm->nphysics_mat_1d);
    }
    if(sm->elem2d_physics_mat!=NULL){
        smat_physics_free_array(sm->elem2d_physics_mat,sm->nphysics_mat_2d);
    }
    if(sm->elem3d_physics_mat!=NULL){
        smat_physics_free_array(sm->elem3d_physics_mat,sm->nphysics_mat_3d);
    }
    //if(sm->node_physics_mat!=NULL){
    //    smat_physics_free_array(sm->node_physics_mat,sm->nphysics_mat_node);
    //}
    if(sm->node_physics_mat!=NULL){
        sm->node_physics_mat = (SMAT_PHYSICS **) tl_free(sizeof(SMAT_PHYSICS *), sm->grid->nnodes, sm->node_physics_mat);
    }

    //free mat ids
    if(sm->elem1d_physics_mat_id!=NULL){
        sm->elem1d_physics_mat_id = (int *) tl_free(sizeof(int), sm->grid->nelems1d, sm->elem1d_physics_mat_id);
    }
    if(sm->elem2d_physics_mat_id!=NULL){
        sm->elem2d_physics_mat_id = (int *) tl_free(sizeof(int), sm->grid->nelems2d, sm->elem2d_physics_mat_id);
    }
    if(sm->elem3d_physics_mat_id!=NULL){
        sm->elem3d_physics_mat_id = (int *) tl_free(sizeof(int), sm->grid->nelems3d, sm->elem3d_physics_mat_id);
    }
//    if(sm->node_physics_mat_id!=NULL){
//        sm->node_physics_mat_id = (int *) tl_free(sizeof(int), sm->grid->nnodes, sm->node_physics_mat_id);
//    }

    //free bc masks and dirichlet data
    if(sm->bc_mask!=NULL){
        sm->bc_mask = (int *) tl_free(sizeof(int), *(sm->ndofs), sm->bc_mask);
    }
    if(sm->dirichlet_data!=NULL){
        sm->dirichlet_data = (double *) tl_free(sizeof(double), *(sm->ndofs), sm->dirichlet_data);
    }

}
