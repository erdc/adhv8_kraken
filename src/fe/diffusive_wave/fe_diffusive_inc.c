/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Increments the diffusive solutions after a Newton iterate.
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Gary Brown, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \author    Gajanan Choudhary, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout]  mod (SMODEL *) pointer to the model struct
 *
 * \note Also calculates area-weighted nodal velocities from elemental velocities.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_diffusive_inc(SSUPER_MODEL *sm, int isubModel) {

    SMODEL *mod = &(sm->submodel[isubModel]);

    int i;
#ifdef _PETSC
    PetscScalar const *values;
    PetscInt ierr;

    ierr = VecGetArrayRead(sm->sol,&(values));
    for (i = 0; i < mod->grid->my_nnodes; i++) {
        mod->sw->d2->head[i] += values[mod->fmap[i]*mod->max_nsys];
    }
    ierr = VecRestoreArrayRead(sm->sol,&(values));
#else
    for (i = 0; i < mod->grid->nnodes; i++) {
        mod->sw->d2->head[i] += sm->sol[mod->fmap[i]];
    }
#endif

    sarray_init_dbl(mod->sw->d2->darray,mod->grid->nnodes);
    svect2d_init_array(mod->sw->d2->vel,mod->grid->nnodes);
    SELEM_2D *elem2d;
    SNODE nodes2d[NDONQUAD];
    SVECT2D elem2d_vel[NDONQUAD];
    SVECT elem2d_nds[NDONQUAD];
    double node2d_z[NDONQUAD];
    double elem2d_head[NDONQUAD];
    SVECT2D vv;
    double elem2d_avg_depth = 0., area = 0., djac = 0.;
    int isTriangle;
    for (int j=0; j<mod->grid->nelems2d; j++){
        elem2d = &(mod->grid->elem2d[j]);
        if (elem2d->bflag==0 || elem2d->bflag==UNSET_INT){

            djac = fabs(elem2d->djac);
    
            global_to_local_dbl(mod->sw->d2->head, elem2d_head, elem2d->nodes, elem2d->nnodes);
            for (i=0; i<elem2d->nnodes; i++) {
                snode_copy(&(nodes2d[i]), mod->grid->node[elem2d->nodes[i]]);
                node2d_z[i] = nodes2d[i].z;
            }
            snode2svect(nodes2d, elem2d_nds, elem2d->nnodes);
    
            elem2d_avg_depth = 0., area = 0.;
            if (elem2d->nnodes == NDONTRI) isTriangle = YES; else isTriangle = NO;
            if (isTriangle == TRUE) {
                area = elem2d->djac;
                elem2d_avg_depth = integrate_triangle_f(1.,1.,elem2d_head); // cjt :: djac = area cancel here
            } else {
                area = integrate_quadrilateral_area(elem2d_nds, 1.);
                elem2d_avg_depth = integrate_quadrilateral_f(elem2d_nds,1./area,elem2d_head);
            }
    
            vv = getDiffusiveWaveVelocities(elem2d->nnodes, elem2d_head, elem2d_avg_depth, node2d_z, elem2d->grad_shp, mod->str_values[elem2d->string].fterms.manningsn, mod->manning_units_constant, elem2d_vel);
#ifdef _DEBUG
        //printf("\nElem2d[%i]: vavg.x=%f, vavg.y=%f", i, vv.x, vv.y);
#endif
        }
        else {
            djac = 0.0;
            vv.x=0.0;
            vv.y=0.0;
            svect2d_init_array(elem2d_vel,NDONQUAD);
        }

        for (int k=0; k<elem2d->nnodes; k++){
            mod->sw->d2->vel[elem2d->nodes[k]].x += djac*vv.x; //elem2d_vel[k].x; //
            mod->sw->d2->vel[elem2d->nodes[k]].y += djac*vv.y; //elem2d_vel[k].y; //
            mod->sw->d2->darray[elem2d->nodes[k]] += djac;
        }
    }

    for (i = 0; i < mod->grid->nnodes; i++) {
        if (mod->sw->d2->darray[i]>SMALL){
            mod->sw->d2->vel[i].x /= mod->sw->d2->darray[i];
            mod->sw->d2->vel[i].y /= mod->sw->d2->darray[i];
        }
        else {
            assert(fabs(mod->sw->d2->vel[i].x)<=NOT_QUITE_SMALL);
            assert(fabs(mod->sw->d2->vel[i].y)<=NOT_QUITE_SMALL);
            mod->sw->d2->vel[i].x = 0.;
            mod->sw->d2->vel[i].y = 0.;
        }
        if (mod->sw->d2->head[i]<=0.0){
            mod->sw->d2->vel[i].x=0.0;
            mod->sw->d2->vel[i].y=0.0;
        }
    }
}
