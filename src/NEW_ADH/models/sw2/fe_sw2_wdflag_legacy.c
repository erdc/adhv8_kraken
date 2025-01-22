/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_sw2_wdflag_legacy.c This file is the routine that
 *  updates the w/d flag on each element       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Assigns a wet/dry flag to all element.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in,out] grid (SGRID *) a grid.  On return, assigned grid->wd_flag[ie]
 * @param[in]     sw2d (SSW_2D *) a shallow water 2D model
 *
 * \note CJT\:: Label wetting and drying AND neighboring elements as wet/dry elements
 * \note CJT\:: flag = 0 :: fully wet and no surrounding elements contain dry node
 *              || flag = 1 :: at least one is dry on the element or neighboring element
 *              || flag = 2 :: all nodes dry on the element at the last 2 time steps
 * \note : Still need to write routine to clear wd_tmp at end of time loop
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "adh.h"
static int isize = 0;           /* the size of the array */
static int *wd_tmp = NULL;          /* nodal flags just used in this routine */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//needs to be called inside an init routine
void fe_sw2_wdflag_legacy(SMODEL_SUPER *sm) {

    int i, ie, ii;                  /* loop counter */
    int gsize;
    int e_wet;
    int isize_prev;
    SGRID *grid = sm->grid;
    SSW *sw2d = sm->sw;
    int nnode_sw2 = sw2d->dvar.nnode_dvar;
    //if nnode!=nnode_sw2, do something more complicated
    int nnode = sm->grid->nnodes;
    int *wd_flag = sw2d->dvar.elem_flags[sw2d->WD_FLAG];
    /////////OR JUST CHECK sm->isSimple flag if we use that/////////
    //Works for simple and complex case
    //reallocate if grid changes
    if (isize < nnode_sw2) {
        isize_prev = isize;
        isize = nnode_sw2;
        wd_tmp = (int *) tl_realloc(sizeof(int), isize, isize_prev, wd_tmp);
    }
    sarray_init_int(wd_tmp,isize);
//    for (i = 0; i < nnode_sw2; i++){
//        printf("initallizing wd_tmp\n");
//        wd_tmp[i] = 0;
//    }
    //simple case only:
    //again, maybe better to check isSimple or some flag we save
    if (nnode == nnode_sw2){
        for (ie = 0; ie < grid->nelems2d; ie++) {
            e_wet = 0;
            for (i = 0; i < NDONTRI; i++) {
                //water head at node i of elem ie
                if (sm->sol[ grid->elem2d[ie].nodes[i]*3 ] > 0.0)
                    e_wet++;
            }
            if (e_wet != NDONTRI) {
                for (i = 0; i < NDPRFC; i++)
                    wd_tmp[ grid->elem2d[ie].nodes[i] ] = 1;
            }
        }


        for (ie = 0; ie < grid->nelems2d; ie++) {
            wd_flag[ie] = 0;
        
            for (i = 0; i < NDONTRI; i++) {
                if (wd_tmp[ grid->elem2d[ie].nodes[i] ] != 0) {
                    wd_flag[ie] = 1;
                }
            }
        
            if (sm->sol[ grid->elem2d[ie].nodes[0]*3 ] < 0 &&
                sm->sol[ grid->elem2d[ie].nodes[1]*3 ] < 0 &&
                sm->sol[ grid->elem2d[ie].nodes[2]*3 ] < 0 &&
                sm->sol_old[ grid->elem2d[ie].nodes[0]*3 ] < 0 &&
                sm->sol_old[ grid->elem2d[ie].nodes[1]*3 ] < 0 &&
                sm->sol_old[ grid->elem2d[ie].nodes[2]*3 ] < 0)
                wd_flag[ie] = 2;
        }
        //sarray_printScreen_dbl(sm->sol,grid->nnodes*3, "sol");
        //sarray_printScreen_int(wd_flag,grid->nelems2d, "wd flag");
    }else{
        //LEAVING INCOMPLETE FOR NOW UNTIL WE FINALIZE MAPS
        assert(1==0);
        //Complex case; will do later, need to get something working first
        //need to use all of our ivar and dvar maps
        //NOTE: SHOULD LOOP OVER SUBGRID THAT SW2 IS ACTIVE ON
        //DO WE WANT BETTER WAY TO HANDLE THIS???
//      for (ie = 0; ie < grid->nelems2d; ie++) {
//            e_wet = 0;
//            for (i = 0; i < NDONTRI; i++) {
//                if (sw2d->head[ grid->elem2d[ie].nodes[i] ] > 0.0)
//                    e_wet++;
//            }
//            if (e_wet != NDONTRI) {
//                for (i = 0; i < NDPRFC; i++)
//                    wd_tmp[ grid->elem2d[ie].nodes[i] ] = 1;
//            }
//        }
//        for (ie = 0; ie < grid->nelems2d; ie++) {
//            grid->wd_flag[ie] = 0;
//        
//            for (i = 0; i < NDONTRI; i++) {
//                if (wd_tmp[ grid->elem2d[ie].nodes[i] ] != 0) {
//                    grid->wd_flag[ie] = 1;
//                }
//            }
//        
//            if (sw2d->head[ grid->elem2d[ie].nodes[0] ] < 0 &&
//                sw2d->head[ grid->elem2d[ie].nodes[1] ] < 0 &&
//                sw2d->head[ grid->elem2d[ie].nodes[2] ] < 0 &&
//                sw2d->old_head[ grid->elem2d[ie].nodes[0] ] < 0 &&
//                sw2d->old_head[ grid->elem2d[ie].nodes[1] ] < 0 &&
//                sw2d->old_head[ grid->elem2d[ie].nodes[2] ] < 0)
//                grid->wd_flag[ie] = 2;
//        }


    }   
#ifdef _MESSG
    comm_update_int(wd_tmp, 1, grid->smpi);
#endif  
}
