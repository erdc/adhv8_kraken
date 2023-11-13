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
 * \note CJT\:: flag = 0 :: fully wet || flag = 1 :: some nodes are dry || flag = 2 :: all nodes dry
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_sw2_wdflag(SSW_2D *sw2d, SGRID *grid) {

    int i, ie, ii;                  /* loop counter */
    int gsize;
    int *wd_tmp;
    int e_wet;
    
    gsize = sizeof(int) * grid->nnodes;
    if (gsize > sw2d->vwork_size) {
        sw2d->vwork = (void *) tl_realloc(1, gsize, sw2d->vwork_size, sw2d->vwork);
        sw2d->vwork_size = gsize;
    }
    wd_tmp = (int *) sw2d->vwork;
    
    for (i = 0; i < grid->nnodes; i++)
        wd_tmp[i] = 0;
    
    for (ie = 0; ie < grid->nelems2d; ie++) {
        e_wet = 0;
        for (i = 0; i < NDONTRI; i++) {
            if (sw2d->head[ grid->elem2d[ie].nodes[i] ] > 0.0)
                e_wet++;
        }
        if (e_wet != NDONTRI) {
            for (i = 0; i < NDPRFC; i++)
                wd_tmp[ grid->elem2d[ie].nodes[i] ] = 1;
        }
    }
    
#ifdef _MESSG
    comm_update_int(wd_tmp, 1, grid->smpi);
#endif
    
    for (ie = 0; ie < grid->nelems2d; ie++) {
        grid->wd_flag[ie] = 0;
        
        for (i = 0; i < NDONTRI; i++) {
            if (wd_tmp[ grid->elem2d[ie].nodes[i] ] != 0) {
                grid->wd_flag[ie] = 1;
            }
        }
        
        if (sw2d->head[ grid->elem2d[ie].nodes[0] ] < 0 &&
            sw2d->head[ grid->elem2d[ie].nodes[1] ] < 0 &&
            sw2d->head[ grid->elem2d[ie].nodes[2] ] < 0 &&
            sw2d->old_head[ grid->elem2d[ie].nodes[0] ] < 0 &&
            sw2d->old_head[ grid->elem2d[ie].nodes[1] ] < 0 &&
            sw2d->old_head[ grid->elem2d[ie].nodes[2] ] < 0)
            grid->wd_flag[ie] = 2;
    }
    
}
