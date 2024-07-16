/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file fe_supermodel_init.c contains the routine that increments all solution variables in
 * a supermodel  */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
static int FILE_DEBUG = OFF;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file intializes a SuperModel's solution structures
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in,out] SSUPER_MODEL *sm :: ONE  Super Model
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_supermodel_init(SSUPER_MODEL *sm) {
    int i;
    //this will update sm solution structs that depend on time
    //is there a better way so if we add a variable we dont have to add bunches of things
    //everywhere?
    //for now any time sm is updated, this needs to be updated so maybe put this in sm folde
    for (i = 0; i < grid->nhead ; i ++){
        sm->older_head[i]  = sm->old_head[i];
        sm->old_head[i] = sm->head[i];
    }
    for (i = 0; i < grid->nvel2d; i++) {
        sm->older_vel2d[i].x = sm->old_vel2d[i].x;
        sm->older_vel2d[i].y = sm->old_vel2d[i].y;
        sm->old_vel2d[i].x = sm->vel2d[i].x;
        sm->old_vel2d[i].y = sm->vel2d[i].y;
    }
    for (i = 0; i < grid->nvel3d; i++) {
        sm->older_vel3d[i].x = sm->old_vel3d[i].x;
        sm->older_vel3d[i].y = sm->old_vel3d[i].y;
        sm->older_vel3d[i].z = sm->old_vel3d[i].z;
        sm->old_vel3d[i].x = sm->vel3d[i].x;
        sm->old_vel3d[i].y = sm->vel3d[i].y;
        sm->old_vel3d[i].z = sm->vel3d[i].x;
    }
    //maybe loop over number of constituents
    //need some sort of convention
    for (i = 0; i < sm->nconcentration ; i ++){
        sm->older_concentration[i]  = sm->old_concentration[i];
        sm->old_concentration[i] = sm->concentration[i];
    }
    for (i=0;i>sm->ndisplacement;i++){
        sm->older_displacement[i] = sm->old_displacement[i];
        sm->old_displacement[i] = sm->displacement[i];
    }
    for (i=0;i>sm->nprs;i++){
        sm->older_prs[i] = sm->old_prs[i];
        sm->old_prs[i] = sm->prs[i];
    }    
    for (i=0;i>sm->nc;i++){
        sm->older_c[i] = sm->old_c[i];
        sm->old_c[i] = sm->c[i];
    } 
}
