/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     The FE engine driver for a given super model
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gary Brown
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gajanan Choudhary, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  \returns YES for a good calculation and NO for a bad calculation
 *
 * @param[in] SSUPER_MODEL *sm :: ONE  Super Model
 * \note CJT \:: matrix size will grow with refinement. It never shrinks with unrefinement!
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
extern double DEBUG_TIME;
static double total_time_mass_flux = 0.;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int fe_main(SAPP *app) {
    
    SSUPER_MODEL *sm = NULL;
    int isuperModel;
    
    //****************************************************************************
    //****************************************************************************
    // Loop over all SuperModels -------------------------------------------------
    for (isuperModel=0; isuperModel<app->nsuperModels; isuperModel++) {
        sm =  &(app->supermodel[isuperModel]);

        //Mark Note: Each super model = one linear system of equations
        
        //****************************************************************************
        // allocates FE matrix memory if necessary -----------------------------------
        fe_allocate_initialize_supermodel_system(sm);
        
        //****************************************************************************
        // initialize submodels of the SuperModel ------------------------------------
        fe_supermodel_init(sm);
        
        //****************************************************************************
        // solve the SuperModel for one dt -------------------------------------------
        if (fe_newton(sm, isuperModel, mod->grid->my_nnodes, mod->grid->nnodes, mod->grid->macro_nnodes,
#ifdef _MESSG
            sm->supersmpi,
#endif
            fe_sw2_init, fe_sw2_update, fe_sw2_resid, fe_sw2_load, fe_sw2_inc) == NO) {
            return (NO);
        }
        
        //Flux coupled models go somewhere here

        //****************************************************************************
        // calculate mass errors for all submodels -----------------------------------
        //fe_calculate_mass_errors(sm);
        
        //****************************************************************************
        // calculate approx continuity errors for adaption  --------------------------
        //fe_calculate_adaption_errors(sm);
    }
    
#ifdef _DEBUG
    if (DEBUG==ON) tl_check_all_pickets(__FILE__, __LINE__);
#endif
    
    return (SUCCESS);
}

