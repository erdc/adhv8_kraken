/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     The FE engine driver for a given Design model, advanves the model one time step
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gary Brown
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gajanan Choudhary, Ph.D.
 *  \author    Mark Loveland, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  \returns YES for a good calculation and NO for a bad calculation
 *
 * @param[in,out] DESIGN_MODEL *dm :: ONE  Design Model
 * \note CJT \:: matrix size will grow with refinement. It never shrinks with unrefinement!
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
extern double DEBUG_TIME;
static double total_time_mass_flux = 0.;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int fe_main(DESIGN_MODEL *dm) {
    
    SSUPER_MODEL *sm = NULL;
    int isuperModel;
    
    //initialize grid stuff?
    SGRID grid = dm->grid;
    SMAT mat = dm->mat;


    //****************************************************************************
    //****************************************************************************
    // Loop over all SuperModels -------------------------------------------------
    for (isuperModel=0; isuperModel<dm->nsuperModels; isuperModel++) {
        sm =  &(dm->supermodel[isuperModel]);


        //****************************************************************************
        // initialize submodels of the SuperModel ------------------------------------
        //importantly this will be used to create the dofmap, ndof
        //this will occur outside of fe main if we stick with how it was done last version
        //update_elem_physics(sm,grid);

        //after maps are set up then reallocate/initialize
        //and allocate any model specific structures for convenient i/o
        fe_supermodel_init(sm);   

        //****************************************************************************
        // allocates memory to solve linear system if necessary ----------------------
        fe_allocate_initialize_linear_system(sm);

        //potentially pass info for coupling between super models?
        //how to we want to pass info from one sm to another?
        //we restrict this to full mesh <-> full surface of mesh


        //****************************************************************************
        // solve the SuperModel for one dt -------------------------------------------

        //need to put grid as surface grid if only surface is defined maybe
        if (fe_newton(sm, isuperModel, grid, mat) == NO) {
            return (NO);
        }
        
        //****************************************************************************
        // calculate mass errors for all submodels -----------------------------------
        //fe_calculate_mass_errors(sm);
        
        //****************************************************************************
        // calculate approx continuity errors for adaption  --------------------------
        //fe_calculate_adaption_errors(sm);

    }
    //free grid
    
#ifdef _DEBUG
    if (DEBUG==ON) tl_check_all_pickets(__FILE__, __LINE__);
#endif
    
    return (SUCCESS);
}

