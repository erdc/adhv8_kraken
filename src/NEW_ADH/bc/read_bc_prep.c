/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file read_bc_prep.c This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for allocation counts only
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] *sm (SUPER_MODEL *)  a double pointer to an array of AdH supermodels
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_bc_prep(SMODEL_SUPER *mod) {
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *token;
    int idum;
    double fdum;
    
    SIO info = *(mod->io);   // alias
    
    assert(mod->io->bc.fp);
    
    // ++++++++++++++++++++++++++++++++++++++++++++++
    // read for counting and flagging first
    // ++++++++++++++++++++++++++++++++++++++++++++++
    while ((read = getline(&line, &len, mod->io->bc.fp)) != -1) {
        get_token(line,&token);
        if (strcmp(token, "OP") == 0) {
            get_next_token(&token);
            if (strcmp(token,"NS2") == 0) {
                mod->flag.NS_FLOW = TRUE;
                mod->flag.NS2_FLOW = TRUE;
            }
            if (strcmp(token,"NS3") == 0) {
                mod->flag.NS_FLOW = TRUE;
                mod->flag.NS3_FLOW = TRUE;
            }
            if (strcmp(token,"SW2") == 0) {
                mod->flag.SW_FLOW = TRUE;
                mod->flag.SW2_FLOW = TRUE;
            }
            if (strcmp(token,"SW3") == 0) {
                mod->flag.SW_FLOW = TRUE;
                mod->flag.SW3_FLOW = TRUE;
            }
            if (strcmp(token,"DW") == 0) {
                mod->flag.DIFFUSIVE_WAVE = TRUE;
            }
#ifdef _ADH_GROUNDWATER
            if (strcmp(token,"GW") == 0) {
                mod->flag.GW_FLOW = TRUE;
            }
#endif
            if (strcmp(token,"TRN") == 0) {
                mod->flag.DIFFUSIVE_WAVE = TRUE;
                mod->ntransport = get_next_token_int(&token);
                if (mod->ntransport > 0) {
                    mod->flag.TRANSPORT = TRUE;
                    if (mod->flag.SW2_FLOW) mod->flag.SW2_TRANSPORT = TRUE;
                    if (mod->flag.SW3_FLOW) mod->flag.SW3_TRANSPORT = TRUE;
                    if (mod->flag.NS2_FLOW) mod->flag.NS2_TRANSPORT = TRUE;
                    if (mod->flag.NS3_FLOW) mod->flag.NS3_TRANSPORT = TRUE;
#ifdef _ADH_GROUNDWATER
                    if (mod->flag.GW_FLOW)  mod->flag.GW_TRANSPORT  = TRUE;
#endif
                }
            }
            if (strcmp(token,"WAVES") == 0) {
                mod->flag.WAVE = TRUE;
            }
            if (strcmp(token,"WINDS") == 0) {
                mod->flag.WIND = TRUE;
            }
        }
        
#ifdef _ADH_STRUCTURES
        if (strcmp(token,"WEIR") == 0) {
            mod->nweir = get_next_token_int(&token);
        }
        if (strcmp(token,"FLAP") == 0) {
            mod->nflap = get_next_token_int(&token);
        }
        if (strcmp(token,"SLUICE") == 0) {
            mod->nsluice = get_next_token_int(&token);
        }
#endif
        if (strcmp(token,"NDS") == 0) {
            idum = get_next_token_int(&token);
            idum = get_next_token_int(&token);
            if (idum > mod->nstring) {
                mod->nstring = idum;
            }
        }
        if (strcmp(token,"EGS") == 0) {
            idum = get_next_token_int(&token);
            idum = get_next_token_int(&token);
            idum = get_next_token_int(&token);
            if (idum > mod->nstring) {
                mod->nstring = idum;
            }
        }
        if (strcmp(token,"MDS") == 0) {
            idum = get_next_token_int(&token);
            idum = get_next_token_int(&token);
            idum = get_next_token_int(&token);
            if (idum > mod->nstring) {
                mod->nstring = idum;
            }
        }
        if (strcmp(token,"SERIES") == 0) {
            idum = get_next_token_int(&token);
            if (idum > mod->nseries) {
                mod->nseries = idum;
            }
        }
        if (strcmp(token,"DB") == 0) {
            get_next_token(&token);
            if (strcmp(token,"DB") == 0) mod->nstring_rt = 1;
        }
        if (strcmp(token,"FOUT") == 0) {
            get_next_token(&token);
            if (strcmp(token,"GRID") == 0) mod->file_output.grid2dm = ON;
        }
    }
    rewind(mod->io->bc.fp);
}

