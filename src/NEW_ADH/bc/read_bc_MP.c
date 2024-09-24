/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_MP.c  This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for reading material properties
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] *sm (SUPER_MODEL *)  a double pointer to an array of AdH supermodels
 * @param[in] **token (CHAR) a BC file line string token
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_bc_MP(SMODEL_SUPER *mod, char **token) {
    
    int imat = 0, itrn = 0;
    int nmat = mod->nmat;                 // alias
    int ntrn = mod->ntransport; // alias
    SMAT mat;
    
    if (strcmp(*token, "ML") == 0) {
        imat = get_id(token,mod->nstring,"ML String Error\n");
        mat = mod->mat[imat];
        if (mod->flag.SW_FLOW) {
            mat.sw->max_lev = get_next_token_int(token);
            if (mat.sw->max_lev > 0) mod->flag.GRID_ADAPTION = ON;
        } else if (mod->flag.NS_FLOW) {
            mat.ns->max_lev = get_next_token_int(token);
            if (mat.ns->max_lev > 0) mod->flag.GRID_ADAPTION = ON;
        }
        
    } else if (strcmp(*token, "COR") == 0) {
        mod->flag.CORIOLIS = ON;
        imat = get_id(token,mod->nstring,"COR String Error\n");
        mat = mod->mat[imat];
        if (mod->flag.SW_FLOW) {mat.sw->coriolis = get_next_token_dbl(token);}
        else if (mod->flag.NS_FLOW) {mat.ns->coriolis = get_next_token_dbl(token);}
        
    } else if (strcmp(*token, "EVS") == 0) {
        imat = get_id(token,mod->nstring,"EVS String Error\n");
        mat = mod->mat[imat];
        if (mod->flag.SW_FLOW) {
            mat.sw->EVSF = YES;
            mat.sw->EEVF = NO;
            mat.sw->ev.xx = get_next_token_dbl(token);
            mat.sw->ev.yy = get_next_token_dbl(token);
            mat.sw->ev.xy = get_next_token_dbl(token);
            if (mod->flag.SW3_FLOW || mod->flag.NS3_FLOW) {
                mat.sw->ev.zz = get_next_token_dbl(token);
                mat.sw->ev.xz = get_next_token_dbl(token);
                mat.sw->ev.yz = get_next_token_dbl(token);
            } else {
                mat.sw->ev.zz = 0.;
                mat.sw->ev.xz = 0.;
                mat.sw->ev.yz = 0.;
            }
            mat.sw->eev_coef = UNSET_FLT;
        } else if (mod->flag.NS_FLOW) {
            mat.ns->EVSF = YES;
            mat.ns->EEVF = NO;
            mat.ns->ev.xx = get_next_token_dbl(token);
            mat.ns->ev.yy = get_next_token_dbl(token);
            mat.ns->ev.xy = get_next_token_dbl(token);
            if (mod->flag.SW3_FLOW) {
                mat.ns->ev.zz = get_next_token_dbl(token);
                mat.ns->ev.xz = get_next_token_dbl(token);
                mat.ns->ev.yz = get_next_token_dbl(token);
            } else {
                mat.ns->ev.zz = 0.;
                mat.ns->ev.xz = 0.;
                mat.ns->ev.yz = 0.;
            }
            mat.ns->eev_coef = UNSET_FLT;
        }
        
    } else if (strcmp(*token, "TUR") == 0) {
        imat = get_id(token,mod->nstring,"TUR String Error\n");
        mat = mod->mat[imat];
        if (mod->flag.SW_FLOW) {
            mat.sw->turbulence_model_xy = get_next_token_int(token);
            /* 0 is Smag, 1 is level 2 MY, 2 is level 2.5 MY, 3 is k-e GSavant */
            mat.sw->turbulence_model_z = get_next_token_int(token);
            /* 0 is Smag, 1 is level 2 MY, 2 is level 2.5 MY, 3 is k-e GSavant */
            if (mat.sw->turbulence_model_z < 0) printf("For Material %d Turbulence Model is OFF\n", imat+1);
            printf("For Material %d Hor is %d Vert is %d\n", imat+1, mat.sw->turbulence_model_xy, mat.sw->turbulence_model_z);
            if (mat.sw->turbulence_model_xy == 0) {
                mat.sw->smag_coeff = get_next_token_dbl(token);
                if (mat.sw->smag_coeff == 0) printf("For Material %d Horizontal Turbulence Model is OFF\n", imat+1);
                printf("For Material %d Smag coeff is %lf\n", imat+1,mat.sw->smag_coeff);
                if (mat.sw->smag_coeff < 0) tl_error("Smagorinski coefficient must be between 0.1 and 0.23: This is the third numeric value on the MP TUR card");
            }
            if (mat.sw->turbulence_model_z > 0) {
                mat.sw->supression_func = get_next_token_int(token);
                if (mat.sw->supression_func == 0) printf("Buoyancy Supression in the Vertical  Turbulence Model is OFF\n");
                printf("The Chosen buoyancy supression function is %d (1: Henderson-Sellers, 2: Munk-Anderson, 3: Kent-Pritchard, 4: Pritchard and 5: French-McCutcheon)\n", mat.sw->supression_func);
                if (mat.sw->turbulence_model_z == 2) {
                    mat.sw->wall_func = get_next_token_int(token);
                    printf("The wall function is %d for material %d (1: Mellor-Yamada 1982, 2: Burchard 1998, 3: Burchard 2001, 4: Blumberg et al. 1992)\n",mat.sw->wall_func, imat+1);
                }
                if (mat.sw->turbulence_model_z > 1) {
                    mat.sw->min_tke = get_next_token_dbl(token);
                    mat.sw->min_tds = get_next_token_dbl(token);
                    printf("For Material %d Minimum TKE is %G Minimum TDS is %G\n", imat+1, mat.sw->min_tke, mat.sw->min_tds);
                }
            }
        } else if (mod->flag.NS_FLOW) {
            mat.ns->turbulence_model_xy = get_next_token_int(token);
            /* 0 is Smag, 1 is level 2 MY, 2 is level 2.5 MY, 3 is k-e GSavant */
            mat.ns->turbulence_model_z = get_next_token_int(token);
            /* 0 is Smag, 1 is level 2 MY, 2 is level 2.5 MY, 3 is k-e GSavant */
            if (mat.ns->turbulence_model_z < 0) printf("For Material %d Turbulence Model is OFF\n", imat+1);
            printf("For Material %d Hor is %d Vert is %d\n", imat+1, mat.ns->turbulence_model_xy, mat.ns->turbulence_model_z);
            if (mat.ns->turbulence_model_xy == 0) {
                mat.ns->smag_coeff = get_next_token_dbl(token);
                if (mat.ns->smag_coeff == 0) printf("For Material %d Horizontal Turbulence Model is OFF\n", imat+1);
                printf("For Material %d Smag coeff is %lf\n", imat+1,mat.ns->smag_coeff);
                if (mat.ns->smag_coeff < 0) tl_error("Smagorinski coefficient must be between 0.1 and 0.23: This is the third numeric value on the MP TUR card");
            }
            if (mat.ns->turbulence_model_z > 0) {
                mat.ns->supression_func = get_next_token_int(token);
                if (mat.ns->supression_func == 0) printf("Buoyancy Supression in the Vertical  Turbulence Model is OFF\n");
                printf("The Chosen buoyancy supression function is %d (1: Henderson-Sellers, 2: Munk-Anderson, 3: Kent-Pritchard, 4: Pritchard and 5: French-McCutcheon)\n", mat.ns->supression_func);
                if (mat.ns->turbulence_model_z == 2) {
                    mat.ns->wall_func = get_next_token_int(token);
                    printf("The wall function is %d for material %d (1: Mellor-Yamada 1982, 2: Burchard 1998, 3: Burchard 2001, 4: Blumberg et al. 1992)\n",mat.ns->wall_func, imat+1);
                }
                if (mat.ns->turbulence_model_z > 1) {
                    mat.ns->min_tke = get_next_token_dbl(token);
                    mat.ns->min_tds = get_next_token_dbl(token);
                    printf("For Material %d Minimum TKE is %G Minimum TDS is %G\n", imat+1, mat.ns->min_tke, mat.ns->min_tds);
                }
            }
            
        } else if (strcmp(*token, "EEV") == 0) {
            imat = get_id(token,mod->nstring,"EEV String Error\n");
            mat = mod->mat[imat];
            if (mod->flag.SW_FLOW) {
                mat.sw->EEVF = YES;
                mat.sw->eev_coef = get_next_token_dbl(token);
                mat.sw->ev.xx = UNSET_FLT;
                mat.sw->ev.xy = UNSET_FLT;
                mat.sw->ev.yy = UNSET_FLT;
                mat.sw->eev_mode = get_next_token_int(token);
                if (mat.sw->eev_mode != 1 && mat.sw->eev_mode != 2 && mat.sw->eev_mode != 3) {
                    printf("Selected EEV mode should be 1,2 or 3.  \n");
                    exit(0);
                }
            } else if (mod->flag.NS_FLOW) {
                mat.ns->EEVF = YES;
                mat.ns->eev_coef = get_next_token_dbl(token);
                mat.ns->ev.xx = UNSET_FLT;
                mat.ns->ev.xy = UNSET_FLT;
                mat.ns->ev.yy = UNSET_FLT;
                mat.ns->eev_mode = get_next_token_int(token);
                if (mat.ns->eev_mode != 1 && mat.ns->eev_mode != 2 && mat.ns->eev_mode != 3) {
                    printf("Selected EEV mode should be 1,2 or 3.  \n");
                    exit(0);
                }
            }
            
        } else if (strcmp(*token, "HYDCON") == 0) {
            imat = get_id(token,mod->nstring,"HYDCON String Error\n");
            mat = mod->mat[imat];
            mod->green_ampt = TRUE;
            if (mod->flag.SW_FLOW) {mat.sw->hyd_conductivity =  get_next_token_dbl(token);}
            else if (mod->flag.NS_FLOW) {mat.ns->hyd_conductivity =  get_next_token_dbl(token);}
            
        } else if (strcmp(*token, "SUCK") == 0) {
            imat = get_id(token,mod->nstring,"SUCK String Error\n");
            mat = mod->mat[imat];
            mod->green_ampt = TRUE;
            if (mod->flag.SW_FLOW) {mat.sw->psi =  get_next_token_dbl(token);}
            else if (mod->flag.NS_FLOW) {mat.ns->psi  =  get_next_token_dbl(token);}
            
        } else if (strcmp(*token, "SATDEP") == 0) {
            imat = get_id(token,mod->nstring,"SATDEP String Error\n");
            mat = mod->mat[imat];
            mod->green_ampt = TRUE;
            if (mod->flag.SW_FLOW) {mat.sw->rooting_depth =  get_next_token_dbl(token);}
            else if (mod->flag.NS_FLOW) {mat.sw->rooting_depth  =  get_next_token_dbl(token);}
            
        } else if (strcmp(*token, "DF") == 0) {
            imat = get_id(token,mod->nstring,"DF String Error\n");
            itrn = get_id(token,mod->nstring,"DF TRN ID String Error \n");
            mat = mod->mat[imat];
            mat.trn[itrn].d_m = get_next_token_dbl(token);
            
        } else if (strcmp(*token, "SRT") == 0) {
            imat = get_id(token,mod->nstring,"SRT String Error\n");
            if (mod->flag.SW_FLOW) {
                mat.sw->refine_tolerance = get_next_token_dbl(token);
                if (mat.sw->refine_tolerance < SMALL)
                    tl_error("Refinement tolerance is too small for the precision of the machine.");
                if (mat.sw->refine_tolerance == 0.0)
                    tl_error("Refinement tolerance must be greater than zero.");
            } else if (mod->flag.NS_FLOW) {
                mat.ns->refine_tolerance = get_next_token_dbl(token);
                if (mat.ns->refine_tolerance < SMALL)
                    tl_error("Refinement tolerance is too small for the precision of the machine.");
                if (mat.ns->refine_tolerance == 0.0)
                    tl_error("Refinement tolerance must be greater than zero.");
            }
            
        } else if (strcmp(*token, "TRT") == 0) {
            imat = get_id(token,mod->nstring,"TRT String Error\n");
            itrn = get_id(token,mod->nstring,"TRT TRN ID String Error \n");
            mat = mod->mat[imat];
            mat.trn[itrn].refine_tolerance = get_next_token_dbl(token);
            if (mat.trn[itrn].refine_tolerance < SMALL)
                tl_error("Refinement tolerance is too small for the precision of the machine.");
            if (mat.trn[itrn].refine_tolerance == 0.0)
                tl_error("Refinement tolerance must be greater than zero.");
            
        } else if (strcmp(*token, "MU") == 0) {
            mod->viscosity = get_next_token_dbl(token);
            
        } else if (strcmp(*token, "MUC") == 0) {
            mod->flag.MUC = ON;
            mod->manning_units_constant = get_next_token_dbl(token);
            
        } else if (strcmp(*token, "G") == 0) {
            mod->gravity = get_next_token_dbl(token);
            
        } else if (strcmp(*token, "RHO") == 0) {
            mod->density = get_next_token_dbl(token);
            
        } else if (strcmp(*token, "DTL") == 0) {
            mod->drying_lower_limit = get_next_token_dbl(token);
            mod->drying_upper_limit = mod->drying_lower_limit;
            if (mod->drying_lower_limit < 0.) {
                printf("  drying limit is negative, %15.6e\n", mod->drying_lower_limit);
                printf("  it has been reset to 0.\n");
                mod->drying_lower_limit = 0.;
                mod->drying_upper_limit = 0.;
            }
            
        } else if (strcmp(*token, "WND") == 0) {
            mod->flag.WIND = ON;
            printf("Winds are active in this simulation\n");
            get_next_token(token);
            if (strcmp(*token, "STR") == 0) {
                imat = get_id(token,mod->nstring,"WND STR String Error\n");
                mat = mod->mat[imat];
                if (mod->flag.SW_FLOW) {
                    mat.sw->wind_flag = get_next_token_int(token);
                    if (mat.sw->wind_flag < 0 || mat.sw->wind_flag > 2)
                        tl_error("Wind stress calculation flag only supports options 0-2.");
                } else if (mod->flag.NS_FLOW) {
                    mat.ns->wind_flag = get_next_token_int(token);
                    if (mat.ns->wind_flag < 0 || mat.ns->wind_flag > 2)
                        tl_error("Wind stress calculation flag only supports options 0-2.");
                }
            } else if (strcmp(*token, "ATT") == 0) {
                imat = get_id(token,mod->nstring,"WND STR String Error\n");
                mat = mod->mat[imat];
                if (mod->flag.SW_FLOW) {
                    mat.sw->windatt = get_next_token_dbl(token);
                    if (mat.sw->windatt < 0.0 || mat.sw->windatt > 1.0)
                        tl_error("Wind attenuation can only be a value from 0.0 to 1.0.");
                } else if (mod->flag.NS_FLOW) {
                    mat.ns->windatt = get_next_token_dbl(token);
                    if (mat.ns->windatt < 0.0 || mat.ns->windatt > 1.0)
                        tl_error("Wind attenuation can only be a value from 0.0 to 1.0.");
                }
            }
        }
        
        
        
    }
}
