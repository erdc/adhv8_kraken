#include "global_header.h"

void read_bc_MP(SMODEL *mod, char *data) {
    
    char line[MAXLINE];           /* the input line */
    char *subdata = NULL;         /* the data after the second card   is read */
    char *subsubdata = NULL;
    
    int imat = 0, itrn = 0;
    SIO info = *(mod->io);    // alias
    int nmat = mod->nmat;                 // alias
    int ntrn = mod->ntransport; // alias
    SMAT mat;
    
    switch (parse_card(data, &subdata)) {
            
        case CARD_ML:
            imat = get_material_id(info, &subdata, nmat);
            mat = mod->mat[imat];
            if (mod->flag.SW_FLOW) {
                mat.sw->max_lev = read_int_field(info, &subdata);
                if (mat.sw->max_lev > 0) mod->flag.GRID_ADAPTION = ON;
            } else if (mod->flag.NS_FLOW) {
                mat.ns->max_lev = read_int_field(info, &subdata);
                if (mat.ns->max_lev > 0) mod->flag.GRID_ADAPTION = ON;
            }
            break;
            
        case CARD_COR:
            mod->flag.CORIOLIS = ON;
            imat = get_material_id(info, &subdata, nmat);
            mat = mod->mat[imat];
            if (mod->flag.SW_FLOW) {mat.sw->coriolis = read_dbl_field(info, &subdata);}
            else if (mod->flag.NS_FLOW) {mat.ns->coriolis = read_dbl_field(info, &subdata);}
            break;
            
        case CARD_EVS:
            imat = get_material_id(info, &subdata, nmat);
            mat = mod->mat[imat];
            if (mod->flag.SW_FLOW) {
                mat.sw->EVSF = YES;
                mat.sw->EEVF = NO;
                mat.sw->ev.xx = read_dbl_field(info, &subdata);
                mat.sw->ev.yy = read_dbl_field(info, &subdata);
                mat.sw->ev.xy = read_dbl_field(info, &subdata);
                if (mod->flag.SW3_FLOW || mod->flag.NS3_FLOW) {
                    mat.sw->ev.zz = read_dbl_field(info, &subdata);
                    mat.sw->ev.xz = read_dbl_field(info, &subdata);
                    mat.sw->ev.yz = read_dbl_field(info, &subdata);
                } else {
                    mat.sw->ev.zz = 0.;
                    mat.sw->ev.xz = 0.;
                    mat.sw->ev.yz = 0.;
                }
                mat.sw->eev_coef = UNSET_FLT;
            } else if (mod->flag.NS_FLOW) {
                mat.ns->EVSF = YES;
                mat.ns->EEVF = NO;
                mat.ns->ev.xx = read_dbl_field(info, &subdata);
                mat.ns->ev.yy = read_dbl_field(info, &subdata);
                mat.ns->ev.xy = read_dbl_field(info, &subdata);
                if (mod->flag.SW3_FLOW) {
                    mat.ns->ev.zz = read_dbl_field(info, &subdata);
                    mat.ns->ev.xz = read_dbl_field(info, &subdata);
                    mat.ns->ev.yz = read_dbl_field(info, &subdata);
                } else {
                    mat.ns->ev.zz = 0.;
                    mat.ns->ev.xz = 0.;
                    mat.ns->ev.yz = 0.;
                }
                mat.ns->eev_coef = UNSET_FLT;
            }
            break;
            
        case CARD_TUR:
            imat = get_material_id(info, &subdata, nmat);
            mat = mod->mat[imat];
            if (mod->flag.SW_FLOW) {
                mat.sw->turbulence_model_xy = read_int_field(info, &subdata);
                /* 0 is Smag, 1 is level 2 MY, 2 is level 2.5 MY, 3 is k-e GSavant */
                mat.sw->turbulence_model_z = read_int_field(info, &subdata);
                /* 0 is Smag, 1 is level 2 MY, 2 is level 2.5 MY, 3 is k-e GSavant */
                if (mat.sw->turbulence_model_z < 0) printf("For Material %d Turbulence Model is OFF\n", imat+1);
                printf("For Material %d Hor is %d Vert is %d\n", imat+1, mat.sw->turbulence_model_xy, mat.sw->turbulence_model_z);
                if (mat.sw->turbulence_model_xy == 0) {
                    mat.sw->smag_coeff = read_dbl_field(info, &subdata);
                    if (mat.sw->smag_coeff == 0) printf("For Material %d Horizontal Turbulence Model is OFF\n", imat+1);
                    printf("For Material %d Smag coeff is %lf\n", imat+1,mat.sw->smag_coeff);
                    if (mat.sw->smag_coeff < 0) tl_error("Smagorinski coefficient must be between 0.1 and 0.23: This is the third numeric value on the MP TUR card");
                }
                if (mat.sw->turbulence_model_z > 0) {
                    mat.sw->supression_func = read_int_field(info, &subdata);
                    if (mat.sw->supression_func == 0) printf("Buoyancy Supression in the Vertical  Turbulence Model is OFF\n");
                    printf("The Chosen buoyancy supression function is %d (1: Henderson-Sellers, 2: Munk-Anderson, 3: Kent-Pritchard, 4: Pritchard and 5: French-McCutcheon)\n", mat.sw->supression_func);
                    if (mat.sw->turbulence_model_z == 2) {
                        mat.sw->wall_func = read_int_field(info, &subdata);
                        printf("The wall function is %d for material %d (1: Mellor-Yamada 1982, 2: Burchard 1998, 3: Burchard 2001, 4: Blumberg et al. 1992)\n",mat.sw->wall_func, imat+1);
                    }
                    if (mat.sw->turbulence_model_z > 1) {
                        mat.sw->min_tke = read_dbl_field(info, &subdata);
                        mat.sw->min_tds = read_dbl_field(info, &subdata);
                        printf("For Material %d Minimum TKE is %G Minimum TDS is %G\n", imat+1, mat.sw->min_tke, mat.sw->min_tds);
                    }
                }
            } else if (mod->flag.NS_FLOW) {
                mat.ns->turbulence_model_xy = read_int_field(info, &subdata);
                /* 0 is Smag, 1 is level 2 MY, 2 is level 2.5 MY, 3 is k-e GSavant */
                mat.ns->turbulence_model_z = read_int_field(info, &subdata);
                /* 0 is Smag, 1 is level 2 MY, 2 is level 2.5 MY, 3 is k-e GSavant */
                if (mat.ns->turbulence_model_z < 0) printf("For Material %d Turbulence Model is OFF\n", imat+1);
                printf("For Material %d Hor is %d Vert is %d\n", imat+1, mat.ns->turbulence_model_xy, mat.ns->turbulence_model_z);
                if (mat.ns->turbulence_model_xy == 0) {
                    mat.ns->smag_coeff = read_dbl_field(info, &subdata);
                    if (mat.ns->smag_coeff == 0) printf("For Material %d Horizontal Turbulence Model is OFF\n", imat+1);
                    printf("For Material %d Smag coeff is %lf\n", imat+1,mat.ns->smag_coeff);
                    if (mat.ns->smag_coeff < 0) tl_error("Smagorinski coefficient must be between 0.1 and 0.23: This is the third numeric value on the MP TUR card");
                }
                if (mat.ns->turbulence_model_z > 0) {
                    mat.ns->supression_func = read_int_field(info, &subdata);
                    if (mat.ns->supression_func == 0) printf("Buoyancy Supression in the Vertical  Turbulence Model is OFF\n");
                    printf("The Chosen buoyancy supression function is %d (1: Henderson-Sellers, 2: Munk-Anderson, 3: Kent-Pritchard, 4: Pritchard and 5: French-McCutcheon)\n", mat.ns->supression_func);
                    if (mat.ns->turbulence_model_z == 2) {
                        mat.ns->wall_func = read_int_field(info, &subdata);
                        printf("The wall function is %d for material %d (1: Mellor-Yamada 1982, 2: Burchard 1998, 3: Burchard 2001, 4: Blumberg et al. 1992)\n",mat.ns->wall_func, imat+1);
                    }
                    if (mat.ns->turbulence_model_z > 1) {
                        mat.ns->min_tke = read_dbl_field(info, &subdata);
                        mat.ns->min_tds = read_dbl_field(info, &subdata);
                        printf("For Material %d Minimum TKE is %G Minimum TDS is %G\n", imat+1, mat.ns->min_tke, mat.ns->min_tds);
                    }
                }
            }
            break;
            
        case CARD_EEV:
            imat = get_material_id(info, &subdata, nmat);
            mat = mod->mat[imat];
            if (mod->flag.SW_FLOW) {
                mat.sw->EEVF = YES;
                mat.sw->eev_coef = read_dbl_field(info, &subdata);
                mat.sw->ev.xx = UNSET_FLT;
                mat.sw->ev.xy = UNSET_FLT;
                mat.sw->ev.yy = UNSET_FLT;
                mat.sw->eev_mode = read_int_field(info, &subdata);
                if (mat.sw->eev_mode != 1 && mat.sw->eev_mode != 2 && mat.sw->eev_mode != 3) {
                    printf("Selected EEV mode should be 1,2 or 3.  \n");
                    exit(0);
                }
            } else if (mod->flag.NS_FLOW) {
                mat.ns->EEVF = YES;
                mat.ns->eev_coef = read_dbl_field(info, &subdata);
                mat.ns->ev.xx = UNSET_FLT;
                mat.ns->ev.xy = UNSET_FLT;
                mat.ns->ev.yy = UNSET_FLT;
                mat.ns->eev_mode = read_int_field(info, &subdata);
                if (mat.ns->eev_mode != 1 && mat.ns->eev_mode != 2 && mat.ns->eev_mode != 3) {
                    printf("Selected EEV mode should be 1,2 or 3.  \n");
                    exit(0);
                }
            }
            break;
        case CARD_HYDCON:
            mod->green_ampt = TRUE;
            imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
            if (mod->flag.SW_FLOW) {mat.sw->hyd_conductivity =  read_dbl_field(info, &subdata);}
            else if (mod->flag.NS_FLOW) {mat.ns->hyd_conductivity =  read_dbl_field(info, &subdata);}
            break;
        case CARD_SUCK:
            mod->green_ampt = TRUE;
            imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
            if (mod->flag.SW_FLOW) {mat.sw->psi =  read_dbl_field(info, &subdata);}
            else if (mod->flag.NS_FLOW) {mat.ns->psi =  read_dbl_field(info, &subdata);}
            break;
        case CARD_SATDEP:
            mod->green_ampt = TRUE;
            imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
            if (mod->flag.SW_FLOW) {mat.sw->rooting_depth =  read_dbl_field(info, &subdata);}
            else if (mod->flag.NS_FLOW) {mat.ns->rooting_depth =  read_dbl_field(info, &subdata);}
            break;
        case CARD_DF:
            imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
            itrn = get_transport_id(info, &subdata, mod->ntransport);
            mat.trn[itrn].d_m = read_dbl_field(info, &subdata);
            break;
            
        case CARD_SRT:
            imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
            if (mod->flag.SW_FLOW) {
                mat.sw->refine_tolerance = read_dbl_field(info, &subdata);
                if (mat.sw->refine_tolerance < SMALL)
                    tl_error("Refinement tolerance is too small for the precision of the machine.");
                if (mat.sw->refine_tolerance == 0.0)
                    tl_error("Refinement tolerance must be greater than zero.");
            } else if (mod->flag.NS_FLOW) {
                mat.ns->refine_tolerance = read_dbl_field(info, &subdata);
                if (mat.ns->refine_tolerance < SMALL)
                    tl_error("Refinement tolerance is too small for the precision of the machine.");
                if (mat.ns->refine_tolerance == 0.0)
                    tl_error("Refinement tolerance must be greater than zero.");
            }
            break;
            
        case CARD_TRT:
            imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
            itrn = get_transport_id(info, &subdata, mod->ntransport);
            mat.trn[itrn].refine_tolerance = read_dbl_field(info, &subdata);
            if (mat.trn[itrn].refine_tolerance < SMALL)
                tl_error("Refinement tolerance is too small for the precision of the machine.");
            if (mat.trn[itrn].refine_tolerance == 0.0)
                tl_error("Refinement tolerance must be greater than zero.");
            break;
            
        case CARD_MU:
            mod->viscosity = read_dbl_field(info, &subdata);
            break;
            
        case CARD_MUC:
            mod->flag.MUC = ON;
            mod->manning_units_constant = read_dbl_field(info, &subdata);
            break;
            
        case CARD_G:
            mod->gravity = read_dbl_field(info, &subdata);
            break;
            
        case CARD_RHO:
            mod->density = read_dbl_field(info, &subdata);
            break;
            
        case CARD_DTL:
            mod->drying_lower_limit = read_dbl_field(info, &subdata);
            mod->drying_upper_limit = mod->drying_lower_limit;
            if (mod->drying_lower_limit < 0.) {
                printf("  drying limit is negative, %15.6e\n", mod->drying_lower_limit);
                printf("  it has been reset to 0.\n");
                mod->drying_lower_limit = 0.;
                mod->drying_upper_limit = 0.;
            }
            break;
            
        case CARD_WND:
            mod->flag.WIND = ON;
            printf("Winds are active in this simulation\n");
            switch (parse_card(subdata, &subsubdata)) {
                case CARD_STR: /* read type of wave stress transform (cjt) */
                    imat = get_material_id(info, &subsubdata, nmat); mat = mod->mat[imat];
                    if (mod->flag.SW_FLOW) {
                        mat.sw->wind_flag = read_int_field(info, &subsubdata);
                        if (mat.sw->wind_flag < 0 || mat.sw->wind_flag > 2)
                            tl_error("Wind stress calculation flag only supports options 0-2.");
                    } else if (mod->flag.NS_FLOW) {
                        mat.ns->wind_flag = read_int_field(info, &subsubdata);
                        if (mat.ns->wind_flag < 0 || mat.ns->wind_flag > 2)
                            tl_error("Wind stress calculation flag only supports options 0-2.");
                    }
                    break;
                    
                case CARD_ATT: /* read wind attenuation factor (cjt) */
                    imat = get_material_id(info, &subsubdata, nmat); mat = mod->mat[imat];
                    if (mod->flag.SW_FLOW) {
                        mat.sw->windatt = read_dbl_field(info, &subsubdata);
                        if (mat.sw->windatt < 0.0 || mat.sw->windatt > 1.0)
                            tl_error("Wind attenuation can only be a value from 0.0 to 1.0.");
                    } else if (mod->flag.NS_FLOW) {
                        mat.ns->windatt = read_dbl_field(info, &subsubdata);
                        if (mat.ns->windatt < 0.0 || mat.ns->windatt > 1.0)
                            tl_error("Wind attenuation can only be a value from 0.0 to 1.0.");
                    }
                    break;
                    
                default:
                    tl_error("\n BC ERROR :: MP WND card must be followed by a wind property card. \n");
                    break;
            }
            break;
        default:
            break;
            
    }
}
