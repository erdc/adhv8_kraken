#include "global_header.h"

#ifdef _DEBUG
static int DEBUG = OFF;
#endif

//*******************************************************************************************//
//*******************************************************************************************//
//*******************************************************************************************//
// This reads the file to count the number of sediment grains for smodel initialization
void read_bc_SEDLIB_prep(SMODEL *mod, char *data) {
    
    char line[MAXLINE];           /* the input line */
    char *subdata = NULL;         /* the data after the second card   is read */
    char *subsubdata = NULL;
    
    int i, j, k, imat, istr;
    int nsed = 0, nclay = 0, nsand = 0, nsfssi = 0, ncti = 0, nlayers = 0, nsilt = 0;
 
    SIO *io = mod->io;  // alias
    int super = strlen(io->sup.filename);

    // was code compiled with SEDLIB?
#ifndef _SEDIMENT
    printf("**** WARNING: SEDLIB card used in bc file, but code was not built with sedlib ****\n");
    return;
#else
    printf("-- reading sedlib input file ---\n");
#endif

    mod->flag.SEDIMENT = ON;
    mod->flag.SEDLIB = ON;

    // open and initialize sedlib input file
    init_adh_in_file(&(io->sedfile));
    sprintf(io->sedfile.filename, "%s.sedlib", io->proj_name);
    open_input_file(&(io->sedfile), "sedlib file", super);
    assert(io->sedfile.fp);

    //**************************************************************************************//
    //**************************************************************************************//
    // read sedlib file to get counts
    while (fgets(line, MAXLINE, io->sedfile.fp) != NULL) {
        io_save_line(io, io->sedfile.fp, io->sedfile.filename, line);
        if (strip_comments(line) <= 1) {    /* ignore empty ('\n') or comment lines */
            continue;
        }
        switch (parse_card(line, &data)) {
            case CARD_GRAIN:
                switch (parse_card(data, &subdata)) {
                    case CARD_CLA:
                        nsed++;
                        nclay++;
                        break;
                        
                    case CARD_SND:
                        nsed++;
                        nsand++;
                        break;
                }
                break;
                
            case CARD_MP:
                switch (parse_card(data, &subdata)) {
                    case CARD_NBL:
                        nlayers = read_int_field(*io, &subdata);
                        if(nlayers < 1) {
                            tl_error("Maximum number of bed layers must be greater than 0.");
                        }
                        break;
                
                    case CARD_NSF:
                        nsfssi = read_int_field(*io, &subdata);
                        if (nsfssi < 2) {
                            tl_error("Maximum number of shear stress intervals must be greater than 1.");
                        }
                        break;
                
                    case CARD_NCP:
                        ncti = read_int_field(*io, &subdata);
                        if (ncti < 2) {
                            tl_error("Maximum number of consolidation times must be greater than 1.");
                        }
                        break;
                }
                break;
                
            default:
                break;
        }
    }
    rewind(io->sedfile.fp);
    io_save_line(io, NULL, "", "");


#ifdef _DEBUG
    if (DEBUG == ON) {
        printf("\n nsed: %d nclay: %d nsand: %d nlayers: %d nsfssi: %d ncti: %d nsilt: %d\n",nsed, nclay,nsand,nlayers,nsfssi,ncti,nsilt);
    }
#endif

    mod->nsed = nsed;
    mod->nclay = nclay;
    mod->nsand = nsand;
    mod->nlayers = nlayers;
    mod->nsfssi = nsfssi;
    mod->nconti = ncti;
    mod->nsilt = nsilt;

}


//*******************************************************************************************//
//*******************************************************************************************//
//*******************************************************************************************//

void read_bc_SEDLIB(SMODEL *mod, char *data) {
    
    char line[MAXLINE];           /* the input line */
    char *subdata = NULL;         /* the data after the second card   is read */
    char *subsubdata = NULL;
    
    int i, j, k, imat, istr, ibc = UNSET_INT, iseries = UNSET_INT;
    int ised = 0, iclay = 0, isand = 0, ilayer = 0;
    double thick = 0., dist_total = 0., dist = 0., tempval = 0.;

    // alias
    int nsed = mod->nsed;
    int nclay = mod->nclay;
    int nsand = mod->nsand;
    int nsfssi = mod->nsfssi;
    int ncti = mod->nconti;
    int nlayers = mod->nlayers;
    int nsilt = mod->nsilt;
    int nsedfl = 0;
    SIO *io = mod->io;  
    int super = strlen(io->sup.filename);
    SIO info = *io;
    int nnodes2d = 0, nnodes3d = 0;
    SGRID *grid = mod->grid;

    int ndim = UNSET_INT;
    
    if (mod->flag.SW3_FLOW) {
        ndim = 3;
        nnodes3d = mod->grid->nnodes;
        nnodes2d = mod->grid->nnodes_bed;
    } else {
        ndim = 2;
        nnodes3d = mod->grid->nnodes;
        nnodes2d = mod->grid->nnodes; 
    }
    
    // was code compiled with SEDLIB?
#ifndef _SEDLIB
    printf("**** WARNING: SEDLIB card used in bc file, but code was not built with sedlib ****\n");
    return;
#else
    printf("-- reading sedlib input file ---\n");
#endif
#ifdef _SEDIMENT
    double rho = mod->density;
    if (mod->density < 1.) {
        rho = 1000.;
    }
    
    //**************************************************************************************//
    //**************************************************************************************//
    // allocate adh and sedlib sediment variables
    ssediment_alloc_init(ndim, nnodes3d, nnodes2d, nlayers, nsand, nclay, nsilt, nsed, nsfssi, ncti, &(mod->sed));

    // sanity check
    if (nsed <= 0) {
        tl_error("Sedlib turned on but no grain classes given.  Check that you are using GRAIN instead of CN.");
    }
	for (imat = 0; imat < mod->nmat; imat++)
	for (ised = 0; ised < nsed; ised ++)
		mod->mat[imat].sed[ised].refine_tolerance = 1;


    //**************************************************************************************//
    //**************************************************************************************//
    
    // variables for sanity checking
    int bed_flag[nlayers];
    for (ilayer=0; ilayer<nlayers; ilayer++) {
        bed_flag[ilayer] = NO;
    }

    //**************************************************************************************//
    //**************************************************************************************//
    
    // read sedlib bc file
    while (fgets(line, MAXLINE, mod->io->sedfile.fp) != NULL) {
        io_save_line(io, io->sedfile.fp, io->sedfile.filename, line);
        if (strip_comments(line) <= 1) {    /* ignore empty ('\n') or comment lines */
            continue;
        }
        switch (parse_card(line, &data)) {
               
            // sediment boundary conditions
            case CARD_NB:
                switch (parse_card(data, &subdata)) {
                    case CARD_SED:
                        ibc = get_string_id(info, &subdata, mod->nstring); // string is defined in bc
                        ised = get_transport_id(info, &subdata, mod->sed->nsed);
                        mod->str_values[ibc].sed[ised].bc_flag = BCT_NEU;
                        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
                        mod->str_values[ibc].sed[ised].isigma = iseries;
                        break;
                }
                break;

            case CARD_DB:
                switch (parse_card(data, &subdata)) {
                    case CARD_SED:
                        // add later
                        break;
                }
                break;
        
            // sediment grain definitions
            case CARD_GRAIN:
                switch (parse_card(data, &subdata)) {
                        
                    case CARD_CLA:
#ifdef _DEBUG
                        if (DEBUG) {
                            printf("CLA CARD FOUND\n");
                        }
#endif
                        mod->flag.SEDIMENT = ON;
                        ised = get_transport_id(info, &subdata, mod->nsed);
                        mod->sed->grain[ised].type = CLA;

                        mod->sed->grain[ised].reference_c = read_dbl_field(info, &subdata) * 1.E-6;
                        if (mod->sed->grain[ised].reference_c <= 0) {
                            printf("WARNING: Your reference concentration is < 0, setting it to 1\n");
                        }
                        mod->sed->grain[ised].diameter = read_dbl_field(info, &subdata);
                        mod->sed->grain[ised].specific_gravity = read_dbl_field(info, &subdata);
                        tempval = read_dbl_field(info, &subdata);
                        mod->sed->grain[ised].porosity = (tempval - mod->density*mod->sed->grain[ised].specific_gravity) / (mod->density * (1 - mod->sed->grain[ised].specific_gravity));
                        mod->sed->grain[ised].clay.tau_ce = read_dbl_field(info, &subdata);
                        mod->sed->grain[ised].clay.erode_const = read_dbl_field(info, &subdata);
                        mod->sed->grain[ised].clay.tau_cd = read_dbl_field(info, &subdata);
                        mod->sed->grain[ised].clay.settling_velocity = read_dbl_field(info, &subdata);
                        iclay++;
                        break;
                        
                    case CARD_SND:
#ifdef _DEBUG
                        if (DEBUG) {
                            printf("SND CARD FOUND\n");
                        }
#endif
                        mod->flag.SEDIMENT = ON;
                        ised = get_transport_id(info, &subdata, mod->nsed);
                        mod->sed->grain[ised].type = SND;
                        mod->sed->grain[ised].reference_c = read_dbl_field(info, &subdata) * 1.E-6;
                        if (mod->sed->grain[ised].reference_c <= 0) {
                            printf("WARNING: Your reference concentration is < 0, setting it to 1\n");
                            mod->sed->grain[ised].reference_c = 1.;
                        } 
                        mod->sed->grain[ised].diameter = read_dbl_field(info, &subdata);
                        mod->sed->grain[ised].specific_gravity = read_dbl_field(info, &subdata);
                        mod->sed->grain[ised].porosity = read_dbl_field(info, &subdata);
                        isand++;
                        break;
                        
                }
                break;

            // sediment material properties
            case CARD_MP:
                switch (parse_card(data, &subdata)) {
                        
                    case CARD_NBL:
#ifdef _DEBUG
                        if (DEBUG) {
                            printf("NBL CARD FOUND\n");
                        }
#endif
                        nlayers = read_int_field(info, &subdata);
                        mod->sed->bed_type_flag = read_int_field(info, &subdata);
                        if(mod->sed->bed_type_flag < 0 || mod->sed->bed_type_flag > 1)
                            tl_error("Bed type must be 0 or 1.");
                        break;
                        
                    case CARD_SBA:
#ifdef _DEBUG
                        if (DEBUG) {
                            printf("SBA CARD FOUND\n");
                        }
#endif
                        ilayer = get_bed_layer(info, &subdata, mod->nlayers);
                        thick = read_dbl_field(info, &subdata);
                        if (mod->sed->bed_type_flag == 0 && thick < 0.) {
                            tl_error("Thickness too small on SBA card.");
                        }
                        for (i = 0; i < nnodes2d; i++) {
                            mod->sed->bed_layer[ilayer].thickness[i] = thick;
                        }
                        dist_total = 0.;
                        for (i = 0; i < nsed; i++) {
                            dist = read_dbl_field(info, &subdata);
                            if (dist < 0. || dist > 1.) {
                                tl_error("Distribution invalid");
                            }
                            dist_total += dist;
                            for (j = 0; j < nnodes2d; j++) {
                                mod->sed->bed_layer[ilayer].distribution[i][j] = dist;
                            }
                        }
                        if (dist_total > (1. + NOT_QUITE_SMALL) || dist_total < (1. - NOT_QUITE_SMALL)) {
                            tl_error("Total bed sediment distribution invalid on SBA card.");
                        }
                        bed_flag[ilayer] = YES;
                        break;
                        
/* THis is the case to define materials by material number */
/* Need to rethink this "SBM" for 3D SW Gsavant */						

					case CARD_SBM:

#ifdef _DEBUG
                        if (DEBUG) {
                            printf("SBM CARD FOUND\n");
                        }
#endif
                        ilayer = get_bed_layer(info, &subdata, mod->nlayers);
                        imat = get_material_id(info, &subdata, mod->nmat);
						thick = read_dbl_field(info, &subdata);
						for (i = 0; i < mod->grid->nelems2d;  i++)
							if (grid->elem2d[i].mat == imat)
							{
								mod->sed->bed_layer[ilayer].thickness[mod->grid->elem2d[i].nodes[0]] = thick; 
								mod->sed->bed_layer[ilayer].thickness[mod->grid->elem2d[i].nodes[1]] = thick; 
								mod->sed->bed_layer[ilayer].thickness[mod->grid->elem2d[i].nodes[2]] = thick; 
							}
                        
                        if (mod->sed->bed_type_flag == 0 && thick < 0.)
						{
							tl_error("Thickness is too small on bed for layer\n");
						}



						dist_total = 0.;
						for (i = 0; i < nsed; i ++) {
							dist = read_dbl_field(info, &subdata);
						     if (dist < 0. || dist > 1.) {
                                tl_error("Distribution invalid");
                            }
							 dist_total += dist;
                         for (j = 0; j < mod->grid->nelems2d;  j++)
							if (grid->elem2d[j].mat == imat)
							{
								mod->sed->bed_layer[ilayer].distribution[i][mod->grid->elem2d[j].nodes[0]] = dist; 
								mod->sed->bed_layer[ilayer].distribution[i][mod->grid->elem2d[j].nodes[1]] = dist; 
								mod->sed->bed_layer[ilayer].distribution[i][mod->grid->elem2d[j].nodes[2]] = dist; 
							}
							
						}
						if (dist_total > (1. + NOT_QUITE_SMALL) || dist_total < (1. - NOT_QUITE_SMALL)) {
                            tl_error("Total bed sediment distribution invalid on SBM card.");
                        }
					bed_flag[ilayer] = YES;	
						break;

/* This is the case where there is no bed displacement allowed in AdH Hydro, sedlib doesn't care about NDM. It is purely for bed update */
                   case CARD_NDM:
                      imat = get_material_id(info, &subdata, mod->nmat);
                      mod->mat[imat].sw->bed_disp_flag = NO;
                      printf("No Bed Displacement Specified for %d\n",imat+1);
                   break;


                    case CARD_CBA:
#ifdef _DEBUG
                        if (DEBUG) {
                            printf("CBA CARD FOUND\n");
                        }
#endif
                        ilayer = get_bed_layer(info, &subdata, mod->nlayers);
                        tempval = read_dbl_field(info, &subdata);
                        for (i = 0; i < nnodes2d; i++) {
                            mod->sed->bed_layer[ilayer].bulk_density[i] = tempval;
                        }
                        tempval = read_dbl_field(info, &subdata);
                        for (i = 0; i < nnodes2d; i++) {
                            mod->sed->bed_layer[ilayer].critical_erosion_shear[i] = tempval;
                        }
                        tempval = read_dbl_field(info, &subdata);
                        for (i = 0; i < nnodes2d; i++) {
                            mod->sed->bed_layer[ilayer].erosion_rate_constant[i] = tempval;
                        }
                        tempval = read_dbl_field(info, &subdata);
                        for (i = 0; i < nnodes2d; i++) {
                            mod->sed->bed_layer[ilayer].erosion_rate_exponent[i] = tempval;
                        }
                        bed_flag[ilayer] = YES;
                        mod->sed->bed_layer[ilayer].cbed_flag = YES;
                        break;
                        
                    case CARD_TRT:
#ifdef _DEBUG
                        if (DEBUG) {
                            printf("TRT CARD FOUND\n");
                        }
#endif
                        imat = get_material_id(info, &subdata, mod->nmat);
                        ised = get_transport_id(info, &subdata, mod->nsed);
                        mod->mat[imat].sed[ised].refine_tolerance = read_dbl_field(info, &subdata);
                        if (mod->mat[imat].sed[ised].refine_tolerance < SMALL)
                            tl_error("Refinement tolerance is too small for the precision of the machine.");
                        if (mod->mat[imat].sed[ised].refine_tolerance == 0.0)
                            tl_error("Refinement tolerance must be greater than zero.");
                        break;
                        
                    case CARD_DF:
#ifdef _DEBUG
                        if (DEBUG) {
                            printf("DF CARD FOUND\n");
                        }
#endif
                        imat = get_material_id(info, &subdata, mod->nmat);
                        ised = get_transport_id(info, &subdata, mod->nsed);
                        mod->mat[imat].sed[ised].d_m = read_dbl_field(info, &subdata);
                        break;
                }
                break;
          
            case CARD_SP:
                switch (parse_card(data, &subdata)) {
                    case CARD_HID:
                        mod->sed->hid_fact_flag = read_int_field(info, &subdata);
                        break;
					case CARD_NSE:
						mod->sed->noncoh_ent_flag = read_int_field(info, &subdata);
						break;
					case CARD_NBE:
						mod->sed->noncoh_ble_flag = read_int_field(info, &subdata);
						break;
					case CARD_CSV:
						mod->sed->cohesive_settling_flag = read_int_field(info, &subdata);
						break;
                }
                break;

            default:
                break;
        }
    }
    io_save_line(io, NULL, "", "");
    fclose(io->sedfile.fp);

#ifdef _DEBUG
    if (DEBUG) {
        tl_check_all_pickets(__FILE__, __LINE__);
    }
#endif
#endif
}
