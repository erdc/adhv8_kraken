#include "global_header.h"

#ifdef _DEBUG
static int DEBUG = OFF;
#endif

//*******************************************************************************************//
//*******************************************************************************************//
//*******************************************************************************************//
// This reads the file to count the number of sediment grains for smodel initialization
void read_bc_NSM(SMODEL *mod, char *data) {
    
    char line[MAXLINE];           /* the input line */
    char *subdata = NULL;         /* the data after the second card   is read */
    char *subsubdata = NULL;
    
    int i, j, k, imat, istr, itrn, ibc, iseries;
     
    SIO *io = mod->io;  // alias
	SIO info = *io;
    int super = strlen(io->sup.filename);
	int nstring = mod->nstring; 

    // was code compiled with SEDLIB?
#ifndef _NSM
    printf("**** WARNING: NSM card used in bc file, but code was not built with NSM ****\n");
    return;
#else
    printf("-- reading NSM input file ---\n");
#endif

    mod->flag.NSM = ON;
	mod->nnsm = 0;


    // open and initialize sedlib input file
    init_adh_in_file(&(io->nsmfile));
    sprintf(io->nsmfile.filename, "%s.nsm", io->proj_name);
    open_input_file(&(io->nsmfile), "nsm file", super);
    assert(io->nsmfile.fp);

    //**************************************************************************************//
    //**************************************************************************************//
    // read sedlib file to get counts
    while (fgets(line, MAXLINE, io->nsmfile.fp) != NULL) {
        io_save_line(io, io->nsmfile.fp, io->nsmfile.filename, line);
        if (strip_comments(line) <= 1) {    /* ignore empty ('\n') or comment lines */
            continue;
        }
        switch (parse_card(line, &data)) {
            case CARD_NSM:
			printf("Reading NSM Cards\n");
                switch (parse_card(data, &subdata)) {
                    case CARD_NO2:
                        itrn = get_transport_id(info, &subdata, mod->ntransport);
                        mod->con[itrn].type = NOO;
						mod->NIO2 = 1;
                        mod->con[itrn].property[0] = read_dbl_field(info, &subdata);
						mod->nnsm += 1;
                        break;
					case CARD_NO3:
                        itrn = get_transport_id(info, &subdata, mod->ntransport);
                        mod->con[itrn].type = NOOO;
						mod->NIO3 = 1;
                        mod->con[itrn].property[0] = read_dbl_field(info, &subdata);
						mod->nnsm += 1;
                        break;
                    case CARD_NH4:
                        itrn = get_transport_id(info, &subdata, mod->ntransport);
                        mod->con[itrn].type = NHHHH;
						mod->NIH4 = 1;
                        mod->con[itrn].property[0] = read_dbl_field(info, &subdata);
						mod->nnsm += 1;
                        break;
					case CARD_ON:
                        itrn = get_transport_id(info, &subdata, mod->ntransport);
                        mod->con[itrn].type = OrN;
						mod->ORN = 1;
                        mod->con[itrn].property[0] = read_dbl_field(info, &subdata);
						mod->nnsm += 1;
                        break;
					case CARD_OP:
                        itrn = get_transport_id(info, &subdata, mod->ntransport);
                        mod->con[itrn].type = OP;
						mod->ORP = 1;
                        mod->con[itrn].property[0] = read_dbl_field(info, &subdata);
						mod->nnsm += 1;
                        break;
					case CARD_DP:
                        itrn = get_transport_id(info, &subdata, mod->ntransport);
                        mod->con[itrn].type = DP;
						mod->DIP = 1;
                        mod->con[itrn].property[0] = read_dbl_field(info, &subdata);
						mod->nnsm += 1;
                        break;
					case CARD_PO4:
                        itrn = get_transport_id(info, &subdata, mod->ntransport);
                        mod->con[itrn].type = POOOO;
						mod->PHO4 = 1;
                        mod->con[itrn].property[0] = read_dbl_field(info, &subdata);
						mod->nnsm += 1;
                        break;
					case CARD_ALG:
                        itrn = get_transport_id(info, &subdata, mod->ntransport);
                        mod->con[itrn].type = ALG;
						mod->Alg = 1;
                        mod->con[itrn].property[0] = read_dbl_field(info, &subdata);
						mod->nnsm += 1;
                        break;
					case CARD_CBOD:
                        itrn = get_transport_id(info, &subdata, mod->ntransport);
                        mod->con[itrn].type = CBOD;
						mod->CABOD = 1;
                        mod->con[itrn].property[0] = read_dbl_field(info, &subdata);
						mod->nnsm += 1;
                        break;
					case CARD_DO:
                        itrn = get_transport_id(info, &subdata, mod->ntransport);
						printf("ITRN %d\n",itrn+1);
                        mod->con[itrn].type = DO;
						mod->DIO = 1;
                        mod->con[itrn].property[0] = read_dbl_field(info, &subdata);
						mod->nnsm += 1;
                        break;
				printf("Total NSM Constituents %d\n",mod->nnsm);
                }
                break;
                
            case CARD_MP:
                switch (parse_card(data, &subdata)) {
                    case CARD_TRT:
                       imat = get_material_id(info, &subdata, mod->nmat); 
                       itrn = get_transport_id(info, &subdata, mod->ntransport);
                       mod->mat[imat].trn[itrn].refine_tolerance = read_dbl_field(info, &subdata);
                       if (mod->mat[imat].trn[itrn].refine_tolerance < SMALL)
                       tl_error("Refinement tolerance is too small for the precision of the machine.");
                       if (mod->mat[imat].trn[itrn].refine_tolerance == 0.0)
                       tl_error("Refinement tolerance must be greater than zero.");
                       break;
                        
                    case CARD_DF:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       itrn = get_transport_id(info, &subdata, mod->ntransport);
                       mod->mat[imat].trn[itrn].d_m = read_dbl_field(info, &subdata);
                        break;

					case CARD_ALP0:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->alpha0 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_ALP1:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->alpha1 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_ALP2:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->alpha2 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_ALP3:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->alpha3 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_ALP4:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->alpha4 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_ALP5:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->alpha5 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_ALP6:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->alpha6 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_BET1:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->beta1 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_BET2:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->beta2 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_BET3:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->beta3 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_BET4:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->beta4 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_K1:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->k1 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_K2:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->k2 = read_dbl_field(info, &subdata);
					     break;

				    case CARD_K3:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->k3 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_K4:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->k4 = read_dbl_field(info, &subdata);
                         break;

				    case CARD_GRO:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->mu = read_dbl_field(info, &subdata);
                         break;

                    case CARD_RES:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->rho = read_dbl_field(info, &subdata);
                         break;

                    case CARD_SIG1:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->sigma1 = read_dbl_field(info, &subdata);
                         break;

                    case CARD_SIG2:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->sigma2 = read_dbl_field(info, &subdata);
                         break;

                    case CARD_SIG3:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->sigma3 = read_dbl_field(info, &subdata);
                         break;

                    case CARD_SIG4:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->sigma4 = read_dbl_field(info, &subdata);
                         break;

                    case CARD_SIG5:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->sigma5 = read_dbl_field(info, &subdata);
                         break;

                    case CARD_KL:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->kl = read_dbl_field(info, &subdata);
                         break;

                    case CARD_KN:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->kn = read_dbl_field(info, &subdata);
                         break;

                    case CARD_KP:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->kp = read_dbl_field(info, &subdata);
                         break;

                    case CARD_PN:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->pn = read_dbl_field(info, &subdata);
                         break;

                    case CARD_LAM0:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->lambda0 = read_dbl_field(info, &subdata);
                         break;

                    case CARD_LAM1:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->lambda1 = read_dbl_field(info, &subdata);
                         break;

                    case CARD_LAM2:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->lambda2 = read_dbl_field(info, &subdata);
                         break;

                    case CARD_MDO:
                       imat = get_material_id(info, &subdata, mod->nmat);
                       mod->mat[imat].wnsm->min_do_conc = read_dbl_field(info, &subdata);
                         break;
                }
                break;

				case CARD_NB:
                switch (parse_card(data, &subdata)) {
                    case CARD_TRN:
                        ibc = get_string_id(info, &subdata, mod->nstring); // string is defined in bc
                        itrn = get_transport_id(info, &subdata, mod->ntransport);
                        mod->str_values[ibc].trans[itrn].bc_flag = BCT_NEU;
                        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
                        mod->str_values[ibc].trans[itrn].isigma = iseries;
                        break;
                }

				case CARD_DB:
                switch (parse_card(data, &subdata)) {
                    case CARD_TRN:
                        ibc = get_string_id(info, &subdata, mod->nstring); // string is defined in bc
                        itrn = get_transport_id(info, &subdata, mod->ntransport);
                        mod->str_values[ibc].trans[itrn].bc_flag = BCT_DIR;
                        iseries = sseries_set_type(mod, &subdata, TIME_SERIES);
                        mod->str_values[ibc].trans[itrn].isigma = iseries;
                        break;
                }
                break;

              
            default:
                break;
        }
    }
    rewind(io->nsmfile.fp);
    io_save_line(io, NULL, "", "");

}


