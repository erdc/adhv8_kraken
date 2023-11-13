/* reads the sizes from the boundary condition file and allocates the memory */

#include "global_header.h"

void read_bc_prep(SMODEL * mod)
{
    char line[MAXLINE];           /* the input line */
    char *data = NULL;            /* the data after the first card is read */
    char *subdata = NULL;         /* the data after the second card is read */
    int i, j, k;                  /* loop counters */
    int ie;                       /* loop counter over the elements */
    int itmp;                     /* tmp variable */
    int ierr = 0;
    int count = 0;
    int idum1, idum2, idum3, idum4;
    CARD card;
    
    SIO info = *(mod->io);   // alias
    
    assert(mod->io->bc.fp);
    mod->nweir = 0;
    mod->nflap = 0;
    mod->ntransport = 0;
    while (fgets(line, MAXLINE, mod->io->bc.fp) != NULL) {
        io_save_line(mod->io, mod->io->bc.fp, mod->io->bc.filename, line);
        if (strip_comments(line) <= 1) {    /* ignore empty ('\n') or comment lines */
            continue;
        }
        switch (parse_card(line, &data)) {
                
            case CARD_OP:
                switch (parse_card(data, &subdata)) {
                    case CARD_NS2:
		        mod->flag.NS_FLOW = TRUE;
                        mod->flag.NS2_FLOW = TRUE;
                        if (mod->ntransport > 0) {
                            mod->flag.NS2_TRANSPORT = TRUE;
                        }
                        break;
                    case CARD_NS3:
		        mod->flag.NS_FLOW = TRUE;
                        mod->flag.NS3_FLOW = TRUE;
                        if (mod->ntransport > 0) {
                            mod->flag.NS3_TRANSPORT = TRUE;
                        }
                        break;
                    case CARD_SW2:
                        mod->flag.SW_FLOW = TRUE;
                        mod->flag.SW2_FLOW = TRUE;
                        if (mod->ntransport > 0) {
                            mod->flag.SW2_TRANSPORT = TRUE;
                        }
                        break;
                    case CARD_DIF:
                        mod->flag.SW_FLOW = TRUE;
                        mod->flag.SW2_FLOW = TRUE;
                        mod->flag.DIFFUSIVE_WAVE = TRUE;
                        break;
                        
                    case CARD_SW3:
                        mod->flag.SW_FLOW = TRUE;
                        mod->flag.SW3_FLOW = TRUE;
                        if (mod->ntransport > 0) {
                            mod->flag.SW3_TRANSPORT = TRUE;
                        }
                        break;
                    case CARD_TRN:
                        mod->ntransport = read_int_field(info, &subdata);
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
                        break;
#ifdef _ADH_ICM
                    case CARD_ICM:
                        mod->ICM = TRUE;
                        break;
#endif
                    case CARD_WAV:
                        mod->flag.WAVE = TRUE;
                        break;
                    case CARD_WND:
                        mod->flag.WIND = TRUE;
                        break;
#ifdef _SEDLIB
                    case CARD_SEDLIB:
                        read_bc_SEDLIB_prep(mod, data);
                        break;
#endif
                    case CARD_WQNSM:
                        mod->flag.NSM = TRUE;
                        break;
#ifdef _ADH_GROUNDWATER
		    case CARD_GW:
		        mod->flag.GW_FLOW = TRUE;
#endif			
                    default:
                        break;
                }
                break;
#ifdef _ADH_STRUCTURES
            case CARD_WER:
                mod->nweir = read_int_field(info, &data);
                break;
            case CARD_FLP:
                mod->nflap = read_int_field(info, &data);
                break;
            case CARD_SLUICE:
                mod->nsluice = read_int_field(info, &data);
                break;
#endif
            case CARD_MP:
                switch (parse_card(data, &subdata)) {
                    case CARD_WND:
                        mod->flag.WIND = ON;
                        break;
                    default:
                    break;}
                break;
                
            case CARD_NDS:
                itmp = read_int_field(info, &data);
                itmp = read_int_field(info, &data);
                if (itmp > mod->nstring) {
                    mod->nstring = itmp;
                }
                break;
            case CARD_EGS:
                itmp = read_int_field(info, &data);
                itmp = read_int_field(info, &data);
                itmp = read_int_field(info, &data);
                if (itmp > mod->nstring) {
                    mod->nstring = itmp;
                }
                break;
            case CARD_MDS:
                itmp = read_int_field(info, &data);
                itmp = read_int_field(info, &data);
                itmp = read_int_field(info, &data);
                if (itmp > mod->nstring) {
                    mod->nstring = itmp;
                }
                break;
            case CARD_TRI:
                idum1 = read_int_field(info, &data);    // 3d elem id
                idum1 = read_int_field(info, &data);    // node 1
                idum1 = read_int_field(info, &data);    // node 2
                idum1 = read_int_field(info, &data);    // node 3
                idum1 = read_int_field(info, &data);    // bflag
                itmp = read_int_field(info, &data); // string is for bc read
                if (itmp > mod->nstring) {
                    mod->nstring = itmp;
                }
                break;
            case CARD_QUAD:
                idum1 = read_int_field(info, &data);    // 3d elem id
                idum1 = read_int_field(info, &data);    // node 1
                idum1 = read_int_field(info, &data);    // node 2
                idum1 = read_int_field(info, &data);    // node 3
                idum1 = read_int_field(info, &data);    // node 4
                idum1 = read_int_field(info, &data);    // bflag
                itmp = read_int_field(info, &data); // string is for bc read
                if (itmp > mod->nstring) {
                    mod->nstring = itmp;
                }
                break;
                
            case CARD_MTS:
                itmp = read_int_field(info, &data);
                itmp = read_int_field(info, &data);
                if (itmp > mod->nstring) {
                    mod->nstring = itmp;
                }
                break;
#ifdef _HEC
                //        case CARD_DSS:
                //            dss.nsers++;
#endif
            case CARD_SERIES:
                card = parse_card(data, &subdata);
                itmp = read_int_field(info, &subdata);
                if (itmp > mod->nseries) {
                    mod->nseries = itmp;
                }
                break;
                /* added to alloc for OB OF edge string - dss */
            case CARD_DB:
                switch (parse_card(data, &subdata)) {
                    case CARD_OVL:
                        // this increments nstrings to add a new edge string when NB OVL is used
                        mod->nstring_rt = 1;
                        break;
                    default:
                        break;
                }
                break;
            case CARD_FR:
                switch (parse_card(data, &subdata)) {
                    default:
                        break;
                }
                break;
            case CARD_NB:  /* Added for alloc of tidal series (cjt) */
                switch (parse_card(data, &subdata)) {
                    default:
                        break;
                }
                break;
            case CARD_FOUT:
                switch (parse_card(data, &subdata)) {
                    case CARD_GRID:
                        mod->file_output.grid2dm = ON;
                        break;
                }
            case CARD_END:
                mod->nstring += mod->nstring_rt;
                break;
            default:
                break;
        }
    }
    rewind(mod->io->bc.fp);
    io_save_line(mod->io, NULL, "", "");
    
    /* If tidal BCs are used, re-read bc file and add to time-series (cjt) */
    if (mod->ntides > 0) {
        /* reads and allocates tidal variables */
        //itide = init_read_tides(nseries);
        mod->nseries += mod->ntides;
    }
    
    /* if 3d, read faces file */
    if (mod->flag.SW3_FLOW || mod->flag.NS3_FLOW
#ifdef _ADH_GROUNDWATER
	|| mod->flag.GW_FLOW
#endif 
	) {
//       int super = strlen(mod->io->sup.filename);  
//       open_input_file( &(mod->io->face), "boundary face file", super);
//       if (mod->io->face.fp == NULL) {
//          printf("WARNING :: NO FACE FILE IS FOUND!\n");
//          return;
//        } 
        
        while (fgets(line, MAXLINE, mod->io->face.fp) != NULL) {
            io_save_line(mod->io, mod->io->face.fp, mod->io->face.filename, line);
            if (strip_comments(line) <= 1) {  /* ignore empty ('\n') or comment lines */
                continue;
            }
            switch (parse_card(line, &data)) {
	            case CARD_FCS:
		      itmp = read_int_field(info, &data);
		      itmp = read_int_field(info, &data);
		      itmp = read_int_field(info, &data);
		      itmp = read_int_field(info, &data);
		      if (itmp > mod->nstring) {
			  mod->nstring = itmp;
		      }
		      break;             
                case CARD_TRI:
                    idum1 = read_int_field(info, &data);    // 3d elem id
                    idum1 = read_int_field(info, &data);    // node 1
                    idum1 = read_int_field(info, &data);    // node 2
                    idum1 = read_int_field(info, &data);    // node 3
                    idum1 = read_int_field(info, &data);    // bflag
                    itmp = read_int_field(info, &data); // string is for bc read
                    if (itmp > mod->nstring) {
                        mod->nstring = itmp;
                    }
                    break;
                case CARD_QUAD:
                    idum1 = read_int_field(info, &data);    // 3d elem id
                    idum1 = read_int_field(info, &data);    // node 1
                    idum1 = read_int_field(info, &data);    // node 2
                    idum1 = read_int_field(info, &data);    // node 3
                    idum1 = read_int_field(info, &data);    // node 4
                    idum1 = read_int_field(info, &data);    // bflag
                    itmp = read_int_field(info, &data); // string is for bc read
                    if (itmp > mod->nstring) {
                        mod->nstring = itmp;
                    }
                    break;
                    
                default:
                    break;
            }
        }
        io_save_line(mod->io, NULL, "", "");
	    rewind(mod->io->face.fp);		
    }
    
#ifdef _ADH_GROUNDWATER
#ifdef _DWGW_COUPLING
    if (mod->flag.GW_FLOW == ON && mod->flag.DIFFUSIVE_WAVE == ON) mod->amICoupled = YES;
#endif
#endif
    //printf("nstring: %d\n",mod->nstring);

    /* checks that the variables have been read */
    if (mod->ntransport == UNSET_INT) {
        tl_error("No OP TRN line.");
    }
    if (mod->nseries <= 0) {
        printf("WARNING: No time series were specified.\n");
    }
    if (mod->nstring <= 0) {
        printf("WARNING: No boundary condition strings were specificed.\n");
    }
    
}
