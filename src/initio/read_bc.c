/* reads the boundary condition file first */

#include "global_header.h"

#if defined _ADH_SW2 || _ADH_SW3
#include "fr_defs.h"
#endif

#ifdef _ADH_BREACH
#include "breach_type.h"
#endif

static int DEBUG = OFF;

/**************************************************/

struct MPElem
{
    int flag;
    char name[CHAR_SIZE];
};

/**************************************************/

void read_bc(SMODEL * mod)
{
    
    char line[MAXLINE];           /* the input line */
    char *data = NULL;            /* the data */
    char *subdata = NULL;         /* the data after the second card is read */
    char *subsubdata = NULL;
    CARD card, sub_card;
    
    int ielem1d=0, ielem2d=0;     /* string counts */
    int iend=OFF;                 /* flag for end of file */
    int istring, type;
    SIO *io = mod->io;  // alias
    /* Defaults */
    mod->solver_info.prec_value = 1;
    mod->nblock = 1;
    mod->nalloc_inc = 50;
    /* check that boundary condition file is ok */
    assert(io && io->bc.fp);
#ifdef _DEBUG
    if (DEBUG) debug.readBC = ON;
    printf("-- MYID %d reading boundary condition file: %s \n\n",mod->myid,io->bc.filename);
#endif
    
    /* read for debug cards first *******************************/
    rewind(io->bc.fp);
    while (fgets(line, MAXLINE, io->bc.fp) != NULL) {
        io_save_line(io, io->bc.fp, io->bc.filename, line);
        if (strip_comments(line) <= 1) {
            continue;
        }
        card = parse_card(line, &data);
        switch (card) {
            case CARD_NOTERM:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_NOTERM(mod,data);
                break;
            case CARD_DEBUG:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_DEBUG(mod,data);
                break;
            case CARD_SOUT:     // optional screen output
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_SCREEN_OUTPUT(mod,data);
                break;
            case CARD_FOUT:     // optional file output
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_FILE_OUTPUT(mod,data);
                break;
            case CARD_TEST:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_TEST(mod,data);
                break;
                //case CARD_STOP:
                //    read_bc_STOP(mod,data);
            case CARD_END:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                iend = ON;
                break;
            default:
                break;
        }
        if (iend) break;
    }
    iend = OFF;
    rewind(io->bc.fp);
    io_save_line(io, NULL, "", "");
    
    /* read bc for series ***************************************/
#ifdef _DEBUG
    if (debug.readBC) {
        printf("******************************************************\n");
        printf("DEBUGGING THE BC READ :: WRITING OUT EACH LINE AS READ\n");
        printf("------------------------------------------------------\n");
        printf("reading bc file for series extraction\n");
        printf("------------------------------------------------------\n");
    }
#endif
    
    while (fgets(line, MAXLINE, io->bc.fp) != NULL) {
        io_save_line(io, io->bc.fp, io->bc.filename, line);
        if (strip_comments(line) <= 1) {
            continue;
        }
        card = parse_card(line, &data);
#ifdef _DEBUG
        if (debug.readBC) {
            printf("  %s", line);
        }
#endif
        
        /* CARD LISTING */
        switch (card) {
            case CARD_SERIES:
                sub_card = parse_card(data, &subdata);
                /* what type of series */
                switch (sub_card) {
                    case CARD_AWRITE:  // auto create output
#if _DEBUG
                        if (DEBUG) {
#if _MESSG
                            tag(cstorm_comm);
#else
                            tag();
#endif
                        }
#endif
                        type = OUTPUT_SERIES;
                        read_bc_AWRITE(mod,subdata);
                        break;
                    case CARD_WIND:
#if _DEBUG
                        if (DEBUG) {
#if _MESSG
                            tag(cstorm_comm);
#else
                            tag();
#endif
                        }
#endif
                        type = WIND_SERIES;
                        break;
                    case CARD_WAVE:
#if _DEBUG
                        if (DEBUG) {
#if _MESSG
                            tag(cstorm_comm);
#else
                            tag();
#endif
                        }
#endif
                        type = WAVE_SERIES;
                        break;
                    case CARD_DT:
#if _DEBUG
                        if (DEBUG) {
#if _MESSG
                            tag(cstorm_comm);
#else
                            tag();
#endif
                        }
#endif
                        type = DT_SERIES;
                        break;
                    case CARD_WRITE:
#if _DEBUG
                        if (DEBUG) {
#if _MESSG
                            tag(cstorm_comm);
#else
                            tag();
#endif
                        }
#endif
                        type = OUTPUT_SERIES;
                        break;
#ifdef _ADH_GROUNDWATER
		    case CARD_CONSTI:
		        type = CONSTITUITIVE_SERIES;
			break;
#endif
                    default:
#if _DEBUG
                        if (DEBUG) {
#if _MESSG
                            tag(cstorm_comm);
#else
                            tag();
#endif
                        }
#endif
                        type = ANY_SERIES;
                        break;
                }
                if (sub_card != CARD_AWRITE) read_bc_SERIES(mod, type, subdata);
                break;
            case CARD_END:
                iend = ON;
                break;
            default:
                break;
        }
        
        if (iend) break;
    }
    iend = OFF;
    rewind(io->bc.fp);
    io_save_line(io, NULL, "", "");
    
    
    /* read bc for strings *****************************************/
#ifdef _DEBUG
    if (debug.readBC) {
        printf("------------------------------------------------------\n");
        printf("reading bc file for string extraction\n");
        printf("------------------------------------------------------\n");
    }
#endif
    
    while (fgets(line, MAXLINE, io->bc.fp) != NULL) {
        io_save_line(io, io->bc.fp, io->bc.filename, line);
        if (strip_comments(line) <= 1) {
            continue;
        }
        card = parse_card(line, &data);
        
#ifdef _DEBUG
        if (debug.readBC) {
            printf("  %s", line);
        }
#endif
        
        /* CARD LISTING */
        switch (card) {
            case CARD_NDS:
            case CARD_EGS:
            case CARD_MDS:
            case CARD_MTS:
                //case CARD_FCS:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_STRINGS(mod, card, data, &ielem1d, &ielem2d);
                break;
            case CARD_END:
                iend = ON;
                break;
            default:
                break;
        }
        
        if (iend) break;
    }
    iend = OFF;
    rewind(io->bc.fp);
    io_save_line(io, NULL, "", "");
    
    /* read the rest of bc *****************************************/
#ifdef _DEBUG
    if (debug.readBC) {
        printf("------------------------------------------------------\n");
        printf("now reading rest of bc file\n");
        printf("------------------------------------------------------\n");
    }
#endif
    while (fgets(line, MAXLINE, io->bc.fp) != NULL) {
#ifdef _DEBUG
        if (debug.readBC) {
            printf("  %s", line);
        }
#endif
        io_save_line(io, io->bc.fp, io->bc.filename, line);
        if (strip_comments(line) <= 1) {
            continue;
        }
        card = parse_card(line, &data);
        
        /* CARD LISTING */
        switch (card) {
            case CARD_OP:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_OP(mod,data);
                break;
            case CARD_IP:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_IP(mod,data);
                break;
            case CARD_PC:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_PC(mod,data);
                break;
            case CARD_TC:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_TC(mod,data);
                break;
            case CARD_OFF:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_OFF(mod,data);
                break;
            case CARD_MP:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_MP(mod,data);
                break;
            case CARD_DB:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_DB(mod,data);
                break;
            case CARD_NB:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_NB(mod,data);
                break;
            case CARD_HY:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                switch (parse_card(data, &subdata)){
                    case CARD_INT:
                        istring = get_string_id(*(mod->io),&subdata, mod->nstring);
                        mod->str_values[istring].ol_flow.bc_flag = BCT_HYBRID_INTERNAL;
                        break;
                }
                break;
            case CARD_FR:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_FR(mod,data);
                break;
            case CARD_CN:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_CN(mod,data);
                break;
#ifdef _SEDLIB
            case CARD_SEDLIB:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                read_bc_SEDLIB(mod, data);
                break;
#endif
                /*            case CARD_NSM:
                 printf("WQ Simulations with NSM active\n");
                 read_bc_NSM(mod, data);
                 break;
                 */
            case CARD_END:
#if _DEBUG
                if (DEBUG) {
#if _MESSG
                    tag(cstorm_comm);
#else
                    tag();
#endif
                }
#endif
                iend = ON;
                break;
            default:
                //tl_error("ERROR: No BC END card");
                break;
        }
        
        if (iend) break;
    }
    iend = OFF;

#ifdef _ADH_GROUNDWATER
    if (mod->flag.GW_FLOW) {
      read_bc_GW(mod,data);
    }
#endif
    
    
    iend = OFF;
    io_save_line(io, NULL, "", "");
    
    /* for DB OVL, we must creat an edge string from the node string to apply mass conservation */
    // cjt :: this is bugged
    if (mod->nstring_rt > 0) {
        create_dbovl_elem1d(mod);
    }
    
    /* since face strings for SW3D are read with geo in separate file, assign here */
    /* this is NOT a good way to do this */
    if (mod->flag.SW3_FLOW || mod->flag.NS3_FLOW || mod->flag.GW_FLOW) {
        for (istring = 0; istring < mod->nstring; istring++) {
            if (mod->str_values[istring].string_type == UNSET_INT) {
                mod->str_values[istring].string_type = STR_FACE;
            }
        }
    }
    
    /* re/calculate 1d element jacobians and 2d element connections after bc read */
    if (mod->grid->nelems1d > 0) {
        elem1d_find_elem2d_init(mod->grid);
        elem1d_outward_nrml(mod->grid);
#ifdef _ADH_GROUNDWATER
#ifdef _DEBUG
        // Gajanan gkc May 2020. Must make sure that elem1d[*]->elem2d is the surface
        // 2D element and not the 2D element of a sidewall. This is because we are 
        // defining 3d, 2d, as well as 1d elements in case of GW-DW coupling in a single
        // grid.
        //for(int i=0; i<mod->grid->nelems1d; i++){
        //    int myelem2d = mod->grid->elem1d[i].elem2d;
        //    printf("\nelem1d[%4i].elem2d=%4i",i, myelem2d);
        //    printf("(%4i, %4i, %4i)", mod->grid->elem2d[myelem2d].nodes[0], mod->grid->elem2d[myelem2d].nodes[1], mod->grid->elem2d[myelem2d].nodes[2]);
        //}
        //exit(-1);
#endif
#endif
    }
    
#ifdef _ADH_GROUNDWATER
#ifdef _DEBUG
//    if (mod->grid->nelems2d > 0) {
//        for(int i=0; i<mod->grid->nelems2d; i++){
//            int myelem3d = mod->grid->elem2d[i].id_3d;
//            printf("\nelem2d[%4i](%4i, %4i, %4i).id_3d=%4i ",i, mod->grid->elem2d[i].nodes[0], mod->grid->elem2d[i].nodes[1], mod->grid->elem2d[i].nodes[2], myelem3d);
//            printf("(%4i, %4i, %4i, %4i)", mod->grid->elem3d[myelem3d].nodes[0], mod->grid->elem3d[myelem3d].nodes[1], mod->grid->elem3d[myelem3d].nodes[2], mod->grid->elem3d[myelem3d].nodes[3]);
//        }
//        exit(-1);
//    }
#endif
#endif
    
    if ((mod->flag.SW2_FLOW || mod->flag.SW3_FLOW) && mod->flag.MUC == NO)
        printf("NOTE: MP MUC card not included, manning units constant set = 1.0\n");
    if (mod->flag.CORIOLIS == NO)
        printf("NOTE:  No COR flag. Latitude for all elements assumed 0.0 \n");
    
    /* check material inputs */
    smat_checkall(mod);
    
    /* check series */
    if (mod->series_dt == NULL) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> please include time-series for dt\n");
    }
    sseries_checklist(mod->nseries, mod->series_head);
    
    /* check strings */
    sstr_value_checkall(mod->str_values, mod->nstring);
    
#ifdef _DEBUG 
    printf("-- MYID %d finished boundary condition file: %s \n\n",mod->grid->smpi->myid,io->bc.filename);
#endif
}
